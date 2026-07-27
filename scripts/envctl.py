#!/usr/bin/env python3
"""Maintain the repository-owned hic_basic development environment."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from urllib.parse import urlparse

import yaml


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "envs" / "dependencies.yml"
ENVIRONMENT_PATH = ROOT / "environment.yml"
LOCK_DIR = ROOT / "envs" / "locks"
CONDA_LOCK_PATH = LOCK_DIR / "hic-basic-linux-64.lock"
PIP_LOCK_PATH = LOCK_DIR / "pip-linux-64.txt"
OFFICIAL_NAME_PATH = ROOT / "envs" / "official-name.txt"
RUNNER = ROOT / "run_in_env.sh"
DEFAULT_CONDA_LOCK = "/share/home/ychi/mambaforge/bin/conda-lock"


class EnvironmentError(RuntimeError):
    """Raised when an environment invariant is not satisfied."""


def load_registry(path: Path = REGISTRY_PATH) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise EnvironmentError("dependency registry must have schema_version=1")
    dependencies = data.get("dependencies")
    if not isinstance(dependencies, list):
        raise EnvironmentError("dependency registry must contain a dependencies list")
    return data


def included(dependency: dict) -> bool:
    return dependency.get("default", True) is not False


def conda_dependencies(registry: dict) -> list[dict]:
    return [
        dependency
        for dependency in registry["dependencies"]
        if included(dependency) and dependency.get("manager") == "conda"
    ]


def pip_dependencies(registry: dict) -> list[dict]:
    return [
        dependency
        for dependency in registry["dependencies"]
        if included(dependency) and dependency.get("manager") == "pip"
    ]


def render_environment(registry: dict) -> str:
    lines = [
        "# Generated from envs/dependencies.yml by scripts/envctl.py.",
        "# Edit the registry, then run `make env-render`; do not hand-edit this file.",
        f"name: {registry['environment_name']}",
        "channels:",
    ]
    lines.extend(f"  - {channel}" for channel in registry["channels"])
    lines.append("dependencies:")
    lines.extend(f"  - {dependency['spec']}" for dependency in conda_dependencies(registry))
    return "\n".join(lines) + "\n"


def validate_registry(registry: dict, require_locks: bool = True) -> list[str]:
    errors: list[str] = []
    identifiers: set[str] = set()
    conda_specs: set[str] = set()
    for dependency in registry["dependencies"]:
        identifier = dependency.get("id")
        manager = dependency.get("manager")
        if not isinstance(identifier, str) or not identifier:
            errors.append("dependency is missing a non-empty id")
            continue
        if identifier in identifiers:
            errors.append(f"duplicate dependency id: {identifier}")
        identifiers.add(identifier)
        if manager not in {"conda", "pip", "external"}:
            errors.append(f"{identifier}: unsupported manager {manager!r}")
        if manager == "conda":
            spec = dependency.get("spec")
            if not isinstance(spec, str) or not spec:
                errors.append(f"{identifier}: Conda dependency is missing spec")
            elif included(dependency) and not re.search(r"(?:[<>=!~]|\s)", spec):
                errors.append(f"{identifier}: default Conda dependency must have a version constraint")
            elif spec in conda_specs:
                errors.append(f"duplicate Conda spec: {spec}")
            else:
                conda_specs.add(spec)
        if manager == "pip" and not dependency.get("requirement"):
            errors.append(f"{identifier}: pip dependency is missing requirement")
        elif manager == "pip" and included(dependency):
            requirement = dependency["requirement"]
            is_local = requirement.startswith("./") or requirement.startswith("/")
            if not is_local and "==" not in requirement and " @ " not in requirement:
                errors.append(f"{identifier}: default pip dependency must be exactly pinned")
            if not dependency.get("distribution") or not dependency.get("expected_version"):
                errors.append(
                    f"{identifier}: default pip dependency needs distribution and expected_version"
                )
        if not dependency.get("reason"):
            errors.append(f"{identifier}: dependency is missing reason")

    expected = render_environment(registry)
    if not ENVIRONMENT_PATH.is_file() or ENVIRONMENT_PATH.read_text(encoding="utf-8") != expected:
        errors.append("environment.yml is not synchronized; run `make env-render`")

    if require_locks:
        if not CONDA_LOCK_PATH.is_file():
            errors.append(f"missing Conda lock: {CONDA_LOCK_PATH.relative_to(ROOT)}")
        if pip_dependencies(registry):
            if not PIP_LOCK_PATH.is_file():
                errors.append(f"missing pip lock: {PIP_LOCK_PATH.relative_to(ROOT)}")
            else:
                pip_lock = PIP_LOCK_PATH.read_text(encoding="utf-8")
                if pip_lock.count("--hash=sha256:") != len(pip_dependencies(registry)):
                    errors.append("pip lock does not contain one sha256 hash per dependency")
                for dependency in pip_dependencies(registry):
                    requirement = dependency["requirement"]
                    if requirement not in pip_lock:
                        errors.append(f"pip lock is missing {dependency['id']}")
                    elif requirement.startswith("./"):
                        artifact = (ROOT / requirement).resolve()
                        if not artifact.is_file():
                            errors.append(f"vendored pip artifact is missing: {requirement}")
                        else:
                            artifact_hash = hashlib.sha256(artifact.read_bytes()).hexdigest()
                            expected_entry = (
                                f"{requirement} \\\n    --hash=sha256:{artifact_hash}"
                            )
                            if expected_entry not in pip_lock:
                                errors.append(f"vendored pip artifact hash drifted: {requirement}")
        if not OFFICIAL_NAME_PATH.is_file():
            errors.append(f"missing official environment name: {OFFICIAL_NAME_PATH.relative_to(ROOT)}")
        else:
            recorded_name = OFFICIAL_NAME_PATH.read_text(encoding="utf-8").strip()
            if not re.fullmatch(r"hic_basic-env-[0-9a-f]{12}", recorded_name):
                errors.append(f"invalid official environment name: {recorded_name!r}")
    return errors


def lock_digest() -> str:
    digest = hashlib.sha256()
    for path in (CONDA_LOCK_PATH, PIP_LOCK_PATH):
        if not path.is_file():
            raise EnvironmentError(f"required lock input is missing: {path.relative_to(ROOT)}")
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def environment_name(kind: str = "candidate") -> str:
    prefix = "hic_basic-candidate" if kind == "candidate" else "hic_basic-env"
    return f"{prefix}-{lock_digest()[:12]}"


def run(command: list[str], *, cwd: Path = ROOT, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, check=check, text=True)


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(content)
        temporary_path.chmod(0o644)
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def render_pip_lock(registry: dict) -> str:
    """Download no dependencies and hash exactly one artifact per pip exception."""
    cache_root = ROOT / ".env-cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    lines = ["# Non-Conda exceptions. Install only with --require-hashes --no-deps."]
    for dependency in pip_dependencies(registry):
        requirement = dependency["requirement"]
        if requirement.startswith("./") or requirement.startswith("/"):
            artifact = (ROOT / requirement).resolve() if requirement.startswith("./") else Path(requirement).resolve()
            if not artifact.is_relative_to(ROOT.resolve()):
                raise EnvironmentError(f"vendored pip artifact escapes repository: {requirement}")
            if not artifact.is_file():
                raise EnvironmentError(f"vendored pip artifact does not exist: {requirement}")
        else:
            with tempfile.TemporaryDirectory(prefix="pip-lock-", dir=cache_root) as temp_dir:
                run(
                    [
                        sys.executable,
                        "-m",
                        "pip",
                        "download",
                        "--disable-pip-version-check",
                        "--no-deps",
                        "--only-binary=:all:",
                        "--dest",
                        temp_dir,
                        requirement,
                    ]
                )
                artifacts = [path for path in Path(temp_dir).iterdir() if path.is_file()]
                if len(artifacts) != 1:
                    raise EnvironmentError(
                        f"pip requirement resolved to {len(artifacts)} artifacts: {requirement}"
                    )
                artifact_hash = hashlib.sha256(artifacts[0].read_bytes()).hexdigest()
            lines.extend([f"{requirement} \\", f"    --hash=sha256:{artifact_hash}"])
            continue

        artifact_hash = hashlib.sha256(artifact.read_bytes()).hexdigest()
        lines.extend([f"{requirement} \\", f"    --hash=sha256:{artifact_hash}"])
    return "\n".join(lines) + "\n"


def explicit_lock_records(path: Path) -> set[tuple[str, str]]:
    records: set[tuple[str, str]] = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or line == "@EXPLICIT":
            continue
        url, separator, package_hash = line.partition("#")
        if not separator or not package_hash:
            raise EnvironmentError(f"explicit lock entry has no package hash: {line}")
        filename = Path(urlparse(url).path).name
        records.add((filename, package_hash))
    return records


def installed_conda_records(prefix: Path) -> tuple[set[tuple[str, str]], list[dict]]:
    records: set[tuple[str, str]] = set()
    payloads: list[dict] = []
    for record_path in sorted((prefix / "conda-meta").glob("*.json")):
        payload = json.loads(record_path.read_text(encoding="utf-8"))
        filename = payload.get("fn")
        package_hash = payload.get("md5") or payload.get("sha256")
        if filename and package_hash:
            records.add((filename, package_hash))
        payloads.append(payload)
    return records, payloads


def installed_distributions() -> dict[str, dict]:
    result: dict[str, dict] = {}
    for distribution in metadata.distributions():
        name = distribution.metadata.get("Name")
        if not name:
            continue
        normalized = re.sub(r"[-_.]+", "-", name).lower()
        installer = (distribution.read_text("INSTALLER") or "").strip() or None
        direct_url_text = distribution.read_text("direct_url.json")
        direct_url = json.loads(direct_url_text) if direct_url_text else None
        result[normalized] = {
            "name": name,
            "version": distribution.version,
            "installer": installer,
            "direct_url": direct_url,
        }
    return result


def current_snapshot() -> dict:
    prefix = Path(sys.prefix).resolve()
    _, conda_payloads = installed_conda_records(prefix)
    pip_check = subprocess.run(
        [sys.executable, "-m", "pip", "check"],
        text=True,
        capture_output=True,
        check=False,
    )
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=False
    )
    return {
        "schema_version": 1,
        "prefix": str(prefix),
        "python": sys.version,
        "python_executable": sys.executable,
        "git_head": git_head.stdout.strip() if git_head.returncode == 0 else None,
        "conda": {
            payload["name"]: {
                "version": payload.get("version"),
                "build": payload.get("build"),
                "channel": payload.get("channel"),
                "filename": payload.get("fn"),
                "md5": payload.get("md5"),
            }
            for payload in conda_payloads
            if payload.get("name")
        },
        "distributions": installed_distributions(),
        "pip_check": {
            "returncode": pip_check.returncode,
            "stdout": pip_check.stdout.strip(),
            "stderr": pip_check.stderr.strip(),
        },
    }


def verify_current(registry: dict) -> dict:
    errors = validate_registry(registry)
    prefix = Path(sys.prefix).resolve()
    locked = explicit_lock_records(CONDA_LOCK_PATH)
    installed, _ = installed_conda_records(prefix)
    missing = sorted(locked - installed)
    unexpected = sorted(installed - locked)
    if missing:
        errors.append(f"{len(missing)} Conda lock records are missing")
    if unexpected:
        errors.append(f"{len(unexpected)} unexpected Conda records are installed")

    distributions = installed_distributions()
    allowed_pip = {
        re.sub(r"[-_.]+", "-", dependency["distribution"]).lower()
        for dependency in pip_dependencies(registry)
    }
    allowed_pip.add("hic-basic")
    pip_owned = {
        name for name, payload in distributions.items() if payload.get("installer") == "pip"
    }
    unexpected_pip = sorted(pip_owned - allowed_pip)
    if unexpected_pip:
        errors.append("unexpected pip-owned distributions: " + ", ".join(unexpected_pip))

    pip_dependency_issues: list[str] = []
    for dependency in pip_dependencies(registry):
        normalized = re.sub(r"[-_.]+", "-", dependency["distribution"]).lower()
        payload = distributions.get(normalized)
        if payload is None:
            pip_dependency_issues.append(f"{dependency['distribution']}: not installed")
        elif payload.get("installer") != "pip":
            pip_dependency_issues.append(
                f"{dependency['distribution']}: owned by {payload.get('installer')!r}, expected pip"
            )
        elif payload.get("version") != dependency["expected_version"]:
            pip_dependency_issues.append(
                f"{dependency['distribution']}: {payload.get('version')} != {dependency['expected_version']}"
            )
    if pip_dependency_issues:
        errors.append("pip dependency verification failed")

    for dependency in conda_dependencies(registry):
        distribution_name = dependency.get("distribution")
        if not distribution_name:
            continue
        normalized = re.sub(r"[-_.]+", "-", distribution_name).lower()
        payload = distributions.get(normalized)
        if payload and payload.get("installer") == "pip":
            errors.append(f"Conda-owned dependency is shadowed by pip: {distribution_name}")

    editable = distributions.get("hic-basic")
    direct_url = editable.get("direct_url") if editable else None
    editable_path = None
    if direct_url and direct_url.get("dir_info", {}).get("editable"):
        editable_path = Path(urlparse(direct_url["url"]).path).resolve()
    if editable_path != ROOT.resolve():
        errors.append(f"hic_basic is not editable from this checkout: {editable_path}")

    import_failures: dict[str, str] = {}
    for dependency in registry["dependencies"]:
        if not included(dependency):
            continue
        for import_name in dependency.get("imports", []):
            try:
                importlib.import_module(import_name)
            except Exception as exc:  # pragma: no cover - environment-dependent detail
                import_failures[import_name] = f"{type(exc).__name__}: {exc}"
    if import_failures:
        errors.append(f"{len(import_failures)} registered imports failed")

    pip_check = subprocess.run(
        [sys.executable, "-m", "pip", "check"],
        text=True,
        capture_output=True,
        check=False,
    )
    if pip_check.returncode != 0:
        errors.append("pip check failed")

    return {
        "ok": not errors,
        "prefix": str(prefix),
        "lock_digest": lock_digest(),
        "errors": errors,
        "missing_conda_records": missing,
        "unexpected_conda_records": unexpected,
        "unexpected_pip_distributions": unexpected_pip,
        "pip_dependency_issues": pip_dependency_issues,
        "import_failures": import_failures,
        "pip_check": (pip_check.stdout + pip_check.stderr).strip(),
        "editable_project": str(editable_path) if editable_path else None,
    }


def command_render(args: argparse.Namespace) -> int:
    rendered = render_environment(load_registry())
    if args.check:
        if not ENVIRONMENT_PATH.is_file() or ENVIRONMENT_PATH.read_text(encoding="utf-8") != rendered:
            print("environment.yml is not synchronized", file=sys.stderr)
            return 1
        print("environment.yml is synchronized")
        return 0
    ENVIRONMENT_PATH.write_text(rendered, encoding="utf-8")
    print(ENVIRONMENT_PATH)
    return 0


def command_validate(_args: argparse.Namespace) -> int:
    errors = validate_registry(load_registry())
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("dependency registry and locks are consistent")
    return 0


def command_lock(args: argparse.Namespace) -> int:
    registry = load_registry()
    pre_lock_errors = validate_registry(registry, require_locks=False)
    if pre_lock_errors:
        raise EnvironmentError("; ".join(pre_lock_errors))
    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    pip_lock = render_pip_lock(registry)
    conda_lock = Path(args.conda_lock)
    conda_lock_template = "envs/locks/.hic-basic-{platform}.lock.tmp"
    temporary_conda_lock = LOCK_DIR / ".hic-basic-linux-64.lock.tmp"
    command = [
        str(conda_lock),
        "lock",
        "--micromamba",
        "--kind",
        "explicit",
        "--platform",
        "linux-64",
        "--filename-template",
        conda_lock_template,
        "--file",
        str(ENVIRONMENT_PATH),
    ]
    try:
        run(command)
        os.replace(temporary_conda_lock, CONDA_LOCK_PATH)
    finally:
        temporary_conda_lock.unlink(missing_ok=True)
    atomic_write_text(PIP_LOCK_PATH, pip_lock)
    print(CONDA_LOCK_PATH)
    print(PIP_LOCK_PATH)
    return 0


def command_lock_pip(_args: argparse.Namespace) -> int:
    registry = load_registry()
    errors = validate_registry(registry, require_locks=False)
    if errors:
        raise EnvironmentError("; ".join(errors))
    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    atomic_write_text(PIP_LOCK_PATH, render_pip_lock(registry))
    print(PIP_LOCK_PATH)
    return 0


def command_name(args: argparse.Namespace) -> int:
    print(environment_name(args.kind))
    return 0


def command_create(args: argparse.Namespace) -> int:
    errors = validate_registry(load_registry())
    if errors:
        raise EnvironmentError("; ".join(errors))
    name = args.env or environment_name("candidate")
    run([str(RUNNER), "--create-env", "--env", name, "--file", str(CONDA_LOCK_PATH)])
    run(
        [
            str(RUNNER),
            "--env",
            name,
            "--network",
            "internet",
            "python",
            str(Path(__file__).resolve()),
            "install-editable",
        ]
    )
    print(name)
    return 0


def command_install_editable(_args: argparse.Namespace) -> int:
    run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-build-isolation",
            "--require-hashes",
            "--no-deps",
            "--requirement",
            str(PIP_LOCK_PATH),
        ]
    )
    run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-build-isolation",
            "--no-deps",
            "--editable",
            str(ROOT),
        ]
    )
    return 0


def command_snapshot(args: argparse.Namespace) -> int:
    payload = json.dumps(current_snapshot(), indent=2, sort_keys=True) + "\n"
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload, encoding="utf-8")
        print(output)
    else:
        print(payload, end="")
    return 0


def command_verify(_args: argparse.Namespace) -> int:
    result = verify_current(load_registry())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["ok"] else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    render_parser = subparsers.add_parser("render", help="render environment.yml from registry")
    render_parser.add_argument("--check", action="store_true")
    render_parser.set_defaults(func=command_render)

    validate_parser = subparsers.add_parser("validate", help="validate registry and lock assets")
    validate_parser.set_defaults(func=command_validate)

    lock_parser = subparsers.add_parser("lock", help="resolve an explicit linux-64 Conda lock")
    lock_parser.add_argument("--conda-lock", default=DEFAULT_CONDA_LOCK)
    lock_parser.set_defaults(func=command_lock)

    pip_lock_parser = subparsers.add_parser(
        "lock-pip", help="download and hash pip exceptions without dependencies"
    )
    pip_lock_parser.set_defaults(func=command_lock_pip)

    name_parser = subparsers.add_parser("name", help="print the content-addressed environment name")
    name_parser.add_argument("--kind", choices=("candidate", "official"), default="candidate")
    name_parser.set_defaults(func=command_name)

    create_parser = subparsers.add_parser("create", help="create and bootstrap an exact candidate env")
    create_parser.add_argument("--env")
    create_parser.set_defaults(func=command_create)

    install_parser = subparsers.add_parser(
        "install-editable", help="install locked pip exceptions and editable hic_basic"
    )
    install_parser.set_defaults(func=command_install_editable)

    snapshot_parser = subparsers.add_parser("snapshot", help="capture the current environment")
    snapshot_parser.add_argument("--output")
    snapshot_parser.set_defaults(func=command_snapshot)

    verify_parser = subparsers.add_parser("verify", help="verify current env against exact locks")
    verify_parser.set_defaults(func=command_verify)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return args.func(args)
    except (EnvironmentError, OSError, subprocess.CalledProcessError) as exc:
        print(f"envctl: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
