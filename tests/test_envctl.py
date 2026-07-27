from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hic_basic_envctl", ROOT / "scripts" / "envctl.py")
assert SPEC is not None and SPEC.loader is not None
envctl = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(envctl)


def test_environment_intent_is_generated_from_registry():
    registry = envctl.load_registry()

    assert envctl.render_environment(registry) == envctl.ENVIRONMENT_PATH.read_text()
    assert envctl.validate_registry(registry, require_locks=False) == []


def test_default_direct_dependencies_are_constrained_and_documented():
    registry = envctl.load_registry()
    default_dependencies = [dependency for dependency in registry["dependencies"] if envctl.included(dependency)]

    assert default_dependencies
    assert all(dependency.get("reason") for dependency in default_dependencies)
    assert all(
        any(operator in dependency["spec"] for operator in "<>=!~")
        for dependency in default_dependencies
        if dependency["manager"] == "conda"
    )
    pip_lock = envctl.PIP_LOCK_PATH.read_text()
    assert pip_lock.count("--hash=sha256:") == len(envctl.pip_dependencies(registry))
    assert all(
        dependency["requirement"] in pip_lock
        for dependency in envctl.pip_dependencies(registry)
    )
    assert all(
        dependency.get("distribution") and dependency.get("expected_version")
        for dependency in envctl.pip_dependencies(registry)
    )


def test_explicit_lock_records_uses_filename_and_hash(tmp_path):
    lock = tmp_path / "linux-64.lock"
    lock.write_text(
        "# platform: linux-64\n"
        "@EXPLICIT\n"
        "https://conda.example/linux-64/numpy-1.0-py_0.conda#abc123\n"
        "https://conda.example/noarch/six-1.0-py_0.tar.bz2#def456\n"
    )

    assert envctl.explicit_lock_records(lock) == {
        ("numpy-1.0-py_0.conda", "abc123"),
        ("six-1.0-py_0.tar.bz2", "def456"),
    }


def test_installed_conda_records_reads_exact_artifacts(tmp_path):
    metadata_dir = tmp_path / "conda-meta"
    metadata_dir.mkdir()
    (metadata_dir / "numpy.json").write_text(
        json.dumps({"name": "numpy", "fn": "numpy-1.0-py_0.conda", "md5": "abc123"})
    )

    records, payloads = envctl.installed_conda_records(tmp_path)

    assert records == {("numpy-1.0-py_0.conda", "abc123")}
    assert payloads[0]["name"] == "numpy"


def test_atomic_write_text_replaces_complete_file(tmp_path):
    output = tmp_path / "lock.txt"
    output.write_text("old")

    envctl.atomic_write_text(output, "complete\n")

    assert output.read_text() == "complete\n"
    assert not list(tmp_path.glob(".lock.txt.*"))


def test_content_addressed_names_change_with_lock_inputs(tmp_path, monkeypatch):
    registry = tmp_path / "dependencies.yml"
    conda_lock = tmp_path / "conda.lock"
    pip_lock = tmp_path / "pip.txt"
    registry.write_text("registry")
    conda_lock.write_text("conda")
    pip_lock.write_text("pip")
    monkeypatch.setattr(envctl, "REGISTRY_PATH", registry)
    monkeypatch.setattr(envctl, "CONDA_LOCK_PATH", conda_lock)
    monkeypatch.setattr(envctl, "PIP_LOCK_PATH", pip_lock)

    candidate = envctl.environment_name("candidate")
    official = envctl.environment_name("official")
    pip_lock.write_text("changed")

    assert candidate.startswith("hic_basic-candidate-")
    assert official.replace("hic_basic-env-", "") == candidate.replace("hic_basic-candidate-", "")
    assert envctl.environment_name("candidate") != candidate


def test_promoted_environment_name_matches_lock_digest():
    assert envctl.OFFICIAL_NAME_PATH.read_text().strip() == envctl.environment_name("official")
