#!/usr/bin/env python3
"""Internal implementation for run_in_env.sh --canonical operations."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import secrets
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from durable_policy import QUEUE_POLICIES, classify_observation


JOB_ACTIONS = {"inspect_only", "adopt", "cancel_before_submit"}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{secrets.token_hex(3)}.tmp")
    with tmp.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def allocate_bundle(project_root: Path) -> tuple[str, Path]:
    root = project_root / ".run_in_env" / "executions"
    root.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S")
    for _ in range(10):
        ref = f"riex-{stamp}-{os.getpid()}-{secrets.token_hex(3)}"
        bundle = root / ref
        try:
            bundle.mkdir(mode=0o700)
        except FileExistsError:
            continue
        return ref, bundle
    raise RuntimeError("could not allocate a unique canonical execution bundle")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_file_hashes(project_root: Path, command: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in command:
        candidate = Path(value)
        candidate = candidate if candidate.is_absolute() else Path.cwd() / candidate
        try:
            resolved = candidate.resolve(strict=True)
            relative = resolved.relative_to(project_root)
        except (FileNotFoundError, ValueError, OSError):
            continue
        if resolved.is_file():
            result[relative.as_posix()] = sha256_file(resolved)
    return result


def load_control(path: Path | None, project_root: Path) -> dict[str, Any]:
    if path is None:
        return {
            "schema_version": 1,
            "owner_ref": "manual",
            "queue_policy": "require_running",
            "previous_job_action": "inspect_only",
            "previous_job_ref": None,
        }
    resolved = path if path.is_absolute() else project_root / path
    resolved = resolved.resolve(strict=True)
    try:
        resolved.relative_to(project_root)
    except ValueError as exc:
        raise ValueError("canonical control file must be inside the project root") from exc
    data = json.loads(resolved.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise ValueError("canonical control must be a schema_version=1 object")
    if data.get("queue_policy") not in QUEUE_POLICIES:
        raise ValueError("canonical control has an invalid queue_policy")
    if data.get("previous_job_action") not in JOB_ACTIONS:
        raise ValueError("canonical control has an invalid previous_job_action")
    if not isinstance(data.get("owner_ref"), str) or not data["owner_ref"]:
        raise ValueError("canonical control has an invalid owner_ref")
    previous_ref = data.get("previous_job_ref")
    if previous_ref is not None and not isinstance(previous_ref, str):
        raise ValueError("canonical control has an invalid previous_job_ref")
    data["control_path"] = resolved.relative_to(project_root).as_posix()
    return data


def extract_command(wrapper_args: list[str]) -> list[str]:
    if "--" in wrapper_args:
        return wrapper_args[wrapper_args.index("--") + 1 :]
    options_with_values = {"--env", "--cpu", "--mem", "--gpu", "--time", "--network"}
    command: list[str] = []
    skip = False
    for value in wrapper_args:
        if skip:
            skip = False
            continue
        if value in options_with_values:
            skip = True
            continue
        if value.startswith("-") and not command:
            continue
        command.append(value)
    return command


def parse_last_json(text: str) -> dict[str, Any]:
    for line in reversed(text.splitlines()):
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    raise ValueError("run_in_env child did not return a JSON object")


def run_child(command: list[str], timeout: float) -> tuple[int, str, str, bool]:
    if timeout <= 0:
        return 124, "", "canonical deadline reached before command start\n", True
    process = subprocess.Popen(
        command,
        cwd=Path.cwd(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=(os.name == "posix"),
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        return process.returncode, stdout, stderr, False
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGTERM)
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.terminate()
            if process.poll() is None:
                process.kill()
        stdout, stderr = process.communicate()
        return 124, stdout, stderr, True


def remaining(deadline: float) -> float:
    return max(0.0, deadline - time.monotonic())


def record_observation(bundle: Path, index: int, observation: dict[str, Any]) -> str:
    path = bundle / "observations" / f"{index:04d}.json"
    write_json(path, observation)
    return str(path)


def prior_submission(project_root: Path, owner_ref: str) -> dict[str, Any] | None:
    candidates: list[tuple[str, dict[str, Any]]] = []
    for path in (project_root / ".run_in_env" / "executions").glob("riex-*/result.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if (
            payload.get("owner_ref") == owner_ref
            and payload.get("action") == "canonical_submit"
            and isinstance(payload.get("job_ref"), str)
        ):
            candidates.append((str(payload.get("finished_at", "")), payload))
    return max(candidates, key=lambda item: item[0])[1] if candidates else None


def finish(
    *,
    bundle: Path,
    execution_ref: str,
    action: str,
    started_at: str,
    started_monotonic: float,
    deadline_seconds: int,
    control: dict[str, Any],
    classification: dict[str, Any],
    child_returncode: int | None,
    job_ref: str | None = None,
    receipt_path: str | None = None,
    observation: dict[str, Any] | None = None,
    supersedes: str | None = None,
    adopted_from: str | None = None,
) -> tuple[dict[str, Any], int]:
    result_path = bundle / "result.json"
    payload: dict[str, Any] = {
        "schema_version": 1,
        "action": action,
        "execution_ref": execution_ref,
        "execution_status": classification["execution_status"],
        "failure_kind": classification.get("failure_kind"),
        "started_at": started_at,
        "finished_at": utc_now(),
        "elapsed_seconds": round(time.monotonic() - started_monotonic, 3),
        "deadline_seconds": deadline_seconds,
        "bundle_path": str(bundle),
        "request_path": str(bundle / "request.json"),
        "result_path": str(result_path),
        "stdout_path": str(bundle / "stdout.log"),
        "stderr_path": str(bundle / "stderr.log"),
        "owner_ref": control["owner_ref"],
        "control_path": control.get("control_path"),
        "queue_policy": control["queue_policy"],
        "previous_job_action": control["previous_job_action"],
        "child_returncode": child_returncode,
        "job_ref": job_ref,
        "receipt_path": receipt_path,
        "observation": observation,
        "supersedes": supersedes,
        "adopted_from": adopted_from,
    }
    write_json(result_path, payload)
    status = payload["execution_status"]
    exit_code = 0 if status in {"success", "submitted_pending"} else (124 if status == "timeout" else 1)
    return payload, exit_code


def execute_run(args: argparse.Namespace, wrapper_args: list[str], control: dict[str, Any]) -> int:
    project_root = args.project_root.resolve()
    execution_ref, bundle = allocate_bundle(project_root)
    started_at = utc_now()
    started_monotonic = time.monotonic()
    request = {
        "schema_version": 1,
        "action": "canonical_run",
        "execution_ref": execution_ref,
        "owner_ref": control["owner_ref"],
        "started_at": started_at,
        "working_directory": str(Path.cwd()),
        "wrapper_args": wrapper_args,
        "command": extract_command(wrapper_args),
        "command_file_sha256": command_file_hashes(project_root, extract_command(wrapper_args)),
        "control_path": control.get("control_path"),
    }
    write_json(bundle / "request.json", request)
    child = [str(args.wrapper), *wrapper_args]
    code, stdout, stderr, timed_out = run_child(child, args.deadline)
    (bundle / "stdout.log").write_text(stdout, encoding="utf-8")
    (bundle / "stderr.log").write_text(stderr, encoding="utf-8")
    if stdout:
        sys.stdout.write(stdout)
    if stderr:
        sys.stderr.write(stderr)
    classification = (
        {"execution_status": "timeout", "failure_kind": "canonical_timeout"}
        if timed_out
        else {
            "execution_status": "success" if code == 0 else "failure",
            "failure_kind": None if code == 0 else "execution_failed",
        }
    )
    payload, exit_code = finish(
        bundle=bundle,
        execution_ref=execution_ref,
        action="canonical_run",
        started_at=started_at,
        started_monotonic=started_monotonic,
        deadline_seconds=args.deadline,
        control=control,
        classification=classification,
        child_returncode=code,
    )
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    return exit_code


def execute_submit(args: argparse.Namespace, wrapper_args: list[str], control: dict[str, Any]) -> int:
    project_root = args.project_root.resolve()
    execution_ref, bundle = allocate_bundle(project_root)
    started_at = utc_now()
    started_monotonic = time.monotonic()
    request = {
        "schema_version": 1,
        "action": "canonical_submit",
        "execution_ref": execution_ref,
        "owner_ref": control["owner_ref"],
        "started_at": started_at,
        "working_directory": str(Path.cwd()),
        "wrapper_args": wrapper_args,
        "command": extract_command(wrapper_args),
        "command_file_sha256": command_file_hashes(project_root, extract_command(wrapper_args)),
        "control_path": control.get("control_path"),
    }
    write_json(bundle / "request.json", request)

    existing_submission = (
        prior_submission(project_root, control["owner_ref"])
        if control.get("control_path")
        else None
    )
    if existing_submission is not None:
        classification = {
            "execution_status": "failure",
            "failure_kind": "duplicate_submission_blocked",
        }
        payload, exit_code = finish(
            bundle=bundle,
            execution_ref=execution_ref,
            action="canonical_submit",
            started_at=started_at,
            started_monotonic=started_monotonic,
            deadline_seconds=args.deadline,
            control=control,
            classification=classification,
            child_returncode=None,
            job_ref=existing_submission["job_ref"],
            receipt_path=existing_submission.get("receipt_path"),
            observation=existing_submission.get("observation"),
        )
        print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
        return exit_code

    supersedes = None
    previous_ref = control.get("previous_job_ref")
    previous_action = control["previous_job_action"]
    if previous_ref and previous_action in {"inspect_only", "adopt"}:
        classification = {"execution_status": "failure", "failure_kind": "active_job_requires_explicit_action"}
        payload, exit_code = finish(
            bundle=bundle, execution_ref=execution_ref, action="canonical_submit",
            started_at=started_at, started_monotonic=started_monotonic,
            deadline_seconds=args.deadline, control=control, classification=classification,
            child_returncode=None,
        )
        print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
        return exit_code
    if previous_ref and previous_action == "cancel_before_submit":
        cancel_cmd = [str(args.wrapper), "--job-cancel", previous_ref]
        code, stdout, stderr, timed_out = run_child(cancel_cmd, remaining(started_monotonic + args.deadline))
        (bundle / "cancel.stdout.log").write_text(stdout, encoding="utf-8")
        (bundle / "cancel.stderr.log").write_text(stderr, encoding="utf-8")
        try:
            cancel = parse_last_json(stdout)
        except ValueError:
            cancel = {}
        if timed_out or code != 0 or cancel.get("final_state") not in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"}:
            classification = {"execution_status": "failure", "failure_kind": "cancel_not_confirmed"}
            payload, exit_code = finish(
                bundle=bundle, execution_ref=execution_ref, action="canonical_submit",
                started_at=started_at, started_monotonic=started_monotonic,
                deadline_seconds=args.deadline, control=control, classification=classification,
                child_returncode=code, observation=cancel,
            )
            print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
            return exit_code
        supersedes = previous_ref

    code, stdout, stderr, timed_out = run_child(
        [str(args.wrapper), *wrapper_args], remaining(started_monotonic + args.deadline)
    )
    (bundle / "stdout.log").write_text(stdout, encoding="utf-8")
    (bundle / "stderr.log").write_text(stderr, encoding="utf-8")
    if timed_out or code != 0:
        classification = {
            "execution_status": "timeout" if timed_out else "failure",
            "failure_kind": "submission_timeout" if timed_out else "submission_failed",
        }
        payload, exit_code = finish(
            bundle=bundle, execution_ref=execution_ref, action="canonical_submit",
            started_at=started_at, started_monotonic=started_monotonic,
            deadline_seconds=args.deadline, control=control, classification=classification,
            child_returncode=code, supersedes=supersedes,
        )
        if stderr:
            sys.stderr.write(stderr)
        print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
        return exit_code

    try:
        receipt = parse_last_json(stdout)
        job_ref = str(receipt["job_ref"])
        receipt_path = str(Path(receipt["bundle_path"]) / "receipt.json")
    except (KeyError, ValueError) as exc:
        raise RuntimeError(f"canonical submit returned an invalid receipt: {exc}") from exc

    observation: dict[str, Any] | None = None
    classification: dict[str, Any] | None = None
    index = 0
    deadline = started_monotonic + args.deadline
    while remaining(deadline) > 0:
        index += 1
        status_code, status_stdout, status_stderr, status_timeout = run_child(
            [str(args.wrapper), "--job-status", job_ref], min(30.0, remaining(deadline))
        )
        if status_timeout:
            break
        try:
            observation = parse_last_json(status_stdout)
        except ValueError:
            observation = {
                "state": "UNKNOWN", "phase": "unknown", "health": "SCHEDULER_UNAVAILABLE",
                "detail": status_stderr, "returncode": status_code,
            }
        record_observation(bundle, index, observation)
        classification = classify_observation(
            observation, control["queue_policy"], deadline_reached=False
        )
        if classification is not None:
            break
        time.sleep(min(args.poll_interval, remaining(deadline)))

    if classification is None:
        observation = observation or {
            "state": "UNKNOWN", "phase": "unknown", "health": "SCHEDULER_UNAVAILABLE"
        }
        classification = classify_observation(
            observation, control["queue_policy"], deadline_reached=True
        )
    payload, exit_code = finish(
        bundle=bundle, execution_ref=execution_ref, action="canonical_submit",
        started_at=started_at, started_monotonic=started_monotonic,
        deadline_seconds=args.deadline, control=control, classification=classification,
        child_returncode=code, job_ref=job_ref, receipt_path=receipt_path,
        observation=observation, supersedes=supersedes,
    )
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    return exit_code


def execute_observe(args: argparse.Namespace, wrapper_args: list[str], control: dict[str, Any]) -> int:
    project_root = args.project_root.resolve()
    execution_ref, bundle = allocate_bundle(project_root)
    started_at = utc_now()
    started_monotonic = time.monotonic()
    previous_ref = control.get("previous_job_ref")
    command_ref = wrapper_args[-1] if wrapper_args else None
    request = {
        "schema_version": 1, "action": f"canonical_{args.action}",
        "execution_ref": execution_ref, "owner_ref": control["owner_ref"],
        "started_at": started_at, "working_directory": str(Path.cwd()),
        "wrapper_args": wrapper_args, "command": [], "command_file_sha256": {},
        "control_path": control.get("control_path"),
    }
    write_json(bundle / "request.json", request)
    if previous_ref and command_ref != previous_ref:
        raise ValueError("canonical job action does not match the control's previous_job_ref")
    code, stdout, stderr, timed_out = run_child([str(args.wrapper), *wrapper_args], args.deadline)
    (bundle / "stdout.log").write_text(stdout, encoding="utf-8")
    (bundle / "stderr.log").write_text(stderr, encoding="utf-8")
    try:
        observation = parse_last_json(stdout)
    except ValueError:
        observation = {"state": "UNKNOWN", "phase": "unknown", "health": "SCHEDULER_UNAVAILABLE"}
    if args.action == "job_cancel":
        ok = not timed_out and code == 0 and observation.get("final_state") in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"}
        classification = {"execution_status": "success" if ok else "failure", "failure_kind": None if ok else "cancel_not_confirmed"}
    else:
        classification = classify_observation(observation, control["queue_policy"], deadline_reached=True)
    payload, exit_code = finish(
        bundle=bundle, execution_ref=execution_ref, action=f"canonical_{args.action}",
        started_at=started_at, started_monotonic=started_monotonic,
        deadline_seconds=args.deadline, control=control, classification=classification,
        child_returncode=code, job_ref=command_ref,
        receipt_path=(str(project_root / ".run_in_env" / "jobs" / str(command_ref) / "receipt.json") if command_ref else None),
        observation=observation,
        adopted_from=(previous_ref if control["previous_job_action"] == "adopt" else None),
    )
    if stderr:
        sys.stderr.write(stderr)
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    return exit_code


def parse_args() -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--wrapper", type=Path, required=True)
    parser.add_argument("--action", choices=["run", "submit", "job_status", "job_collect", "job_cancel"], required=True)
    parser.add_argument("--deadline", type=int, default=240)
    parser.add_argument("--poll-interval", type=float, default=2.0)
    parser.add_argument("--control", type=Path)
    parser.add_argument("wrapper_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    wrapper_args = args.wrapper_args[1:] if args.wrapper_args[:1] == ["--"] else args.wrapper_args
    if not 1 <= args.deadline <= 240:
        parser.error("canonical deadline must be 1-240 seconds")
    if args.poll_interval < 0:
        parser.error("poll interval must be non-negative")
    return args, wrapper_args


def main() -> int:
    try:
        args, wrapper_args = parse_args()
        control = load_control(args.control, args.project_root.resolve())
        if args.action == "run":
            return execute_run(args, wrapper_args, control)
        if args.action == "submit":
            return execute_submit(args, wrapper_args, control)
        return execute_observe(args, wrapper_args, control)
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        print(f"canonical_driver.py: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
