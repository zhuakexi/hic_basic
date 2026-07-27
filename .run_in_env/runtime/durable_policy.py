#!/usr/bin/env python3
"""Pure durable handoff classification shared by runtime and evals."""

from __future__ import annotations

from typing import Any


QUEUE_POLICIES = {"require_running", "accept_pending"}
TERMINAL_FAILURES = {"FAILED", "CANCELLED", "TIMEOUT"}


def classify_observation(
    observation: dict[str, Any],
    queue_policy: str,
    *,
    deadline_reached: bool,
) -> dict[str, Any] | None:
    """Return a terminal handoff classification, or None while polling may continue."""
    if queue_policy not in QUEUE_POLICIES:
        raise ValueError(f"unsupported queue policy: {queue_policy}")

    state = str(observation.get("state", "UNKNOWN"))
    phase = str(observation.get("phase", "unknown"))
    health = str(observation.get("health", "UNKNOWN"))

    if state == "COMPLETED" and health == "OK":
        return {"execution_status": "success", "failure_kind": None}
    if state == "RUNNING" and phase == "command_started" and health == "OK":
        return {"execution_status": "success", "failure_kind": None}
    if state == "PENDING" and queue_policy == "accept_pending":
        return {"execution_status": "submitted_pending", "failure_kind": None}
    if state in TERMINAL_FAILURES:
        return {"execution_status": "failure", "failure_kind": "execution_failed"}
    if state == "UNKNOWN" or health in {
        "BOOTSTRAP_FAILED",
        "COMMAND_FAILED",
        "SCHEDULER_UNAVAILABLE",
        "UNKNOWN",
    }:
        return {"execution_status": "failure", "failure_kind": "scheduler_unknown"}
    if deadline_reached:
        return {"execution_status": "failure", "failure_kind": "start_not_verified"}
    return None
