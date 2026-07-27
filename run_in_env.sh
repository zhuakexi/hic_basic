#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_VERSION="1.5.2"
JOB_SCHEMA_VERSION=1
ENV_ROOT="${RUN_IN_ENV_ENV_ROOT:-/path/to/conda/envs}"
DEFAULT_ENV_NAME="${RUN_IN_ENV_DEFAULT_ENV_NAME:-default}"

# Deploy this script with a reviewed site profile and project-local override.
IMAGE="${RUN_IN_ENV_IMAGE:-/path/to/light_base.sif}"
ENV_PREFIX="${RUN_IN_ENV_PREFIX:-${ENV_ROOT}/${DEFAULT_ENV_NAME}}"
INITIAL_ENV_PREFIX="${ENV_PREFIX}"
ENV_PREFIX_WAS_DEFAULT=0
if [[ -z "${RUN_IN_ENV_PREFIX:-}" ]]; then
  ENV_PREFIX_WAS_DEFAULT=1
fi
RUNTIME="${RUN_IN_ENV_RUNTIME:-singularity}"
ENV_MANAGER_BIN="${RUN_IN_ENV_MANAGER_BIN:-}"
BINDS=()
RUNTIME_ARGS=()
NODE03_HOST="${RUN_IN_ENV_NODE03_HOST:-}"
SSH_ARGS=()
INTERNET_CAPABLE_STRATEGIES=()
SLURM_PARTITION="${RUN_IN_ENV_SLURM_PARTITION:-}"
SLURM_OUTPUT_TEMPLATE="${RUN_IN_ENV_SLURM_OUTPUT:-slurm-%j.out}"
JOBS_ROOT="${RUN_IN_ENV_JOBS_DIR:-${PROJECT_ROOT}/.run_in_env/jobs}"
CANONICAL_DRIVER="${RUN_IN_ENV_CANONICAL_DRIVER:-${PROJECT_ROOT}/.run_in_env/runtime/canonical_driver.py}"
CLUSTER_POLICY_PATH="${RUN_IN_ENV_CLUSTER_POLICY:-${PROJECT_ROOT}/run_in_env.cluster.json}"
CLUSTER_HEALTH_ROOT="${RUN_IN_ENV_CLUSTER_HEALTH_DIR:-${PROJECT_ROOT}/.run_in_env/cluster-health}"
SITE_PROFILE_PATH="${RUN_IN_ENV_SITE_PROFILE:-${PROJECT_ROOT}/run_in_env.site.sh}"

GIB_MIB=$((1024))
TIB_MIB=$((1024 * 1024))
NODE03_CPU_THRESHOLD=16
NODE03_MEM_THRESHOLD_MIB=$((48 * GIB_MIB))
SLURM_CPU_THRESHOLD=96
SLURM_MEM_THRESHOLD_MIB=$((256 * GIB_MIB))
SLURM_TIME_THRESHOLD_SECONDS=3600

if [[ -f "${SITE_PROFILE_PATH}" ]]; then
  # shellcheck disable=SC1090
  source "${SITE_PROFILE_PATH}"
fi

if [[ -f "${PROJECT_ROOT}/run_in_env.local.sh" ]]; then
  # shellcheck disable=SC1091
  source "${PROJECT_ROOT}/run_in_env.local.sh"
fi

if [[ "${ENV_PREFIX_WAS_DEFAULT}" == "1" && "${ENV_PREFIX}" == "${INITIAL_ENV_PREFIX}" ]]; then
  ENV_PREFIX="${ENV_ROOT}/${DEFAULT_ENV_NAME}"
fi

ACTION="run"
CANONICAL_MODE=0
CANONICAL_CONTROL=""
ENV_REQUEST=""
ENV_REQUEST_EXPLICIT=0
SPEC_FILE=""
JOB_REF_REQUEST=""
POSITIONAL_ARGS=()
EXTRA_BINDS=()

CPU_REQUEST=4
GPU_REQUEST=0
MEM_REQUEST_RAW="16G"
MEM_REQUEST_MIB=16384
TIME_REQUEST_RAW="-1"
TIME_REQUEST_SECONDS=-1
NETWORK_REQUEST="none"

SELECTED_STRATEGY=""
SELECTED_REASON=""
LOCAL_EXEC_MODE=""

CLUSTER_POLICY_LOADED=0
CLUSTER_POLICY_DIGEST="none"
CLUSTER_PREFLIGHT_TIMEOUT_SECONDS=2
CLUSTER_BOOTSTRAP_EXIT_CODE=78
CLUSTER_QUARANTINE_SECONDS=3600
CLUSTER_HARD_EXCLUDES=()
CLUSTER_DYNAMIC_EXCLUDES=()
CLUSTER_EXCLUDES=()
CLUSTER_REQUIRED_PATHS=()
CLUSTER_REQUIRED_KINDS=()

usage() {
  cat <<'EOF'
Usage:
  ./run_in_env.sh --env <env-name-or-prefix> --self-check [resource-options]
  ./run_in_env.sh --self-check [resource-options]
  ./run_in_env.sh --env <env-name-or-prefix> [resource-options] <command> [args...]
  ./run_in_env.sh [resource-options] <command> [args...]
  ./run_in_env.sh --canonical [--control <attempt-control.json>] [resource-options] -- <command> [args...]
  ./run_in_env.sh --create-env --env <env-name-or-prefix> --file <spec.yaml|lock>
  ./run_in_env.sh --submit [resource-options] -- <command> [args...]
  ./run_in_env.sh --canonical --submit [--control <attempt-control.json>] [resource-options] -- <command> [args...]
  ./run_in_env.sh --job-status <job-ref>
  ./run_in_env.sh --job-collect <job-ref>
  ./run_in_env.sh --job-diagnose <job-ref>
  ./run_in_env.sh --job-cancel <job-ref>

This is the only project-supported compute entrypoint.
The environment, image, binds, remote host, and scheduler partition are deployment configuration.

Resource options:
  --cpu <N>                 Requested CPU count (default: 4)
  --mem <SIZE>              Requested memory (default: 16G; accepts 16, 16G, 16384M, 1T)
  --gpu <N>                 Requested GPU count (default: 0)
  --time <VALUE>            Requested time limit (default: -1; accepts -1, 1h, 90m, 01:00:00)
  --network <MODE>          Network requirement: none (default) or internet

Canonical delivery mode:
  --canonical               Bound this canonical run, or submit-and-start verification, to 240 seconds
  --control <FILE>          Immutable attempt control emitted by the workflow helper

Exploratory commands do not use --canonical and are not covered by this delivery deadline.
Each later canonical retry receives its own independent 240-second deadline.

Execution strategy (automatic):
  - gpu > 0                     -> ssh to the configured remote host
  - cpu >= 96 and mem > 256G
    and time > 1h               -> SLURM run-and-wait on the configured partition
  - cpu >= 16 or mem >= 48G     -> ssh to the configured remote host
  - otherwise                   -> local execution

Examples:
  ./run_in_env.sh python train.py
  ./run_in_env.sh --network internet python fetch_inputs.py
  ./run_in_env.sh --canonical -- python analysis.py
  ./run_in_env.sh --cpu 32 --mem 64G python heavy_job.py
  ./run_in_env.sh --gpu 1 --time 4h python train_gpu.py
  ./run_in_env.sh --cpu 96 --mem 300G --time 2h python full_pipeline.py
  ./run_in_env.sh --submit --cpu 96 --mem 300G --time 4h -- python full_pipeline.py
  ./run_in_env.sh --canonical --submit --control analysis/1.full_pipeline/1.full_pipeline.attempt1.control.json --cpu 96 --mem 300G --time 4h -- python analysis/1.full_pipeline/1.full_pipeline.py
  ./run_in_env.sh --job-status rie-20260717T120000-12345-ab12cd
  ./run_in_env.sh --job-collect rie-20260717T120000-12345-ab12cd
  ./run_in_env.sh --job-diagnose rie-20260717T120000-12345-ab12cd
  ./run_in_env.sh --job-cancel rie-20260717T120000-12345-ab12cd
EOF
}

usage_error() {
  local message="$1"
  echo "run_in_env.sh: ${message}" >&2
  usage >&2
  exit 64
}

legacy_flag_error() {
  local flag="$1"
  echo "run_in_env.sh: legacy scheduling flag is no longer supported: ${flag}" >&2
  echo "run_in_env.sh: declare resources with --cpu/--mem/--gpu/--time and let the wrapper choose the backend" >&2
  exit 64
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "run_in_env.sh: required command not found: ${cmd}" >&2
    return 127
  fi
}

current_runner() {
  if command -v micromamba >/dev/null 2>&1; then
    printf '%s\n' "micromamba"
    return 0
  fi
  if command -v mamba >/dev/null 2>&1; then
    printf '%s\n' "mamba"
    return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    printf '%s\n' "conda"
    return 0
  fi
  return 1
}

absolute_existing_path() {
  local raw_path="$1"
  local abs_path
  local dir_path
  local base_name

  if [[ "${raw_path}" == /* ]]; then
    abs_path="${raw_path}"
  else
    abs_path="${PWD}/${raw_path}"
  fi

  dir_path="$(dirname "${abs_path}")"
  base_name="$(basename "${abs_path}")"

  if [[ ! -e "${dir_path}/${base_name}" ]]; then
    echo "run_in_env.sh: path not found: ${raw_path}" >&2
    return 127
  fi

  (
    cd "${dir_path}" &&
    printf '%s/%s\n' "$(pwd -P)" "${base_name}"
  )
}

resolve_env_prefix() {
  local raw_env="$1"
  if [[ "${raw_env}" == /* ]]; then
    printf '%s\n' "${raw_env}"
    return 0
  fi

  printf '%s/%s\n' "${ENV_ROOT}" "${raw_env}"
}

print_missing_env_hint() {
  echo "run_in_env.sh: create it first with: ./run_in_env.sh --create-env --env ${ENV_PREFIX} --file /path/to/spec.yaml" >&2
}

ensure_env_prefix_visible() {
  local context_label="$1"
  if [[ ! -d "${ENV_PREFIX}" ]]; then
    echo "run_in_env.sh: target env prefix is not visible in ${context_label}: ${ENV_PREFIX}" >&2
    print_missing_env_hint
    return 127
  fi
  if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
    echo "run_in_env.sh: target env python is not executable in ${context_label}: ${ENV_PREFIX}/bin/python" >&2
    return 127
  fi
}

add_bind_if_dir() {
  local dir_path="$1"
  if [[ -n "${dir_path}" && -d "${dir_path}" ]]; then
    EXTRA_BINDS+=("${dir_path}:${dir_path}")
  fi
}

prepare_extra_binds() {
  EXTRA_BINDS=()
  add_bind_if_dir "$(dirname "${ENV_PREFIX}")"
  if [[ -n "${SPEC_FILE}" ]]; then
    add_bind_if_dir "$(dirname "${SPEC_FILE}")"
  fi
}

build_bind_args() {
  local bind
  local -a args=()
  for bind in "${BINDS[@]}"; do
    args+=(-B "${bind}")
  done
  for bind in "${EXTRA_BINDS[@]}"; do
    args+=(-B "${bind}")
  done
  if ((${#args[@]} > 0)); then
    printf '%s\0' "${args[@]}"
  fi
}

runtime_supports_nv() {
  local runtime_base="${RUNTIME##*/}"
  [[ "${runtime_base}" == singularity* || "${runtime_base}" == apptainer* ]]
}

runtime_args_contain() {
  local needle="$1"
  shift || true
  local arg
  for arg in "$@"; do
    if [[ "${arg}" == "${needle}" ]]; then
      return 0
    fi
  done
  return 1
}

build_runtime_args() {
  local add_nv="${1:-0}"
  local -a args=()

  if declare -p RUNTIME_ARGS >/dev/null 2>&1 && ((${#RUNTIME_ARGS[@]} > 0)); then
    args=("${RUNTIME_ARGS[@]}")
  fi

  if [[ "${add_nv}" == "1" ]] && runtime_supports_nv; then
    if ((${#args[@]} == 0)) || ! runtime_args_contain "--nv" "${args[@]}"; then
      args+=(--nv)
    fi
  fi

  if ((${#args[@]} > 0)); then
    printf '%s\0' "${args[@]}"
  fi
}

normalize_mem_to_mib() {
  local raw_value="$1"
  local normalized
  normalized="$(echo "${raw_value}" | tr '[:upper:]' '[:lower:]')"
  normalized="${normalized// /}"
  if [[ -z "${normalized}" ]]; then
    return 1
  fi

  python3 - <<'PY' "${normalized}"
import math
import re
import sys

value = sys.argv[1]
m = re.fullmatch(r'([0-9]+(?:\.[0-9]+)?)([a-z]*)', value)
if not m:
    sys.exit(1)

amount = float(m.group(1))
unit = m.group(2) or "g"
factors = {
    "m": 1,
    "mb": 1,
    "mib": 1,
    "g": 1024,
    "gb": 1024,
    "gib": 1024,
    "t": 1024 * 1024,
    "tb": 1024 * 1024,
    "tib": 1024 * 1024,
}
if unit not in factors:
    sys.exit(1)

result = math.ceil(amount * factors[unit])
if result <= 0:
    sys.exit(1)
print(result)
PY
}

format_mem_mib() {
  local mib="$1"
  if ((mib % TIB_MIB == 0)); then
    printf '%sT\n' "$((mib / TIB_MIB))"
    return 0
  fi
  if ((mib % GIB_MIB == 0)); then
    printf '%sG\n' "$((mib / GIB_MIB))"
    return 0
  fi
  printf '%sM\n' "${mib}"
}

normalize_time_to_seconds() {
  local raw_value="$1"
  local normalized
  normalized="$(echo "${raw_value}" | tr '[:upper:]' '[:lower:]')"
  normalized="${normalized// /}"

  if [[ "${normalized}" == "-1" ]]; then
    printf '%s\n' "-1"
    return 0
  fi

  python3 - <<'PY' "${normalized}"
import re
import sys

value = sys.argv[1]

if re.fullmatch(r'[0-9]+', value):
    seconds = int(value)
    if seconds <= 0:
        sys.exit(1)
    print(seconds)
    sys.exit(0)

m = re.fullmatch(r'([0-9]+(?:\.[0-9]+)?)([smhd])', value)
if m:
    amount = float(m.group(1))
    unit = m.group(2)
    factors = {"s": 1, "m": 60, "h": 3600, "d": 86400}
    seconds = int(amount * factors[unit])
    if amount * factors[unit] != seconds:
        seconds += 1
    if seconds <= 0:
        sys.exit(1)
    print(seconds)
    sys.exit(0)

m = re.fullmatch(r'([0-9]+):([0-9]{2})', value)
if m:
    hours = int(m.group(1))
    minutes = int(m.group(2))
    if minutes >= 60:
        sys.exit(1)
    seconds = hours * 3600 + minutes * 60
    if seconds <= 0:
        sys.exit(1)
    print(seconds)
    sys.exit(0)

m = re.fullmatch(r'([0-9]+):([0-9]{2}):([0-9]{2})', value)
if m:
    hours = int(m.group(1))
    minutes = int(m.group(2))
    seconds_part = int(m.group(3))
    if minutes >= 60 or seconds_part >= 60:
        sys.exit(1)
    seconds = hours * 3600 + minutes * 60 + seconds_part
    if seconds <= 0:
        sys.exit(1)
    print(seconds)
    sys.exit(0)

m = re.fullmatch(r'([0-9]+)-([0-9]{2}):([0-9]{2}):([0-9]{2})', value)
if m:
    days = int(m.group(1))
    hours = int(m.group(2))
    minutes = int(m.group(3))
    seconds_part = int(m.group(4))
    if hours >= 24 or minutes >= 60 or seconds_part >= 60:
        sys.exit(1)
    seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds_part
    if seconds <= 0:
        sys.exit(1)
    print(seconds)
    sys.exit(0)

sys.exit(1)
PY
}

format_time_seconds() {
  local seconds="$1"
  if [[ "${seconds}" == "-1" ]]; then
    printf '%s\n' "-1"
    return 0
  fi
  if ((seconds % 86400 == 0)); then
    printf '%sd\n' "$((seconds / 86400))"
    return 0
  fi
  if ((seconds % 3600 == 0)); then
    printf '%sh\n' "$((seconds / 3600))"
    return 0
  fi
  if ((seconds % 60 == 0)); then
    printf '%sm\n' "$((seconds / 60))"
    return 0
  fi
  printf '%ss\n' "${seconds}"
}

format_slurm_time() {
  local seconds="$1"
  local days hours minutes secs

  days=$((seconds / 86400))
  seconds=$((seconds % 86400))
  hours=$((seconds / 3600))
  seconds=$((seconds % 3600))
  minutes=$((seconds / 60))
  secs=$((seconds % 60))

  if ((days > 0)); then
    printf '%d-%02d:%02d:%02d\n' "${days}" "${hours}" "${minutes}" "${secs}"
    return 0
  fi
  printf '%02d:%02d:%02d\n' "${hours}" "${minutes}" "${secs}"
}

validate_positive_int() {
  local label="$1"
  local raw_value="$2"
  if [[ ! "${raw_value}" =~ ^[0-9]+$ ]] || ((raw_value <= 0)); then
    usage_error "${label} must be a positive integer"
  fi
}

validate_nonnegative_int() {
  local label="$1"
  local raw_value="$2"
  if [[ ! "${raw_value}" =~ ^[0-9]+$ ]]; then
    usage_error "${label} must be a non-negative integer"
  fi
}

normalize_resource_requests() {
  local parsed_mem
  local parsed_time
  local network_strategy

  validate_positive_int "--cpu" "${CPU_REQUEST}"
  validate_nonnegative_int "--gpu" "${GPU_REQUEST}"

  if ! parsed_mem="$(normalize_mem_to_mib "${MEM_REQUEST_RAW}")"; then
    usage_error "--mem has unsupported format: ${MEM_REQUEST_RAW}"
  fi
  MEM_REQUEST_MIB="${parsed_mem}"
  MEM_REQUEST_RAW="$(format_mem_mib "${MEM_REQUEST_MIB}")"

  if ! parsed_time="$(normalize_time_to_seconds "${TIME_REQUEST_RAW}")"; then
    usage_error "--time has unsupported format: ${TIME_REQUEST_RAW}"
  fi
  TIME_REQUEST_SECONDS="${parsed_time}"
  TIME_REQUEST_RAW="$(format_time_seconds "${TIME_REQUEST_SECONDS}")"

  if [[ "${NETWORK_REQUEST}" != "none" && "${NETWORK_REQUEST}" != "internet" ]]; then
    usage_error "--network must be one of: none, internet"
  fi
  if ((${#INTERNET_CAPABLE_STRATEGIES[@]} > 0)); then
    for network_strategy in "${INTERNET_CAPABLE_STRATEGIES[@]}"; do
      case "${network_strategy}" in
        local|node03|slurm) ;;
        *)
          echo "run_in_env.sh: invalid INTERNET_CAPABLE_STRATEGIES entry: ${network_strategy}" >&2
          return 78
          ;;
      esac
    done
  fi
}

select_strategy() {
  SELECTED_REASON="default-local"
  SELECTED_STRATEGY="local"

  if ((GPU_REQUEST > 0)); then
    SELECTED_STRATEGY="node03"
    SELECTED_REASON="gpu-requested"
    return 0
  fi

  if ((CPU_REQUEST >= SLURM_CPU_THRESHOLD)) && ((MEM_REQUEST_MIB > SLURM_MEM_THRESHOLD_MIB)) && ((TIME_REQUEST_SECONDS > SLURM_TIME_THRESHOLD_SECONDS)); then
    SELECTED_STRATEGY="slurm"
    SELECTED_REASON="large-cpu-mem-long-run"
    return 0
  fi

  if ((CPU_REQUEST >= NODE03_CPU_THRESHOLD)) && ((MEM_REQUEST_MIB >= NODE03_MEM_THRESHOLD_MIB)); then
    SELECTED_STRATEGY="node03"
    SELECTED_REASON="cpu+mem-threshold"
    return 0
  fi
  if ((CPU_REQUEST >= NODE03_CPU_THRESHOLD)); then
    SELECTED_STRATEGY="node03"
    SELECTED_REASON="cpu-threshold"
    return 0
  fi
  if ((MEM_REQUEST_MIB >= NODE03_MEM_THRESHOLD_MIB)); then
    SELECTED_STRATEGY="node03"
    SELECTED_REASON="mem-threshold"
    return 0
  fi
}

strategy_supports_internet() {
  local requested_strategy="$1" candidate
  if ((${#INTERNET_CAPABLE_STRATEGIES[@]} > 0)); then
    for candidate in "${INTERNET_CAPABLE_STRATEGIES[@]}"; do
      if [[ "${candidate}" == "${requested_strategy}" ]]; then
        return 0
      fi
    done
  fi
  return 1
}

validate_network_capability() {
  local strategy="$1"
  if [[ "${NETWORK_REQUEST}" != "internet" ]]; then
    return 0
  fi
  if strategy_supports_internet "${strategy}"; then
    return 0
  fi
  echo "run_in_env.sh: strategy ${strategy} has no configured direct internet access" >&2
  echo "run_in_env.sh: fetch network inputs with a small local --network internet command, store them on shared storage, then run compute with --network none" >&2
  return 69
}

print_execution_plan() {
  echo "run_in_env.sh: strategy=${SELECTED_STRATEGY} reason=${SELECTED_REASON} cpu=${CPU_REQUEST} mem=${MEM_REQUEST_RAW} gpu=${GPU_REQUEST} time=${TIME_REQUEST_RAW} network=${NETWORK_REQUEST}"
}

resolve_local_context() {
  if [[ -n "${SINGULARITY_NAME:-}" || -n "${APPTAINER_NAME:-}" ]]; then
    ensure_env_prefix_visible "the current context"
    if ! EXEC_RUNNER="$(current_runner)"; then
      echo "run_in_env.sh: inside the container, but no supported runner is available (micromamba/mamba/conda)" >&2
      return 127
    fi
    LOCAL_EXEC_MODE="in_container"
    return 0
  fi

  require_command "${RUNTIME}"
  if [[ ! -f "${IMAGE}" ]]; then
    echo "run_in_env.sh: container image not found: ${IMAGE}" >&2
    return 127
  fi
  ensure_env_prefix_visible "the host"
  prepare_extra_binds
  EXEC_RUNNER="micromamba@target_image"
  LOCAL_EXEC_MODE="host_via_target_image"
}

resolve_target_image_context() {
  local require_existing_env="${1:-1}"
  require_command "${RUNTIME}"
  if [[ ! -f "${IMAGE}" ]]; then
    echo "run_in_env.sh: container image not found: ${IMAGE}" >&2
    return 127
  fi
  if [[ "${require_existing_env}" == "1" ]]; then
    ensure_env_prefix_visible "the host"
  fi
  prepare_extra_binds
}

build_target_image_exec_cmd() {
  local mode="$1"
  local add_nv="$2"
  shift 2

  local -a bind_args=()
  local -a runtime_args=()
  local -a cmd=()
  local arg

  while IFS= read -r -d '' arg; do
    bind_args+=("${arg}")
  done < <(build_bind_args)

  while IFS= read -r -d '' arg; do
    runtime_args+=("${arg}")
  done < <(build_runtime_args "${add_nv}")

  cmd=("${RUNTIME}" exec)
  if ((${#runtime_args[@]} > 0)); then
    cmd+=("${runtime_args[@]}")
  fi
  if ((${#bind_args[@]} > 0)); then
    cmd+=("${bind_args[@]}")
  fi

  if [[ "${mode}" == "run" ]]; then
    cmd+=("${IMAGE}" micromamba run -p "${ENV_PREFIX}" "$@")
  elif [[ "${mode}" == "create" ]]; then
    cmd+=("${IMAGE}" micromamba create -p "${ENV_PREFIX}" --yes --file "${SPEC_FILE}")
  else
    echo "run_in_env.sh: unsupported target-image mode: ${mode}" >&2
    return 64
  fi

  printf '%s\0' "${cmd[@]}"
}

exec_in_current_container() {
  exec "${EXEC_RUNNER}" run -p "${ENV_PREFIX}" "$@"
}

create_in_current_container() {
  exec "${EXEC_RUNNER}" create -p "${ENV_PREFIX}" --yes --file "${SPEC_FILE}"
}

exec_via_target_image() {
  local -a cmd=()
  local arg

  while IFS= read -r -d '' arg; do
    cmd+=("${arg}")
  done < <(build_target_image_exec_cmd "run" 0 "$@")

  exec "${cmd[@]}"
}

create_via_target_image() {
  local -a cmd=()
  local arg

  while IFS= read -r -d '' arg; do
    cmd+=("${arg}")
  done < <(build_target_image_exec_cmd "create" 0)

  exec "${cmd[@]}"
}

exec_via_node03() {
  local -a remote_cmd=()
  local arg
  local remote_payload
  local add_nv=0

  if [[ -z "${NODE03_HOST}" ]]; then
    echo "run_in_env.sh: remote execution selected but NODE03_HOST is not configured" >&2
    return 78
  fi
  require_command ssh
  resolve_target_image_context 1

  if ((GPU_REQUEST > 0)); then
    add_nv=1
  fi

  while IFS= read -r -d '' arg; do
    remote_cmd+=("${arg}")
  done < <(build_target_image_exec_cmd "run" "${add_nv}" "$@")

  remote_payload="$(printf '%q ' "${remote_cmd[@]}")"
  remote_shell_cmd="$(printf 'bash -lc %q' "${remote_payload}")"

  if ((${#SSH_ARGS[@]} > 0)); then
    exec ssh "${SSH_ARGS[@]}" "${NODE03_HOST}" "${remote_shell_cmd}"
  fi
  exec ssh "${NODE03_HOST}" "${remote_shell_cmd}"
}

load_cluster_policy() {
  local policy_output line key value kind path
  local -A seen_excludes=()

  if ((CLUSTER_POLICY_LOADED)); then
    return 0
  fi

  if ! policy_output="$(python3 - "${CLUSTER_POLICY_PATH}" "${CLUSTER_HEALTH_ROOT}" <<'PY'
import hashlib
import json
import re
import sys
import time
from pathlib import Path

policy_path = Path(sys.argv[1])
health_root = Path(sys.argv[2])
policy = {
    "schema_version": 1,
    "hard_exclude_nodes": [],
    "required_paths": [],
    "preflight_timeout_seconds": 2,
    "bootstrap_exit_code": 78,
    "dynamic_quarantine_seconds": 3600,
}
digest = "none"

if policy_path.exists():
    raw = policy_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    loaded = json.loads(raw)
    if not isinstance(loaded, dict) or loaded.get("schema_version") != 1:
        raise SystemExit("cluster policy must be an object with schema_version=1")
    policy.update(loaded)

node_re = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
hard = policy.get("hard_exclude_nodes", [])
if not isinstance(hard, list) or any(not isinstance(node, str) or not node_re.fullmatch(node) for node in hard):
    raise SystemExit("cluster policy hard_exclude_nodes is invalid")

required = policy.get("required_paths", [])
if not isinstance(required, list):
    raise SystemExit("cluster policy required_paths must be a list")
normalized_required = []
for item in required:
    if not isinstance(item, dict):
        raise SystemExit("cluster policy required path must be an object")
    path = item.get("path")
    kind = item.get("kind", "path")
    if kind not in {"path", "file", "directory"}:
        raise SystemExit("cluster policy required path kind is invalid")
    if not isinstance(path, str) or not path.startswith("/") or any(char in path for char in "\n\r\t"):
        raise SystemExit("cluster policy required path is invalid")
    normalized_required.append((kind, path))

timeout = policy.get("preflight_timeout_seconds", 2)
bootstrap_exit = policy.get("bootstrap_exit_code", 78)
quarantine_seconds = policy.get("dynamic_quarantine_seconds", 3600)
if not isinstance(timeout, int) or not 1 <= timeout <= 10:
    raise SystemExit("cluster policy preflight timeout must be 1-10 seconds")
if not isinstance(bootstrap_exit, int) or not 1 <= bootstrap_exit <= 255:
    raise SystemExit("cluster policy bootstrap exit code must be 1-255")
if not isinstance(quarantine_seconds, int) or not 60 <= quarantine_seconds <= 604800:
    raise SystemExit("cluster policy quarantine duration must be 60-604800 seconds")

dynamic = []
nodes_dir = health_root / "nodes"
if nodes_dir.is_dir():
    now = time.time()
    for record_path in sorted(nodes_dir.glob("*.json")):
        try:
            record = json.loads(record_path.read_text(encoding="utf-8"))
            node = record.get("node")
            if (
                record.get("schema_version") == 1
                and record.get("status") == "quarantined"
                and isinstance(node, str)
                and node_re.fullmatch(node)
                and float(record.get("expires_at_epoch", 0)) > now
            ):
                dynamic.append(node)
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            continue

print("DIGEST\t" + digest)
print("TIMEOUT\t" + str(timeout))
print("BOOTSTRAP_EXIT\t" + str(bootstrap_exit))
print("QUARANTINE_SECONDS\t" + str(quarantine_seconds))
for node in hard:
    print("HARD_EXCLUDE\t" + node)
for node in dynamic:
    print("DYNAMIC_EXCLUDE\t" + node)
for kind, path in normalized_required:
    print("REQUIRED\t" + kind + "\t" + path)
PY
  )"; then
    echo "run_in_env.sh: invalid cluster policy: ${CLUSTER_POLICY_PATH}" >&2
    return 78
  fi

  while IFS= read -r line; do
    IFS=$'\t' read -r key value path <<< "${line}"
    case "${key}" in
      DIGEST) CLUSTER_POLICY_DIGEST="${value}" ;;
      TIMEOUT) CLUSTER_PREFLIGHT_TIMEOUT_SECONDS="${value}" ;;
      BOOTSTRAP_EXIT) CLUSTER_BOOTSTRAP_EXIT_CODE="${value}" ;;
      QUARANTINE_SECONDS) CLUSTER_QUARANTINE_SECONDS="${value}" ;;
      HARD_EXCLUDE) CLUSTER_HARD_EXCLUDES+=("${value}") ;;
      DYNAMIC_EXCLUDE) CLUSTER_DYNAMIC_EXCLUDES+=("${value}") ;;
      REQUIRED)
        kind="${value}"
        CLUSTER_REQUIRED_KINDS+=("${kind}")
        CLUSTER_REQUIRED_PATHS+=("${path}")
        ;;
    esac
  done <<< "${policy_output}"

  if ((${#CLUSTER_HARD_EXCLUDES[@]} > 0)); then
    for value in "${CLUSTER_HARD_EXCLUDES[@]}"; do
      if [[ -n "${value}" && -z "${seen_excludes[${value}]:-}" ]]; then
        seen_excludes["${value}"]=1
        CLUSTER_EXCLUDES+=("${value}")
      fi
    done
  fi
  if ((${#CLUSTER_DYNAMIC_EXCLUDES[@]} > 0)); then
    for value in "${CLUSTER_DYNAMIC_EXCLUDES[@]}"; do
      if [[ -n "${value}" && -z "${seen_excludes[${value}]:-}" ]]; then
        seen_excludes["${value}"]=1
        CLUSTER_EXCLUDES+=("${value}")
      fi
    done
  fi
  CLUSTER_POLICY_LOADED=1
}

cluster_exclude_csv() {
  local IFS=,
  if ((${#CLUSTER_EXCLUDES[@]} > 0)); then
    printf '%s' "${CLUSTER_EXCLUDES[*]}"
  fi
}

slurm_preflight() {
  local topology excludes index kind path

  load_cluster_policy
  require_command timeout
  require_command scontrol
  require_command sinfo

  if [[ ! -d "${PWD}" || ! -r "${PWD}" || ! -x "${PWD}" ]]; then
    echo "run_in_env.sh: SLURM preflight failed: working directory is unavailable: ${PWD}" >&2
    return 78
  fi
  for ((index = 0; index < ${#CLUSTER_REQUIRED_PATHS[@]}; index++)); do
    kind="${CLUSTER_REQUIRED_KINDS[index]}"
    path="${CLUSTER_REQUIRED_PATHS[index]}"
    case "${kind}" in
      directory)
        if [[ ! -d "${path}" || ! -r "${path}" || ! -x "${path}" ]]; then
          echo "run_in_env.sh: SLURM preflight failed: required directory unavailable: ${path}" >&2
          return 78
        fi
        ;;
      file)
        if [[ ! -f "${path}" || ! -r "${path}" ]]; then
          echo "run_in_env.sh: SLURM preflight failed: required file unavailable: ${path}" >&2
          return 78
        fi
        ;;
      *)
        if [[ ! -e "${path}" ]]; then
          echo "run_in_env.sh: SLURM preflight failed: required path unavailable: ${path}" >&2
          return 78
        fi
        ;;
    esac
  done

  if ! timeout "${CLUSTER_PREFLIGHT_TIMEOUT_SECONDS}s" scontrol ping >/dev/null 2>&1; then
    echo "run_in_env.sh: SLURM preflight failed: controller did not respond within ${CLUSTER_PREFLIGHT_TIMEOUT_SECONDS}s" >&2
    return 69
  fi

  if ! topology="$(timeout "${CLUSTER_PREFLIGHT_TIMEOUT_SECONDS}s" sinfo -h -N -p "${SLURM_PARTITION}" -o '%N|%T|%c|%m' 2>/dev/null)"; then
    echo "run_in_env.sh: SLURM preflight failed: could not inspect partitions ${SLURM_PARTITION}" >&2
    return 69
  fi
  excludes="$(cluster_exclude_csv)"

  if ! RUN_IN_ENV_TOPOLOGY="${topology}" \
    RUN_IN_ENV_EXCLUDES="${excludes}" \
    RUN_IN_ENV_CPU="${CPU_REQUEST}" \
    RUN_IN_ENV_MEM_MIB="${MEM_REQUEST_MIB}" \
    python3 - <<'PY'
import os
import re

excluded = {node for node in os.environ["RUN_IN_ENV_EXCLUDES"].split(",") if node}
cpu = int(os.environ["RUN_IN_ENV_CPU"])
memory = int(os.environ["RUN_IN_ENV_MEM_MIB"])
eligible = []
fitting = []
bad_states = ("down", "drain", "fail", "maint", "power", "not_responding", "unknown")

for raw in os.environ["RUN_IN_ENV_TOPOLOGY"].splitlines():
    fields = raw.split("|")
    if len(fields) != 4:
        continue
    node, state, node_cpu, node_memory = fields
    state = re.sub(r"[^a-z_]", "", state.lower())
    if not node or node in excluded or any(state.startswith(prefix) for prefix in bad_states):
        continue
    try:
        node_cpu = int(node_cpu)
        node_memory = int(node_memory)
    except ValueError:
        continue
    eligible.append(node)
    if node_cpu >= cpu and node_memory >= memory:
        fitting.append(node)

if not eligible:
    raise SystemExit("no scheduler-eligible nodes remain after cluster exclusions")
if not fitting:
    raise SystemExit(
        f"no eligible node can satisfy cpu={cpu} mem_mib={memory}; "
        f"eligible={','.join(eligible)}"
    )
PY
  then
    echo "run_in_env.sh: SLURM preflight failed for partitions ${SLURM_PARTITION}" >&2
    return 69
  fi
}

require_slurm() {
  if [[ -z "${SLURM_PARTITION}" ]]; then
    echo "run_in_env.sh: SLURM execution selected but SLURM_PARTITION is not configured" >&2
    return 78
  fi
  require_command sbatch
  require_command squeue
  require_command scontrol
  require_command sinfo
  ensure_env_prefix_visible "the current context"
  slurm_preflight
}

slurm_build_sbatch_args() {
  local -a args=()

  args+=("--partition=${SLURM_PARTITION}")
  if ((${#CLUSTER_EXCLUDES[@]} > 0)); then
    args+=("--exclude=$(cluster_exclude_csv)")
  fi
  args+=("--cpus-per-task=${CPU_REQUEST}")
  args+=("--mem=${MEM_REQUEST_RAW}")
  args+=("--output=${SLURM_OUTPUT_TEMPLATE}")
  if ((TIME_REQUEST_SECONDS > 0)); then
    args+=("--time=$(format_slurm_time "${TIME_REQUEST_SECONDS}")")
  fi

  printf '%s\n' "${args[@]}"
}

slurm_build_payload() {
  local -a cmd=("$@")
  printf '%q ' "${cmd[@]}"
}

slurm_generate_header() {
  local index kind path
  echo "#!/bin/bash"
  echo "# Generated by run_in_env.sh on $(date)"
  echo "# Environment: ${ENV_PREFIX}"
  echo ""

  local line
  while IFS= read -r line; do
    if [[ -n "${line}" ]]; then
      echo "#SBATCH ${line}"
    fi
  done < <(slurm_build_sbatch_args)

  echo ""
  echo "set -euo pipefail"
  if [[ -n "${ENV_MANAGER_BIN}" ]]; then
    printf 'export PATH=%q:"${PATH}"\n' "${ENV_MANAGER_BIN}"
  fi
  echo 'if command -v micromamba >/dev/null 2>&1; then'
  echo '  _ENV_RUNNER=micromamba'
  echo 'elif command -v mamba >/dev/null 2>&1; then'
  echo '  _ENV_RUNNER=mamba'
  echo 'elif command -v conda >/dev/null 2>&1; then'
  echo '  _ENV_RUNNER=conda'
  echo 'else'
  echo '  echo "run_in_env.sh: no conda/mamba/micromamba found on compute node" >&2'
  echo '  exit 127'
  echo 'fi'
  for ((index = 0; index < ${#CLUSTER_REQUIRED_PATHS[@]}; index++)); do
    kind="${CLUSTER_REQUIRED_KINDS[index]}"
    path="${CLUSTER_REQUIRED_PATHS[index]}"
    case "${kind}" in
      directory)
        printf 'if [[ ! -d %q || ! -r %q || ! -x %q ]]; then echo %q >&2; exit %q; fi\n' \
          "${path}" "${path}" "${path}" "run_in_env.sh: bootstrap failed: required directory unavailable: ${path}" "${CLUSTER_BOOTSTRAP_EXIT_CODE}"
        ;;
      file)
        printf 'if [[ ! -f %q || ! -r %q ]]; then echo %q >&2; exit %q; fi\n' \
          "${path}" "${path}" "run_in_env.sh: bootstrap failed: required file unavailable: ${path}" "${CLUSTER_BOOTSTRAP_EXIT_CODE}"
        ;;
      *)
        printf 'if [[ ! -e %q ]]; then echo %q >&2; exit %q; fi\n' \
          "${path}" "run_in_env.sh: bootstrap failed: required path unavailable: ${path}" "${CLUSTER_BOOTSTRAP_EXIT_CODE}"
        ;;
    esac
  done
  printf 'if [[ ! -x %q ]]; then echo %q >&2; exit %q; fi\n' \
    "${ENV_PREFIX}/bin/python" "run_in_env.sh: bootstrap failed: environment unavailable: ${ENV_PREFIX}/bin/python" "${CLUSTER_BOOTSTRAP_EXIT_CODE}"
  echo ""
}

slurm_run_cmd() {
  local tmp_script
  local job_id
  local default_out

  require_slurm

  tmp_script=$(mktemp "${PROJECT_ROOT}/.run_in_env_slurm_XXXXXX.slurm")
  {
    slurm_generate_header
    printf '%s run -p %s %s\n' '"${_ENV_RUNNER}"' "${ENV_PREFIX}" "$(slurm_build_payload "$@")"
    echo ""
  } > "${tmp_script}"

  echo "run_in_env.sh: submitting SLURM job..."
  job_id="$(sbatch "${tmp_script}" 2>&1 | grep -o '[0-9]*$' || true)"
  rm -f "${tmp_script}"

  if [[ -z "${job_id}" ]]; then
    echo "run_in_env.sh: failed to submit SLURM job" >&2
    return 1
  fi

  echo "run_in_env.sh: job submitted (ID: ${job_id}), waiting for completion..."
  while squeue -j "${job_id}" -h 2>/dev/null | grep -q .; do
    sleep 5
  done

  echo "run_in_env.sh: job ${job_id} completed"

  default_out="${SLURM_OUTPUT_TEMPLATE//%j/${job_id}}"
  if [[ -f "${PROJECT_ROOT}/${default_out}" && -s "${PROJECT_ROOT}/${default_out}" ]]; then
    echo "--- SLURM job output ---"
    cat "${PROJECT_ROOT}/${default_out}"
    echo "--- end output ---"
  fi
}

utc_now() {
  date -u +'%Y-%m-%dT%H:%M:%SZ'
}

validate_job_ref() {
  local job_ref="$1"
  if [[ ! "${job_ref}" =~ ^rie-[0-9]{8}T[0-9]{6}-[0-9]+-[0-9a-f]{6}$ ]]; then
    echo "run_in_env.sh: invalid job reference: ${job_ref}" >&2
    return 64
  fi
}

create_job_bundle() {
  local stamp random job_ref bundle attempt
  mkdir -p "${JOBS_ROOT}"
  chmod 700 "${JOBS_ROOT}" 2>/dev/null || true
  stamp="$(date -u +'%Y%m%dT%H%M%S')"

  for attempt in 1 2 3 4 5; do
    random="$(printf '%06x' "$(((RANDOM << 1 ^ RANDOM ^ attempt) & 0xffffff))")"
    job_ref="rie-${stamp}-$$-${random}"
    bundle="${JOBS_ROOT}/${job_ref}"
    if (umask 077 && mkdir "${bundle}") 2>/dev/null; then
      JOB_REF="${job_ref}"
      JOB_BUNDLE="${bundle}"
      return 0
    fi
  done

  echo "run_in_env.sh: could not allocate a unique durable job record" >&2
  return 1
}

write_preflight_json() {
  local target="${JOB_BUNDLE}/preflight.json"
  local tmp="${target}.tmp.$$"
  local required="" index

  for ((index = 0; index < ${#CLUSTER_REQUIRED_PATHS[@]}; index++)); do
    required+="${CLUSTER_REQUIRED_KINDS[index]}|${CLUSTER_REQUIRED_PATHS[index]}"$'\n'
  done

  RUN_IN_ENV_PREFLIGHT_AT="$(utc_now)" \
  RUN_IN_ENV_POLICY_PATH="${CLUSTER_POLICY_PATH}" \
  RUN_IN_ENV_POLICY_DIGEST="${CLUSTER_POLICY_DIGEST}" \
  RUN_IN_ENV_EXCLUDES="$(cluster_exclude_csv)" \
  RUN_IN_ENV_REQUIRED="${required}" \
  RUN_IN_ENV_PARTITION="${SLURM_PARTITION}" \
  RUN_IN_ENV_CPU="${CPU_REQUEST}" \
  RUN_IN_ENV_MEM="${MEM_REQUEST_RAW}" \
  python3 - "${tmp}" <<'PY'
import json
import os
import sys

required = []
for raw in os.environ["RUN_IN_ENV_REQUIRED"].splitlines():
    if not raw:
        continue
    kind, path = raw.split("|", 1)
    required.append({"kind": kind, "path": path})

data = {
    "schema_version": 1,
    "checked_at": os.environ["RUN_IN_ENV_PREFLIGHT_AT"],
    "policy_path": os.environ["RUN_IN_ENV_POLICY_PATH"],
    "policy_digest": os.environ["RUN_IN_ENV_POLICY_DIGEST"],
    "excluded_nodes": [node for node in os.environ["RUN_IN_ENV_EXCLUDES"].split(",") if node],
    "required_paths": required,
    "partition": os.environ["RUN_IN_ENV_PARTITION"],
    "cpu": int(os.environ["RUN_IN_ENV_CPU"]),
    "mem": os.environ["RUN_IN_ENV_MEM"],
    "result": "passed",
}
with open(sys.argv[1], "x", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
PY
  if ! ln "${tmp}" "${target}"; then
    rm -f "${tmp}"
    echo "run_in_env.sh: refusing to overwrite preflight record: ${target}" >&2
    return 1
  fi
  rm -f "${tmp}"
}

write_job_json() {
  local target="$1"
  local mode="$2"
  shift 2
  local tmp="${target}.tmp.$$"

  RUN_IN_ENV_JSON_MODE="${mode}" \
  RUN_IN_ENV_JSON_REF="${JOB_REF}" \
  RUN_IN_ENV_JSON_BUNDLE="${JOB_BUNDLE}" \
  RUN_IN_ENV_JSON_SCRIPT="${JOB_SCRIPT}" \
  RUN_IN_ENV_JSON_STDOUT="${JOB_STDOUT}" \
  RUN_IN_ENV_JSON_STDERR="${JOB_STDERR}" \
  RUN_IN_ENV_JSON_CWD="${PWD}" \
  RUN_IN_ENV_JSON_ENV_PREFIX="${ENV_PREFIX}" \
  RUN_IN_ENV_JSON_CPU="${CPU_REQUEST}" \
  RUN_IN_ENV_JSON_MEM="${MEM_REQUEST_RAW}" \
  RUN_IN_ENV_JSON_GPU="${GPU_REQUEST}" \
  RUN_IN_ENV_JSON_TIME="${TIME_REQUEST_RAW}" \
  RUN_IN_ENV_JSON_NETWORK="${NETWORK_REQUEST}" \
  RUN_IN_ENV_JSON_PARTITION="${SLURM_PARTITION}" \
  RUN_IN_ENV_JSON_WRAPPER_VERSION="${WRAPPER_VERSION}" \
  RUN_IN_ENV_JSON_SCHEMA_VERSION="${JOB_SCHEMA_VERSION}" \
  python3 - "${tmp}" "$@" <<'PY'
import json
import os
import sys

target = sys.argv[1]
mode = os.environ["RUN_IN_ENV_JSON_MODE"]
data = {
    "schema_version": int(os.environ["RUN_IN_ENV_JSON_SCHEMA_VERSION"]),
    "wrapper_version": os.environ["RUN_IN_ENV_JSON_WRAPPER_VERSION"],
    "job_ref": os.environ["RUN_IN_ENV_JSON_REF"],
    "backend": "slurm",
    "bundle_path": os.environ["RUN_IN_ENV_JSON_BUNDLE"],
    "submission_script": os.environ["RUN_IN_ENV_JSON_SCRIPT"],
    "working_directory": os.environ["RUN_IN_ENV_JSON_CWD"],
    "env_prefix": os.environ["RUN_IN_ENV_JSON_ENV_PREFIX"],
    "resources": {
        "cpu": int(os.environ["RUN_IN_ENV_JSON_CPU"]),
        "mem": os.environ["RUN_IN_ENV_JSON_MEM"],
        "gpu": int(os.environ["RUN_IN_ENV_JSON_GPU"]),
        "time": os.environ["RUN_IN_ENV_JSON_TIME"],
        "network": os.environ["RUN_IN_ENV_JSON_NETWORK"],
        "partition": os.environ["RUN_IN_ENV_JSON_PARTITION"],
    },
    "stdout_path": os.environ["RUN_IN_ENV_JSON_STDOUT"],
    "stderr_path": os.environ["RUN_IN_ENV_JSON_STDERR"],
}
if mode == "request":
    data.update({"created_at": sys.argv[2], "argv": sys.argv[3:]})
elif mode == "receipt":
    data.update({"submitted_at": sys.argv[2], "slurm_job_id": sys.argv[3]})
elif mode == "submission_error":
    data.update({"failed_at": sys.argv[2], "error": sys.argv[3]})
else:
    raise SystemExit("unsupported JSON mode")

with open(target, "x", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
PY
  if ! ln "${tmp}" "${target}"; then
    rm -f "${tmp}"
    echo "run_in_env.sh: refusing to overwrite durable record: ${target}" >&2
    return 1
  fi
  rm -f "${tmp}"
}

generate_durable_script() {
  local target="$1"
  shift
  local tmp="${target}.tmp.$$"
  local index kind path
  (umask 077; set -o noclobber; {
    echo "#!/bin/bash"
    echo "# Generated by run_in_env.sh ${WRAPPER_VERSION}"
    echo "set -euo pipefail"
    printf '_RIE_EVENT_FILE=%q\n' "${JOB_BUNDLE}/lifecycle.jsonl"
    printf '_RIE_HEALTH_ROOT=%q\n' "${CLUSTER_HEALTH_ROOT}"
    printf '_RIE_JOB_REF=%q\n' "${JOB_REF}"
    printf '_RIE_QUARANTINE_SECONDS=%q\n' "${CLUSTER_QUARANTINE_SECONDS}"
    printf '_RIE_BOOTSTRAP_EXIT=%q\n' "${CLUSTER_BOOTSTRAP_EXIT_CODE}"
    cat <<'SH'
_rie_event() {
  local event="$1" detail="${2:-}" exit_code="${3:-}"
  python3 - "${_RIE_EVENT_FILE}" "${event}" "${detail}" "${exit_code}" <<'PY'
import datetime
import json
import os
import socket
import sys

record = {
    "at": datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
    "event": sys.argv[2],
    "node": socket.gethostname(),
}
if sys.argv[3]:
    record["detail"] = sys.argv[3]
if sys.argv[4]:
    record["exit_code"] = int(sys.argv[4])
line = json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
fd = os.open(sys.argv[1], os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
try:
    os.write(fd, line.encode())
    os.fsync(fd)
finally:
    os.close(fd)
PY
}

_rie_quarantine_node() {
  local failure_class="$1" detail="$2"
  python3 - "${_RIE_HEALTH_ROOT}" "${_RIE_JOB_REF}" "${failure_class}" "${detail}" "${_RIE_QUARANTINE_SECONDS}" <<'PY'
import datetime
import json
import os
import socket
import sys
import time
from pathlib import Path

root = Path(sys.argv[1]) / "nodes"
root.mkdir(parents=True, exist_ok=True)
node = socket.gethostname()
now = time.time()
record = {
    "schema_version": 1,
    "status": "quarantined",
    "node": node,
    "failure_class": sys.argv[3],
    "detail": sys.argv[4],
    "job_ref": sys.argv[2],
    "observed_at": datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
    "expires_at_epoch": now + int(sys.argv[5]),
}
target = root / f"{node}.json"
tmp = root / f".{node}.{os.getpid()}.tmp"
with open(tmp, "x", encoding="utf-8") as handle:
    json.dump(record, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
os.replace(tmp, target)
PY
}

_rie_bootstrap_fail() {
  local failure_class="$1" detail="$2"
  _rie_event bootstrap_failed "${failure_class}:${detail}" "${_RIE_BOOTSTRAP_EXIT}"
  _rie_quarantine_node "${failure_class}" "${detail}" || true
  echo "run_in_env.sh: bootstrap failed on $(hostname): ${failure_class}: ${detail}" >&2
  exit "${_RIE_BOOTSTRAP_EXIT}"
}

_rie_event allocated
_rie_event bootstrap_started
SH
    for ((index = 0; index < ${#CLUSTER_REQUIRED_PATHS[@]}; index++)); do
      kind="${CLUSTER_REQUIRED_KINDS[index]}"
      path="${CLUSTER_REQUIRED_PATHS[index]}"
      case "${kind}" in
        directory)
          printf 'if [[ ! -d %q || ! -r %q || ! -x %q ]]; then _rie_bootstrap_fail required_directory_unavailable %q; fi\n' "${path}" "${path}" "${path}" "${path}"
          ;;
        file)
          printf 'if [[ ! -f %q || ! -r %q ]]; then _rie_bootstrap_fail required_file_unavailable %q; fi\n' "${path}" "${path}" "${path}"
          ;;
        *)
          printf 'if [[ ! -e %q ]]; then _rie_bootstrap_fail required_path_unavailable %q; fi\n' "${path}" "${path}"
          ;;
      esac
    done
    printf 'if [[ ! -x %q ]]; then _rie_bootstrap_fail environment_unavailable %q; fi\n' "${ENV_PREFIX}/bin/python" "${ENV_PREFIX}/bin/python"
    if [[ -n "${ENV_MANAGER_BIN}" ]]; then
      printf 'export PATH=%q:"${PATH}"\n' "${ENV_MANAGER_BIN}"
    fi
    echo 'if command -v micromamba >/dev/null 2>&1; then'
    echo '  _ENV_RUNNER=micromamba'
    echo 'elif command -v mamba >/dev/null 2>&1; then'
    echo '  _ENV_RUNNER=mamba'
    echo 'elif command -v conda >/dev/null 2>&1; then'
    echo '  _ENV_RUNNER=conda'
    echo 'else'
    echo '  _rie_bootstrap_fail environment_runner_missing micromamba_mamba_conda'
    echo 'fi'
    echo '_rie_event bootstrap_ok'
    echo '_rie_event command_started'
    echo 'set +e'
    printf '%s run -p %q %s\n' '"${_ENV_RUNNER}"' "${ENV_PREFIX}" "$(slurm_build_payload "$@")"
    echo '_RIE_COMMAND_EXIT=$?'
    echo 'set -e'
    echo '_rie_event command_exited "" "${_RIE_COMMAND_EXIT}"'
    echo 'exit "${_RIE_COMMAND_EXIT}"'
  } > "${tmp}")
  chmod 700 "${tmp}"
  if ! ln "${tmp}" "${target}"; then
    rm -f "${tmp}"
    echo "run_in_env.sh: refusing to overwrite submission script: ${target}" >&2
    return 1
  fi
  rm -f "${tmp}"
}

submit_durable_job() {
  local -a sbatch_args=()
  local submitted_at submit_output job_id error_at

  if [[ $# -eq 0 ]]; then
    usage_error "--submit requires a command"
  fi
  normalize_resource_requests
  if ((GPU_REQUEST > 0)); then
    echo "run_in_env.sh: --submit v1 does not support GPU requests on this SLURM backend" >&2
    return 64
  fi
  validate_network_capability slurm
  require_slurm

  create_job_bundle
  JOB_SCRIPT="${JOB_BUNDLE}/submit.slurm"
  JOB_STDOUT="${JOB_BUNDLE}/stdout.log"
  JOB_STDERR="${JOB_BUNDLE}/stderr.log"
  write_preflight_json
  generate_durable_script "${JOB_SCRIPT}" "$@"
  write_job_json "${JOB_BUNDLE}/request.json" request "$(utc_now)" "$@"

  sbatch_args=(
    --parsable
    "--partition=${SLURM_PARTITION}"
    "--cpus-per-task=${CPU_REQUEST}"
    "--mem=${MEM_REQUEST_RAW}"
    "--chdir=${PWD}"
    "--output=${JOB_STDOUT}"
    "--error=${JOB_STDERR}"
  )
  if ((${#CLUSTER_EXCLUDES[@]} > 0)); then
    sbatch_args+=("--exclude=$(cluster_exclude_csv)")
  fi
  if ((TIME_REQUEST_SECONDS > 0)); then
    sbatch_args+=("--time=$(format_slurm_time "${TIME_REQUEST_SECONDS}")")
  fi

  if ! submit_output="$(sbatch "${sbatch_args[@]}" "${JOB_SCRIPT}" 2>&1)"; then
    error_at="$(utc_now)"
    write_job_json "${JOB_BUNDLE}/submission_error.json" submission_error "${error_at}" "${submit_output}"
    echo "run_in_env.sh: SLURM submission failed; record=${JOB_BUNDLE}" >&2
    echo "${submit_output}" >&2
    return 1
  fi
  if [[ ! "${submit_output}" =~ ^([0-9]+)(\;[^[:space:];]+)?$ ]]; then
    error_at="$(utc_now)"
    write_job_json "${JOB_BUNDLE}/submission_error.json" submission_error "${error_at}" "unexpected sbatch --parsable output: ${submit_output}"
    echo "run_in_env.sh: unexpected sbatch --parsable output: ${submit_output}" >&2
    return 1
  fi

  job_id="${BASH_REMATCH[1]}"
  submitted_at="$(utc_now)"
  write_job_json "${JOB_BUNDLE}/receipt.json" receipt "${submitted_at}" "${job_id}"
  python3 - "${JOB_BUNDLE}/receipt.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    print(json.dumps(json.load(handle), sort_keys=True, separators=(",", ":")))
PY
}

load_job_receipt() {
  local job_ref="$1"
  local receipt
  validate_job_ref "${job_ref}"
  JOB_REF="${job_ref}"
  JOB_BUNDLE="${JOBS_ROOT}/${JOB_REF}"
  receipt="${JOB_BUNDLE}/receipt.json"
  if [[ ! -f "${receipt}" ]]; then
    echo "run_in_env.sh: submitted job receipt not found: ${receipt}" >&2
    return 66
  fi

  mapfile -t JOB_FIELDS < <(python3 - "${receipt}" "${JOB_REF}" "${JOB_SCHEMA_VERSION}" "${JOB_BUNDLE}" <<'PY'
import json
import os
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    data = json.load(handle)
bundle = sys.argv[4]
if data.get("schema_version") != int(sys.argv[3]) or data.get("job_ref") != sys.argv[2]:
    raise SystemExit("malformed or unsupported receipt")
if data.get("backend") != "slurm" or data.get("bundle_path") != bundle:
    raise SystemExit("malformed receipt identity")
expected = {
    "stdout_path": os.path.join(bundle, "stdout.log"),
    "stderr_path": os.path.join(bundle, "stderr.log"),
    "submission_script": os.path.join(bundle, "submit.slurm"),
}
for key in ("slurm_job_id", "stdout_path", "stderr_path", "submission_script"):
    value = data.get(key)
    if not isinstance(value, str) or not value or "\n" in value:
        raise SystemExit("malformed receipt field: " + key)
    if key in expected and value != expected[key]:
        raise SystemExit("unexpected receipt path: " + key)
    print(value)
PY
  ) || {
    echo "run_in_env.sh: malformed job receipt: ${receipt}" >&2
    return 65
  }
  if [[ ${#JOB_FIELDS[@]} -ne 4 || ! "${JOB_FIELDS[0]}" =~ ^[0-9]+$ ]]; then
    echo "run_in_env.sh: malformed job receipt: ${receipt}" >&2
    return 65
  fi
  JOB_ID="${JOB_FIELDS[0]}"
  JOB_STDOUT="${JOB_FIELDS[1]}"
  JOB_STDERR="${JOB_FIELDS[2]}"
  JOB_SCRIPT="${JOB_FIELDS[3]}"
}

canonical_job_state() {
  local state="${1%%+*}"
  state="${state%% *}"
  case "${state}" in
    PENDING|CONFIGURING|REQUEUED|REQUEUE_FED|REQUEUE_HOLD|RESIZING) echo PENDING ;;
    RUNNING|COMPLETING|SUSPENDED|SIGNALING|STAGE_OUT) echo RUNNING ;;
    COMPLETED) echo COMPLETED ;;
    CANCELLED) echo CANCELLED ;;
    TIMEOUT) echo TIMEOUT ;;
    FAILED|BOOT_FAIL|DEADLINE|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|REVOKED|SPECIAL_EXIT) echo FAILED ;;
    *) echo UNKNOWN ;;
  esac
}

state_is_terminal() {
  [[ "$1" == COMPLETED || "$1" == FAILED || "$1" == CANCELLED || "$1" == TIMEOUT ]]
}

write_observation() {
  local observed_at="$1" source="$2" raw_state="$3" state="$4" exit_code="$5"
  local submit_time="$6" start_time="$7" end_time="$8" node_list="$9"
  local reason="${10}" run_time="${11}" raw="${12}"
  local target="${JOB_BUNDLE}/scheduler-observation.json" tmp="${JOB_BUNDLE}/scheduler-observation.json.tmp.$$"
  python3 - "${tmp}" "${JOB_REF}" "${JOB_ID}" "${observed_at}" "${source}" "${raw_state}" "${state}" "${exit_code}" "${submit_time}" "${start_time}" "${end_time}" "${node_list}" "${reason}" "${run_time}" "${raw}" <<'PY'
import json
import os
import sys

keys = ("job_ref", "slurm_job_id", "observed_at", "source", "raw_state", "state", "exit_code", "submit_time", "start_time", "end_time", "node_list", "reason", "run_time", "raw")
data = dict(zip(keys, sys.argv[2:]))
data["schema_version"] = 1
with open(sys.argv[1], "x", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
PY
  mv "${tmp}" "${target}"
}

load_terminal_observation() {
  local observation="${JOB_BUNDLE}/scheduler-observation.json"
  [[ -f "${observation}" ]] || return 1
  mapfile -t OBS_FIELDS < <(python3 - "${observation}" "${JOB_REF}" "${JOB_ID}" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    data = json.load(handle)
if data.get("job_ref") != sys.argv[2] or data.get("slurm_job_id") != sys.argv[3]:
    raise SystemExit(1)
if data.get("state") not in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT"}:
    raise SystemExit(1)
for key in ("observed_at", "source", "raw_state", "state", "exit_code", "submit_time", "start_time", "end_time", "node_list", "reason", "run_time"):
    print(data.get(key, ""))
PY
  ) || return 1
  [[ ${#OBS_FIELDS[@]} -eq 11 ]] || return 1
  OBSERVED_AT="${OBS_FIELDS[0]}"
  case "${OBS_FIELDS[1]}" in
    scontrol|sacct) OBS_SOURCE="scheduler-cache:${OBS_FIELDS[1]}" ;;
    *) OBS_SOURCE="observation-cache:${OBS_FIELDS[1]}" ;;
  esac
  OBS_RAW_STATE="${OBS_FIELDS[2]}"
  OBS_STATE="${OBS_FIELDS[3]}"
  OBS_EXIT_CODE="${OBS_FIELDS[4]}"
  OBS_SUBMIT_TIME="${OBS_FIELDS[5]}"
  OBS_START_TIME="${OBS_FIELDS[6]}"
  OBS_END_TIME="${OBS_FIELDS[7]}"
  OBS_NODE_LIST="${OBS_FIELDS[8]}"
  OBS_REASON="${OBS_FIELDS[9]}"
  OBS_RUN_TIME="${OBS_FIELDS[10]}"
}

load_terminal_runner_lifecycle() {
  local lifecycle="${JOB_BUNDLE}/lifecycle.jsonl"
  [[ -f "${lifecycle}" ]] || return 1
  mapfile -t RUNNER_EXIT_FIELDS < <(python3 - "${lifecycle}" <<'PY'
import json
import sys

try:
    with open(sys.argv[1], encoding="utf-8") as handle:
        events = [json.loads(line) for line in handle if line.strip()]
except (OSError, json.JSONDecodeError):
    raise SystemExit(1)
if not events or any(not isinstance(event, dict) for event in events):
    raise SystemExit(1)
last = events[-1]
exit_code = last.get("exit_code")
if (
    last.get("event") != "command_exited"
    or not isinstance(exit_code, int)
    or isinstance(exit_code, bool)
    or not 0 <= exit_code <= 255
    or not any(event.get("event") == "command_started" for event in events[:-1])
):
    raise SystemExit(1)
at = last.get("at")
node = last.get("node")
if not isinstance(at, str) or not at or not isinstance(node, str) or not node:
    raise SystemExit(1)
print(exit_code)
print(node)
print(json.dumps(last, sort_keys=True, separators=(",", ":")))
PY
  ) || return 1
  [[ ${#RUNNER_EXIT_FIELDS[@]} -eq 3 ]] || return 1
}

observe_job() {
  local active raw raw_state exit_code submit_time start_time end_time node_list reason run_time state sacct_raw
  node_list=""
  reason=""
  run_time=""
  OBSERVED_AT="$(utc_now)"
  active="$(squeue -h -j "${JOB_ID}" -o '%T' 2>/dev/null || true)"

  if raw="$(scontrol show job "${JOB_ID}" -o 2>/dev/null)" && [[ -n "${raw}" ]]; then
    raw_state="" exit_code="" submit_time="" start_time="" end_time="" node_list="" reason="" run_time=""
    [[ "${raw}" =~ JobState=([^[:space:]]+) ]] && raw_state="${BASH_REMATCH[1]}"
    [[ "${raw}" =~ ExitCode=([^[:space:]]+) ]] && exit_code="${BASH_REMATCH[1]}"
    [[ "${raw}" =~ SubmitTime=([^[:space:]]+) ]] && submit_time="${BASH_REMATCH[1]}"
    [[ "${raw}" =~ StartTime=([^[:space:]]+) ]] && start_time="${BASH_REMATCH[1]}"
    [[ "${raw}" =~ EndTime=([^[:space:]]+) ]] && end_time="${BASH_REMATCH[1]}"
    [[ "${raw}" =~ (^|[[:space:]])NodeList=([^[:space:]]+) ]] && node_list="${BASH_REMATCH[2]}"
    [[ "${raw}" =~ Reason=([^[:space:]]+) ]] && reason="${BASH_REMATCH[1]}"
    [[ "${raw}" =~ RunTime=([^[:space:]]+) ]] && run_time="${BASH_REMATCH[1]}"
    state="$(canonical_job_state "${raw_state}")"
    OBS_SOURCE="scontrol"
  else
    sacct_raw="$(sacct -n -X -P -j "${JOB_ID}" --format=JobIDRaw,State,ExitCode,Submit,Start,End 2>/dev/null | awk -F'|' -v id="${JOB_ID}" '$1 == id {print; exit}' || true)"
    if [[ -n "${sacct_raw}" ]]; then
      IFS='|' read -r _ raw_state exit_code submit_time start_time end_time <<< "${sacct_raw}"
      raw="${sacct_raw}"
      state="$(canonical_job_state "${raw_state}")"
      OBS_SOURCE="sacct"
    elif load_terminal_observation; then
      return 0
    elif [[ -z "${active}" ]] && load_terminal_runner_lifecycle; then
      raw_state="COMMAND_EXITED"
      exit_code="${RUNNER_EXIT_FIELDS[0]}"
      if [[ "${exit_code}" == "0" ]]; then
        state="COMPLETED"
        reason="command_exit_zero"
      else
        state="FAILED"
        reason="command_exit_nonzero"
      fi
      submit_time=""
      start_time=""
      end_time=""
      node_list="${RUNNER_EXIT_FIELDS[1]}"
      run_time=""
      raw="${RUNNER_EXIT_FIELDS[2]}"
      OBS_SOURCE="runner-lifecycle"
    else
      raw_state="${active:-NOT_FOUND}"
      state="UNKNOWN"
      exit_code=""
      submit_time=""
      start_time=""
      end_time=""
      node_list=""
      reason=""
      run_time=""
      raw="scheduler record unavailable; squeue=${active:-empty}"
      OBS_SOURCE="unavailable"
    fi
  fi

  OBS_RAW_STATE="${raw_state}"
  OBS_STATE="${state}"
  OBS_EXIT_CODE="${exit_code}"
  OBS_SUBMIT_TIME="${submit_time}"
  OBS_START_TIME="${start_time}"
  OBS_END_TIME="${end_time}"
  OBS_NODE_LIST="${node_list}"
  OBS_REASON="${reason}"
  OBS_RUN_TIME="${run_time}"
  write_observation "${OBSERVED_AT}" "${OBS_SOURCE}" "${OBS_RAW_STATE}" "${OBS_STATE}" "${OBS_EXIT_CODE}" "${OBS_SUBMIT_TIME}" "${OBS_START_TIME}" "${OBS_END_TIME}" "${OBS_NODE_LIST}" "${OBS_REASON}" "${OBS_RUN_TIME}" "${raw}"
}

print_job_result_json() {
  local action="$1"
  python3 - "${action}" "${JOB_REF}" "${JOB_ID}" "${OBS_STATE}" "${OBS_RAW_STATE}" "${OBS_EXIT_CODE}" "${OBS_SOURCE}" "${OBSERVED_AT}" "${JOB_STDOUT}" "${JOB_STDERR}" "${OBS_NODE_LIST}" "${OBS_REASON}" "${OBS_RUN_TIME}" "${JOB_BUNDLE}" <<'PY'
import datetime
import json
import os
import re
import sys

keys = ("action", "job_ref", "slurm_job_id", "state", "raw_state", "exit_code", "source", "observed_at", "stdout_path", "stderr_path", "node_list", "reason", "run_time", "bundle_path")
data = dict(zip(keys, sys.argv[1:]))

events = []
lifecycle_path = os.path.join(data["bundle_path"], "lifecycle.jsonl")
try:
    with open(lifecycle_path, encoding="utf-8") as handle:
        for line in handle:
            try:
                event = json.loads(line)
                if isinstance(event, dict) and isinstance(event.get("event"), str):
                    events.append(event)
            except json.JSONDecodeError:
                continue
except FileNotFoundError:
    pass

last = events[-1] if events else {}
event_name = last.get("event", "")
phase_map = {
    "allocated": "allocated",
    "bootstrap_started": "bootstrap",
    "bootstrap_ok": "bootstrap_ok",
    "bootstrap_failed": "bootstrap_failed",
    "command_started": "command_started",
    "command_exited": "command_exited",
}
phase = phase_map.get(event_name, "submitted" if data["state"] == "PENDING" else "unknown")
diagnostic = False

if data["state"] == "PENDING":
    health = "STARTING"
elif data["state"] == "RUNNING":
    if event_name == "bootstrap_failed":
        health = "BOOTSTRAP_FAILED"
        diagnostic = True
    elif event_name == "command_exited" and int(last.get("exit_code", 0)) != 0:
        health = "COMMAND_FAILED"
        diagnostic = True
    elif event_name in {"command_started", "command_exited"}:
        health = "OK"
    elif event_name in {"allocated", "bootstrap_started", "bootstrap_ok"}:
        health = "STARTING"
    else:
        health = "UNKNOWN"
        diagnostic = True
elif data["state"] == "COMPLETED":
    health = "OK"
    if phase == "unknown":
        phase = "completed"
elif data["state"] in {"FAILED", "CANCELLED", "TIMEOUT"}:
    health = "BOOTSTRAP_FAILED" if event_name == "bootstrap_failed" else "COMMAND_FAILED"
    diagnostic = True
elif data["source"] == "unavailable":
    health = "SCHEDULER_UNAVAILABLE"
    diagnostic = True
else:
    health = "UNKNOWN"
    diagnostic = True

def file_info(path):
    try:
        stat = os.stat(path)
        return stat.st_size, datetime.datetime.fromtimestamp(stat.st_mtime, datetime.timezone.utc).isoformat().replace("+00:00", "Z")
    except FileNotFoundError:
        return 0, None

def duration_seconds(value):
    match = re.fullmatch(r"(?:(\d+)-)?(\d+):(\d{2}):(\d{2})", value or "")
    if not match:
        return None
    days, hours, minutes, seconds = match.groups()
    return int(days or 0) * 86400 + int(hours) * 3600 + int(minutes) * 60 + int(seconds)

stdout_bytes, stdout_mtime = file_info(data["stdout_path"])
stderr_bytes, stderr_mtime = file_info(data["stderr_path"])
data.update({
    "phase": phase,
    "health": health,
    "last_event_at": last.get("at"),
    "lifecycle_path": lifecycle_path,
    "run_seconds": duration_seconds(data.pop("run_time")),
    "stdout_bytes": stdout_bytes,
    "stderr_bytes": stderr_bytes,
    "stdout_mtime": stdout_mtime,
    "stderr_mtime": stderr_mtime,
    "diagnostic_recommended": diagnostic,
    "suggested_action": f"./run_in_env.sh --job-diagnose {data['job_ref']}" if diagnostic else None,
})
if data["node_list"] in {"", "(null)", "None"}:
    data["node_list"] = None
if data["reason"] in {"", "(null)", "None"}:
    data["reason"] = None
print(json.dumps(data, sort_keys=True, separators=(",", ":")))
PY
}

job_status() {
  load_job_receipt "$1"
  require_command squeue
  require_command scontrol
  observe_job
  print_job_result_json status
}

job_collect() {
  load_job_receipt "$1"
  require_command squeue
  require_command scontrol
  observe_job
  print_job_result_json collect
  case "${OBS_STATE}" in
    COMPLETED)
      [[ "${OBS_EXIT_CODE}" == "0:0" || "${OBS_EXIT_CODE}" == "0" ]] || return 1
      return 0
      ;;
    FAILED|CANCELLED|TIMEOUT) return 1 ;;
    PENDING|RUNNING)
      echo "run_in_env.sh: job ${JOB_REF} is not terminal (${OBS_STATE})" >&2
      return 75
      ;;
    *)
      echo "run_in_env.sh: scheduler has no authoritative terminal record for ${JOB_REF}" >&2
      return 69
      ;;
  esac
}

write_cancellation_record() {
  local requested_at="$1" prior_state="$2" command_result="$3" outcome="$4"
  local stamp target tmp
  stamp="$(date -u +'%Y%m%dT%H%M%SZ')-$$"
  target="${JOB_BUNDLE}/cancellation-${stamp}.json"
  tmp="${target}.tmp"
  python3 - "${tmp}" "${WRAPPER_VERSION}" "${JOB_REF}" "${JOB_ID}" "${requested_at}" "${prior_state}" "${command_result}" "${OBSERVED_AT}" "${OBS_STATE}" "${outcome}" <<'PY'
import json
import os
import sys

keys = (
    "wrapper_version", "job_ref", "slurm_job_id", "requested_at",
    "prior_state", "command_result", "observed_at", "final_state", "outcome",
)
data = dict(zip(keys, sys.argv[2:]))
data["schema_version"] = 1
data["action"] = "cancel"
data["record_path"] = os.path.abspath(sys.argv[1][:-4])
with open(sys.argv[1], "x", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
os.replace(sys.argv[1], sys.argv[1][:-4])
print(json.dumps(data, sort_keys=True, separators=(",", ":")))
PY
}

job_cancel() {
  local requested_at prior_state command_result outcome
  local wait_seconds="${RUN_IN_ENV_CANCEL_WAIT_SECONDS:-30}" deadline

  if [[ ! "${wait_seconds}" =~ ^[0-9]+$ ]] || (( wait_seconds > 60 )); then
    echo "run_in_env.sh: RUN_IN_ENV_CANCEL_WAIT_SECONDS must be an integer from 0 to 60" >&2
    return 64
  fi

  load_job_receipt "$1"
  require_command squeue
  require_command scontrol
  require_command scancel
  observe_job
  requested_at="$(utc_now)"
  prior_state="${OBS_STATE}"

  if state_is_terminal "${OBS_STATE}"; then
    write_cancellation_record "${requested_at}" "${prior_state}" "not_requested" "already_terminal"
    return 0
  fi
  if [[ "${OBS_STATE}" == "UNKNOWN" ]]; then
    write_cancellation_record "${requested_at}" "${prior_state}" "not_requested" "scheduler_unavailable"
    echo "run_in_env.sh: scheduler has no authoritative record for ${JOB_REF}; cancellation not attempted" >&2
    return 69
  fi

  if ! scancel "${JOB_ID}"; then
    command_result="failed"
    write_cancellation_record "${requested_at}" "${prior_state}" "${command_result}" "cancel_failed"
    echo "run_in_env.sh: scheduler rejected cancellation for ${JOB_REF}" >&2
    return 1
  fi
  command_result="accepted"
  deadline=$((SECONDS + wait_seconds))
  while true; do
    observe_job
    if state_is_terminal "${OBS_STATE}"; then
      if [[ "${OBS_STATE}" == "CANCELLED" ]]; then
        outcome="cancelled"
      else
        outcome="terminal_after_request"
      fi
      write_cancellation_record "${requested_at}" "${prior_state}" "${command_result}" "${outcome}"
      return 0
    fi
    if [[ "${OBS_STATE}" == "UNKNOWN" || ${SECONDS} -ge ${deadline} ]]; then
      break
    fi
    sleep 1
  done

  if [[ "${OBS_STATE}" == "UNKNOWN" ]]; then
    outcome="confirmation_unavailable"
    write_cancellation_record "${requested_at}" "${prior_state}" "${command_result}" "${outcome}"
    echo "run_in_env.sh: cancellation was requested but scheduler confirmation is unavailable for ${JOB_REF}" >&2
    return 69
  fi
  outcome="confirmation_timeout"
  write_cancellation_record "${requested_at}" "${prior_state}" "${command_result}" "${outcome}"
  echo "run_in_env.sh: cancellation was requested but not confirmed within ${wait_seconds}s for ${JOB_REF}" >&2
  return 75
}

job_diagnose() {
  local status_json diagnostic_path stamp node_raw="" sstat_raw="" runtime_raw=""
  local stdout_tail="" stderr_tail="" runtime_rc=0 runtime_attempted=0 batch_pids=""
  local diagnostic_script

  load_job_receipt "$1"
  require_command squeue
  require_command scontrol
  observe_job
  status_json="$(print_job_result_json diagnose)"

  if [[ -n "${OBS_NODE_LIST}" && "${OBS_NODE_LIST}" != "(null)" ]]; then
    node_raw="$(timeout 3s scontrol show node "${OBS_NODE_LIST}" -o 2>&1 || true)"
  fi
  if command -v sstat >/dev/null 2>&1; then
    sstat_raw="$(timeout 3s sstat -n -P -j "${JOB_ID}.batch" --format=JobID,NTasks,Pids,AveCPU,AveDiskRead,AveDiskWrite,MaxRSS 2>&1 || true)"
    batch_pids="$(awk -F'|' 'NF >= 3 {print $3; exit}' <<< "${sstat_raw}")"
  fi

  if [[ -f "${JOB_STDOUT}" ]]; then
    stdout_tail="$(timeout 2s tail -c 8192 "${JOB_STDOUT}" 2>&1 || true)"
  fi
  if [[ -f "${JOB_STDERR}" ]]; then
    stderr_tail="$(timeout 2s tail -c 8192 "${JOB_STDERR}" 2>&1 || true)"
  fi

  if [[ "${OBS_STATE}" == "RUNNING" ]] && command -v srun >/dev/null 2>&1; then
    diagnostic_script='printf "host=%s at=%s\n" "$(hostname)" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"; printf "%s\n" "-- mounts --"; for path in /share /shareb /sharec /shared; do findmnt -rn -T "$path" -o TARGET,SOURCE,FSTYPE,OPTIONS 2>/dev/null || true; done; printf "%s\n" "-- processes --"; if [[ -n "${RUN_IN_ENV_DIAG_PIDS:-}" ]]; then echo "process_scope=batch"; ps -p "${RUN_IN_ENV_DIAG_PIDS}" -o pid,ppid,state,etimes,wchan:32,comm --sort=pid 2>/dev/null | head -n 200; else echo "process_scope=user-fallback"; ps -u "$(id -un)" -o pid,ppid,state,etimes,wchan:32,comm --sort=pid 2>/dev/null | head -n 200; fi; printf "%s\n" "-- nfs counters --"; awk '\''/ mounted on \/(share|shareb|sharec|shared) / {show=1; print; next} show && ($1=="bytes:" || $1=="xprt:") {print} show && /^$/ {show=0}'\'' /proc/self/mountstats 2>/dev/null | head -n 100'
    set +e
    runtime_attempted=1
    runtime_raw="$(RUN_IN_ENV_DIAG_PIDS="${batch_pids}" timeout 15s srun --jobid="${JOB_ID}" --overlap --exact --nodes=1 --ntasks=1 --time=00:00:10 bash -lc "${diagnostic_script}" 2>&1)"
    runtime_rc=$?
    set -e
  fi

  stamp="$(date -u +'%Y%m%dT%H%M%SZ')"
  diagnostic_path="${JOB_BUNDLE}/diagnostic-${stamp}-$$.json"
  RUN_IN_ENV_DIAG_STATUS="${status_json}" \
  RUN_IN_ENV_DIAG_NODE="${node_raw:0:32768}" \
  RUN_IN_ENV_DIAG_SSTAT="${sstat_raw:0:32768}" \
  RUN_IN_ENV_DIAG_RUNTIME="${runtime_raw:0:65536}" \
  RUN_IN_ENV_DIAG_RUNTIME_RC="${runtime_rc}" \
  RUN_IN_ENV_DIAG_RUNTIME_ATTEMPTED="${runtime_attempted}" \
  RUN_IN_ENV_DIAG_STDOUT="${stdout_tail:0:8192}" \
  RUN_IN_ENV_DIAG_STDERR="${stderr_tail:0:8192}" \
  python3 - "${diagnostic_path}" "${JOB_BUNDLE}/scheduler-observation.json" "${JOB_BUNDLE}/preflight.json" <<'PY'
import json
import os
import re
import sys
import datetime

def load_optional(path):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return None

status = json.loads(os.environ["RUN_IN_ENV_DIAG_STATUS"])
runtime_output = os.environ["RUN_IN_ENV_DIAG_RUNTIME"]
classification = status["health"]
if (
    "process_scope=batch" in runtime_output
    and re.search(r"(?m)^\s*\d+\s+\d+\s+D\s+.*(?:rpc_|nfs|wait_on_page|io_schedule)", runtime_output, re.IGNORECASE)
):
    classification = "SUSPECT_STALLED"
report = {
    "schema_version": 1,
    "action": "diagnose",
    "job_ref": status["job_ref"],
    "slurm_job_id": status["slurm_job_id"],
    "created_at": datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
    "classification": classification,
    "status": status,
    "scheduler_observation": load_optional(sys.argv[2]),
    "preflight": load_optional(sys.argv[3]),
    "node_raw": os.environ["RUN_IN_ENV_DIAG_NODE"],
    "sstat_raw": os.environ["RUN_IN_ENV_DIAG_SSTAT"],
    "runtime_probe": {
        "attempted": os.environ["RUN_IN_ENV_DIAG_RUNTIME_ATTEMPTED"] == "1",
        "exit_code": int(os.environ["RUN_IN_ENV_DIAG_RUNTIME_RC"]),
        "output": os.environ["RUN_IN_ENV_DIAG_RUNTIME"],
    },
    "stdout_tail": os.environ["RUN_IN_ENV_DIAG_STDOUT"],
    "stderr_tail": os.environ["RUN_IN_ENV_DIAG_STDERR"],
    "read_only": True,
}
with open(sys.argv[1], "x", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")
    handle.flush()
    os.fsync(handle.fileno())
print(json.dumps({
    "action": "diagnose",
    "job_ref": status["job_ref"],
    "slurm_job_id": status["slurm_job_id"],
    "state": status["state"],
    "health": classification,
    "diagnostic_path": sys.argv[1],
}, sort_keys=True, separators=(",", ":")))
PY
}

run_cmd() {
  if [[ $# -eq 0 ]]; then
    usage >&2
    return 64
  fi

  normalize_resource_requests
  select_strategy
  validate_network_capability "${SELECTED_STRATEGY}"
  print_execution_plan

  case "${SELECTED_STRATEGY}" in
    local)
      resolve_local_context
      if [[ "${LOCAL_EXEC_MODE}" == "in_container" ]]; then
        exec_in_current_container "$@"
      fi
      exec_via_target_image "$@"
      ;;
    node03)
      exec_via_node03 "$@"
      ;;
    slurm)
      slurm_run_cmd "$@"
      ;;
    *)
      echo "run_in_env.sh: unsupported strategy: ${SELECTED_STRATEGY}" >&2
      return 64
      ;;
  esac
}

self_check() {
  local -a self_check_cmd=(
    python -c "import os, sys; print('python_executable=' + sys.executable); print('cwd=' + os.getcwd())"
  )

  normalize_resource_requests
  select_strategy
  validate_network_capability "${SELECTED_STRATEGY}"

  echo "run_in_env.sh: self-check"
  echo "strategy=${SELECTED_STRATEGY}"
  echo "reason=${SELECTED_REASON}"
  echo "runtime=${RUNTIME}"
  echo "image=${IMAGE}"
  echo "site_profile=${SITE_PROFILE_PATH}"
  echo "env_root=${ENV_ROOT}"
  echo "default_env_name=${DEFAULT_ENV_NAME}"
  echo "env_prefix=${ENV_PREFIX}"
  echo "cpu=${CPU_REQUEST}"
  echo "mem=${MEM_REQUEST_RAW}"
  echo "gpu=${GPU_REQUEST}"
  echo "time=${TIME_REQUEST_RAW}"
  echo "network=${NETWORK_REQUEST}"
  run_cmd "${self_check_cmd[@]}"
}

create_env() {
  if [[ "${ENV_REQUEST_EXPLICIT}" != "1" ]]; then
    usage_error "--create-env requires --env"
  fi
  if [[ -z "${SPEC_FILE}" ]]; then
    usage_error "--create-env requires --file"
  fi
  if [[ -e "${ENV_PREFIX}" ]]; then
    echo "run_in_env.sh: target env already exists: ${ENV_PREFIX}" >&2
    return 1
  fi
  if [[ "${CPU_REQUEST}" != "4" || "${GPU_REQUEST}" != "0" || "${MEM_REQUEST_RAW}" != "16G" || "${TIME_REQUEST_RAW}" != "-1" || "${NETWORK_REQUEST}" != "none" ]]; then
    usage_error "--create-env does not accept --cpu/--mem/--gpu/--time/--network"
  fi

  if [[ -n "${SINGULARITY_NAME:-}" || -n "${APPTAINER_NAME:-}" ]]; then
    if ! EXEC_RUNNER="$(current_runner)"; then
      echo "run_in_env.sh: inside the container, but no supported runner is available (micromamba/mamba/conda)" >&2
      return 127
    fi
    create_in_current_container
  fi

  resolve_target_image_context 0
  create_via_target_image
}

canonical_dispatch() {
  local deadline="${RUN_IN_ENV_TEST_CANONICAL_DEADLINE_SECONDS:-240}"
  local poll_interval="${RUN_IN_ENV_TEST_CANONICAL_POLL_SECONDS:-2}"
  local -a driver_args wrapper_args

  if [[ ! "${deadline}" =~ ^[0-9]+$ ]] || ((deadline < 1 || deadline > 240)); then
    echo "run_in_env.sh: invalid internal canonical deadline" >&2
    return 78
  fi
  if [[ ! -f "${CANONICAL_DRIVER}" ]]; then
    echo "run_in_env.sh: canonical runtime component not found: ${CANONICAL_DRIVER}" >&2
    return 78
  fi
  require_command python3

  driver_args=(
    "${CANONICAL_DRIVER}"
    --project-root "${PROJECT_ROOT}"
    --wrapper "${PROJECT_ROOT}/run_in_env.sh"
    --action "${ACTION}"
    --deadline "${deadline}"
    --poll-interval "${poll_interval}"
  )
  if [[ -n "${CANONICAL_CONTROL}" ]]; then
    driver_args+=(--control "${CANONICAL_CONTROL}")
  fi

  case "${ACTION}" in
    run|submit)
      if [[ "${ENV_REQUEST_EXPLICIT}" == "1" ]]; then
        wrapper_args+=(--env "${ENV_REQUEST}")
      fi
      wrapper_args+=(
        --cpu "${CPU_REQUEST}"
        --mem "${MEM_REQUEST_RAW}"
        --gpu "${GPU_REQUEST}"
        --time "${TIME_REQUEST_RAW}"
        --network "${NETWORK_REQUEST}"
      )
      if [[ "${ACTION}" == "submit" ]]; then
        wrapper_args+=(--submit)
      fi
      wrapper_args+=(-- "${POSITIONAL_ARGS[@]}")
      ;;
    job_status) wrapper_args=(--job-status "${JOB_REF_REQUEST}") ;;
    job_collect) wrapper_args=(--job-collect "${JOB_REF_REQUEST}") ;;
    job_cancel) wrapper_args=(--job-cancel "${JOB_REF_REQUEST}") ;;
    *)
      usage_error "--canonical does not support --${ACTION/_/-}"
      ;;
  esac

  python3 "${driver_args[@]}" -- "${wrapper_args[@]}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --self-check)
      if [[ "${ACTION}" != "run" ]]; then
        usage_error "choose only one action"
      fi
      ACTION="self_check"
      shift
      ;;
    --create-env)
      if [[ "${ACTION}" != "run" ]]; then
        usage_error "choose only one action"
      fi
      ACTION="create_env"
      shift
      ;;
    --submit)
      if [[ "${ACTION}" != "run" ]]; then
        usage_error "choose only one action"
      fi
      ACTION="submit"
      shift
      ;;
    --canonical)
      if [[ "${CANONICAL_MODE}" == "1" ]]; then
        usage_error "--canonical may be specified only once"
      fi
      CANONICAL_MODE=1
      shift
      ;;
    --control)
      if [[ $# -lt 2 ]]; then
        usage_error "--control requires a value"
      fi
      CANONICAL_CONTROL="$2"
      shift 2
      ;;
    --job-status|--job-collect|--job-diagnose|--job-cancel)
      if [[ "${ACTION}" != "run" ]]; then
        usage_error "choose only one action"
      fi
      if [[ $# -lt 2 ]]; then
        usage_error "$1 requires a job reference"
      fi
      if [[ "$1" == "--job-status" ]]; then
        ACTION="job_status"
      elif [[ "$1" == "--job-diagnose" ]]; then
        ACTION="job_diagnose"
      elif [[ "$1" == "--job-cancel" ]]; then
        ACTION="job_cancel"
      else
        ACTION="job_collect"
      fi
      JOB_REF_REQUEST="$2"
      shift 2
      ;;
    --env)
      if [[ $# -lt 2 ]]; then
        usage_error "--env requires a value"
      fi
      ENV_REQUEST="$2"
      ENV_REQUEST_EXPLICIT=1
      shift 2
      ;;
    --file)
      if [[ $# -lt 2 ]]; then
        usage_error "--file requires a value"
      fi
      SPEC_FILE="$2"
      shift 2
      ;;
    --cpu)
      if [[ $# -lt 2 ]]; then
        usage_error "--cpu requires a value"
      fi
      CPU_REQUEST="$2"
      shift 2
      ;;
    --mem)
      if [[ $# -lt 2 ]]; then
        usage_error "--mem requires a value"
      fi
      MEM_REQUEST_RAW="$2"
      shift 2
      ;;
    --gpu)
      if [[ $# -lt 2 ]]; then
        usage_error "--gpu requires a value"
      fi
      GPU_REQUEST="$2"
      shift 2
      ;;
    --time)
      if [[ $# -lt 2 ]]; then
        usage_error "--time requires a value"
      fi
      TIME_REQUEST_RAW="$2"
      shift 2
      ;;
    --network)
      if [[ $# -lt 2 ]]; then
        usage_error "--network requires a value"
      fi
      NETWORK_REQUEST="$2"
      shift 2
      ;;
    --slurm-*)
      legacy_flag_error "$1"
      ;;
    --)
      shift
      POSITIONAL_ARGS+=("$@")
      break
      ;;
    -*)
      usage_error "unknown option: $1"
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      while [[ $# -gt 0 ]]; do
        POSITIONAL_ARGS+=("$1")
        shift
      done
      ;;
  esac
done

ENV_PREFIX="$(resolve_env_prefix "${ENV_REQUEST:-${ENV_PREFIX}}")"

if [[ -n "${SPEC_FILE}" ]]; then
  SPEC_FILE="$(absolute_existing_path "${SPEC_FILE}")"
fi

if [[ "${CANONICAL_MODE}" == "0" && -n "${CANONICAL_CONTROL}" ]]; then
  usage_error "--control requires --canonical"
fi
if [[ "${CANONICAL_MODE}" == "1" ]]; then
  if [[ -n "${SPEC_FILE}" ]]; then
    usage_error "--canonical does not accept --file"
  fi
  canonical_dispatch
  exit $?
fi

case "${ACTION}" in
  self_check)
    if [[ ${#POSITIONAL_ARGS[@]} -ne 0 ]]; then
      usage_error "--self-check does not accept a command"
    fi
    if [[ -n "${SPEC_FILE}" ]]; then
      usage_error "--self-check does not accept --file"
    fi
    self_check
    ;;
  create_env)
    if [[ ${#POSITIONAL_ARGS[@]} -ne 0 ]]; then
      usage_error "--create-env does not accept a command"
    fi
    create_env
    ;;
  submit)
    if [[ -n "${SPEC_FILE}" ]]; then
      usage_error "--file is only valid with --create-env"
    fi
    submit_durable_job "${POSITIONAL_ARGS[@]}"
    ;;
  job_status|job_collect|job_diagnose|job_cancel)
    if [[ ${#POSITIONAL_ARGS[@]} -ne 0 ]]; then
      usage_error "--${ACTION/_/-} does not accept a command"
    fi
    if [[ -n "${SPEC_FILE}" || "${ENV_REQUEST_EXPLICIT}" == "1" || "${CPU_REQUEST}" != "4" || "${GPU_REQUEST}" != "0" || "${MEM_REQUEST_RAW}" != "16G" || "${TIME_REQUEST_RAW}" != "-1" || "${NETWORK_REQUEST}" != "none" ]]; then
      usage_error "--${ACTION/_/-} accepts only a job reference"
    fi
    case "${ACTION}" in
      job_status) job_status "${JOB_REF_REQUEST}" ;;
      job_collect) job_collect "${JOB_REF_REQUEST}" ;;
      job_diagnose) job_diagnose "${JOB_REF_REQUEST}" ;;
      job_cancel) job_cancel "${JOB_REF_REQUEST}" ;;
    esac
    ;;
  run)
    if [[ -n "${SPEC_FILE}" ]]; then
      usage_error "--file is only valid with --create-env"
    fi
    run_cmd "${POSITIONAL_ARGS[@]}"
    ;;
  *)
    usage_error "unsupported action: ${ACTION}"
    ;;
esac
