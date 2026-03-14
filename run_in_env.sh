#!/usr/bin/env bash
set -euo pipefail

IMAGE="/shareb/ychi/ana/envs/light_base_1_1.sif"
ENV_PREFIX="/share/home/ychi/mambaforge/envs/hic_basic_v096"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINDS=(
  "/shared/mwang:/shared/mwang"
  "/shared/hxie:/shared/hxie"
  "/shareb/mliu:/shareb/mliu"
  "/sharec/ychi:/sharec/ychi"
  "/shareb/ychi:/shareb/ychi"
  "/shared/ychi:/shared/ychi"
  "/share/home/ychi:/share/home/ychi"
)

EXEC_MODE=""
EXEC_RUNNER=""

usage() {
  cat <<'EOF'
Usage:
  ./run_in_env.sh --self-check
  ./run_in_env.sh <command> [args...]

This is the only project-supported compute entrypoint.
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "run_in_env.sh: required command not found: $cmd" >&2
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

build_bind_args() {
  local bind
  local -a args=()
  for bind in "${BINDS[@]}"; do
    args+=(-B "$bind")
  done
  printf '%s\0' "${args[@]}"
}

ensure_env_prefix_visible() {
  if [[ ! -d "${ENV_PREFIX}" ]]; then
    echo "run_in_env.sh: target env prefix is not visible in the current context: ${ENV_PREFIX}" >&2
    return 127
  fi
  if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
    echo "run_in_env.sh: target env python is not executable in the current context: ${ENV_PREFIX}/bin/python" >&2
    return 127
  fi
}

resolve_execution_context() {
  if [[ -n "${SINGULARITY_NAME:-}" ]]; then
    ensure_env_prefix_visible
    if ! EXEC_RUNNER="$(current_runner)"; then
      echo "run_in_env.sh: inside Singularity, but no supported runner is available (micromamba/mamba/conda)" >&2
      return 127
    fi
    EXEC_MODE="in_singularity_current_container"
    return 0
  fi

  require_command singularity
  if [[ ! -f "${IMAGE}" ]]; then
    echo "run_in_env.sh: container image not found: ${IMAGE}" >&2
    return 127
  fi
  if [[ ! -d "${ENV_PREFIX}" ]]; then
    echo "run_in_env.sh: target env prefix is not visible on the host: ${ENV_PREFIX}" >&2
    return 127
  fi

  EXEC_MODE="host_via_target_image"
  EXEC_RUNNER="micromamba@target_image"
}

exec_in_current_container() {
  exec "${EXEC_RUNNER}" run -p "${ENV_PREFIX}" "$@"
}

exec_via_target_image() {
  local -a bind_args=()
  local arg
  while IFS= read -r -d '' arg; do
    bind_args+=("$arg")
  done < <(build_bind_args)

  exec singularity exec --nv "${bind_args[@]}" "${IMAGE}" \
    micromamba run -p "${ENV_PREFIX}" "$@"
}

run_cmd() {
  if [[ $# -eq 0 ]]; then
    usage >&2
    return 64
  fi

  resolve_execution_context

  if [[ "${EXEC_MODE}" == "in_singularity_current_container" ]]; then
    exec_in_current_container "$@"
  fi

  exec_via_target_image "$@"
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "--self-check" ]]; then
  shift
  if [[ $# -ne 0 ]]; then
    usage >&2
    exit 64
  fi

  resolve_execution_context

  echo "run_in_env.sh: self-check"
  echo "mode=${EXEC_MODE}"
  echo "runner=${EXEC_RUNNER}"
  echo "env_prefix=${ENV_PREFIX}"

  if [[ "${EXEC_MODE}" == "in_singularity_current_container" ]]; then
    exec_in_current_container python -c \
      "import os, sys; print('python_executable=' + sys.executable); print('cwd=' + os.getcwd())"
  fi

  exec_via_target_image python -c \
    "import os, sys; print('python_executable=' + sys.executable); print('cwd=' + os.getcwd())"
fi

run_cmd "$@"
