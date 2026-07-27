# Versioned mgt deployment example from docker_files@5ea0fd7 (2026-07-17).
# Project-specific IMAGE and extra binds belong in run_in_env.local.sh.

ENV_ROOT="${RUN_IN_ENV_ENV_ROOT:-/share/home/ychi/mambaforge/envs}"
DEFAULT_ENV_NAME="${RUN_IN_ENV_DEFAULT_ENV_NAME:-hic_basic_v096}"
ENV_MANAGER_BIN="${RUN_IN_ENV_MANAGER_BIN:-/share/home/ychi/mambaforge/bin}"
RUNTIME="${RUN_IN_ENV_RUNTIME:-singularity}"

BINDS=(
  "/share/home/ychi:/share/home/ychi"
  "/shareb/ychi:/shareb/ychi"
  "/shared:/shared"
)
RUNTIME_ARGS=()

NODE03_HOST="${RUN_IN_ENV_NODE03_HOST:-node03}"
SSH_ARGS=()
# mgt local has direct internet; node03 does not. SLURM stays unlisted until validated.
INTERNET_CAPABLE_STRATEGIES=(local)
SLURM_PARTITION="${RUN_IN_ENV_SLURM_PARTITION:-fatcomp,comp}"
SLURM_OUTPUT_TEMPLATE="${RUN_IN_ENV_SLURM_OUTPUT:-slurm-%j.out}"

NODE03_CPU_THRESHOLD=16
NODE03_MEM_THRESHOLD_MIB=$((48 * GIB_MIB))
SLURM_CPU_THRESHOLD=96
SLURM_MEM_THRESHOLD_MIB=$((256 * GIB_MIB))
SLURM_TIME_THRESHOLD_SECONDS=3600
