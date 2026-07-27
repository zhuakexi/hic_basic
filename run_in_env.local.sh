# Project-owned overrides for the maintained run_in_env runtime.
IMAGE="/shareb/ychi/ana/envs/light_base_1_1.sif"
DEFAULT_ENV_NAME="hic_basic_v096"

# The mgt site profile already binds /share/home/ychi, /shareb/ychi, and /shared.
BINDS+=(
  "/shareb/mliu:/shareb/mliu"
  "/sharec/ychi:/sharec/ychi"
)
