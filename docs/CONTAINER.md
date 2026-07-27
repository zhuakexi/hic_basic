# Container requirements

`hic_basic` keeps Python and Conda dependencies in this repository, but it does
not own the Singularity image. `run_in_env.sh` combines a site-managed image
with a repository-managed Conda prefix.

## Minimal contract

The supported Linux container must provide:

- a glibc-based x86_64 userspace compatible with the linux-64 lock;
- `micromamba` for creating and running prefix environments;
- `libgl1` and `libgomp1`, required by the Open3D feature;
- standard POSIX shell and core utilities.

CUDA is not required by the core library. Git, Node.js, LibreOffice,
build-essential, SSH, and SLURM clients are useful site/runtime capabilities,
but they are not minimum `hic_basic` package dependencies.

## Known-working implementation

The maintained, broader implementation is:

```text
/share/home/ychi/dev/docker_files/light_base/Dockerfile
```

Its publication chain is documented in
`/share/home/ychi/dev/docker_files/light_base/README.md` and mirrored at
<https://github.com/zhuakexi/light_base>. The image intentionally contains
operational tools beyond this repository's minimum contract.

The current mgt deployment selects its SIF in `run_in_env.local.sh`. Changing
the image and changing the Conda lock are separate changes and must be tested
separately.

## Smoke checks

Run these through the repository runtime:

```bash
./run_in_env.sh --self-check
./run_in_env.sh python -c "import open3d; print(open3d.__version__)"
./run_in_env.sh python -c "from pymol import cmd; print(cmd is not None)"
```
