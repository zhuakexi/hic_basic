# hic_basic environment files

The repository-owned environment starts at [`../environment.yml`](../environment.yml).
It is generated from `dependencies.yml`, resolved to exact locks under `locks/`,
and managed only through `../run_in_env.sh` plus `../scripts/envctl.py`.

For normal maintenance:

```bash
make env-create
make env-check
make env-test
```

See [`../docs/ENVIRONMENT.md`](../docs/ENVIRONMENT.md) before changing a
dependency. It documents candidate creation, lock review, verification, and
promotion.

`hic_basic.yaml` and `pymol.yaml` are retained as historical compatibility
evidence. They are not current install inputs and must not be used as the base
for a new official environment. In particular, `hic_basic.yaml` was derived
from the mixed `hic_basic_v096` prefix and has no exact transitive lock.
