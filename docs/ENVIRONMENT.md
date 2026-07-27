# Environment maintenance

The official development environment belongs to this repository. Historical
`hic_basic_v09x` prefixes mixed dependencies from several analysis projects;
they are compatibility references, not environment sources or rollback bases.

## Quick start

For a first checkout, create the official prefix directly from the committed
locks, install the pip exceptions without dependency resolution, install this
checkout editable, and verify it:

```bash
make env-bootstrap
./run_in_env.sh --self-check
```

All later operations also go through `run_in_env.sh`:

```bash
make env-lock
make env-create
make env-check
make env-test
```

`make env-create` uses a content-addressed candidate name. Override it with
`make env-create ENV=<name>` when needed.

## Owned artifacts

- `envs/dependencies.yml`: direct dependency registry and rationale.
- `environment.yml`: generated, human-readable Conda intent.
- `envs/locks/hic-basic-linux-64.lock`: exact Conda URLs and hashes.
- `envs/locks/pip-linux-64.txt`: hashed exceptions unavailable from Conda.
- `envs/official-name.txt`: promoted content-addressed prefix used for bootstrap.
- `scripts/envctl.py`: render, lock, create, install, snapshot, and verify logic.

Do not add a transitive dependency to the registry. Do not hand-edit generated
locks or `environment.yml`.

## Adding a dependency

1. Register it in `envs/dependencies.yml` with manager, scope, feature, reason,
   imports, and an explicit constraint.
2. Prefer Conda. A pip entry requires evidence that the configured Conda
   channels provide no usable artifact; document missing packages or broken
   dependency metadata in the registry and vendor notes. Remote pip entries
   must be exact versions with wheels; source exceptions must be vendored.
3. Run `make env-render` and review the intent diff.
4. Run `make env-lock` and review every changed Conda URL and pip artifact
   hash. Resolution downloads no pip dependencies.
5. Create a new candidate; never install into an official prefix.
6. Run `make env-check` and `make env-test` against the candidate.
7. Run feature-specific tests and compare the full suite with the prior
   official environment. No new failure is acceptable.

Pip exceptions install with `--require-hashes --no-deps`. The repository itself
installs last with `python -m pip install --no-deps -e .`. These rules prevent
pip from replacing Conda-owned NumPy, pandas, xarray, or other compiled stacks.

## Verification

`envctl verify` checks:

- every installed Conda record against the explicit lock;
- unexpected Conda packages;
- unexpected pip-owned distributions;
- missing or wrong-version pip exceptions;
- pip shadowing of Conda-owned packages;
- all registered imports;
- `pip check`;
- editable `hic_basic` provenance.

Snapshots under `.env-backups/` record actual Conda and Python distribution
metadata. They are local audit/rollback evidence and are not committed.

## Browser renderer asset

The `structure-render` feature registers and locks Playwright plus Chromium's
runtime libraries as normal Conda dependencies. The matching headless Chromium
binary is maintained by Playwright in its user cache rather than inside a Conda
prefix. Install that binary explicitly after environment creation:

```bash
make renderer-bootstrap
```

Direct rendering never performs this download implicitly. A missing or
incompatible browser therefore fails with the bootstrap command instead of
changing the environment during a scientific run. The renderer reports the
Playwright, browser, protocol, Three.js, and WebGL versions with every result.

## Promotion and rollback

Candidates use `hic_basic-candidate-<lock-sha>`. A promoted environment is
rebuilt from the same locks as `hic_basic-env-<lock-sha>`, verified again, then
recorded in `envs/official-name.txt` and selected in `run_in_env.local.sh`.
Refresh `run_in_env.runtime.json` whenever that project-local runtime
configuration changes.

The name digest covers the exact Conda and pip locks, not registry prose or the
editable source tree. The same install artifacts therefore retain the same
name while source edits remain dynamically testable.

Keep at least the previous repository-owned official prefix. The mixed
`hic_basic_v096` prefix may remain as an emergency compatibility reference, but
it is not part of the official rollback chain.

Container requirements and the maintained broader image are documented in
`docs/CONTAINER.md`.
