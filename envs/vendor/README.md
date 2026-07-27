# Vendored environment sources

`hires_utils-0.1+hicbasic.45eab4f.tar.gz` is a source-distribution backup of
`zhuakexi/hires_utils` commit
`45eab4f32f61ea9d895d5196744917355bafec4d`. That commit was the clean base of
the previously deployed editable checkout but was not available from its
private upstream when this environment was locked.

The upstream checkout contains no license file. This artifact is retained for
the repository owner's internal reproducibility; resolve its licensing before
public redistribution.

The source is unchanged except for `hires-utils-package-data.patch`. The patch
gives the artifact a provenance-bearing local version and includes
`reference/chrom_alias.csv`, which `hires_utils.hires_io` reads during import
but upstream `setup.py` omitted from non-editable installations.

The artifact is installed with `--no-deps` and a SHA-256 requirement hash.
Never replace it in place. Export a clean source revision, reapply and review
the packaging patch, use a new filename/version, update `dependencies.yml`,
and regenerate the pip and Conda locks.

`python_tsp-0.4.2+hicbasic.1-py3-none-any.whl` is the upstream MIT-licensed
0.4.2 wheel with the small change recorded in `python-tsp-no-tsplib.patch`.
`hic_basic` uses only Euclidean distances and simulated annealing, but upstream
eagerly imported the optional TSPLIB reader and required `tsplib95`. That
reader pins obsolete NetworkX and tabulate versions. The patched wheel loads
TSPLIB support only if called and removes its unconditional metadata
requirement; NumPy and requests requirements are unchanged.
