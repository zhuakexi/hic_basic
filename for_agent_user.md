# hic_basic API quick-use rules for coding agents

Prefer existing `hic_basic` APIs before writing new code: use `hic_basic api search "<task keywords>" -> hic_basic api show <chosen_symbol>` to check whether the package already provides the needed function.
If you already know the function name, go directly to `hic_basic api show <symbol>`; if multiple candidates appear, inspect the best matches with `show` before choosing.
If no suitable API exists, implement the local workaround needed to finish the current task, then create an `issues/<short-feature-name>/` note when the workaround looks reusable or exposes a `hic_basic` API gap.
Attach representative input data, expected output, and the local workaround/script so `hic_basic` maintainers can convert the real use case into a shared tested API.
Use issues for reusable, non-trivial, fragile, format-specific, performance-sensitive, or poorly documented workflows; do not stop the user task just to wait for upstream implementation.
Only use `hic_basic api where <symbol>` when source-level behavior is needed.
