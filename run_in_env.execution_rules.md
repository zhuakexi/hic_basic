## Execution Rules

Run all compute through `./run_in_env.sh`; do not call raw language/Conda/container runtimes, SSH, or SLURM commands. Use the default environment unless another is required, and declare realistic `--cpu`, `--mem`, `--gpu`, and `--time` for non-trivial work. Use `--self-check` only to diagnose the execution setup. If the wrapper fails, report it rather than bypassing it.

Put multi-line exploratory code in a descriptively named `.tmp/` script, run it through the wrapper, then delete it or move reusable code into the canonical analysis.

Normal execution:

```bash
./run_in_env.sh --cpu 8 --mem 32G --time 2h python analysis.py
```

Exploratory commands, coding, and diagnosis between retries use normal
execution and have no special four-minute deadline. Run an authoritative
canonical script with `--canonical --control <attempt-control.json>`. Each
canonical invocation has its own non-renewable 240-second deadline; after a
timeout, improve the script and deliberately run it again with a fresh
per-invocation deadline.

Declare `--network internet` for downloads, APIs, remote databases, or package
resolution. Network capability is site-specific. If the selected backend lacks
it, fetch with a small local wrapper command into shared storage, then run the
large or GPU computation offline; do not wait for node03 network calls to hang.

Long canonical work uses the same entrypoint with durable submission:

```bash
./run_in_env.sh --canonical --submit --control <attempt-control.json> --cpu 8 --mem 32G --time 4h -- python analysis.py
./run_in_env.sh --job-status <job-ref>
./run_in_env.sh --job-diagnose <job-ref>
./run_in_env.sh --job-collect <job-ref>
./run_in_env.sh --job-cancel <job-ref>
```

The canonical submit form performs bounded submit/start verification and emits
one execution result. Keep its `execution_ref` for handoff validation. Status
and collect never resubmit. Scheduler evidence remains authoritative for the
scheduler outcome. After scheduler expiry, validated wrapper lifecycle evidence
may establish canonical command completion or failure; its source remains
`runner-lifecycle`, and normal output validation is still required.

Treat shared-storage I/O separately from CPU: avoid recursive scans and unbounded fan-out, use node-local scratch for reusable indexes and sort temporaries, and benchmark one sample before parallel I/O-heavy work. The wrapper classifies durable observations from `state`, `phase`, `health`, and the immutable attempt control; do not recreate that state table in prompts. `PENDING` is accepted only by an `accept_pending` control. For a prior active job, obey its `inspect_only`, `adopt`, or `cancel_before_submit` action. Cancellation must be scheduler-confirmed before replacement submission.
