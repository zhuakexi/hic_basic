# tests/data

This directory holds **minimal, reproducible** test fixtures for hic_basic.

Guidelines
- Keep files as small as possible while still reproducing the behavior.
- Prefer truncated or synthetic data; avoid full datasets.
- Do not add project-specific names, sample IDs, or sensitive labels.
- Use formats already supported by the library (pairs, cool/scool, 3dg, tsv, csv, fastq).

When adding new data
- Add a short note in the test describing **how the file was created** (command or script).
- Use stable filenames that describe content, not provenance.
- Update or add tests that demonstrate why the file exists.

Tip: if you need to shrink large sequencing files, use domain-aware tools
(`samtools`, `seqtk`, `cooler`, etc.) and record the exact command used.
