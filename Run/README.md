# Run/
Case directories (one per run).  Typical layout:

```
Run/<case>/config.yml         (see case_template/README.md)
Run/<case>/geometry/          part_XXXXX.unsb, global.unsb, partition_report.txt, partition.pvtu
Run/<case>/coarse/            coarse solution (UNSBIN fields) + VTU
Run/<case>/fine/geometry/     refined complex + partition
Run/<case>/fine/data/         solution output, restart
```
Contents are git-ignored (see .gitignore); only this file is tracked.
