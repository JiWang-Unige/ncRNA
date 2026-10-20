# ncrna-annot

Ab initio annotation of lncRNA (+ small ncRNA) gene structures jointly with protein-coding genes.

- `docs/PHASE0_PROTOCOL.md`, `docs/PHASE0_STATUS.md`, `docs/phase0_results/`: Phase 0 (benchmarks, ruler, probes, unknown-handling test).
- `docs/PHASE1_PROTOCOL.md`: Phase 1 protocol v1 (frozen 2026-10-20): coding ⊗ lncRNA HMM-CRF decoder.
- `scripts/phase0/`: Phase 0 code (ruler.py = evaluation ruler v1.3).
- `env/`: environment installers (GENATATOR inference env, small-RNA tools).
- `kaggle/`: Kaggle kernel builders.

Runs: Baobab (this directory, source of truth) -> mirror on lab server `ont:~/ncrna-annot` via `~/bin/ont_sync.sh`
(`PROJECT=<this dir> REMOTE_ROOT=ncrna-annot`). data/, outputs/, logs/ are not tracked.
