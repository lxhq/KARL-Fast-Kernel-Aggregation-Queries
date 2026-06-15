# CPU-KARL Branch Policy

Use `main` for clean CPU-KARL baseline timing. This branch contains the
published `Publish_codes` source extracted from `Publish_codes.zip`, without
profiling counters, per-node timers, generated result logs, or build outputs.

Use `cpu-karl-profiling` for diagnostic instrumentation only. This branch should
contain no bound redesigns and no experimental method IDs. Build with
`-DKARL_PROFILE` to enable counters and detailed timing output. Results
collected from that build should be stored as profiling or diagnostic data, not
baseline timing data.

Use `cpu-karl-bound-experiments` for CPU-KARL bound redesign attempts such as
anchor-factorized bounds, adaptive anchor variants, residual-budget traversal,
and related diagnostic scripts.
