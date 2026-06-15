# CPU-KARL Branch Policy

Use `main` for clean CPU-KARL baseline timing. This branch contains the
published `Publish_codes` source extracted from `Publish_codes.zip`, without
profiling counters, per-node timers, generated result logs, or build outputs.

Use `cpu-karl-profiling` for diagnostic instrumentation. Results collected from
that branch should be stored as profiling or diagnostic data, not baseline
timing data.
