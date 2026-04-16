# Square-Heisenberg PBC PEPS+VMC Status Note

**Date:** 2026-03-25

This is the merged `0325` version. It absorbs the stable `8x8` conclusions from the `2026-03-21` reorganization and the new `16x16` benchmark results confirmed on `2026-03-25`. From this point on, `0325` is the only maintained version.

## Conventions

- **CU**: Yu-Bin's `2x2` cluster-update tensor used as the initial PEPS state.
- **SU**: our simple-update initial tensor.
- **Measured**: energy from a separate measurement stage. These are the safest numbers to compare across methods.
- **Trace**: energy taken from the optimization log or `run_history.json`. These are useful for optimizer comparison, but they are not yet standalone measurement results.

## 1. Main Messages

- For **8x8**, the core physics conclusion is unchanged: **VMC optimization matters much more than increasing `D` alone**. At `J2=0`, `D=8` fixed-lr SR already reaches `-0.673316`, very close to QMC `-0.673487`.
- For **8x8, J2=0.5**, the canonical `D=8` benchmark story is now:  
  `Yu-Bin GO < our fixed-lr SR < our StepSelector / MinSR`,  
  with one cleanup to the presentation: the canonical StepSelector comparison now uses the higher-sample continuation job `560212`, which extends the green trace and replaces the old dark-blue low-sample continuation in `fig8_yubin_comparison.png`.
- For **16x16, J2=0**, the story has materially changed. SR and MinSR both complete `60/60` iterations with zero crashes and tie at `-0.66967`. This is no longer a production-status teaser; it is a completed optimizer benchmark against Yu-Bin and QMC.

## 2. Stable 8x8 Results

### 2.1 Measured endpoints

These are the cleanest numbers because they come from explicit measurement stages.

| Case | Init | Optimizer / stage | Evidence | E/site |
| --- | --- | --- | --- | --- |
| 8x8, `J2=0`, `D=5` | CU | verify | measured | `-0.671811(33)` |
| 8x8, `J2=0`, `D=5` | CU | SR, 30 iters | measured | `-0.673108(23)` |
| 8x8, `J2=0`, `D=8` | CU | verify | measured | `-0.671875(256)` |
| 8x8, `J2=0`, `D=8` | CU | SR, 30 iters | measured | `-0.673316(20)` |
| 8x8, `J2=0.5`, `D=5` | CU | verify | measured | `-0.494904(56)` |
| 8x8, `J2=0.5`, `D=5` | CU | SR, 30 iters | measured | `-0.497786(47)` |
| 8x8, `J2=0.5`, `D=8` | CU | verify | measured | `-0.496481(278)` |

Two conclusions are already stable:

- Increasing `D` from `5` to `8` changes the **CU baseline** only weakly at `J2=0`, so **larger D alone is not the main gain**.
- Once VMC is turned on, the picture changes immediately. For `8x8, J2=0`, `D=8` SR reaches `-0.673316`, only about `0.03%` above QMC.

### 2.2 Optimizer comparisons on 8x8

These remain trace-level optimizer comparisons, but the canonical presentation is now cleaner than before.

| Case | What is being compared | Value type | Best / final E/site | Takeaway |
| --- | --- | --- | --- | --- |
| 8x8, `J2=0`, `D=8` | StepSelector, 60 iters | trace | `-0.673489` | Reaches the QMC level in the optimization trace. |
| 8x8, `J2=0`, `D=8` | MinSR, lr=`0.1`, 60 iters | trace | `~ -0.67340` | Comparable to fixed-lr SR; validates MinSR on a clean benchmark. |
| 8x8, `J2=0.5`, `D=8` | CU -> SR, 30 iters | trace | `~ -0.49814` | Already below Yu-Bin's `D=8` GO value `-0.4978`. |
| 8x8, `J2=0.5`, `D=8` | SU -> SR, 30 iters | trace | `~ -0.49752` | Worse starting point than CU, but same general basin. |
| 8x8, `J2=0.5`, `D=8` | StepSelector, canonical continuation = job `560212` | trace | `-0.49850` | Higher-sample continuation still beats Yu-Bin `D=8` GO; this is now the canonical green trace in the merged figures. |
| 8x8, `J2=0.5`, `D=8` | MinSR, lr=`0.1`, 60 iters | trace | `-0.49855` | Roughly tied with the canonical StepSelector continuation and slightly lower in this trace-level comparison. |

The practical reading is:

- For `J2=0`, fixed-lr SR is already enough to establish the physics point. StepSelector and MinSR matter mainly as optimizer refinements.
- For `J2=0.5`, the canonical comparison should no longer lean on the old blue low-sample continuation. The clearer benchmark plot is the one that keeps the original green StepSelector segment and extends it with job `560212`, the `12824`-sample continuation.

## 3. 8x8, J2 = 0.5: canonical figure update

The old `fig8_yubin_comparison.png` right panel split our StepSelector story into:

- a green `lr = 0.3` segment for iterations `0-29`,
- and a dark-blue `lr = 0.1` continuation for iterations `30-59`.

That is no longer the right canonical comparison for the merged `0325` version.

For the merged deck and report, the right way to present this benchmark is:

- keep the original green StepSelector trajectory for the first `30` iterations,
- extend that same green trace with job `560212`,
- drop the old dark-blue low-sample continuation from the main comparison.

This makes the visual comparison match the actual benchmark point we want to emphasize. Job `560212` is the `12824`-sample continuation and still lands at `-0.49850`, below Yu-Bin's `D=8` GO benchmark `-0.4978`. The point of that figure is therefore not "the lower-sample continuation once went a bit deeper," but "the higher-sample continuation still confirms the gain."

## 4. 16x16, J2 = 0: SR and MinSR now tie

The local run records are:

- `run/16x16J2=0D8/20260319_icpx_sr_fresh/run_history.json`
- `run/16x16J2=0D8/20260319_icpx_minsr_fresh/run_history.json`

The run-history summaries already contain the main message:

- Job `554726` reports `60/60` iterations completed, best `E/site = -0.66967`, better than Yu-Bin's `-0.6696`, with gap to QMC `4.6e-4`.
- Job `554923` reports `60/60` iterations completed, best `E/site = -0.66967`, and explicitly notes that MinSR matches SR at this scale.

The optimization logs give the more precise trace-level best values:

| Job | Best iter | Best trace `E/site` |
| --- | --- | --- |
| `554726` SR | `45` | `-0.669671(47)` |
| `554923` MinSR | `34` | `-0.669670(47)` |

At the present noise level, the SR-MinSR difference is negligible. The useful conclusion is simply that both optimizers now reach the same `16x16` benchmark energy, and MinSR no longer needs to be presented as only an `8x8` validation.

For the merged slide, the Yu-Bin `16x16` trajectory is still taken from the digitized draft curve in `report/data/yubin_heis_grad_L16_auto.csv`, but the curve is aligned to the exact draft-table benchmark value `-0.6696` so that the plotted trajectory and the cited endpoint stay consistent.

## 5. What this merged `0325` version establishes

Three short conclusions are now safe.

1. The `8x8` physics story is already stable and remains the foundation of the report.
2. The canonical `8x8, J2=0.5` StepSelector benchmark should use job `560212` as the continuation shown in the main comparison figure.
3. `16x16, J2=0` is no longer an in-progress status item. We have two completed `60/60` optimizer runs, and SR and MinSR tie at `-0.66967`.

The remaining caveat is unchanged: these are still **optimization-trace best values**, not standalone final measurements of the best saved states.

## 6. Files used for the merged version

| Topic | Primary source |
| --- | --- |
| 8x8 measured CU / SR results | `report/data/measure_results.csv` |
| 8x8 `J2=0` StepSelector trace | `report/data/8x8_J2=0_D8_stepselector_vmc.csv` |
| 8x8 `J2=0` MinSR trace | `report/data/8x8_J2=0_D8_minsr_lr01.csv` |
| 8x8 `J2=0.5` fixed-lr SR and SU comparison | `report/data/8x8_J2=0.5_D8_yubin_vmc_partial.csv`, `report/data/8x8_J2=0.5_D8_su_vmc_partial.csv` |
| 8x8 `J2=0.5` StepSelector canonical continuation | `run/8x8J2=0.5D8/20260307_yubin_vmc_stepselector/run_history.json`, `run/8x8J2=0.5D8/20260307_yubin_vmc_stepselector/vmc/energy/optimization_log.jsonl` |
| 8x8 `J2=0.5` MinSR | `run/8x8J2=0.5D8/20260319_icpx_minsr_60iter/run_history.json`, `run/8x8J2=0.5D8/20260319_icpx_minsr_60iter/vmc/energy/optimization_log.jsonl` |
| 16x16 `J2=0` SR summary and trace | `run/16x16J2=0D8/20260319_icpx_sr_fresh/run_history.json`, `run/16x16J2=0D8/20260319_icpx_sr_fresh/vmc/energy/optimization_log.jsonl` |
| 16x16 `J2=0` MinSR summary and trace | `run/16x16J2=0D8/20260319_icpx_minsr_fresh/run_history.json`, `run/16x16J2=0D8/20260319_icpx_minsr_fresh/vmc/energy/optimization_log.jsonl` |
| Yu-Bin draft benchmarks | `report/data/yubin_draft_extracted.md`, `report/data/yubin_heis_grad_L16_auto.csv`, `report/data/yubin_j05_D8_trajectory.csv`, `report/data/yubin_j05_D9_trajectory.csv` |
