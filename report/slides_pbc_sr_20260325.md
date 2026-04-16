---
marp: true
theme: default
paginate: true
math: mathjax
size: 16:9
style: |
  :root {
    --bg1: #fbf8f1;
    --bg2: #f1ebdf;
    --ink: #172433;
    --muted: #5c6773;
    --line: #d8cfbf;
    --blue: #0f5fa8;
    --green: #20744f;
    --gold: #c88400;
    --red: #b54b5a;
  }

  section {
    font-family: "Avenir Next", "Helvetica Neue", "Segoe UI", sans-serif;
    color: var(--ink);
    background: linear-gradient(180deg, var(--bg1) 0%, var(--bg2) 100%);
    padding: 46px 58px;
    font-size: 21px;
    line-height: 1.28;
  }

  section.title {
    justify-content: center;
  }

  section.dark {
    background: linear-gradient(180deg, #203242 0%, #13202d 100%);
    color: #f8f5ef;
  }

  section.dark h1,
  section.dark h2,
  section.dark h3,
  section.dark strong {
    color: #f8f5ef;
  }

  section.figure-slide {
    padding: 34px 44px 28px;
  }

  section.figure-slide h2 {
    margin-bottom: 4px;
  }

  section.figure-slide h3 {
    font-size: 22px;
    margin-bottom: 8px;
    color: var(--muted);
  }

  h1 {
    font-size: 42px;
    margin: 0 0 10px 0;
    color: var(--ink);
  }

  h2 {
    font-size: 30px;
    margin: 0 0 14px 0;
    color: var(--ink);
    border-bottom: 3px solid #cdbb97;
    padding-bottom: 8px;
  }

  h3 {
    font-size: 22px;
    margin: 0 0 10px 0;
    color: var(--blue);
  }

  p, li {
    color: inherit;
  }

  strong {
    color: var(--blue);
  }

  img {
    background: transparent;
  }

  .cards {
    display: grid;
    grid-template-columns: 1fr 1fr 1fr;
    gap: 15px;
  }

  .cards.two {
    grid-template-columns: 1fr 1fr;
  }

  .card {
    background: rgba(255,255,255,0.82);
    border: 1px solid var(--line);
    border-left: 8px solid var(--blue);
    border-radius: 18px;
    padding: 15px 18px;
    box-shadow: 0 10px 26px rgba(23, 36, 51, 0.06);
  }

  .dark .card,
  .dark .band {
    background: rgba(255,255,255,0.09);
    border-color: rgba(255,255,255,0.16);
    box-shadow: none;
  }

  .card.green { border-left-color: var(--green); }
  .card.gold  { border-left-color: var(--gold); }
  .card.red   { border-left-color: var(--red); }

  .kpi {
    font-size: 31px;
    font-weight: 800;
    line-height: 1.02;
    margin: 6px 0 10px 0;
    color: inherit;
  }

  .muted {
    color: var(--muted);
    font-size: 15px;
  }

  .dark .muted {
    color: rgba(248,245,239,0.74);
  }

  .small {
    font-size: 16px;
  }

  .band {
    margin-top: 12px;
    padding: 13px 17px;
    border-radius: 16px;
    border: 1px solid var(--line);
    background: rgba(255,255,255,0.78);
  }

  .takeaway {
    margin-top: 10px;
    padding: 12px 16px;
    border-radius: 16px;
    border: 1px solid var(--line);
    background: rgba(255,255,255,0.84);
    font-size: 17px;
    line-height: 1.3;
    box-shadow: 0 10px 26px rgba(23, 36, 51, 0.05);
  }

  .takeaway strong {
    color: var(--ink);
  }

  .foot {
    position: absolute;
    left: 58px;
    right: 58px;
    bottom: 18px;
    font-size: 13px;
    color: var(--muted);
  }

  .dark .foot {
    color: rgba(248,245,239,0.68);
  }
---

<!-- _class: title -->

# PBC PEPS+VMC
## Unified `0325` Version

Square Heisenberg, PBC, TRG contraction  
2026-03-25

<div class="cards">
  <div class="card">
    <div class="muted">Absorbed from 2026-03-21</div>
    <div class="kpi">8×8 physics picture</div>
    <div class="small">`J2=0` and `J2=0.5` conclusions remain the stable core of the story.</div>
  </div>
  <div class="card gold">
    <div class="muted">Added on 2026-03-25</div>
    <div class="kpi">16×16 benchmark</div>
    <div class="small">SR and MinSR now both finish `60/60` and tie at <strong>-0.66967</strong>.</div>
  </div>
  <div class="card green">
    <div class="muted">Canonical deck rule</div>
    <div class="kpi">Only maintain `0325`.</div>
    <div class="small">This is the single version to keep updating.</div>
  </div>
</div>

<div class="band small">
  <strong>Figure cleanup in this merge:</strong> the canonical `J2=0.5` StepSelector comparison now uses job `560212` as the continuation of the green trace; the old dark-blue continuation is removed.
</div>

---

<!-- _class: figure-slide -->

## 8×8, J2 = 0
### One page is enough: baseline SR, adaptive-lr `StepSelector`, `MinSR`, Yu-Bin GO, and QMC already tell the whole story

![width:1450px](figures/slide_fig_j2_0_overview.png)

<div class="takeaway"><strong>Read this page in one sentence:</strong> the `CU` baselines at `D=5` and `D=8` are almost the same, while VMC closes the gap to QMC; `StepSelector` is the adaptive-lr version of the same SR story, and `MinSR` lands in the same low-energy window.</div>

---

<!-- _class: figure-slide -->

## 8×8, J2 = 0.5
### Page 1 asks only about bond dimension and starting tensor: `D=5` vs `D=8`, `CU` vs `SU`, with Yu-Bin GO as reference

![width:1450px](figures/slide_fig_j2_05_initial_state.png)

<div class="takeaway"><strong>Answer:</strong> `CU` starts lower than `SU`, `D=8` finishes lower than `D=5`, and the frustrated case gains more from VMC than `J2=0`.</div>

---

<!-- _class: figure-slide -->

## 8×8, J2 = 0.5
### Canonical `D = 8` benchmark page: baseline SR, higher-sample `StepSelector`, `MinSR`, and Yu-Bin GO in one plot

![width:1450px](figures/slide_fig_j2_05_optimizer_merged.png)

<div class="takeaway"><strong>Read this page in one sentence:</strong> fixed-lr SR already beats Yu-Bin `GO D=8`, `MinSR` stays slightly lower still, and the canonical StepSelector trace is now the green line extended by job `560212`, not the old dark-blue low-sample continuation.</div>

---

<!-- _class: figure-slide -->

## 16×16, J2 = 0
### Benchmark page: completed SR and MinSR runs against Yu-Bin GO `D=5` and the QMC reference

![width:1450px](figures/slide_fig_16x16_j2_0_benchmark.png)

<div class="takeaway"><strong>Read this page in one sentence:</strong> SR and MinSR finish the full `60/60` iterations at the same best energy `-0.66967`, beat the Yu-Bin benchmark `-0.6696`, and sit only `4.6e-4` above QMC.</div>

---

<!-- _class: dark -->

## Bottom Line

<div class="cards two">
  <div class="card">
    <div class="muted">J2 = 0</div>
    <div class="kpi">The 8×8 conclusion did not change.</div>
    <div class="small">VMC matters much more than increasing `D` alone, and `D=8` is already close to QMC.</div>
  </div>
  <div class="card gold">
    <div class="muted">J2 = 0.5</div>
    <div class="kpi">Use the merged benchmark figure.</div>
    <div class="small">Job `560212` now carries the canonical StepSelector continuation in the `D=8` comparison.</div>
  </div>
  <div class="card green">
    <div class="muted">16×16 benchmark</div>
    <div class="kpi">Now complete, not provisional.</div>
    <div class="small">SR and MinSR both finish `60/60` and tie at <strong>-0.66967</strong>.</div>
  </div>
  <div class="card red">
    <div class="muted">Caveat</div>
    <div class="kpi">These are still trace-level best values.</div>
    <div class="small">The optimizer benchmark story is clear, but standalone remeasurement of the best states is still useful.</div>
  </div>
</div>

<div class="foot">This `0325` deck is the merged canonical version.</div>
