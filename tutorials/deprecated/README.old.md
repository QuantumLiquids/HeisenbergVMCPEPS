## User Guide and Project Plan

This directory is the single source of truth for how a user runs the code and what the project delivers. Keep it short, practical, and stable.

### What this project is
- PEPS + VMC workflows for 2D Heisenberg-like models
- Square lattice first-class; triangular/kagome supported by binaries where present

### Deliverables (user-facing)
- Simple Update to prepare a PEPS state
- VMC optimization (stochastic reconfiguration via CG)
- MC measurements producing basic observables

### Minimal tutorial set (max 5 files)
- [01-quick-start.md](01-quick-start.md)
- [02-parameter-system.md](02-parameter-system.md)
- [03-parameter-reference.md](03-parameter-reference.md)
- [04-algorithm-workflows.md](04-algorithm-workflows.md)

### Data artifacts (placeholders — verified by small runs)
- `peps/` and `tpsfinal/`: serialized tensors
- `energy/energy_trajectory.csv`: columns = [update_index, energy_per_site, stderr]
- `data/energy.dat`: placeholder columns = [sample_index, energy_per_site]
- `data/correlation.dat`: placeholder schema TBD (to be validated by measurement run)

If an artifact format differs on real runs, update the exact schema in 04-algorithm-workflows.md.
