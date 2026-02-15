## Appendix: `src_dmrg` Benchmark Workflow (Triangle)

This appendix documents DMRG utilities as a benchmark/reference path.

It is not the primary PEPS workflow.

### Scope

Built benchmark binaries (from default CMake in this repo):

- `triangle_vmps`
- `triangle_dmrg`
- `triangle_dmrg_measure`
- `triangle_dmrg_measure_sf`

Main input file:

- `src_dmrg/dmrg_params.json`

### Minimal Run Commands

Run from build directory:

```bash
./triangle_vmps ../src_dmrg/dmrg_params.json
mpirun -n 1 ./triangle_dmrg ../src_dmrg/dmrg_params.json
./triangle_dmrg_measure ../src_dmrg/dmrg_params.json
./triangle_dmrg_measure_sf ../src_dmrg/dmrg_params.json
```

Optional explicit bond-dimension sweep list for DMRG executables:

```bash
mpirun -n 1 ./triangle_dmrg ../src_dmrg/dmrg_params.json --D=50,100,200
```

### What this is used for

- Provide independent baseline energies/observables for narrow-system comparisons.
- Cross-check PEPS/VMC trends on selected geometries and couplings.

### What this is not

- Not a replacement for PEPS workflow in this repository.
- Not integrated into `simple_update -> vmc_optimize -> mc_measure` state flow.

### Comparison Guidance

When comparing DMRG and PEPS results:

- match lattice geometry and couplings exactly
- verify boundary-condition conventions
- compare converged trends, not one-step snapshots
- document bond dimensions and truncation controls on both sides

For the main PEPS workflow, go back to:

- `tutorials/01-quick-start.md`
- `tutorials/03-recipes.md`
