# Simple update MPI

The driver chooses execution from `MPI_COMM_WORLD` size. Direct execution and
`mpirun -n 1` preserve the existing serial algorithm; `mpirun -n 2` or more
ranks enables tiled NN/NNN MPI updates. No extra flag or JSON key is needed.

```sh
# From build/; inputs are read by all ranks.
./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json
mpirun -n 4 ./simple_update ../params/physics_params.json ../params/simple_update_algorithm_params.json
```

Use PEPS commit `9808ed4` or later with the NN/NNN `ExecuteMPI` APIs, and
TensorToolkit containing the split-communicator tag-bound fix `4c6f732`.
Rebuild this executable after upgrading the installed headers.

MPI partitions a square lattice into rectangular tiles. For 16×16 and four
ranks, the process grid is 2×2. Same-layer gates have disjoint sites; only
needed boundary Gamma/lambda tensors are exchanged between layers. Excessive
rank counts produce a warning from PEPS. All ranks must use the same inputs.

MPI uses a different gate order from the serial driver, so finite-tau states
can differ between single-rank serial and multirank MPI runs. Changing the
number of MPI ranks preserves the layered order, up to floating-point effects.

SquareHeisenberg and SquareXY support NN/NNN MPI execution with OBC or PBC.
TriangleHeisenberg remains serial-only; requesting multiple ranks fails with
a diagnostic. Rank zero loads or initializes `peps/` and PEPS distributes the
initial tensors to the other ranks.

Only rank zero writes final PEPS, `tpsfinal/`, stage checkpoints and schedule
JSON/CSV files. Full states are restored before stage checkpointing or starting
the next tau stage. Advanced stopping is collective. Each driver's existing
policy for a stage that fails to converge, and its nonzero exit status, are
preserved. Exceptions are reported with the failing rank and abort a multirank
job, preventing other ranks from waiting indefinitely in a collective.

`ThreadNum` remains the per-rank numerical-thread request. Start with 1 when
using several ranks. Performance depends on bond dimension, lattice size,
communication costs and the numerical backend.

## Local regression

```sh
# From build/; all simulation files are created under a temporary directory.
make simple_update -j4
python3 ../tests/test_simple_update_mpi.py ./simple_update
```

This checks serial/MPI dispatch, NN/NNN on 2 and 4 ranks from a shared saved
input, byte-identical persisted tensors across MPI rank counts, stage dumps,
advanced stopping, checkpoint restart, and nonconverged-stage exit behavior.
The spin test includes PBC and rejection of a triangle MPI request; the t-J
test includes holes, chemical potential and site-dependent pinning.
