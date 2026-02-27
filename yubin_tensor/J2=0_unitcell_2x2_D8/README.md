# Yubin Tensor (D8, Unit Cell 2x2, assumed J2=0)

This directory stores the converted SplitIndexTPS from:

- `/Users/wanghaoxin/GitHub/HeisenbergVMCPEPS/out_tensor-D8.txt`

## Conversion command

```bash
./build/txt_to_tps_sitps \
  --input-txt out_tensor-D8.txt \
  --output-sitps-dir yubin_tensor/J2=0_unitcell_2x2_D8/tpsfinal \
  --ly 2 --lx 2 \
  --phy-dim 2 \
  --dim-r 8 --dim-l 8 --dim-u 8 --dim-d 8 \
  --block-site-map '0,0;0,1;1,0;1,1' \
  --src-index-order 'phy,r,l,u,d' \
  --src-fast-order 'd,u,l,r,phy' \
  --boundary Periodic
```

## Mapping details

- Block to site mapping:
  - block 1 -> `(0,0)` (top-left)
  - block 2 -> `(0,1)` (top-right)
  - block 3 -> `(1,0)` (bottom-left)
  - block 4 -> `(1,1)` (bottom-right)
- Source axis labels: `[phy, r, l, u, d]`
- Source flatten order: `d -> u -> l -> r -> phy` (`d` fastest)
- Target TPS axis order (fixed by qlpeps): `[l, d, r, u, phy]`
- Boundary condition: `Periodic`

## Output

- Canonical converted SITPS directory: `yubin_tensor/J2=0_unitcell_2x2_D8/tpsfinal`
- Metadata file: `yubin_tensor/J2=0_unitcell_2x2_D8/tpsfinal/tps_meta.txt`

## Verification status (2026-02-27)

- Verified by small MC sampling on 8x8 PBC (`56 x 10 = 560` samples).
- Result: `E/site = -0.671874588271331` (consistent with known D5 reference trend).
- Record file:
  - `/Users/wanghaoxin/GitHub/HeisenbergVMCPEPS/run/8x8J2=0D8/20260227_yubin_verify/run_history.json`
