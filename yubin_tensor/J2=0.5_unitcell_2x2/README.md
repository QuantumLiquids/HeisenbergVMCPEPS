# Yubin Tensor (J2=0.5, Unit Cell 2x2)

This directory stores the converted SplitIndexTPS for the older J2=0.5 dataset.

At conversion time, the source file path was:
- `/Users/wanghaoxin/GitHub/HeisenbergVMCPEPS/out_tensor.txt`

## Conversion command

```bash
./build/txt_to_tps_sitps \
  --input-txt out_tensor.txt \
  --output-sitps-dir yubin_tensor/J2=0.5_unitcell_2x2/tpsfinal \
  --ly 2 --lx 2 \
  --phy-dim 2 \
  --dim-r 5 --dim-l 5 --dim-u 5 --dim-d 5 \
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

- Canonical converted SITPS directory: `yubin_tensor/J2=0.5_unitcell_2x2/tpsfinal`
- Metadata file: `yubin_tensor/J2=0.5_unitcell_2x2/tpsfinal/tps_meta.txt`
