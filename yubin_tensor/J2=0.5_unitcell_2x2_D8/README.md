# Yubin Tensor (J2=0.5, D8, Unit Cell 2x2)

This directory stores the converted SplitIndexTPS for J2=0.5, D=8.

Source file:
- `/Users/wanghaoxin/GitHub/HeisenbergVMCPEPS/out_tensor-D8-J2=0.5.txt`

Note:
- The source file contains mixed blocks:
  - blocks 1-4: D5-sized (`1250` values each)
  - blocks 5-8: D8-sized (`8192` values each)
- This conversion uses the D8 subset (blocks 5-8) only.

## Conversion command

```bash
./build/txt_to_tps_sitps \
  --input-txt /tmp/yubin_extract/out_tensor-D8-J2=0.5_D8_only.txt \
  --output-sitps-dir yubin_tensor/J2=0.5_unitcell_2x2_D8/tpsfinal \
  --ly 2 --lx 2 \
  --phy-dim 2 \
  --dim-r 8 --dim-l 8 --dim-u 8 --dim-d 8 \
  --block-site-map '0,0;0,1;1,0;1,1' \
  --src-index-order 'phy,r,l,u,d' \
  --src-fast-order 'd,u,l,r,phy' \
  --boundary Periodic
```

## Mapping details

- Block-to-site mapping:
  - block 1 -> `(0,0)`
  - block 2 -> `(0,1)`
  - block 3 -> `(1,0)`
  - block 4 -> `(1,1)`
- Source axis labels: `[phy, r, l, u, d]`
- Source flatten order: `d -> u -> l -> r -> phy` (`d` fastest)
- Target TPS axis order: `[l, d, r, u, phy]`
- Boundary condition: `Periodic`

## Output

- SITPS directory: `yubin_tensor/J2=0.5_unitcell_2x2_D8/tpsfinal`
- Metadata file: `yubin_tensor/J2=0.5_unitcell_2x2_D8/tpsfinal/tps_meta.txt`

