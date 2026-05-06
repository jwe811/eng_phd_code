# UofS Coordinate Format

Most generated polygon files in this repository use the UofS direction format.
It is compact, line-oriented, and stable across TransferMatrix, MonteCarlo, and
CreatorAll workflows.

## Coordinate Convention

The span direction is `x`, lattice width `L` is `y`, and lattice height `M` is
`z`.

Older files may use `z` as the span direction. The Python post-processing tools
support those files with `--span-axis z`.

## File Structure

```text
UofS
x y z
direction
direction
-111
x y z
direction
...
-111
-999
```

- `UofS` identifies the file format.
- `x y z` starts one polygon at an integer lattice coordinate.
- Direction codes describe unit steps.
- `-111` ends one polygon.
- `-999` ends the file.
- In 2SAP files, two consecutive polygons form one object.

## Direction Codes

| Code | Move |
| :--- | :--- |
| `1` | `+x` |
| `2` | `-x` |
| `3` | `+y` |
| `4` | `-y` |
| `5` | `+z` |
| `6` | `-z` |

## Useful Checks

```bash
python3 scripts/uofs_tool.py summary data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt
python3 scripts/uofs_tool.py validate data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt -L 1 -M 1 -s 2
python3 scripts/uofs_tool.py validate data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt -p 2 -L 2 -M 1 -s 2
```
