# Quasinormal-Mode Collider Spectroscopy

Open **QCS.nb** in Mathematica and evaluate the setup cell. The notebook provides the calculations, plots and tests. We also include the check of our results, which is summarized in [NUMERICAL_CHECKS.md](NUMERICAL_CHECKS.md).

## Quick start

```wolfram
Get[FileNameJoin[{NotebookDirectory[], "Kernel", "QCS.wl"}]]
QCSCheckData[]
pair = QCSPreparePair[];
recovery = QCSReconstructPair[pair];
QCSFigures[True]
```
