# PR 4612 Gadi figure checks

These are outputs from the Python diagnostics in PR 4612, run with ESMValTool and ESMValCore 2.15.0 in Gadi `analysis3-26.08`. Both recipe runs completed all tasks and produced plots and NetCDF outputs. The source was staged under `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_1600deb/`.

| Figure | Gadi PBS job | Input and result | Reference comparison |
| --- | --- | --- | --- |
| 1 | `180167947.gadi-pbs` | CMIP5 CanESM2 and GFDL-ESM2M, `Omon tauuo`, 1986–2005. The polar maps contain the Southern Ocean wind stress band. | The existing documentation image uses `tauu`, while the runnable recipe now uses `tauuo`; compare broad geography, not pixel values. The old `Amon tauu` selection was rejected by ESMValCore on Gadi. |
| 2 | `180167684.gadi-pbs` | CMIP5 CanESM2 and GFDL-ESM2M, `Omon tauuo`, 1986–2005. | The two wind stress curves and peak amplitudes agree visually with the documentation image. Legend placement and rendering differ. |

The `Fig*_Gadi_*.png` files are the unmodified diagnostic plots. The `Fig*_documentation_comparison.jpg` files place each new plot beside the existing documentation image. These visual comparisons do not establish numerical equivalence for every diagnostic.
