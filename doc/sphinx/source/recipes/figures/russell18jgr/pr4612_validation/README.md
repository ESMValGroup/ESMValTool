# PR 4612 Gadi figure checks

These are outputs from the Python diagnostics in PR 4612, run with ESMValTool and ESMValCore 2.15.0 in Gadi `analysis3-26.08`. Both recipe runs completed all tasks and produced plots and NetCDF outputs. The source was staged under `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_1600deb/`.

| Figure | Gadi PBS job | Input and result | Reference comparison |
| --- | --- | --- | --- |
| 1 | `180167947.gadi-pbs` | CMIP5 CanESM2 and GFDL-ESM2M, `Omon tauuo`, 1986–2005. The polar maps contain the Southern Ocean wind stress band. | The existing documentation image uses `tauu`, while the runnable recipe now uses `tauuo`; compare broad geography, not pixel values. The old `Amon tauu` selection was rejected by ESMValCore on Gadi. |
| 2 | `180167684.gadi-pbs` | CMIP5 CanESM2 and GFDL-ESM2M, `Omon tauuo`, 1986–2005. | The two wind stress curves and peak amplitudes agree visually with the documentation image. Legend placement and rendering differ. |

The `Fig*_Gadi_*.png` files are the unmodified diagnostic plots. The `Fig*_documentation_comparison.jpg` files place each new plot beside the existing documentation image. These visual comparisons do not establish numerical equivalence for every diagnostic.

## Representative all-diagnostic run

PBS job `180167943.gadi-pbs` ran a reduced-model copy of the recipe to exercise all 16 diagnostic groups. Figures 1–7i completed and produced 12 plots before Figure 8 failed on the NCL color name `blue4`, which Matplotlib does not recognize. The Python diagnostic now maps that color to `#00008b`. Figure 8 then completed separately in PBS job `180168287.gadi-pbs`, producing three pages of pH maps. The first broad run was cancelled after it stopped advancing following the Figure 8 failure. A clean all-diagnostic rerun completed as PBS job `180170034.gadi-pbs` with one parallel task. Numerical comparison then exposed a masked-land fill-value error in the heat integration for Figures 9a and 9c. The shared integration routine now preserves the input masks; a focused Figure 9a–c rerun completed successfully as PBS job `180171313.gadi-pbs` (85 seconds, exit status 0). The full-model recipe (PBS `180171319.gadi-pbs`) reached Figure 5, then failed when Cartopy inferred latitude edges beyond the valid range on MRI-CGCM3’s global grid. Figure 5 now plots only rows intersecting the Southern Ocean extent, while retaining the full data in its NetCDF output. The corrected diagnostic completed all 27 sea-ice models and seven plot pages in PBS `180180778.gadi-pbs` (exit status 0). [Page 6 from that run](Fig5_Gadi_full_models_page6.png) includes the formerly failing MRI-CGCM3 panel. The resumed complete recipe (PBS `180181010.gadi-pbs`) passed Figure 5 but stopped at Figure 7h: the latitude of a curvilinear model grid was not monotonic, so Iris rejected it as a dimension coordinate in the diagnostic NetCDF file. Figures 7h and 7i now save nonmonotonic row latitude as an auxiliary coordinate. Figure 7h passed with the complete model selection in PBS `180184606.gadi-pbs` (exit status 0). PBS `180184961.gadi-pbs` stopped before diagnostics: resuming from a run that had itself been resumed left metadata pointing to missing copied files. PBS `180185565.gadi-pbs` resumed from the original preprocessing output and completed the full recipe (exit status 0, 27m56s). All 16 figure groups passed, producing 39 PNG plots and 157 diagnostic NetCDF files. [Full-model Figure 7h](Fig7h_Gadi_full_models.png) and [Figure 7i](Fig7i_Gadi_full_models.png) show the corrected latitude handling.

The comparison images below pair the existing documentation figures with outputs from these Gadi runs. The reduced-model recipe uses fewer datasets than the documentation examples, so compare the common model patterns and values rather than the number of panels or lines.

| Figure | Comparison | Observation |
| --- | --- | --- |
| 3b Subantarctic Front | [side by side](Fig3b_documentation_comparison.jpg) | CanESM2 and GFDL-ESM2M front positions follow the documented patterns. |
| 3b Polar Front | [side by side](Fig3b-2_documentation_comparison.jpg) | The common model front positions have similar longitude structure. |
| 4 | [side by side](Fig4_documentation_comparison.jpg) | CanESM2 Drake Passage section and 151.9 Sv net transport match the documentation. |
| 5 | [side by side](Fig5_documentation_comparison.jpg) | CanESM2 sea-ice maximum and minimum extents have the documented geography. |
| 5g | [side by side](Fig5g_documentation_comparison.jpg) | The common-model seasonal cycles and peak magnitudes agree. |
| 6a | [side by side](Fig6a_documentation_comparison.jpg) | CanESM2 transport profile reproduces the documented layer values. |
| 6b | [side by side](Fig6b_documentation_comparison.jpg) | CanESM2 energy profile reproduces the documented layer values. |
| 7 | [side by side](Fig7_documentation_comparison.jpg) | CanESM2 air-sea CO₂ flux has the same broad positive and negative bands. |
| 7h | [side by side](Fig7h_documentation_comparison.jpg) | CanESM2 zonal flux peaks and southern minimum align with the documentation. |
| 7i | [side by side](Fig7i_documentation_comparison.jpg) | The common-model integrated flux curves show the same direction and scale. |
| 8 | [side by side](Fig8_documentation_comparison.jpg) | MRI-ESM1 pH contours reproduce the broad pattern in the documentation. |
| 9a | [side by side](Fig9a_documentation_comparison.jpg) | All four model heat-uptake points and the regression agree with the documentation after preserving masks. |
| 9b | [side by side](Fig9b_documentation_comparison.jpg) | All four model carbon-uptake points and the regression agree. |
| 9c | [side by side](Fig9c_documentation_comparison.jpg) | All three model heat-versus-carbon points and the regression agree after preserving masks. |

The Gadi validation recipe and PBS scripts are retained under the stage path above. The representative recipe differs from the PR recipe only in its smaller model selection; Figure 9a–c used the same model choices as the PR recipe. The complete model recipe passed in PBS `180185565.gadi-pbs`. The reduced-model comparisons above provide a visual check against the existing documentation; the full-model outputs also remain on Gadi in the run directory.
