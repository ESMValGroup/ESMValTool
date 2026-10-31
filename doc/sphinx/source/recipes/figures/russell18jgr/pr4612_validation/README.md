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

## CMIP7 `im55` time-varying `volcello` check

The published ACCESS-ESM1-6 historical monthly files in `/g/data/im55/publications` include time-varying `volcello`. A focused 1986–1987 CMIP7 recipe run using the current PR code passed Figures 4, 6a, and 6b on Gadi `analysis3-26.08` (PBS `182713673.gadi-pbs`, exit 0). The unmodified diagnostic plots are [Figure 4](Fig4_Gadi_CMIP7_ACCESS-ESM1-6_1986-1987.png), [Figure 6a](Fig6a_Gadi_CMIP7_ACCESS-ESM1-6_1986-1987.png), and [Figure 6b](Fig6b_Gadi_CMIP7_ACCESS-ESM1-6_1986-1987.png). These show a different model and period from the documentation comparisons above, so their numerical values are not directly comparable. The long model name also overlaps the top annotations in these raw plots.

For all 240 months of 1986–2005, a section-only calculation compared the PR's product of climatological velocity and volume with the mean of monthly velocity-times-volume, holding its current grid/indexing and area formulas fixed. The Figure 4 difference was `-0.000033 Sv` (155.149515 versus 155.149482 Sv; PBS `182713624.gadi-pbs`, exit 0). At 30°S the Figure 6 **net volume transport** difference was `+0.000106 Sv` (1.050500 versus 1.050605 Sv; PBS `182713719.gadi-pbs`, exit 0). Thus monthly volume variability has a negligible effect on these two aggregate estimates for this model and period. Density-layer membership and heat transport are checked separately below.

The velocity and tracer-grid volume coordinates are staggered (about 37 km at Drake Passage and 74 km at 30°S at the inspected array indices), yet the diagnostic combines them by index. The model also supplies native `umo` and `vmo` mass transports, but these are on different face coordinates from `uo` and `vo`: along the Drake section the `umo` latitude differs from `uo` by 0.13–0.45°, and at 30°S the `vmo` longitude differs from `vo` by 0.5° (PBS `182715369.gadi-pbs`, exit 0). **A previously reported same-index cellwise comparison with `umo`/`vmo` is therefore not a valid estimate of the PR's physical transport error.** Neighboring-cell volume choices and simple velocity interpolation did not establish a conservative face-aligned comparison (PBS `182714724` and `182715607`, both exit 0). The section-integrated native-flux benchmark below avoids the cellwise alignment problem. The Gadi scripts and PBS logs are under `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_1600deb/jobs/`; the CMIP7 recipe output is under `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_cmip7_4920f472/output_smoke/`.

### Model-native absolute transport benchmark

For the same 1986–2005 period, I integrated the published monthly `umo`/`vmo` mass fluxes and the corresponding raw MOM5 monthly `tx_trans`/`ty_trans` at their **own native faces**, rather than pairing their cells with the offset velocity grid. The Drake section is 291°E, 76°S–48°S; the 30°S section spans every longitude and depth. The adjacent Drake end rows are land, so shifting either end by one row leaves the native total unchanged. Mass flux (kg/s) is expressed as approximate volume flux using a fixed reference density of 1035 kg/m³; the raw and published fields agree closely in their section totals. These were checked in PBS `182719669`, `182720421`, and `182720616` (all exit 0).

| 1986–2005 net section transport | PR climatological diagnostic | Raw MOM5 resolved face flux | Published face flux |
| --- | ---: | ---: | ---: |
| Figure 4, Drake volume | +155.1495 Sv | `tx_trans`: +157.9770 Sv | `umo`: +157.9860 Sv |
| Figure 6a, 30°S volume | +1.0505 Sv | `ty_trans`: +0.8618 Sv | `vmo`: +0.8762 Sv |

The Drake estimate is 2.83 Sv (1.8%) below the model-native resolved face flux. At 30°S, the 0.17–0.19 Sv difference is in a small residual of much larger opposing northward and southward flows, so a percentage of the net is misleading. As a source-run check, the raw and published monthly 1986 face fields correlate at 0.99996 (`tx_trans`/`umo`) and 0.99447 (`ty_trans`/`vmo`). This provides an **independent absolute net-volume benchmark** for this model, but does not establish cellwise equivalence or validate each density layer.

The raw MOM5 annual `temp_yflux_adv` diagnostic, described in the file as `cp*rho*dzt*dxt*v*temp` and stored in watts, gives **−0.2927 PW** of resolved advective heat transport at 30°S over 1986–2005 (PBS `182720421`, exit 0). Figure 6b gives **−0.1853 PW** from climatological fields, 0.1074 PW less southward in magnitude; the mean of monthly PR-style classified transports is **−0.2457 PW**. The native flux includes velocity–temperature covariation and uses the model's own face geometry and heat-capacity convention, so these differences cannot be assigned solely to `volcello` averaging. The PR result reproduces the NCL-style diagnostic, but its absolute heat and density-layer values should not be presented as validated native physical transports.

### Figure 6 density-layer and heat check

A separate calculation classified each of the 240 monthly `thetao`/`so` sections during 1986–2005 and averaged the resulting volume and heat transports, using the PR's density equation, layer boundaries, section, and **same-index** velocity/volume geometry (PBS `182714649` and `182714819`, both exit 0). As a cross-check, its 1986–1987 climatological layer values reproduce the recipe's saved diagnostic NetCDF values to within `6e-7 Sv` and `4e-8 PW`.

| 1986–2005 quantity | PR: transport from climatological fields | Mean of monthly classified transport |
| --- | ---: | ---: |
| Figure 6a net volume | +1.050500 Sv | +1.050605 Sv |
| Figure 6a surface layer | −0.257610 Sv | −1.965939 Sv |
| Figure 6a 26.10 σ₀ layer | +7.470736 Sv | +8.913927 Sv |
| Figure 6b net heat | −0.185335 PW | −0.245695 PW |
| Figure 6b surface layer | −0.421713 PW | −0.595747 PW |

The largest difference among the 11 non-net layers is `1.708 Sv` for volume and `0.174 PW` for heat. Reclassifying monthly density while holding climatological transport fixed accounts for differences up to `1.417 Sv` and `0.094 PW`; monthly velocity, volume, and temperature products add further differences. The [comparison plot](Fig6_CMIP7_monthly_vs_climatology_1986-2005.png) and [complete layer values](Fig6_CMIP7_monthly_vs_climatology_1986-2005.csv) make the result reviewable.

The PR faithfully reproduces the NCL-style diagnostic based on climatological fields. It should be described as such: **its Figure 6 layer values are not the 20-year mean of monthly density-layer transports**, and its net heat transport is also sensitive to monthly covariation. The negligible net-volume effect of time-averaging `volcello` does not establish equivalence for the layers or heat. The model-native benchmarks above quantify the absolute net differences for this model and period, while the density-layer transports remain without a face-aligned physical benchmark.
