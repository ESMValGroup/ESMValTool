# PR 4612 Gadi figure checks

These are outputs from the Python diagnostics in PR 4612, run with ESMValTool and ESMValCore 2.15.0 in Gadi `analysis3-26.08`. The focused and full recipe runs completed and produced plots and NetCDF outputs. The source was staged under `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_1600deb/`.

| Figure | Gadi PBS job | Input and result | Reference comparison |
| --- | --- | --- | --- |
| 1 | `182726068.gadi-pbs` | CMIP5 CanESM2 and GFDL-ESM2M, `Amon tauu`, 1986–2005. [New plot](Fig1_Gadi_tauu_1986-2005.png). | The input now matches the [documentation figure](../Fig1_polar-contour_tauu_1986-2005.png); the NCL colour bins are reproduced. The earlier [`tauuo` plot](Fig1_Gadi_tauuo_1986-2005.png) explains the previous field difference. |
| 2 | `180167684.gadi-pbs` | CMIP5 CanESM2 and GFDL-ESM2M, `Omon tauuo`, 1986–2005. | The two wind stress curves and peak amplitudes agree visually with the documentation image. Legend placement and rendering differ. |

The `Fig*_Gadi_*.png` files are the unmodified diagnostic plots. The `Fig*_documentation_comparison.jpg` files place the initial Gadi plots beside the existing documentation images; the Figure 1 comparison JPEG still shows the earlier `tauuo` run. These visual comparisons do not establish numerical equivalence for every diagnostic.

The original Figure 1 documentation image is labelled `tauu`. The first PR validation run used `tauuo`, which has a visibly different wind-stress field. A follow-up recipe using the available CMIP5 `Amon tauu` files passed for both models (PBS `182725881.gadi-pbs`, exit 0). The NCL `ExcludeOuterBoxes` colour setting also uses the six interior palette colours across eight intervals; the Python mapping now reproduces the reference colour bar to within one 8-bit RGB count. The final rerun passed in PBS `182726068.gadi-pbs` (exit 0). Coastlines, graticule labels, and contour rendering still differ between NCL and Cartopy/Matplotlib. The reference PNG has no companion NetCDF file here, so exact data-field equivalence cannot be established from the images alone.

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

## CMIP6 companion including ACCESS-ESM1-5

`recipe_russell18jgr_cmip6.yml` adds historical 1986–2005 simulations from ACCESS-ESM1-5, CanESM5, GFDL-ESM4 and IPSL-CM6A-LR, all `r1i1p1f1`. The original CMIP5 recipe is byte-identical to the previous PR version. The [companion documentation](../../../recipe_russell18jgr_cmip6.rst) describes the model-family choices, grids and field-availability limits. ACCESS-ESM1-5 is used in every group except pH, which was unavailable in the inspected archives. Figures 6a/b use its fixed `Ofx volcello`.

### Completed Gadi run

PBS `183398542.gadi-pbs` completed successfully on 7 November 2026 with `conda/analysis3-26.08` and ESMValCore 2.15.0: **16 diagnostic groups, 47 tasks, 19 PNG figures and 55 diagnostic NetCDF files; exit status 0**. The recipe source was commit `704232e`. An output audit found finite data and provenance records for all outputs. The run resumed completed preprocessing after Figure 2 was fixed to accept nonmonotonic grid-row latitudes.

The complete output is at `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_cmip6/output/recipe_russell18jgr_cmip6_20261107_001156/`. Job scripts and logs are under the same stage's `jobs/` directory. Input-compliance warnings concerned standardized names and missing `positive` attributes for some flux variables. Expected warnings also reported unavailable `volcello` for the two velocity-only Figure 4 models, the default ACCESS plot style, and legacy `history` attribute placement.

### Corrections found during output review

The full run exposed plotting issues that successful execution alone did not detect. Polar maps needed longitude unwrapping and removal of redundant masked halo columns; otherwise filled contours could show the wrong sign or a false gap. Periodic front copies needed separate line segments, and IPSL front extraction needed longitude ordering within each row. **The latter can change derived front positions on affected grids.** Map field values and transport formulas are unchanged. Titles in Figures 4/6 now avoid their annotations. Converted carbon-flux outputs and Figure 6 outputs now carry the corresponding physical units.

The subsequent local checks passed 12 diagnostic regression tests, including row-rotation invariance, periodic front separation and cyclic halo handling, plus the three recipe checks and YAML lint. Corrected carbon maps were rendered from the saved Gadi fields; their saved numerical values agree with the baseline to floating-point precision. The recipe documentation passed a focused Sphinx build with warnings treated as errors.

**Final Gadi validation passed:** PBS `183400209.gadi-pbs` completed with exit 0 using diagnostic source `bce5290`, `analysis3-26.08` and ESMValCore 2.15.0. All 47 tasks completed (16 diagnostic groups, 19 PNG figures and 55 diagnostic NetCDF files). The audit checked finite data, provenance and an output for every selected model in each group; no errors were found. The run is `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_cmip6/output/recipe_russell18jgr_cmip6_20261107_010435/`. It reused 94 verified preprocessed files through `resume_cache_bce5290`. Earlier retries encountered an environment I/O error, an allocation reservation failure and a nested-resume path issue; these are not failures of the final diagnostic run. The existing caveats about climatological transport estimates and differences from the paper still apply.

The final front plots remove the artificial periodic connecting lines and IPSL zigzags. Gaps remain where the legacy longest-contour/nearby-segment selection does not retain a continuous path. These gaps must not be interpreted as evidence that a physical front is absent.

The [numerical comparison](../cmip6/baseline_comparison.json) found all 47 non-front fields unchanged within `rtol=1e-6, atol=1e-8`, with masks preserved. Six output unit labels were corrected. The six non-IPSL front coordinate sets preserve their finite vertices; the two IPSL fronts change after fixing row-dependent longitude ordering. The [output audit](../cmip6/validation_summary.json) records values, units and coverage.

All figures below are unmodified PNGs from the successful final Gadi run. They replace the earlier local carbon-map rerenders.

| Diagnostic | Final Gadi figures |
| --- | --- |
| 1 | [page 1](../cmip6/Fig1_page1.png), [page 2](../cmip6/Fig1_page2.png) |
| 2 | [figure](../cmip6/Fig2.png) |
| 3b | [figure](../cmip6/Fig3b.png) |
| 3b-2 | [figure](../cmip6/Fig3b-2.png) |
| 4 | [page 1](../cmip6/Fig4_page1.png), [page 2](../cmip6/Fig4_page2.png) |
| 5 | [figure](../cmip6/Fig5.png) |
| 5g | [figure](../cmip6/Fig5g.png) |
| 6a | [figure](../cmip6/Fig6a.png) |
| 6b | [figure](../cmip6/Fig6b.png) |
| 7 | [page 1](../cmip6/Fig7_page1.png), [page 2](../cmip6/Fig7_page2.png) |
| 7h | [figure](../cmip6/Fig7h.png) |
| 7i | [figure](../cmip6/Fig7i.png) |
| 8 | [figure](../cmip6/Fig8.png) |
| 9a | [figure](../cmip6/Fig9a.png) |
| 9b | [figure](../cmip6/Fig9b.png) |
| 9c | [figure](../cmip6/Fig9c.png) |


## CMIP7 companion: completed Gadi validation

`recipe_russell18jgr_cmip7.yml` adds a separate historical 1986–2005 companion using ACCESS-ESM1-6, CanESM6-0-MR and UKESM1-3-LL where their fields are available. The CMIP5 and CMIP6 recipe files are unchanged. The [CMIP7 documentation](../../../recipe_russell18jgr_cmip7.rst) records selections, members, published grids and scientific limitations. There are 13 diagnostic groups; Figures 9a–c are omitted because the selected archive does not provide a multi-model set with the required carbon flux and cell areas. CanESM fixed cell areas come from its `1pctCO2` publication on the same `g147` grid; no time-varying fields are borrowed from that experiment.

### Execution and consistency checks

- **PBS `183400498.gadi-pbs`: exit 0**, 7 November 2026, on `analysis3-26.08` with ESMValCore 2.15.0. All **35/35 tasks** completed in 26:25 wall time; peak memory was 81.12 GB.
- **13 diagnostic groups, 17 PNG figures and 30 diagnostic NetCDF files** passed the [output/model/provenance audit](../cmip7/validation_summary.json), with no reported errors. Outputs contain finite data; front-coordinate NaNs intentionally separate disconnected contour segments.
- All **32 preprocessed files retaining time** contain exactly 240 ordered months, January 1986–December 2005: [monthly coverage audit](../cmip7/monthly_coverage_checks.json).
- ACCESS and CanESM sea-ice and cell-area latitude/longitude coordinates match exactly: [grid audit](../cmip7/area_grid_checks.json). This verifies CanESM's reused fixed area grid; it does not validate velocity/tracer staggering in the transport diagnostics.
- The full run used recipe commit `d5786c0` and the diagnostic code tested for CMIP6 at `bce5290`. [Source checksums](../cmip7/sources.sha256) matched the local files. The three recipe-loading/reference/maintainer checks, YAML lint and focused Sphinx build passed locally.
- Nonfatal input warnings remain: missing `positive` attributes for `tauu`, `tauuo` and `fgco2`, and recommendations to supply a land/ocean fraction field for masking. The diagnostics retain their documented variable sign conventions. These warnings are not suppressed or described as CMOR compliance.

Run directory: `/g/data/tm70/rb5533/ESMValTool_Russel/pr4612_cmip7/output/recipe_russell18jgr_cmip7_20261107_011215/`. Job scripts, logs and audits are retained in the stage's `jobs/` directory.

### Scientific interpretation and visual review

ACCESS Figure 6a gives a net climatological estimate of **1.0504988 Sv**, and Figure 6b gives **−0.1853350 PW**. The 11 main density-layer estimates reproduce the earlier same-period climatological check within `4.32e-7 Sv` and `2.63e-8 PW`, respectively: [comparison](../cmip7/climatological_comparison.json). UKESM gives **0.8398934 Sv** and **−0.6031913 PW**. These are consistency checks, **not validation against native monthly transports**. The earlier ACCESS native-flux discrepancies, monthly-covariance effects, density-layer migration and staggered-grid caveats still apply; no independent native-flux validation is claimed for UKESM.

The wind, carbon and pH maps wrap without a spurious radial seam. Front plots retain gaps from the legacy contour selection; a gap does not imply that a physical front is absent. Figure 5 retains the legacy March/September selection. Visual review identified clipping of a UKESM Figure 6b sublayer and the ACCESS Figure 7i curve; their axis ranges were expanded without changing the calculations. The corrected diagnostics (commit `3f61c63`) were rerun on Gadi as PBS `183401838.gadi-pbs`, **exit 0**, in 5:10. Both heat-transport NetCDF arrays and the cumulative carbon-flux array are exactly equal to the full-run outputs. The gallery below uses these corrected plots. The [plotting-rerun audit](../cmip7/plotfix_validation.json) and [tested source hashes](../cmip7/plotfix_sources.sha256) record this final step. Ruff, all 12 Russell diagnostic unit tests and the focused documentation build also passed. An initial plotting job (`183401747`) failed before running a diagnostic because of a batch working-directory typo; the corrected job above supersedes it.

### CMIP7 figures

| Diagnostic | Full figures |
| --- | --- |
| Figure 1, atmospheric wind stress | [UKESM/ACCESS](../cmip7/Fig1_page1.png), [CanESM](../cmip7/Fig1_page2.png) |
| Figure 2, ocean surface wind stress | [Figure 2](../cmip7/Fig2.png) |
| Figure 3b, Subantarctic Front | [Figure 3b](../cmip7/Fig3b.png) |
| Figure 3b-2, Polar Front | [Figure 3b-2](../cmip7/Fig3b-2.png) |
| Figure 4, Drake Passage | [ACCESS/UKESM](../cmip7/Fig4_page1.png), [CanESM velocity only](../cmip7/Fig4_page2.png) |
| Figure 5, sea ice | [Figure 5](../cmip7/Fig5.png) |
| Figure 5g, sea-ice area | [ACCESS/CanESM](../cmip7/Fig5g.png) |
| Figure 6a, climatological volume transport | [ACCESS](../cmip7/Fig6a_ACCESS-ESM1-6.png), [UKESM](../cmip7/Fig6a_UKESM1-3-LL.png) |
| Figure 6b, climatological heat transport | [ACCESS](../cmip7/Fig6b_ACCESS-ESM1-6.png), [UKESM](../cmip7/Fig6b_UKESM1-3-LL.png) |
| Figure 7, carbon flux maps | [ACCESS/UKESM](../cmip7/Fig7.png) |
| Figure 7h, zonal carbon flux | [ACCESS/UKESM](../cmip7/Fig7h.png) |
| Figure 7i, cumulative carbon flux | [ACCESS](../cmip7/Fig7i.png) |
| Figure 8, surface pH | [UKESM](../cmip7/Fig8.png) |
