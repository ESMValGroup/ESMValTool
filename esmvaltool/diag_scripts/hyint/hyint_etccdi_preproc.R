######################################################
#-------------ETCCDI preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
############################################################################
# ABOUT: This function pre-process ETCCDI files obtained with the
#        CRESCENDO_extremeEvents namelist remapping the data from
#        gaussian to lonlat, changing longitude range from 0/360 to -180/180
#        and merging all indices into the HyInt indices file.

hyint_etccdi_preproc <-
  function(work_dir,
             etccdi_dir,
             etccdi_list_import,
             cdo_grid,
             model_idx,
             season,
             yrmon = "yr") {
    year1 <- toString(models_start_year[model_idx])
    year2 <- toString(models_end_year[model_idx])
    print(str(c(year1, year2)))
    hyint_file <-
      getfilename_indices(work_dir, diag_base, model_idx, season)
    etccdi_files <-
      getfilename_etccdi(etccdi_dir, etccdi_list_import, model_idx,
        yrmon = "yr"
      )
    etccdi_files_tmp <- c()
    for (sfile in etccdi_files) {
      sfile_tmp0 <- cdo("delvar", args = "time_bnds", input = sfile)
      if (rgrid != F) {
        sfile_tmp <- cdo("setgrid", args = cdo_grid, input = sfile_tmp0)
      } else {
        sfile_tmp <- cdo("sellonlatbox",
          args = "-180,180,-90,90",
          input = sfile_tmp0
        )
      }
      etccdi_files_tmp <- c(etccdi_files_tmp, sfile_tmp)
      unlink(sfile_tmp0)
    }
    hyint_file_tmp <- tempfile()
    mv_command <- paste("mv -n ", hyint_file, hyint_file_tmp)
    system(mv_command)
    print(paste0("HyInt: merging ", length(etccdi_files), " ETCCDI files"))
    hyint_file_tmp_sel <-
      cdo("sellonlatbox",
        args = "-180,180,-90,90",
        input = hyint_file_tmp
      )
    cdo(
      "merge",
      options = "-O",
      input = c(hyint_file_tmp_sel, etccdi_files_tmp),
      output = hyint_file
    )
    unlink(c(etccdi_files_tmp, hyint_file_tmp, hyint_file_tmp_sel))

    return(0)
  }
