######################################################
#-------------ETCCDI preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
############################################################################
# ABOUT: This function pre-process ETCCDI files obtained with the
#        extreme_events recipe remapping the data from
#        gaussian to lonlat, changing longitude range from 0/360 to -180/180
#        and merging all indices into the HyInt indices file.

hyint_etccdi_preproc <-
  function(work_dir,
             etccdi_dir,
             etccdi_list_import,
             model_idx,
             season,
             prov_info,
             yrmon = "yr") {
    year1 <- toString(models_start_year[model_idx])
    year2 <- toString(models_end_year[model_idx])
    hyint_file <-
      getfilename_indices(work_dir, diag_base, model_idx, season)
    etccdi_files <-
      getfilename_etccdi(etccdi_dir, etccdi_list_import, model_idx,
        yrmon = "yr"
      )
    etccdi_files_tmp <- c()
    for (sfile in etccdi_files) {
      print(paste0("HyInt: pre-processing ", sfile))
      sfile_tmp0 <- cdo("delvar", args = "time_bnds", input = sfile)
      gridf <- tempfile()
      cdo("griddes", input = hyint_file, stdout = gridf)
      sfile_tmp1 <- cdo("remapcon2",
        args = gridf,
        input = sfile_tmp0
      )
      sfile_tmp <- cdo("sellonlatbox",
        args = "-180,180,-90,90",
        input = sfile_tmp1
      )
      etccdi_files_tmp <- c(etccdi_files_tmp, sfile_tmp)
      unlink(c(sfile_tmp0, sfile_tmp1))
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
    cdo("merge",
      options = "-O",
      input = c(hyint_file_tmp_sel, etccdi_files_tmp),
      output = hyint_file
    )
    unlink(c(etccdi_files_tmp, hyint_file_tmp, hyint_file_tmp_sel))

    # Update provenance with etccdi files
    prov_info[[hyint_file]]$ancestors <-
          c(prov_info[[hyint_file]]$ancestors, etccdi_files)

    return(prov_info)
  }
