def obs_filenames():

    """
    List of filenames and paths to obs data required for hydrocycle assessment

    NOTE that some are roots:
        e.g.  dasilva_file = '/project/cma/clim/hydrocycle/lh.' points to a set
              of filenames with seasonal endings (djf, son etc)

    Monthly precip files should contain as many years as required of 12
    months of data.

    """

    # TODO local paths
    gpcp_file = '/project/cma/clim/hydrocycle/gpcp_v2_psg_19792006_monthly.pp'
    cmap_file = '/project/cma/clim/hydrocycle/cmap79_01_monthly.pp'

    dasilva_file_lh = '/project/cma/clim/hydrocycle/dasilva/lh.'
    dasilva_file_sh = '/project/cma/clim/hydrocycle/dasilva/sh.'
    nocs_file = '/project/cma/clim/hydrocycle/nocs2.0_'

    return dict(gpcpm=gpcp_file,
                cmapm=cmap_file,
                dasilvalh=dasilva_file_lh,
                dasilvash=dasilva_file_sh,
                nocss=nocs_file)
