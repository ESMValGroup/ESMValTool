###############################################################################
# namelist_flato13ipcc.yml
---
documentation:
description: |
    Reproducing selected figures from IPCC AR5, chap. 9 (Flato et al., 2013):
    9.42a

authors:
    - schl_ma

references:
    - flato13ipcc

projects:
    - esmval
    - crescendo


preprocessors:

    noop:
        {}


diagnostics:

    fig9-42a:
        description: ECS vs. GMSAT.
        variables:
            tas:
                preprocessor: noop
                mip: Amon
                field: T2Ms
            rtmt:
                preprocessor: noop
                mip: Amon
                field: T2Ms
        additional_models:
            # Day is out of range for month (wait for iris > 2.0)
            # - {model: ACCESS1-0,     project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            # - {model: ACCESS1-0,     project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year:  300, end_year:  449}
            # - {model: ACCESS1-0,     project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year:  300, end_year:  449}
            - {model: bcc-csm1-1,    project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: bcc-csm1-1,    project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year:  160, end_year:  309}
            - {model: bcc-csm1-1,    project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year:  160, end_year:  309}
            - {model: bcc-csm1-1-m,  project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: bcc-csm1-1-m,  project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year:  240, end_year:  389}
            - {model: bcc-csm1-1-m,  project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year:  240, end_year:  389}
            - {model: CanESM2,       project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: CanESM2,       project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 2015, end_year: 2164}
            - {model: CanESM2,       project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: CCSM4,         project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: CCSM4,         project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year:  800, end_year:  949}
            - {model: CCSM4,         project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: CNRM-CM5,      project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: CNRM-CM5,      project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: CNRM-CM5,      project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: CSIRO-Mk3-6-0, project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: CSIRO-Mk3-6-0, project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year:    1, end_year:  150}
            - {model: CSIRO-Mk3-6-0, project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year:    1, end_year:  150}
            - {model: GFDL-CM3,      project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: GFDL-CM3,      project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year:    1, end_year:  150}
            - {model: GFDL-CM3,      project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year:    1, end_year:  150}
            - {model: GISS-E2-H,     project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: GISS-E2-H,     project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 1200, end_year: 1349}
            - {model: GISS-E2-H,     project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: GISS-E2-R,     project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: GISS-E2-R,     project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 3331, end_year: 3480}
            - {model: GISS-E2-R,     project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: inmcm4,        project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: inmcm4,        project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 2090, end_year: 2239}
            - {model: inmcm4,        project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 2090, end_year: 2239}
            - {model: IPSL-CM5A-LR,  project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: IPSL-CM5A-LR,  project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: IPSL-CM5A-LR,  project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: IPSL-CM5B-LR,  project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: IPSL-CM5B-LR,  project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: IPSL-CM5B-LR,  project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}
            - {model: MIROC-ESM,     project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: MIROC-ESM,     project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 1800, end_year: 1949}
            - {model: MIROC-ESM,     project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year:    1, end_year:  150}
            - {model: MIROC5,        project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: MIROC5,        project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 2100, end_year: 2249}
            - {model: MIROC5,        project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 2100, end_year: 2249}
            - {model: MPI-ESM-LR,    project: CMIP5, exp: historical,  ensemble: r1i1p1, start_year: 1961, end_year: 1990}
            - {model: MPI-ESM-LR,    project: CMIP5, exp: piControl,   ensemble: r1i1p1, start_year: 2015, end_year: 2164}
            - {model: MPI-ESM-LR,    project: CMIP5, exp: abrupt4xCO2, ensemble: r1i1p1, start_year: 1850, end_year: 1999}

        scripts:
            fig9-42a:
                script: ipcc_ar5/ch09_fig09-42a.py
                plot_ecs_regression: True
