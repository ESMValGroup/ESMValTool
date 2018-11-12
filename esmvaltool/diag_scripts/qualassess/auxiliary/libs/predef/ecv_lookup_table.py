def ecv_lookup(ecv_name):
    ecv_dict = {
        "albdiffbnd": "diffuse surface albedo for each band",
        "albdirbnd": "direct surface albedo for each band",
        "albisccp": "single scattering cloud albedo",
        "burntArea": "fraction of grid cell that is covered by burnt vegetation",
        "cct": "air pressure at convective cloud top",
        "cfc11": "mole concentration of CFC-11 in sea water",
        "cfc12": "mole concentration of CFC-12 in sea water",
        "clt": "total cloud fraction",
        "cod": "atmosphere optical thickness due to cloud",
        "hfds": "surface downward heat flux in sea water",
        "hur": "relative humidity",
        "intdic": "dissolved inorganic carbon content",
        "lai": "leaf area index",
        "lwp": "cloud liquid water path",
        "mrsos": "moisture in upper portion of soil column",
        "od550aer": "aerosol optical depth at 550nm",
        "phnat": "natural pH",
        "pr": "precipitation",
        "reffclwtop": "cloud-top effective droplet radius",
        "rhs": "near-surface relative humidity",
        "rlds": "surface net downwelling longwave radiation",
        "rldscs": "surface downwelling longwave flux in air assuming clear sky",
        "rlus": "surface upwelling longwave radiation",
        "rluscs": "surface upwelling longwave flux in air assuming clear sky",
        "rlut": "outgoing longwave flux at TOA",
        "rsds": "surface downwelling shortwave radiation",
        "rsdscs": "surface downwelling shortwave flux in air assuming clear sky",
        "rsdt": "incoming shortwave flux at TOA",
        "rsus": "surface upwelling shortwave radiation",
        "rsuscs": "surface upwelling shortwave flux in air assuming clear sky",
        "rsut": "outgoing shortwave flux at TOA",
        "sf6": "mole concentration of SF6 in sea water",
        "sfcWind": "near-surface wind speed",
        "sftgif": "fraction of grid cell covered with glacier",
        "sic": "sea ice concentration",
        "sidivvel": "divergence of sea ice velocity",
        "sidmassdyn": "tendency of sea ice amount due to dynamics",
        "sidmassth": "tendency of sea ice amount due to thermodynamics",
        "siextentn": "sea ice extent north",
        "siextents": "sea ice extent south",
        "sispeed": "sea ice speed",
        "sit": "sea ice thickness",
        "siu": "sea ice x velocity",
        "siv": "sea ice y velocity",
        "sm": "volumetric soil moisture",
        "snc": "area covered by snow",
        "snd": "snow depth",
        "swe": "snow water equivalent",
        "so": "sea water salinity",
        "spco2": "surface partial pressure of carbon dioxide in sea water",
        "ta": "temperature",
        "talknat": "natural total alkalinity",
        "tas": "surface temperature",
        "tasmax": "daily maximum surface temperature",
        "tasmin": "daily minimum surface temperature",
        "tauuo": "surface downward X stress",
        "tauvo": "surface downward Y stress",
        "tos": "sea surface temperature",
        "toz": "total column ozone",
        "tpf": "permafrost layer thickness",
        "ts": "surface temperature (Model - do not use)",
        "ua": "eastward wind",
        "uas": "eastward near-surface wind",
        "va": "northward wind",
        "vas": "northward near-surface wind",
        "vmrch4": "CH4 mixing ratio",
        "vmrco2": "CO2 mixing ratio",
        "vmrstrato3": "Ozone mixing ratio in the stratosphere",
        "vmro3": "Ozone mixing ratio",
        "xco2": "CO2 column",
        "xch4": "CH4 column",
        "zos": "sea surface height above geoid",
    }

    return ecv_dict[ecv_name]
