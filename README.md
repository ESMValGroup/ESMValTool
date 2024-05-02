# Instructions to run FWI calculations

## General information 

The instructions presented here are applicable only for the calculations of Fire Weather Indices (FWI) for attribution of extreme 2023 wildfire season in Canada (Kirchmeier-Young et al., 2024 in prep). 

## Installation

Please follow installation guide from 
https://docs.esmvaltool.org/en/latest/quickstart/installation.html#install-from-source


## ESMValCore Version 

While preparing the data for publication, we have made certain modifications (windspeed derivation, certain model fixes) to ESMValCore, which were sucessully implemeted into the ``main`` branch thus resulting in the original repositories of ESMValCore used for analysis being deleted. 

However, code should work for ``esmvalcore=2.11`` with additional installation of 
``xclim`` package. In case any issues with running of the FWIs will arise, please, contact [Liza](mailto:elizaveta.malinina-rieger@ec.gc.ca). 

The repository used for relative humidity derivation from ERA5 is still not merged and can be found [here](https://github.com/ESMValGroup/ESMValCore/tree/dev_derive_hurs). 

## How to reproduce the calculations

### ERA5 

1. Load ERA5 data using ``era5cli`` software
2. CMORIze ERA5 data using the recipes in **esmvaltool/recipes/cmorizers/recipe_daily_era5_{variable}.yml**.
It is recommended to either use dask settings from ESMValTool or CMORize data year by year not to run into memory issues. 
3. Run the FWI calculations by running **esmvaltool/recipes/recipe_canada_fwi.yml**. It is recommended as well to 
run the recipes year by year, because ``xclim`` FWI calculations are very memory heavy and can't be daskified. 

### CMIP6 

0. (If CMIP6 data is not available in your computer you can use ``search_esgf: True`` in your ``config-user.yml`` to download the data while running the recipes. For more information please see [official ESMValTool documentation](https://docs.esmvaltool.org/en/latest/index.html). But, please, be aware that loading of HighResMIP data in particular, can last very long time. Alternatively any method to download CMIP6 data is fine. We've used ``pyesgf`` for data download.)
1. Run the recipe for DAMIP models **esmvaltool/recipes/recipe_canada_fwi_damip.yml**
2. Run the recipe for HighResMIP models **esmvaltool/recipes/recipe_canada_fwi_highresmip.yml**

**NB:** If running into the memory issues with the FWI calculations, it is recommended to shorten the periods, e.g., running the recipes in 10-year chunks. 

The result of these calculations will be ``.nc`` files with multiple fire weather indices and regression predictors. We do not share the code for attribution itself. 

The code is distributed under [CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/) licence. 