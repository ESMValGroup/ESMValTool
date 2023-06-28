Code to produce figures from Gillett (2023) - 'Is warming proportional to cumulative carbon emissions because heat and carbon are mixed into the ocean by the same processes?'

Code versions used:
ESMValCore: 2.5.0
ESMValTool: 2.5.0

For information installing and running ESMValTool, please see:
https://docs.esmvaltool.org/en/latest/quickstart/index.html
https://tutorial.esmvaltool.org/index.html
Typical install time 30-60 minutes.

The recipes to produce the figures for Gillett (2023) are:
esmvaltool run esmvaltool/recipes/recipe_proportionality.yml
esmvaltool run esmvaltool/recipes/recipe_proportionality_zec.yml

(Note - you will first need to download the CMIP6 data which is used in these recipes if it isn't already available locally).

