# xml2yml Converting Version 1 Namelists to Version 2 Namelists

This converter can turn the old xml namelists into new-style yml namelists.
It is implemented as a xslt stylesheet that needs a processor that is xslt 2.0 capable.
With this, you simply process your old namelist with the stylesheet xml2yml.xsl to produce
a new yml namelist.

After the conversion you need to manually check the mip information in the variables!
Also, check the caveats below!

## Howto
One freely available processor is the Java based [saxon](http://saxon.sourceforge.net/).
You can download the free he edition [here](https://sourceforge.net/projects/saxon/files/latest/download?source=files).
Unpack the zip file into a new directory. Then, provided you have Java installed, you can convert your namelist
simply with:
```
java -jar $SAXONDIR/saxon9he.jar -xsl:xml2yml.xsl -s:namelist.xml -o:namelist.yml
```

## Caveats/Known Limitations
* At the moment, not all model schemes (OBS, CMIP5, CMIP5_ETHZ...) are supported.
  They are, however, relatively easy to add, so if you need help adding a new one,
  please let me know!
* The documentation section (namelist_summary in the old file) is not
  automatically converted.
* In version 1, one could name an exclude, similar to the reference model. This
  is no longer possible and the way to do it is to include the models with
  another `additional_models` tag in the variable section. That conversion is
  not performed by this tool.

## Author
* **Klaus Zimmermann**, direct questions and comments to klaus.zimmermann@smhi.se
