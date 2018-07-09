<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="2.0">
  <xsl:output method="text" omit-xml-declaration="yes" indent="no"/>
  <xsl:strip-space elements="*"/>
  <xsl:template match="namelist_summary">
    <xsl:for-each select="tokenize(., '\r?\n')">
      <xsl:text># </xsl:text>
      <xsl:sequence select="."/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:for-each>
    <xsl:text>&#xa;---&#xa;</xsl:text>
  </xsl:template>

  <xsl:template match="GLOBAL">
    <!-- <xsl:for-each select="*"> -->
    <!--   <xsl:text>  </xsl:text> -->
    <!--   <xsl:value-of select="local-name()"/>: <xsl:value-of select="replace(., '^\s+|\s+$', '')"/> -->
    <!--   <xsl:text>&#xa;</xsl:text> -->
    <!-- </xsl:for-each> -->
  </xsl:template>

  <xsl:template match="MODELS">
    <xsl:if test="count(model)>0">
      <xsl:text>datasets:&#xa;</xsl:text>
    </xsl:if>
  </xsl:template>

  <xsl:template match="DIAGNOSTICS">
    <xsl:text>preprocessors:&#xa;</xsl:text>
    <xsl:text>  prep0: {}&#xa;&#xa;</xsl:text>
    <xsl:text>diagnostics:&#xa;</xsl:text>
    <xsl:apply-templates select="diag"/>
  </xsl:template>

  <xsl:template match="diag">
    <xsl:variable name="diag_name" select="normalize-space(tokenize(diag_script[1], '\.')[1])"/>
    <!-- Put diag_script as temporary name -->
    <xsl:text>  </xsl:text>
    <xsl:value-of select="$diag_name"/>
    <xsl:text>:&#xa;</xsl:text>
    <!-- Add description -->
    <xsl:text>    description: </xsl:text>
    <xsl:value-of select="normalize-space(description)"/>
    <xsl:text>&#xa;</xsl:text>
    <!-- Add variable information -->
    <xsl:text>    variables:&#xa;</xsl:text>
    <xsl:variable name="no_field_types" select="count(field_type)"/>
    <xsl:for-each select="variable">
      <xsl:variable name="ind" select="min((position(), $no_field_types))" />
      <!-- Add variable name -->
      <xsl:text>      </xsl:text>
      <xsl:value-of select="normalize-space()"/>
      <xsl:text>:&#xa;</xsl:text>
      <!-- Add default preprocessor -->
      <xsl:text>        preprocessor: prep0&#xa;</xsl:text>
      <!-- Add references, excludes, and mip-->
      <xsl:apply-templates select="@*"/>
      <xsl:choose>
        <xsl:when test="./@MIP">
          <xsl:text>        mip: </xsl:text>
          <xsl:value-of select="."/>
          <xsl:text>&#xa;</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>        mip: #missing&#xa;</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
      <!-- Add field type -->
      <xsl:text>        field: </xsl:text>
      <xsl:value-of select="normalize-space(../field_type[position()=$ind])"/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:for-each>
    <!-- Add datasets -->
    <xsl:if test="count(model)>0">
      <xsl:text>    additional_datasets:&#xa;</xsl:text>
      <xsl:apply-templates select="model"/>
    </xsl:if>
    <!-- Add scripts -->
    <xsl:text>    scripts:&#xa;</xsl:text>
    <xsl:for-each select="diag_script">
      <xsl:variable name="script" select="normalize-space(.)"/>
      <xsl:text>      </xsl:text>
      <xsl:value-of select="tokenize($script, '\.')[1]"/>
      <xsl:text>:&#xa;        script: </xsl:text>
      <xsl:value-of select="$script"/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:for-each>
    <xsl:text>&#xa;</xsl:text>
  </xsl:template>

  <xsl:template match="@ref_model">
    <xsl:variable name="ref_models" select="tokenize(normalize-space(.), ',')"/>
    <xsl:text>        reference_dataset: </xsl:text>
    <xsl:value-of select="normalize-space($ref_models[1])"/>
    <xsl:text>&#xa;</xsl:text>
    <xsl:if test="count($ref_models)>1">
      <xsl:text>        alternative_dataset: </xsl:text>
      <xsl:value-of select="normalize-space($ref_models[2])"/>
      <xsl:text>&#xa;</xsl:text>
    </xsl:if>
  </xsl:template>

  <xsl:template match="@exclude">
    <xsl:text>        exclude_dataset: </xsl:text>
    <xsl:value-of select="."/>
    <xsl:text>&#xa;</xsl:text>
  </xsl:template>

  <xsl:template match="model">
    <xsl:variable name="model_line" select="tokenize(normalize-space(.))"/>
    <xsl:variable name="dataset_scheme" select="$model_line[1]"/>
    <xsl:choose>
      <!-- <xsl:when test="$dataset_scheme='CMIP5'"> -->
      <!-- </xsl:when> -->
      <xsl:when test="$dataset_scheme='CMIP5_ETHZ'">
        <xsl:text>      - {dataset: </xsl:text>
        <xsl:value-of select="$model_line[2]"/>
        <xsl:text>, project: CMIP5</xsl:text>
        <xsl:text>, exp: </xsl:text>
        <xsl:value-of select="$model_line[4]"/>
        <xsl:text>, ensemble: </xsl:text>
        <xsl:value-of select="$model_line[5]"/>
        <xsl:text>, start_year: </xsl:text>
        <xsl:value-of select="$model_line[6]"/>
        <xsl:text>, end_year: </xsl:text>
        <xsl:value-of select="$model_line[7]"/>
        <xsl:text>}&#xa;</xsl:text>
      </xsl:when>
      <xsl:when test="$dataset_scheme='OBS'">
        <xsl:text>      - {dataset: </xsl:text>
        <xsl:value-of select="$model_line[2]"/>
        <xsl:text>, project: </xsl:text>
        <xsl:value-of select="$model_line[1]"/>
        <xsl:text>, type: </xsl:text>
        <xsl:value-of select="$model_line[3]"/>
        <xsl:text>, version: </xsl:text>
        <xsl:value-of select="$model_line[4]"/>
        <xsl:text>, start_year: </xsl:text>
        <xsl:value-of select="$model_line[5]"/>
        <xsl:text>, end_year: </xsl:text>
        <xsl:value-of select="$model_line[6]"/>
        <xsl:text>, tier: </xsl:text>
        <xsl:value-of select="substring(tokenize($model_line[7], '/')[2], 5, 1)"/>
        <xsl:text>}&#xa;</xsl:text>
      </xsl:when>
      <!-- <xsl:when test="$dataset_scheme='OBS_gridfile'"> -->
      <!-- </xsl:when> -->
      <!-- <xsl:when test="$dataset_scheme='obs4mips'"> -->
      <!-- </xsl:when> -->
      <!-- <xsl:when test="$dataset_scheme='ana4mips'"> -->
      <!-- </xsl:when> -->
      <xsl:otherwise>
        <xsl:message terminate="yes">
          ERROR: unknown dataset scheme <xsl:value-of select="$dataset_scheme"/>
        </xsl:message>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

</xsl:stylesheet>
