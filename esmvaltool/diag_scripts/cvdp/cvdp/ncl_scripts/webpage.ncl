; Create the index.html, methodology.html, and metrics.html web pages.
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: webpage.ncl")

  OUTDIR      = getenv("OUTDIR") 
  VERSION     = getenv("VERSION")
  OT          = getenv("OUTPUT_DATA")
  OPT_CLIMO   = getenv("OPT_CLIMO")
  CLIMO_SYEAR = getenv("CLIMO_SYEAR")
  CLIMO_EYEAR = getenv("CLIMO_EYEAR")
  OBS         = getenv("OBS")

  quote = str_get_dq()
  
  if (OPT_CLIMO.eq."Full") then
     subtxt = "Climatological Period Used: Full"
  end if
  if (OPT_CLIMO.eq."Custom") then
     if (toint(CLIMO_SYEAR).lt.0) then
        subtxt = "Climo. Period (relative to record end): "+CLIMO_SYEAR+"-"+CLIMO_EYEAR
     else
        subtxt = "Climatological Period Used: "+CLIMO_SYEAR+"-"+CLIMO_EYEAR
     end if
  end if
  if (OBS.eq."True") then
     namelist_txt = "Input Namelists: <a href="+quote+"namelist_obs"+quote+">Observations</a> | <a href="+quote+"namelist"+quote+">Models</a>"
  else
     namelist_txt = "Input Namelist: <a href="+quote+"namelist"+quote+">Models</a>"
  end if
  nsim_noobs = numAsciiRow("namelist")
  na = asciiread("namelist",(/nsim_noobs/),"string")
  names = new(nsim_noobs,"string")
  delim = "|"
  do gg = 0,nsim_noobs-1
     names(gg) = str_strip(str_get_field(na(gg),1,delim))
  end do
;  print(names)
  
  txt = new(700,"string")
  quote = str_get_dq()
  
  txt(0) = "<HTML><HEAD><TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
  txt(1) = "body {background-color: WhiteSmoke; font-family: Verdana, Geneva, sans-serif}"
  txt(2) = "img.title {width:990; height:105; align:left}"
  txt(3) = "p.subtitle {font-size:20pt; color:black; text-align:center; font-weight:bold}"
  txt(4) = "p.btext {font-size:11pt; text-align:left}"
  txt(5) = "font.s1 {font-size:13.5pt; text-align:left}"
  txt(7) = "td {text-align:center; font-size:12pt}"
  txt(8) = "td.s1 {text-align:left; font-size:11pt}"
  txt(9) = "td.s2 {border-bottom: 4px double darkgray}"
  txt(10) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
  txt(11) = "</STYLE></HEAD><BODY>"
  txt(12) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"
  txt(13) = "<TABLE width='990' border=0 padding=0 cellpadding=0>"
  txt(14) = "<tr><td width='365' colspan=2><p class=btext><font class=s1><a href="+quote+"methodology.html"+quote+">Methodology and Definitions</a></font>"
  if (isfilepresent2(OUTDIR+"metrics.txt")) then
     if (isfilepresent2(OUTDIR+"metrics.table_1.gif")) then
        txt(15) = "<br>Metrics Tables: <a href="+quote+"metrics.table_0.html"+quote+">Pattern Correlations</a> | <a href="+quote+"metrics.table_14.html"+quote+">RMS</a>"
     else
        txt(15) = "<br><a href="+quote+"metrics.html"+quote+">Metrics Table</a>"
     end if
  end if
  txt(16) = "<br>"+namelist_txt+"</td><td rowspan=3><p class=subtitle>"+webtitle+"</p></td></tr>"
  txt(17) = "<tr><td class=s1 width='90'>Derived<br>Namelists:</td><td class=s1 width='260'><a href="+quote+"namelist_moc"+quote+">MOC</a> | <a href="+quote+"namelist_prect"+quote+">PR</a> | <a href="+quote+"namelist_psl"+quote+">PSL</a> | <a href="+quote+"namelist_aice_nh"+quote+">SIC NH</a><br><a href="+quote+"namelist_aice_sh"+quote+">SIC SH</a> | <a href="+quote+"namelist_snowdp"+quote+">SND</a> | <a href="+quote+"namelist_trefht"+quote+">TAS</a> | <a href="+quote+"namelist_ts"+quote+">TS</a></td></tr>"
  txt(18) = "<tr><td class=s1 width='365' colspan=2>"+subtxt+"<br>Created: "+systemfunc("date")+"<br>CVDP Version "+VERSION+"</td></tr></table><hr width=990 align=left>"
  delete([/subtxt,namelist_txt/])

  txt(26) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(27) = "<CAPTION>Means</CAPTION>"
  txt(28) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.ann.png","Annual")+"</td></tr>"    
  txt(29) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.ann.png","Annual")+"</td></tr>"    
  txt(30) = "<tr><td><b>PSL</b></td><td>"+table_link_setup(OUTDIR,"psl.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.ann.png","Annual")+"</td></tr>"    
  txt(31) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.ann.png","Annual")+"</td></tr>"    
  txt(32) = "<tr><td><b>[PR]</b></td><td>"+table_link_setup(OUTDIR,"pr.za.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pr.za.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"pr.za.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pr.za.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"pr.za.ann.png","Annual")+"</td></tr>"
  txt(33) = "<tr><td><b>SND</b></td><td>"+table_link_setup(OUTDIR,"snd.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"snd.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"snd.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"snd.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"snd.mean.ann.png","Annual")+"</td></tr>"      
  txt(34) = "<tr><td><b>SIC NH</b></td><td>"+table_link_setup(OUTDIR,"aice.mean.nh.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.nh.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.nh.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.nh.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.nh.ann.png","Annual")+"</td></tr>"      
  txt(35) = "<tr><td><b>SIC SH</b></td><td>"+table_link_setup(OUTDIR,"aice.mean.sh.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.sh.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.sh.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.sh.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"aice.mean.sh.ann.png","Annual")+"</td></tr>"      
  txt(36) = "</table>"  
  
  txt(39) = ""
  txt(40) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(41) = "<CAPTION>Standard Deviations</CAPTION>"
  txt(42) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.ann.png","Annual")+"</td></tr>"    
  txt(43) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.ann.png","Annual")+"</td></tr>"    
  txt(44) = "<tr><td><b>PSL</b></td><td>"+table_link_setup(OUTDIR,"psl.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.ann.png","Annual")+"</td></tr>"    
  txt(45) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.ann.png","Annual")+"</td></tr>"    
  txt(46) = "<tr><td><b>SND</b></td><td>"+table_link_setup(OUTDIR,"snd.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"snd.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"snd.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"snd.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"snd.stddev.ann.png","Annual")+"</td></tr>"    
  txt(47) = "<tr><td><b>SIC NH</b></td><td>"+table_link_setup(OUTDIR,"aice.stddev.nh.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.nh.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.nh.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.nh.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.nh.ann.png","Annual")+"</td></tr>"      
  txt(48) = "<tr><td><b>SIC SH</b></td><td>"+table_link_setup(OUTDIR,"aice.stddev.sh.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.sh.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.sh.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.sh.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"aice.stddev.sh.ann.png","Annual")+"</td></tr>"      
  txt(49) = "</table>"  
  
  txt(69) = ""
  txt(70) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(71) = "<CAPTION>Coupled Modes of Variability</CAPTION>"  
  txt(72) = "<tr><td class=s2 rowspan=5><b>ENSO</b></td>"
  txt(73) = "<td rowspan=2>Spatial Composites</td><td>SST/TAS/PSL<br>"+table_link_setup(OUTDIR,"nino34.spatialcomp.jja0.png","JJA<sup>0</sup>")+" "+table_link_setup(OUTDIR,"nino34.spatialcomp.son0.png","SON<sup>0</sup>")+"<br>"+table_link_setup(OUTDIR,"nino34.spatialcomp.djf1.png","DJF<sup>+1</sup>")+" "+table_link_setup(OUTDIR,"nino34.spatialcomp.mam1.png","MAM<sup>+1</sup>")+"</td>"
  txt(74) = "<td>PR<br>"+table_link_setup(OUTDIR,"nino34.spatialcomp.pr.jja0.png","JJA<sup>0</sup>")+" "+table_link_setup(OUTDIR,"nino34.spatialcomp.pr.son0.png","SON<sup>0</sup>")+"<br>"+table_link_setup(OUTDIR,"nino34.spatialcomp.pr.djf1.png","DJF<sup>+1</sup>")+" "+table_link_setup(OUTDIR,"nino34.spatialcomp.pr.mam1.png","MAM<sup>+1</sup>")+"</td></tr>"
  txt(75) = "<tr><td>"+table_link_setup(OUTDIR,"nino34.hov.elnino.png","El Ni&ntilde;o Hovm&ouml;ller")+"</td><td>"+table_link_setup(OUTDIR,"nino34.hov.lanina.png","La Ni&ntilde;a Hovm&ouml;ller")+"</td></tr>"
  txt(76) = "<tr class=s1><td class=s2 rowspan=3>Ni&ntilde;o3.4</td><td>"+table_link_setup(OUTDIR,"nino34.timeseries.png","Timeseries")+"</td><td>"+table_link_setup(OUTDIR,"nino34.powspec.png","Power Spectra")+"</td></tr>"
  txt(77) = "<tr><td>"+table_link_setup(OUTDIR,"nino34.monstddev.png","Monthly  Std. Dev.")+"</td><td>"+table_link_setup(OUTDIR,"nino34.runstddev.png","Running Std. Dev.")+"</td></tr>"
  txt(78) = "<tr><td class=s2>"+table_link_setup(OUTDIR,"nino34.autocor.png","Autocorrelation")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nino34.wavelet.png","Wavelet Analysis")+"</td></tr>"

  txt(90) = "<tr><td><b>PDO</b></td>"  
  txt(91) = "<td>Regr: &nbsp;"+table_link_setup(OUTDIR,"pdo.png","SST")+"&nbsp;"+table_link_setup(OUTDIR,"pdo.tasreg.png","TAS")+"&nbsp;"+table_link_setup(OUTDIR,"pdo.prreg.png","PR")+"</td>"
  txt(92) = "<td>"+table_link_setup(OUTDIR,"pdo.timeseries.png","Timeseries")+"</td>"
  txt(93) = "<td>"+table_link_setup(OUTDIR,"pdo.powspec.png","Power Spectra")+"</td></tr>" 
  txt(94) = "<tr><td class=s2><b>IPO</b></td>"  
  txt(95) = "<td class=s2>Regr: &nbsp;"+table_link_setup(OUTDIR,"ipo.png","SST")+"&nbsp"+table_link_setup(OUTDIR,"ipo.tasreg.png","TAS")+"&nbsp;"+table_link_setup(OUTDIR,"ipo.prreg.png","PR")+"</td>"
  txt(96) = "<td class=s2>"+table_link_setup(OUTDIR,"ipo.timeseries.png","Timeseries")+"</td>"
  txt(97) = "<td class=s2>"+table_link_setup(OUTDIR,"ipo.powspec.png","Power Spectra")+"</td></tr>" 

  txt(101) = "<tr><td rowspan=2><b>AMO</b></td>"  
  txt(102) = "<td>Regr: &nbsp;"+table_link_setup(OUTDIR,"amo.png","SST")+"&nbsp;"+table_link_setup(OUTDIR,"amo.tasreg.png","TAS")+"&nbsp;"+table_link_setup(OUTDIR,"amo.prreg.png","PR")+"</td>"
  txt(103) = "<td rowspan=2>"+table_link_setup(OUTDIR,"amo.timeseries.png","Timeseries")+"</td>"
  txt(104) = "<td rowspan=2>"+table_link_setup(OUTDIR,"amo.powspec.png","Power Spectra")+"</td></tr>" 
  txt(105) = "<tr><td>Regr LP: &nbsp;"+table_link_setup(OUTDIR,"amo.lp.png","SST")+"&nbsp;"+table_link_setup(OUTDIR,"amo.lp.tasreg.png","TAS")+"&nbsp;"+table_link_setup(OUTDIR,"amo.lp.prreg.png","PR")+"</td></tr>"
  txt(106) = "<tr><td class=s2 rowspan=3><b>AMOC</b></td>"  
  txt(107) = "<td>"+table_link_setup(OUTDIR,"amoc.mean.ann.png","Means")+"</td>"
  txt(108) = "<td>"+table_link_setup(OUTDIR,"amoc.stddev.ann.png","Standard Deviations")+"</td>"
  txt(109) = "<td>"+table_link_setup(OUTDIR,"amoc.ann.png","Patterns")+"</td></tr>"
  txt(110) = "<tr><td>"+table_link_setup(OUTDIR,"amoc.timeseries.ann.png","Timeseries")+"</td>"
  txt(111) = "<td>"+table_link_setup(OUTDIR,"amoc.sstreg.ann.png","SST Regressions")+"</td>"
  txt(112) = "<td>"+table_link_setup(OUTDIR,"amoc.tasreg.ann.png","TAS Regressions")+"</td></tr>"
  txt(113) = "<tr><td class=s2>"+table_link_setup(OUTDIR,"amoc.powspec.ann.png","Spectra")+"</td>"
  txt(114) = "<td class=s2 colspan=2>"+table_link_setup(OUTDIR,"amoc_amo.leadlag.ann.png","AMO/AMOC Lag Correlations")+"</td></tr></table>" 



  txt(130) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray><CAPTION>Atmospheric Modes of Variability</CAPTION>"  
  txt(131) = "<tr><td rowspan=4 class=s2><b>NAM</b></td>"
  txt(132) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"nam.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nam.mam.png","MAM")+"</td>"
  txt(133) = "<td>"+table_link_setup(OUTDIR,"nam.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nam.son.png","SON")+"</td>"
  txt(134) = "<td>"+table_link_setup(OUTDIR,"nam.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nam.mon.png","Monthly")+"</td></tr>"
  txt(135) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.mam.png","MAM")+"</td>"
  txt(136) = "<td>"+table_link_setup(OUTDIR,"nam.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.son.png","SON")+"</td>"
  txt(137) = "<td>"+table_link_setup(OUTDIR,"nam.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(138) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"nam.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nam.tempreg.mam.png","MAM")+"</td>"
  txt(139) = "<td>"+table_link_setup(OUTDIR,"nam.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nam.tempreg.son.png","SON")+"</td>"
  txt(140) = "<td>"+table_link_setup(OUTDIR,"nam.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nam.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(141) = "<tr><td class=s2>PR Regressions</td><td class=s2>"+table_link_setup(OUTDIR,"nam.prreg.djf.png","DJF")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nam.prreg.mam.png","MAM")+"</td>"
  txt(142) = "<td class=s2>"+table_link_setup(OUTDIR,"nam.prreg.jja.png","JJA")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nam.prreg.son.png","SON")+"</td>"
  txt(143) = "<td class=s2>"+table_link_setup(OUTDIR,"nam.prreg.ann.png","Annual")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nam.prreg.mon.png","Monthly")+"</td></tr>"

  txt(159) = ""
  txt(160) = "<tr><td rowspan=4 class=s2><b>SAM</b></td>"
  txt(161) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"sam.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sam.mam.png","MAM")+"</td>"
  txt(162) = "<td>"+table_link_setup(OUTDIR,"sam.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sam.son.png","SON")+"</td>"
  txt(163) = "<td>"+table_link_setup(OUTDIR,"sam.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"sam.mon.png","Monthly")+"</td></tr>"
  txt(164) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.mam.png","MAM")+"</td>"
  txt(165) = "<td>"+table_link_setup(OUTDIR,"sam.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.son.png","SON")+"</td>"
  txt(166) = "<td>"+table_link_setup(OUTDIR,"sam.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(167) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"sam.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sam.tempreg.mam.png","MAM")+"</td>"
  txt(168) = "<td>"+table_link_setup(OUTDIR,"sam.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sam.tempreg.son.png","SON")+"</td>"
  txt(169) = "<td>"+table_link_setup(OUTDIR,"sam.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"sam.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(170) = "<tr><td class=s2>PR Regressions</td><td class=s2>"+table_link_setup(OUTDIR,"sam.prreg.djf.png","DJF")+"</td><td class=s2>"+table_link_setup(OUTDIR,"sam.prreg.mam.png","MAM")+"</td>"
  txt(171) = "<td class=s2>"+table_link_setup(OUTDIR,"sam.prreg.jja.png","JJA")+"</td><td class=s2>"+table_link_setup(OUTDIR,"sam.prreg.son.png","SON")+"</td>"
  txt(172) = "<td class=s2>"+table_link_setup(OUTDIR,"sam.prreg.ann.png","Annual")+"</td><td class=s2>"+table_link_setup(OUTDIR,"sam.prreg.mon.png","Monthly")+"</td></tr>"

  txt(189) = ""
  txt(190) = "<tr><td rowspan=4 class=s2><b>NAO</b></td>"
  txt(191) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"nao.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nao.mam.png","MAM")+"</td>"
  txt(192) = "<td>"+table_link_setup(OUTDIR,"nao.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nao.son.png","SON")+"</td>"
  txt(193) = "<td>"+table_link_setup(OUTDIR,"nao.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nao.mon.png","Monthly")+"</td></tr>"
  txt(194) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.mam.png","MAM")+"</td>"
  txt(195) = "<td>"+table_link_setup(OUTDIR,"nao.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.son.png","SON")+"</td>"
  txt(196) = "<td>"+table_link_setup(OUTDIR,"nao.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(197) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"nao.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nao.tempreg.mam.png","MAM")+"</td>"
  txt(198) = "<td>"+table_link_setup(OUTDIR,"nao.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nao.tempreg.son.png","SON")+"</td>"
  txt(199) = "<td>"+table_link_setup(OUTDIR,"nao.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nao.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(200) = "<tr><td class=s2>PR Regressions</td><td class=s2>"+table_link_setup(OUTDIR,"nao.prreg.djf.png","DJF")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nao.prreg.mam.png","MAM")+"</td>"
  txt(201) = "<td class=s2>"+table_link_setup(OUTDIR,"nao.prreg.jja.png","JJA")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nao.prreg.son.png","SON")+"</td>"
  txt(202) = "<td class=s2>"+table_link_setup(OUTDIR,"nao.prreg.ann.png","Annual")+"</td><td class=s2>"+table_link_setup(OUTDIR,"nao.prreg.mon.png","Monthly")+"</td></tr>"

  txt(219) = ""
  txt(220) = "<tr><td rowspan=4 class=s2><b>PNA</b></td>"
  txt(221) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"pna.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pna.mam.png","MAM")+"</td>"
  txt(222) = "<td>"+table_link_setup(OUTDIR,"pna.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pna.son.png","SON")+"</td>"
  txt(223) = "<td>"+table_link_setup(OUTDIR,"pna.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"pna.mon.png","Monthly")+"</td></tr>"
  txt(224) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.mam.png","MAM")+"</td>"
  txt(225) = "<td>"+table_link_setup(OUTDIR,"pna.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.son.png","SON")+"</td>"
  txt(226) = "<td>"+table_link_setup(OUTDIR,"pna.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(227) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"pna.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pna.tempreg.mam.png","MAM")+"</td>"
  txt(228) = "<td>"+table_link_setup(OUTDIR,"pna.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pna.tempreg.son.png","SON")+"</td>"
  txt(229) = "<td>"+table_link_setup(OUTDIR,"pna.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"pna.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(230) = "<tr><td class=s2>PR Regressions</td><td class=s2>"+table_link_setup(OUTDIR,"pna.prreg.djf.png","DJF")+"</td><td class=s2>"+table_link_setup(OUTDIR,"pna.prreg.mam.png","MAM")+"</td>"
  txt(231) = "<td class=s2>"+table_link_setup(OUTDIR,"pna.prreg.jja.png","JJA")+"</td><td class=s2>"+table_link_setup(OUTDIR,"pna.prreg.son.png","SON")+"</td>"
  txt(232) = "<td class=s2>"+table_link_setup(OUTDIR,"pna.prreg.ann.png","Annual")+"</td><td class=s2>"+table_link_setup(OUTDIR,"pna.prreg.mon.png","Monthly")+"</td></tr>"

  txt(249) = ""
  txt(250) = "<tr><td rowspan=4 class=s2><b>NPO</b></td>"
  txt(251) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"npo.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"npo.mam.png","MAM")+"</td>"
  txt(252) = "<td>"+table_link_setup(OUTDIR,"npo.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"npo.son.png","SON")+"</td>"
  txt(253) = "<td>"+table_link_setup(OUTDIR,"npo.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"npo.mon.png","Monthly")+"</td></tr>"
  txt(254) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.mam.png","MAM")+"</td>"
  txt(255) = "<td>"+table_link_setup(OUTDIR,"npo.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.son.png","SON")+"</td>"
  txt(256) = "<td>"+table_link_setup(OUTDIR,"npo.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(257) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"npo.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"npo.tempreg.mam.png","MAM")+"</td>"
  txt(258) = "<td>"+table_link_setup(OUTDIR,"npo.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"npo.tempreg.son.png","SON")+"</td>"
  txt(259) = "<td>"+table_link_setup(OUTDIR,"npo.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"npo.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(260) = "<tr><td class=s2>PR Regressions</td><td class=s2>"+table_link_setup(OUTDIR,"npo.prreg.djf.png","DJF")+"</td><td class=s2>"+table_link_setup(OUTDIR,"npo.prreg.mam.png","MAM")+"</td>"
  txt(261) = "<td class=s2>"+table_link_setup(OUTDIR,"npo.prreg.jja.png","JJA")+"</td><td class=s2>"+table_link_setup(OUTDIR,"npo.prreg.son.png","SON")+"</td>"
  txt(262) = "<td class=s2>"+table_link_setup(OUTDIR,"npo.prreg.ann.png","Annual")+"</td><td class=s2>"+table_link_setup(OUTDIR,"npo.prreg.mon.png","Monthly")+"</td></tr>"

  txt(279) = ""
  txt(280) = "<tr><td rowspan=4 class=s2><b>PSA1</b></td>"
  txt(281) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"psa1.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa1.mam.png","MAM")+"</td>"
  txt(282) = "<td>"+table_link_setup(OUTDIR,"psa1.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa1.son.png","SON")+"</td>"
  txt(283) = "<td>"+table_link_setup(OUTDIR,"psa1.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa1.mon.png","Monthly")+"</td></tr>"
  txt(284) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.mam.png","MAM")+"</td>"
  txt(285) = "<td>"+table_link_setup(OUTDIR,"psa1.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.son.png","SON")+"</td>"
  txt(286) = "<td>"+table_link_setup(OUTDIR,"psa1.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(287) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"psa1.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa1.tempreg.mam.png","MAM")+"</td>"
  txt(288) = "<td>"+table_link_setup(OUTDIR,"psa1.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa1.tempreg.son.png","SON")+"</td>"
  txt(289) = "<td>"+table_link_setup(OUTDIR,"psa1.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa1.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(290) = "<tr><td class=s2>PR Regressions</td><td class=s2>"+table_link_setup(OUTDIR,"psa1.prreg.djf.png","DJF")+"</td><td class=s2>"+table_link_setup(OUTDIR,"psa1.prreg.mam.png","MAM")+"</td>"
  txt(291) = "<td class=s2>"+table_link_setup(OUTDIR,"psa1.prreg.jja.png","JJA")+"</td><td class=s2>"+table_link_setup(OUTDIR,"psa1.prreg.son.png","SON")+"</td>"
  txt(292) = "<td class=s2>"+table_link_setup(OUTDIR,"psa1.prreg.ann.png","Annual")+"</td><td class=s2>"+table_link_setup(OUTDIR,"psa1.prreg.mon.png","Monthly")+"</td></tr>"

  txt(309) = ""
  txt(310) = "<tr><td rowspan=4><b>PSA2</b></td>"
  txt(311) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"psa2.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.mam.png","MAM")+"</td>"
  txt(312) = "<td>"+table_link_setup(OUTDIR,"psa2.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.son.png","SON")+"</td>"
  txt(313) = "<td>"+table_link_setup(OUTDIR,"psa2.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.mon.png","Monthly")+"</td></tr>"
  txt(314) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.mam.png","MAM")+"</td>"
  txt(315) = "<td>"+table_link_setup(OUTDIR,"psa2.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.son.png","SON")+"</td>"
  txt(316) = "<td>"+table_link_setup(OUTDIR,"psa2.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(317) = "<tr><td>SST/TAS Regressions</td><td>"+table_link_setup(OUTDIR,"psa2.tempreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.tempreg.mam.png","MAM")+"</td>"
  txt(318) = "<td>"+table_link_setup(OUTDIR,"psa2.tempreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.tempreg.son.png","SON")+"</td>"
  txt(319) = "<td>"+table_link_setup(OUTDIR,"psa2.tempreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.tempreg.mon.png","Monthly")+"</td></tr>"
  txt(320) = "<tr><td>PR Regressions</td><td>"+table_link_setup(OUTDIR,"psa2.prreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.prreg.mam.png","MAM")+"</td>"
  txt(321) = "<td>"+table_link_setup(OUTDIR,"psa2.prreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.prreg.son.png","SON")+"</td>"
  txt(322) = "<td>"+table_link_setup(OUTDIR,"psa2.prreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.prreg.mon.png","Monthly")+"</td></tr>"
  txt(323) = "</tr></table>"

  txt(330) = ""
  txt(331) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(332) = "<CAPTION>Ice Extent Timeseries</CAPTION>"
  txt(333) = "<tr><td rowspan=2><b>SIC NH</b></td><td>"+table_link_setup(OUTDIR,"aice.extent.nh.djf.png","DJF")+"</td>"
  txt(334) = "<td>"+table_link_setup(OUTDIR,"aice.extent.nh.mam.png","MAM")+"</td>"
  txt(335) = "<td>"+table_link_setup(OUTDIR,"aice.extent.nh.jja.png","JJA")+"</td>"
  txt(336) = "<td>"+table_link_setup(OUTDIR,"aice.extent.nh.son.png","SON")+"</td>"
  txt(337) = "<td>"+table_link_setup(OUTDIR,"aice.extent.nh.ann.png","Annual")+"</td>"
  txt(338) = "<td>"+table_link_setup(OUTDIR,"aice.extent.nh.mon.png","Monthly")+"</td></tr>"
  txt(339) = "<tr><td colspan=2>"+table_link_setup(OUTDIR,"aice.extent.nh.feb.png","Feb")+"&nbsp;&nbsp;"+table_link_setup(OUTDIR,"aice.extent.nh.mar.png","Mar")+"&nbsp;&nbsp;"+table_link_setup(OUTDIR,"aice.extent.nh.sep.png","Sep")+"</td>"
  txt(341) = "<td colspan=2>"+table_link_setup(OUTDIR,"aice.extent.anom.nh.mon.png","Monthly Anomalies")+"</td>"
  txt(342) = "<td colspan=2>"+table_link_setup(OUTDIR,"aice.extent.nh.climo.png","Climatology")+"</td></tr>"
  txt(343) = "<tr><td rowspan=2><b>SIC SH</b></td><td>"+table_link_setup(OUTDIR,"aice.extent.sh.djf.png","DJF")+"</td>"
  txt(344) = "<td>"+table_link_setup(OUTDIR,"aice.extent.sh.mam.png","MAM")+"</td>"
  txt(345) = "<td>"+table_link_setup(OUTDIR,"aice.extent.sh.jja.png","JJA")+"</td>"
  txt(346) = "<td>"+table_link_setup(OUTDIR,"aice.extent.sh.son.png","SON")+"</td>"
  txt(347) = "<td>"+table_link_setup(OUTDIR,"aice.extent.sh.ann.png","Annual")+"</td>"
  txt(348) = "<td>"+table_link_setup(OUTDIR,"aice.extent.sh.mon.png","Monthly")+"</td></tr>"
  txt(349) = "<tr><td colspan=2>"+table_link_setup(OUTDIR,"aice.extent.sh.feb.png","Feb")+"&nbsp;&nbsp;"+table_link_setup(OUTDIR,"aice.extent.sh.mar.png","Mar")+"&nbsp;&nbsp;"+table_link_setup(OUTDIR,"aice.extent.sh.sep.png","Sep")+"</td>"
  txt(351) = "<td colspan=2>"+table_link_setup(OUTDIR,"aice.extent.anom.sh.mon.png","Monthly Anomalies")+"</td>"
  txt(352) = "<td colspan=2>"+table_link_setup(OUTDIR,"aice.extent.sh.climo.png","Climatology")+"</td></tr></table>"


  txt(361) = ""
  txt(362) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(363) = "<CAPTION>Global Trend Maps</CAPTION>"
  txt(364) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.trends.djf.png","DJF")+"</td>"
  txt(365) = "<td>"+table_link_setup(OUTDIR,"sst.trends.mam.png","MAM")+"</td>"
  txt(366) = "<td>"+table_link_setup(OUTDIR,"sst.trends.jja.png","JJA")+"</td>"
  txt(367) = "<td>"+table_link_setup(OUTDIR,"sst.trends.son.png","SON")+"</td>"
  txt(368) = "<td>"+table_link_setup(OUTDIR,"sst.trends.ann.png","Annual")+"</td>"
  txt(369) = "<td>"+table_link_setup(OUTDIR,"sst.trends.mon.png","Monthly")+"</td></tr>"
  txt(370) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.trends.djf.png","DJF")+"</td>"
  txt(371) = "<td>"+table_link_setup(OUTDIR,"tas.trends.mam.png","MAM")+"</td>"
  txt(372) = "<td>"+table_link_setup(OUTDIR,"tas.trends.jja.png","JJA")+"</td>"
  txt(373) = "<td>"+table_link_setup(OUTDIR,"tas.trends.son.png","SON")+"</td>"
  txt(374) = "<td>"+table_link_setup(OUTDIR,"tas.trends.ann.png","Annual")+"</td>"
  txt(375) = "<td>"+table_link_setup(OUTDIR,"tas.trends.mon.png","Monthly")+"</td></tr>"
  txt(376) = "<tr><td><b>PSL</b></td><td>"+table_link_setup(OUTDIR,"psl.trends.djf.png","DJF")+"</td>"
  txt(377) = "<td>"+table_link_setup(OUTDIR,"psl.trends.mam.png","MAM")+"</td>"
  txt(378) = "<td>"+table_link_setup(OUTDIR,"psl.trends.jja.png","JJA")+"</td>"
  txt(379) = "<td>"+table_link_setup(OUTDIR,"psl.trends.son.png","SON")+"</td>"
  txt(380) = "<td>"+table_link_setup(OUTDIR,"psl.trends.ann.png","Annual")+"</td>"
  txt(381) = "<td>"+table_link_setup(OUTDIR,"psl.trends.mon.png","Monthly")+"</td></tr>"
  txt(382) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.trends.djf.png","DJF")+"</td>"
  txt(383) = "<td>"+table_link_setup(OUTDIR,"pr.trends.mam.png","MAM")+"</td>"
  txt(384) = "<td>"+table_link_setup(OUTDIR,"pr.trends.jja.png","JJA")+"</td>"
  txt(385) = "<td>"+table_link_setup(OUTDIR,"pr.trends.son.png","SON")+"</td>"
  txt(386) = "<td>"+table_link_setup(OUTDIR,"pr.trends.ann.png","Annual")+"</td>"
  txt(387) = "<td>"+table_link_setup(OUTDIR,"pr.trends.mon.png","Monthly")+"</td></tr>" 
  txt(388) = "<tr><td><b>SND</b></td><td>"+table_link_setup(OUTDIR,"snd.trends.djf.png","DJF")+"</td>"
  txt(389) = "<td>"+table_link_setup(OUTDIR,"snd.trends.mam.png","MAM")+"</td>"
  txt(390) = "<td>"+table_link_setup(OUTDIR,"snd.trends.jja.png","JJA")+"</td>"
  txt(391) = "<td>"+table_link_setup(OUTDIR,"snd.trends.son.png","SON")+"</td>"
  txt(392) = "<td>"+table_link_setup(OUTDIR,"snd.trends.ann.png","Annual")+"</td>"
  txt(393) = "<td>"+table_link_setup(OUTDIR,"snd.trends.mon.png","Monthly")+"</td></tr>"
  txt(394) = "<tr><td><b>SIC NH</b></td><td>"+table_link_setup(OUTDIR,"aice.trends.nh.djf.png","DJF")+"</td>"
  txt(395) = "<td>"+table_link_setup(OUTDIR,"aice.trends.nh.mam.png","MAM")+"</td>"
  txt(396) = "<td>"+table_link_setup(OUTDIR,"aice.trends.nh.jja.png","JJA")+"</td>"
  txt(397) = "<td>"+table_link_setup(OUTDIR,"aice.trends.nh.son.png","SON")+"</td>"
  txt(398) = "<td>"+table_link_setup(OUTDIR,"aice.trends.nh.ann.png","Annual")+"</td>"
  txt(399) = "<td>"+table_link_setup(OUTDIR,"aice.trends.nh.mon.png","Monthly")+"</td></tr>" 
  txt(400) = "<tr><td><b>SIC SH</b></td><td>"+table_link_setup(OUTDIR,"aice.trends.sh.djf.png","DJF")+"</td>"
  txt(401) = "<td>"+table_link_setup(OUTDIR,"aice.trends.sh.mam.png","MAM")+"</td>"
  txt(402) = "<td>"+table_link_setup(OUTDIR,"aice.trends.sh.jja.png","JJA")+"</td>"
  txt(403) = "<td>"+table_link_setup(OUTDIR,"aice.trends.sh.son.png","SON")+"</td>"
  txt(404) = "<td>"+table_link_setup(OUTDIR,"aice.trends.sh.ann.png","Annual")+"</td>"
  txt(405) = "<td>"+table_link_setup(OUTDIR,"aice.trends.sh.mon.png","Monthly")+"</td></tr></table>"

  txt(429) = ""
  txt(430) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(431) = "<CAPTION>Globally-Averaged Timeseries</CAPTION>"
  txt(432) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.timeseries.djf.png","DJF")+"</td>"
  txt(433) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.mam.png","MAM")+"</td>"
  txt(434) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.jja.png","JJA")+"</td>"
  txt(435) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.son.png","SON")+"</td>"
  txt(436) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.ann.png","Annual")+"</td>"
  txt(437) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(438) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.timeseries.djf.png","DJF")+"</td>"
  txt(439) = "<td>"+table_link_setup(OUTDIR,"tas.timeseries.mam.png","MAM")+"</td>"
  txt(440) = "<td>"+table_link_setup(OUTDIR,"tas.timeseries.jja.png","JJA")+"</td>"
  txt(441) = "<td>"+table_link_setup(OUTDIR,"tas.timeseries.son.png","SON")+"</td>"
  txt(442) = "<td>"+table_link_setup(OUTDIR,"tas.timeseries.ann.png","Annual")+"</td>"
  txt(443) = "<td>"+table_link_setup(OUTDIR,"tas.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(444) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.timeseries.djf.png","DJF")+"</td>"
  txt(445) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.mam.png","MAM")+"</td>"
  txt(446) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.jja.png","JJA")+"</td>"
  txt(447) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.son.png","SON")+"</td>"
  txt(448) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.ann.png","Annual")+"</td>"
  txt(449) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(450) = "</table>"

  txt(574) = ""
  txt(575) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(576) = "<CAPTION>Running Trends of Globally-Averaged Monthly Timeseries</CAPTION>"
  txt(577) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.8yr_runtrend.mon.png","8yr")+"</td>"
  txt(578) = "<td>"+table_link_setup(OUTDIR,"sst.10yr_runtrend.mon.png","10yr")+"</td>"
  txt(579) = "<td>"+table_link_setup(OUTDIR,"sst.12yr_runtrend.mon.png","12yr")+"</td>"
  txt(580) = "<td>"+table_link_setup(OUTDIR,"sst.14yr_runtrend.mon.png","14yr")+"</td>"
  txt(581) = "<td>"+table_link_setup(OUTDIR,"sst.16yr_runtrend.mon.png","16yr")+"</td></tr>"
  txt(582) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.8yr_runtrend.mon.png","8yr")+"</td>"
  txt(583) = "<td>"+table_link_setup(OUTDIR,"tas.10yr_runtrend.mon.png","10yr")+"</td>"
  txt(584) = "<td>"+table_link_setup(OUTDIR,"tas.12yr_runtrend.mon.png","12yr")+"</td>"
  txt(585) = "<td>"+table_link_setup(OUTDIR,"tas.14yr_runtrend.mon.png","14yr")+"</td>"
  txt(586) = "<td>"+table_link_setup(OUTDIR,"tas.16yr_runtrend.mon.png","16yr")+"</td></tr>"
  txt(587) = "</table>"
  
  txt(609) = ""
  txt(610) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(611) = "<CAPTION>Additional Indices</CAPTION>"
  txt(612) = "<tr><td>"+table_link_setup(OUTDIR,"amm.timeseries.png","Atlantic Meridional Mode")+"</td>"
  txt(613) = "<td>"+table_link_setup(OUTDIR,"atl3.timeseries.png","Atlantic Ni&ntilde;o (ATL3)")+"</td>"
  txt(614) = "<td>"+table_link_setup(OUTDIR,"tna.timeseries.png","Tropical North Atlantic SST")+"</td>"
  txt(615) = "<td>"+table_link_setup(OUTDIR,"tsa.timeseries.png","Tropical South Atlantic SST")+"</td></tr>"
  txt(616) = "<tr><td>"+table_link_setup(OUTDIR,"nino12.timeseries.png","ni&ntilde;o1+2")+"</td>"
  txt(617) = "<td>"+table_link_setup(OUTDIR,"nino3.timeseries.png","ni&ntilde;o3")+"</td>"
  txt(618) = "<td>"+table_link_setup(OUTDIR,"nino4.timeseries.png","ni&ntilde;o4")+"</td>"
  txt(619) = "<td>"+table_link_setup(OUTDIR,"npi.timeseries.ndjfm.png","North Pacific PSL")+"</td></tr>"
  txt(620) = "<tr><td>"+table_link_setup(OUTDIR,"iod.timeseries.png","Indian Ocean SST Dipole")+"</td>"
  txt(621) = "<td>"+table_link_setup(OUTDIR,"tio.timeseries.png","Tropical Indian Ocean SST")+"</td>"
  txt(622) = "<td>"+table_link_setup(OUTDIR,"socn.timeseries.png","Southern Ocean SST")+"</td></tr>"
  txt(623) = "</table>"
  
;  if (OT.eq."True") then
;     txt(265) = "<p class=btext>Output timeseries written to "+OUTDIR+"</p>"
;  else
     txt(635) = "<br>"
;  end if  
  txt(636) = "</body></html>"
  
  tt = ind(.not.ismissing(txt))
  txt2 = txt(:tt(dimsizes(tt)-1)+1)
  tt2 = ind(.not.ismissing(txt2))
  txt3 = txt2(tt2)
  asciiwrite(OUTDIR+"index.html",txt3)
  delete([/txt,tt,tt2,txt2,txt3,quote,names,na,nsim_noobs,OT/])
;----------------------------------------------------------------------------
;--   Create calculation description webpage
;----------------------------------------------------------------------------
  txt = new(200,"string")
  quote = str_get_dq()
  
  txt(0) = "<HTML><HEAD><TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
  txt(1) = "body {background-color: WhiteSmoke; font-family: Verdana, Geneva, sans-serif}"
  txt(2) = "img.title {width:990; height:105; align:left}"
  txt(3) = "p.title {font-size:20pt; color:black; font-weight:bold; text-align:center}"
  txt(4) = "p.subtitle1 {font-size:12pt; color:gray}"
  txt(5) = "p.subtitle {font-size:16pt; color:black; font-weight:bold; text-align:left}"
  txt(6) = "p.btext {font-size:9pt}"
  txt(7) = "td {text-align:left; font-size:12pt}"
  txt(8) = "td.s1 {text-align:left; font-size:11pt}"
  txt(9) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
  txt(10) = "</STYLE></HEAD><BODY>"
  txt(11) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"  
  txt(12) = "<TABLE width='990' border=0 cellpadding=6><tr><td class=s1 width='300'><a href="+quote+"index.html"+quote+">Back to Diagnostics Plots</a></td>"
  txt(13) = "<td><p class=title>Methodology and Definitions</p></td></tr></table><hr width=990 align=left>" 

  txt(15) = "<p class=subtitle>General Notes</p>"
  txt(16) = "<TABLE width='800' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(17) = "<tr><td><ul><li>TS is surface ("+quote+"skin"+quote+") temperature and is used in lieu of sea surface temperatures (SSTs), TAS is 2m air temperature (CESM equivalent = TREFHT)."
  txt(18) = "<li>PR is total precipitation (CESM equivalent = PRECC+PRECL), PSL is sea level pressure, SND is snow depth (SNOWDP), and SIC is ice concentration (aice)."
  txt(19) = "<li>The annual cycle is removed prior to every calculation by subtracting the long-term monthly means. Exception: The annual cycle is not removed for mean spatial maps."
  txt(20) = "<li>Area-averages are always based on cosine of latitude weighting."
  txt(21) = "<li>The following calculations use linearly detrended data: standard deviations, nino3.4 spectra, and ENSO spatial composites."
  txt(22) = "<li>For visual clarity, the Y-axis <i>may</i> differ amongst individual panels on a particular plot."
  txt(23) = "<li>Climatological Zonal Averages: Climatological means are zonally averaged over the globe." 
  txt(24) = "<li>Power Spectra: The best-fit first-order Markov red noise spectrum (red curve) and its 95% (blue curve) and 99% (green curve) confidence bounds are shown on each panel. Top X-axis shows the period (in years), and the bottom X-axis shows the frequency (cycles/mo). If calculated, the observational spectrum is overlaid in gray on each model spectrum. The spectra are displayed in variance-conserving form." 
  txt(25) = "<li>Wavelet Analysis:   “A wavelet transform is computed using a Morlet wavelet with a wavenumber of 6. Areas significant at the 95% based on a chi-square test are stippled and the “cone of influence” is hatched. See Torrence, C. and G. P. Compo, 1998: A Practical Guide to Wavelet Analysis. Bull. Amer. Meteor. Soc., 79, 61-78. doi: http://dx.doi.org/10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2."
  txt(26) = "<li>Running Trends: N-year running trends are computed by calculating the linear trend over the N-year interval beginning at each successive timestep.  For instance, for a global timeseries that runs from 1970-2012, the 8yr running trend value for January 1970 is the linear trend during January 1970 - December 1977, and the value for January 2005 is the linear trend during January 2005 - December 2012."
  txt(27) = "<li>Running Standard Deviations: 30-year running standard deviations are computed by calculating the standard deviation centered over the 30 year interval at each successive timestep, with a minimum of 35 years of data required. For example, for a timeseries that runs from 1920-2010, the 30yr running standard deviation value for January 1935 is calculated over the period January 1920-December 1949. The right subtitle shows the minimum / mean / maximum of the computed timeseries." 
  txt(28) = "<li>Metrics Tables: Area-weighted pattern correlations and rms differences are calculated between observations and each model simulation (regridded to match the observational grid) for 11 climate metrics. The Total Score column shows the average of the 11 pattern correlations (Z-transformed) and rms differences. The following domains are used to compute the pattern correlations and rms differences: Means, standard deviations, ENSO, AMO, and PDO: Global; global for means, standard deviations, ENSO, AMO, and PDO; entire longitude/temporal range shown for El Nino and La Nina Hovm&ouml;llers; entire domain shown for NAM (20:90&deg;N) and SAM (20:90&deg;S)."
  txt(29) = "<li>EOF significance: If an eigenvalue is significantly seperated from neighboring values a star is appended to the percent variance explained on the plot. Significance is calculated following North et al. (MWR, 1982)."
  txt(30) = "<li>Ice Extent: Any grid cell defined as having a value at or above 15% is assumed to be 100% ice covered. The area of these grid cells is summed to create ice extent. For data sets with a northern hemisphere pole hole the CVDP looks for an attribute named pole_hole_area that is attached to the AREA variable in the netCDF file. The format of pole_hole_area should be start month (YYYYMM), end month (YYYYMM), area value, (repeat as necessary). If pole_hole_area is detected the CVDP will add the area value to the calculated northern hemisphere ice extent timeseries from the specified start month to the specified end month, and a * is placed after the dataset name in the ice extent plots."    
  txt(31) = "<li>For more information on observational datasets and climate indices, see the <a href="+quote+"http://climatedataguide.ucar.edu/"+quote+" target="+quote+"_blank"+quote+">Climate Data Guide</a>.</ul>"
  txt(32) = "</td></tr></table>"

  txt(40) = "<p class=subtitle>Modes of Variability</p>"
  txt(41) = "<TABLE width='800' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(42) = "<tr><td><b>ENSO Spatial Composites</b></td>"
  txt(43) = "<td>The normalized December nino3.4 timeseries is used to composite all years greater than 1 standard deviation (El Ni&ntilde;o) and all years less that -1 standard deviation (La Ni&ntilde;a). The number of El Ni&ntilde;o/La Ni&ntilde;a events composited is shown in the right subtitle. The December nino3.4 timeseries is based on the December values of the monthly nino3.4 time series smoothed with a 3-point binomial filter. TS/TAS/PSL composites: Temperatures are color shaded and in units of Celsius. Sea level pressure is contoured from -16 to 16hPa by 2hPa; negative contours are dashed. PR composites: Precipitation is color shaded and is in units of mm/day. See Deser, C., A. S. Phillips, R. A. Tomas, Y. Okumura, M. A. Alexander, A. Capotondi, J. D. Scott, Y. -O. Kwon, and M. Ohba, 2012: ENSO and Pacific Decadal Variability in Community Climate System Model Version 4. <i>J. Climate</i>, 25, 2622-2651, doi: 10.1175/JCLI-D-11-00301.1.</td></tr>"   
  txt(44) = "<tr><td><b>ENSO Hovm&ouml;llers</b></td>"
  txt(45) = "<td>A 1-2-1 running average is applied to the monthly nino3.4 timeseries, and then December values are selected and normalized. Meridional averages are calculated by averaging from 3&deg;S:3&deg;N, and spatial patterns are formed by compositing -1yr->+2yrs where the nino3.4 (1-2-1 weighted) December index is greater than 1 (El Ni&ntilde;o) and those years where the index is less than -1 (La Ni&ntilde;a). See Deser, C., A. S. Phillips, R. A. Tomas, Y. Okumura, M. A. Alexander, A. Capotondi, J. D. Scott, Y. -O. Kwon, and M. Ohba, 2012: ENSO and Pacific Decadal Variability in Community Climate System Model Version 4. J. Climate, 25, 2622-2651, doi: 10.1175/JCLI-D-11-00301.1.</td></tr>"  

  txt(46) = "<tr><td><b>PDO</b> (Pacific Decadal Oscillation)</td>"  
  txt(47) = "<td>Monthly index timeseries defined as the leading principal component (PC) of North Pacific (20:70&deg;N, 110&deg;E:100&deg;W) area-weighted SST* anomalies, where SST* denotes that the global mean SST anomaly has been removed at each timestep. Pattern created by regressing SST anomalies (in Celsius) at each grid box onto the normalized PC timeseries. Low pass-filtered timeseries (black curve) is based on a a 61-month running mean. See Deser, C., M. A. Alexander, S. -P. Xie, and A. S. Phillips, 2010: Sea surface temperature variability: patterns and mechanisms. <i>Ann. Rev. Mar. Sci.</i>, 2010.2, 115-143, doi:10.1146/annurev-marine-120408-151453. Also see Mantua, N. J., S. R. Hare, Y. Zhang, J. M. Wallace, and R. Francis, 1997: A Pacific interdecadal climate oscillation with impacts on salmon production. <i>Bull. Amer. Met. Soc.</i>, 1069-1079. For more information on the PDO see the <a href="+quote+"https://climatedataguide.ucar.edu/climate-data/pacific-decadal-oscillation-pdo-definition-and-indices?qt-climatedatasetmaintabs=1"+quote+" target="+quote+"_blank"+quote+">Climate Data Guide</a>.</td></tr>"
  txt(48) = "<tr><td><b>IPO</b> (Interdecadal Pacific Oscillation)</td>"  
  txt(49) = "<td>Monthly index timeseries defined as the leading principal component (PC) of 13yr low pass filtered Pacific (40&deg;S:60&deg;N, 110&deg;E:70&deg;W) area-weighted SST* anomalies, where SST* denotes that the global mean SST anomaly has been removed at each timestep.  Pattern created by regressing SST anomalies (in Celsius) at each grid box onto the normalized PC timeseries. At least 40 years of data are required for the IPO to be calculated. See Meehl, G.A. and A. Hu, 2007: Megadroughts in the Indian Monsoon Region and Southwest North America and a Mechanism for Associated Multidecadal Pacific Sea Surface Temperature Anomalies, <i>J. Clim</i>, 19, 1605-1623, doi: 10.1175/JCLI3675.1.</td></tr>"

  txt(50) = "<tr><td><b>AMO</b> (Atlantic Multidecadal Oscillation)</td>"  
  txt(51) = "<td>Monthly index timeseries defined as area-weighted SST* anomalies averaged over the North Atlantic (0:60&deg;N, 80&deg;W:0&deg;E), where SST* denotes that the global (60&deg;S:60&deg;N) mean SST anomaly has been removed at each timestep. Pattern created by regressing SST* anomalies onto the index timeseries and smoothing with a 9-point spatial filter. Low pass-filtered timeseries (black curve superimposed on the monthly timeseries) is based on a a 61-month running mean. Based on Trenberth, K. E., and D. J. Shea, 2006: Atlantic hurricanes and natural variability in 2005, <i>Geophys. Res. Lett.</i>, 33, L12704, doi:10.1029/2006GL026894. Low-pass filtered regression maps ("+quote+"Regr LP"+quote+") use a 10-year running mean on both the index timeseries and the field being regressed. For more information on the AMO see the <a href="+quote+"https://climatedataguide.ucar.edu/climate-data/atlantic-multi-decadal-oscillation-amo?qt-climatedatasetmaintabs=1"+quote+" target="+quote+"_blank"+quote+">Climate Data Guide</a>.</td></tr>"
  txt(52) = "<tr><td><b>AMOC</b></td>"
  txt(53) = "<td>The Atlantic Meridional Overturning Circulation (AMOC) is defined as the oceanic meridional mass transport (Sv) in the Atlantic sector. To compute AMOC, we follow the methods of Danabasoglu et al. (2012). Here we use annual averages of the AMOC, weighted by the cosine of the latitude and vertical extent of each model layer. Areas in which AMOC variance is low (standard deviation < 1e-6 Sv) are set to missing values for clarity. The leading EOF and associated principal component (PC) timeseries are computed over the Atlantic basin from 33&deg;S to 90&deg;N. The AMOC patterns are created by regressing the AMOC anomalies (in Sv) onto the normalized PC timeseries. The SST/TAS patterns associated with AMOC variations are created by regressing TAS/SST anomalies (in Celsius) at each grid box over the globe onto the normalized AMOC PC timeseries. A 15-point low-pass Lanczos filter is applied to the AMOC PC (and AMO) timeseries prior to computing lead/lag correlations, with a minimum of 90 years of data required. The data are not detrended (unlike Danabasoglu et al., 2012). See Danabasoglu, G., S. G. Yeager, Y. -O. Kwon, J. J. Tribbia, A. S. Phillips, and J. W. Hurrell, 2012. Variability of the Atlantic Meridional Overturning Circulation in CCSM4. J. Climate, 25, 5153-5172, doi: 10.1175/JCLI-D-11-00463.1.<br><br>For CCSM4 and CESM1, the MOC variable is read in, the Eulerian Mean, Eddy-Induced and Submeso components are summed, and the Atlantic Ocean + Mediterranean Sea + Labrador Sea + GIN Sea + Arctic Ocean + Hudson Bay transport region is selected. For CCSM2 and CCSM3 the same transport region is selected but only the Eulerian Mean component is used as that is all that is available. For CMIP5 data the msftmyz variable is read in and the atlantic_arctic_ocean basin is used. For CMIP3 data, the stfmmc variable is read in and the atlantic_ocean geo_region is used.</td></tr>" 

  txt(55) = "<tr><td><b>NAM</b> (Northern Annular Mode)</td>"  
  txt(56) = "<td>Seasonal/annual PSL averages are formed, square root of the cosine of the latitude weighting is applied, and then the leading EOF and associated principal component (PC) timeseries are computed over 20:90&deg;N, 0:360&deg;E. Pattern created by regressing global PSL anomalies (in hPa) onto normalized PC timeseries. Based on Hurrell, J. W., and C. Deser, 2009: North Atlantic climate variability: The role of the North Atlantic Oscillation. <i>J. Mar. Syst.</i>, 78, 28-41, doi:10.1016/j.jmarsys.2008.11.026. Also see Thompson, D. W. J., and J. M. Wallace, 2000: Annular modes in the extratropical circulation. Part I: Month-to-month variability. <i>J. Climate</i>, 13, 1000-1016.</td></tr>"  
  txt(57) = "<tr><td><b>NAO</b> (North Atlantic Oscillation)</td>"  
  txt(58) = "<td>Seasonal/annual PSL averages are formed, square root of the cosine of latitude weighting is applied, and then the leading EOF and associated principal component (PC) timeseries are computed over 20:80&deg;N, 90&deg;W:40&deg;E. Pattern created by regressing global PSL anomalies (in hPa) onto normalized PC timeseries. Based on Hurrell, J. W. and C. Deser, 2009: North Atlantic climate variability: The role of the North Atlantic Oscillation. <i>J. Mar. Syst.</i>, 78, 28-41, doi:10.1016/j.jmarsys.2008.11.026. For more information on the NAO see the <a href="+quote+"https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-pc-based?qt-climatedatasetmaintabs=1"+quote+" target="+quote+"_blank"+quote+">Climate Data Guide</a>.</td></tr>"  
  txt(59) = "<tr><td><b>SAM/PSA1/PSA2</b> (Southern Annular Mode, Pacific South American Patterns 1/2)</td>"  
  txt(60) = "<td>Seasonal/annual PSL averages are formed, square root of the cosine of latitude weighting is applied, and then the 1st (SAM), 2nd (PSA1) and 3rd (PSA2) EOFs and associated principal component (PC) timeseries are computed over 20:90&deg;S, 0:360&deg;E. Patterns created by regressing global PSL anomalies (in hPa) onto normalized PC timeseries. SAM calculation based on Thompson, D. W. J. and J.M. Wallace, 2000: Annular modes in the extratropical circulation. Part I: Month-to-month variability. <i>J. Climate</i>, 13, 1000-1016.</td></tr>"  
  txt(61) = "<tr><td><b>PNA/NPO</b> (Pacific North American Pattern, North Pacific Oscillation)</td>" 
  txt(62) = "<td>Seasonal/annual PSL averages are formed, the square root of the cosine of the latitude weighting is applied, and then the 1st (PNA) and 2nd (NPO) EOFs and associated principal component (PC) timeseries are computed over 20:85&deg;N, 120&deg;E:120&deg;W. Patterns created by regressing global PSL anomalies (in hPa) onto normalized PC timeseries.</td></tr>" 
  txt(63) = "<tr><td><b>SST Regressions</b></td>"
  txt(64) = "<td>SST anomalies (in Celsius) at each grid box are regressed upon the normalized atmospheric mode timeseries.</td></tr>" 
  txt(65) = "<tr><td><b>TAS Regressions</b></td>"
  txt(66) = "<td>TAS anomalies (in Celsius) at each grid box are regressed upon the normalized atmospheric mode timeseries.</td></tr>" 
  txt(67) = "<tr><td><b>PR Regressions</b></td>"
  txt(68) = "<td>PR anomalies (in mm/day) at each grid box are regressed upon the normalized atmospheric mode timeseries.</td></tr>" 
  txt(69) = "</table>"
  
  txt(80) = "<p class=subtitle>Climate Indices</p>"
  txt(81) = "<TABLE width='800' border=1 frame=void rules=all cellpadding=6 bordercolor=darkgray>"
  txt(82) = "<tr><td><b>Atlantic Meridional Mode</b></td>"
  txt(83) = "<td>Defined as the difference between area-averaged SST anomalies computed over 5:15&deg;N, 20:50&deg;W and area-averaged SST anomalies computed over 5:15&deg;S, 20&degW:10&deg;E. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See Doi, T., T. Tozuka and T. Yamagata (2009), Interannual variability of the Guinea Dome and its possible link with the Atlantic Meridional Mode. Climate Dynamics, 33, 985-998, doi:10.1007/s00382-009-0574-z.</td></tr>"
  txt(84) = "<tr><td><b>Atlantic Ni&ntilde;o (ATL3)</b></td>"
  txt(85) = "<td>Area-averaged SST anomalies computed over 3&deg;S:3&deg;N, 20&deg;W:0&deg;E. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See Zebiak, S. E., (1993): Air–sea interaction in the equatorial Atlantic region. Journal of Climate, 6, 1567–1586.</td></tr>"
  txt(86) = "<tr><td><b>Tropical North Atlantic SST</b></td>"
  txt(87) = "<td>Area-averaged SST anomalies computed over 5.5:23.5&deg;N, 15:57.5&deg;W. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See Enfield, D.B., A.M. Mestas, D.A. Mayer, and L. Cid-Serrano (1999), How ubiquitous is the dipole relationship in tropical Atlantic sea surface temperatures?, JGR-O, 104, 7841-7848.</td></tr>"
  txt(88) = "<tr><td><b>Tropical South Atlantic SST</b></td>"
  txt(89) = "<td>Area-averaged SST anomalies computed over 0:20&deg;S, 30&deg;W:10&deg;E. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See Enfield, D.B., A.M. Mestas, D.A. Mayer, and L. Cid-Serrano (1999), How ubiquitous is the dipole relationship in tropical Atlantic sea surface temperatures?, JGR-O, 104, 7841-7848.</td></tr>"
  txt(90) = "<tr><td><b>ni&ntilde;o1+2</b></td>"
  txt(91) = "<td>Area-averaged SST anomalies computed over 0:10&deg;S, 80:90&deg;W. Red/blue shading denotes positive/negative departures from the best-fit linear trend line.</td></tr>"
  txt(92) = "<tr><td><b>ni&ntilde;o3</b></td>"
  txt(93) = "<td>Area-averaged SST anomalies computed over 5&deg;S:5&deg;N, 90:150&deg;W. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See Trenberth, K. E. (1997) The Definition of El Niño. Bulletin of the American Meteorological Society, 78, 2771-2777.</td></tr>"
  txt(94) = "<tr><td><b>ni&ntilde;o3.4</b></td>"
  txt(95) = "<td>Area-averaged SST anomalies computed over 5&deg;S:5&deg;N, 120:170&deg;W. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See Trenberth, K. E. (1997) The Definition of El Niño. Bulletin of the American Meteorological Society, 78, 2771-2777.</td></tr>"
  txt(96) = "<tr><td><b>ni&ntilde;o4</b></td>"
  txt(97) = "<td>Area-averaged SST anomalies computed over 5&deg;S:5&deg;N, 160&deg;E:150&deg;W. Red/blue shading denotes positive/negative departures from the best-fit linear trend line.</td></tr>"
  txt(98) = "<tr><td><b>NPI</b> (North Pacific PSL Index)</td>"
  txt(99) = "<td>Winter (December-March) average PSL anomalies (in hPa) area-averaged over 30&deg;:65&deg;N, 160&deg;E:140&deg;W. Based on Trenberth, K. E. and J. W. Hurrell, 1994: Decadal atmosphere-ocean variations in the Pacific, <i>Climate Dynamics</i>, 9, 303-319.</td></tr>"    

  txt(100) = "<tr><td><b>Indian Ocean SST Dipole</b></td>"
  txt(101) = "<td>Defined as the difference between area-averaged SST anomalies computed over 10&deg;S:10&deg;N, 50:70&deg;E and area-averaged SST anomalies computed over 0:10&deg;S, 90:110&deg;E. Red/blue shading denotes positive/negative departures from the best-fit linear trend line. See: Saji N.H., Goswami B.N., Vinayachandran P.N., Yamagata T., 1999: A dipole mode in the tropical Indian Ocean, Nature, 401, 360-363.</td></tr>"
  txt(102) = "<tr><td><b>Tropical Indian Ocean SST</b></td>"
  txt(103) = "<td>Area-averaged SST anomalies computed over 15&deg;S:15&deg;N, 40:110&deg;E. Red/blue shading denotes positive/negative departures from the best-fit linear trend line.</td></tr>"
  txt(104)= "<tr><td><b>Southern Ocean SST</b></td>"
  txt(105)= "<td>Area-averaged SST anomalies computed over 50&deg;:70&deg;S, 0:360&deg;E. Red/blue shading denotes positive/negative departures from the best-fit linear trend line.</td></tr>"
  txt(106) = "</table>"
  txt(107) = ""
  
  txt(112) = "<br>Created "+systemfunc("date")
  txt(113) = "<p class=btext>CVDP Version "+VERSION+"</p>"
  txt(114) = "</body></html>"
  
  tt = ind(.not.ismissing(txt))
  txt2 = txt(:tt(dimsizes(tt)-1)+1)
  txt2 = where(ismissing(txt2),"",txt2)
  asciiwrite(OUTDIR+"methodology.html",txt2)
  delete([/tt,txt,txt2/])
;----------------------------------------------------------------------------
;--   Create metrics webpage
;----------------------------------------------------------------------------
  if (isfilepresent2(OUTDIR+"metrics.txt")) then
     txt = new(500,"string")  
     quote = str_get_dq()
  
     txt(0) = "<HTML><HEAD><TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
     txt(1) = "body {background-color: WhiteSmoke; font-family: Verdana, Geneva, sans-serif}"
     txt(2) = "img.title {width:990; height:105; align:left}"
     txt(3) = "p.subtitle1 {font-size:15pt; color:gray}"
     txt(4) = "p.subtitle {font-size:20pt; color:black; text-align:center; font-weight:bold}"
     txt(5) = "p.btext {font-size:12pt; text-align:left}"
     txt(6) = "p.btext2 {font-size:9pt; text-align:left}"
     txt(7) = "td {text-align:center; font-size:12pt}"
     txt(8) = "td.s1 {text-align:left; font-size:11pt}"
     txt(9) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
     txt(10) = "</STYLE></HEAD><BODY>"
     txt(11) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"
     txt(12)= "<TABLE width='990' border=0 cellpadding=6>"
     txt(13) = "<tr><td width='300'><p class=btext><a href="+quote+"index.html"+quote+">Back to Diagnostics Plots</a></p></td>"
     txt(14) = "<td><p class=subtitle>"+webtitle+"</p></td></tr></table><hr width=990 align=left><p class=subtitle1>Metrics Table</p>" 
     
     z = asciiread(OUTDIR+"metrics.txt",(/-1/),"string")
     nlines = dimsizes(z)
     txt(24) = "<PRE>"
     do gg = 0,nlines-1
        txt(25+gg) = z(gg)
     end do
     txt(25+nlines) = "</PRE>"
     txt(25+nlines+1) = "<p class=btext><a href="+quote+"namelist_obs"+quote+">Observations Used</a></p>"
     txt(25+nlines+2) = "Created "+systemfunc("date")
     txt(25+nlines+3) = "<p class=btext2>CVDP Version "+VERSION+"</p></body></html>"
  
     tt = ind(.not.ismissing(txt))
     txt2 = txt(:tt(dimsizes(tt)-1)+1)
     txt2 = where(ismissing(txt2),"",txt2)
     asciiwrite(OUTDIR+"metrics.html",txt2)
     delete([/txt,txt2,tt/])
  end if 
  if (isfilepresent2(OUTDIR+"metrics.table_1.gif")) then  ; create sorted metrics table pages
     txt = new(500,"string")  
     quote = str_get_dq()
  
     txt(0) = "<HTML><HEAD><TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
     txt(1) = "body {background-color: WhiteSmoke; font-family: Verdana, Geneva, sans-serif}"
     txt(2) = "img.title {width:990; height:105; align:left}"
     txt(3) = "p.subtitle1 {font-size:20pt; color:darkgray}"
     txt(4) = "p.subtitle {font-size:20pt; color:black; text-align:center; font-weight:bold}"
     txt(5) = "p.btext {font-size:12pt; text-align:left}"
     txt(6) = "p.btext2 {font-size:9pt; text-align:left}"
     txt(7) = "td {text-align:center; font-size:12pt}"
     txt(8) = "td.s1 {text-align:left; font-size:11pt}"
     txt(9) = "td.tablehdr {text-align:center; font-size:17pt}"
     txt(10) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
     txt(11) = "</STYLE></HEAD><BODY>"
     txt(12) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"
     txt(13)= "<TABLE width='990' border=0 cellpadding=6>"


     sort_txt = (/"<a href="+quote+"metrics.table_0.html"+quote+">Namelist (default)</a> | ","<a href="+quote+"metrics.table_1.html"+quote+">Namelist (Alphabetically)</a> | ","<a href="+quote+"metrics.table_2.html"+quote+">ENSO TAS</a> | ","<a href="+quote+"metrics.table_3.html"+quote+">ENSO PSL</a> | ", \
                  "<a href="+quote+"metrics.table_4.html"+quote+">El Ni&ntilde;o Hovm&ouml;ller</a><br>","<a href="+quote+"metrics.table_5.html"+quote+">La Ni&ntilde;a Hovm&ouml;ller</a> | ","<a href="+quote+"metrics.table_6.html"+quote+">AMO</a> | ","<a href="+quote+"metrics.table_7.html"+quote+">PDO</a> | ", \
                  "<a href="+quote+"metrics.table_8.html"+quote+">NAM</a> | ","<a href="+quote+"metrics.table_9.html"+quote+">SAM</a> | ","<a href="+quote+"metrics.table_10.html"+quote+">SST std dev</a> | ","<a href="+quote+"metrics.table_11.html"+quote+">PSL std dev</a> | ", \
                  "<a href="+quote+"metrics.table_12.html"+quote+">PR std dev</a> | ","<a href="+quote+"metrics.table_13.html"+quote+">Mean Score</a>"/)
     sort_txt2 = (/"Namelist (default)","Namelist (Alphabetically)","ENSO TAS","ENSO PSL","El Ni&ntilde;o Hovm&ouml;ller","La Ni&ntilde;a Hovm&ouml;ller","AMO","PDO","NAM","SAM","SST std dev","PSL std dev","PR std dev","Mean Score"/)

     sort_txtA = (/"<a href="+quote+"metrics.table_14.html"+quote+">Namelist (default)</a> | ","<a href="+quote+"metrics.table_15.html"+quote+">Namelist (Alphabetically)</a> | ","<a href="+quote+"metrics.table_16.html"+quote+">ENSO TAS</a> | ","<a href="+quote+"metrics.table_17.html"+quote+">ENSO PSL</a> | ", \
                  "<a href="+quote+"metrics.table_18.html"+quote+">El Ni&ntilde;o Hovm&ouml;ller</a><br>","<a href="+quote+"metrics.table_19.html"+quote+">La Ni&ntilde;a Hovm&ouml;ller</a> | ","<a href="+quote+"metrics.table_20.html"+quote+">AMO</a> | ","<a href="+quote+"metrics.table_21.html"+quote+">PDO</a> | ", \
                  "<a href="+quote+"metrics.table_22.html"+quote+">NAM</a> | ","<a href="+quote+"metrics.table_23.html"+quote+">SAM</a> | ","<a href="+quote+"metrics.table_24.html"+quote+">SST std dev</a> | ","<a href="+quote+"metrics.table_25.html"+quote+">PSL std dev</a> | ", \
                  "<a href="+quote+"metrics.table_26.html"+quote+">PR std dev</a> | ","<a href="+quote+"metrics.table_27.html"+quote+">Mean Score</a>"/)

     strarr = new(14,string)
     strarr = " | "
     strarr(13) = " "
     do gg = 0,13
           txt(14) = "<tr><td width='50'><p class=btext>Go to:</p></td><td class=s1 width='300'><a href="+quote+"methodology.html"+quote+">Methodology and Definitions</a><br><a href="+quote+"index.html"+quote+">Diagnostics Plots</a><br><a href="+quote+"metrics.table_"+(gg+14)+".html"+quote+">RMS Metric Tables</a></td>"
           txt(15) = "<td><p class=subtitle>"+webtitle+"</p></td></tr></table><hr width=990 align=left>"
           sort_txtT = sort_txt
           sort_txtT(gg) = "<b>"+sort_txt2(gg)+"</b>"+strarr(gg)
           txt(24) = "<p class=subtitle1><b>Pattern Correlation Metrics Tables</b></p>"
           txt(25) = "<TABLE width='900' border=0 cellpadding=16><tr><td width=50><b>Sort By:</b></td><td>"+str_concat(sort_txtT)+"</td></tr></table>" 
           txt(26) = "<img width=980 src="+quote+"metrics.table_"+gg+".gif"+quote+" />"
           txt(27) = "<p class=btext><a href="+quote+"namelist_obs"+quote+">Observations Used</a></p>"
           txt(28) = "Created "+systemfunc("date")
           txt(29) = "<p class=btext2>CVDP Version "+VERSION+"</p></body></html>"
           tt = ind(.not.ismissing(txt))
           txt2 = txt(:tt(dimsizes(tt)-1)+1)
           txt2 = where(ismissing(txt2),"",txt2)
           asciiwrite(OUTDIR+"metrics.table_"+gg+".html",txt2)
           delete([/txt2,tt,sort_txtT/])
     end do

     
     do gg = 0,13
           txt(14) = "<tr><td width='50'><p class=btext>Go to:</p></td><td class=s1 width='300'><a href="+quote+"methodology.html"+quote+">Methodology and Definitions</a><br><a href="+quote+"index.html"+quote+">Diagnostics Plots</a><br><a href="+quote+"metrics.table_"+gg+".html"+quote+">Pattern Correlation Metric Tables</a></td>"
           txt(15) = "<td><p class=subtitle>"+webtitle+"</p></td></tr></table><hr width=990 align=left>"
           sort_txtT = sort_txtA
           sort_txtT(gg) = "<b>"+sort_txt2(gg)+"</b>"+strarr(gg)
           txt(24) = "<p class=subtitle1><b>RMS Metrics Tables</b></p>"
           txt(25) = "<TABLE width='900' border=0 cellpadding=16><tr><td width=50><b>Sort By:</b></td><td>"+str_concat(sort_txtT)+"</td></tr></table>"
           txt(26) = "<img width=980 src="+quote+"metrics.table_"+(gg+14)+".gif"+quote+" />"
           txt(27) = "<p class=btext><a href="+quote+"namelist_obs"+quote+">Observations Used</a></p>"
           txt(28) = "Created "+systemfunc("date")
           txt(29) = "<p class=btext2>CVDP Version "+VERSION+"</p></body></html>"
           tt = ind(.not.ismissing(txt))
           txt2 = txt(:tt(dimsizes(tt)-1)+1)
           txt2 = where(ismissing(txt2),"",txt2)
           asciiwrite(OUTDIR+"metrics.table_"+(gg+14)+".html",txt2)
           delete([/txt2,tt/])
     end do
     delete([/txt,strarr/])
  end if
  delete(OUTDIR)
  print("Finished: webpage.ncl")  
end

