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
  WEBTITLE    = getenv("WEBTITLE")
  OPT_CLIMO   = getenv("OPT_CLIMO")
  CLIMO_SYEAR = getenv("CLIMO_SYEAR")
  CLIMO_EYEAR = getenv("CLIMO_EYEAR")
  OBS         = getenv("OBS")

  WEBTITLE = str_sub_str(WEBTITLE,"%"," ")

  quote = str_get_dq()
  
  if (OPT_CLIMO.eq."Full") then
     subtxt = "Climatological Period Used: Full"
  end if
  if (OPT_CLIMO.eq."Common") then
     subtxt = "Climatological Period Used: "
  end if
  if (OPT_CLIMO.eq."Custom") then
     subtxt = "Climatological Period Used: "+CLIMO_SYEAR+"-"+CLIMO_EYEAR
  end if
  if (OBS.eq."True") then
     namelist_txt = "Input Namelists: <a href="+quote+"namelist_obs"+quote+">OBS</a> | <a href="+quote+"namelist"+quote+">Models</a><br>Derived Namelists: <a href="+quote+"namelist_moc"+quote+">MOC</a> | <a href="+quote+"namelist_prect"+quote+">PR</a> | <a href="+quote+"namelist_psl"+quote+">PSL</a><br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="+quote+"namelist_snowdp"+quote+">SND</a> | <a href="+quote+"namelist_trefht"+quote+">TAS</a> | <a href="+quote+"namelist_ts"+quote+">TS</a>"          
  else
     namelist_txt = "Input Namelist: <a href="+quote+"namelist"+quote+">Models</a><br>Derived Namelists: <a href="+quote+"namelist_moc"+quote+">MOC</a> | <a href="+quote+"namelist_prect"+quote+">PR</a> | <a href="+quote+"namelist_psl"+quote+">PSL</a><br>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="+quote+"namelist_snowdp"+quote+">SND</a> | <a href="+quote+"namelist_trefht"+quote+">TAS</a> | <a href="+quote+"namelist_ts"+quote+">TS</a>"
  end if
  nsim_noobs = numAsciiRow("namelist")
  na = asciiread("namelist",(/nsim_noobs/),"string")
  names = new(nsim_noobs,"string")
  delim = "|"
  do gg = 0,nsim_noobs-1
     names(gg) = str_strip(str_get_field(na(gg),1,delim))
  end do
;  print(names)
  
  txt = new(300,"string")
  quote = str_get_dq()
  
  txt(0) = "<HTML><HEAD>"
  txt(1) = "<TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
  txt(2) = "img.title {width:990; height:105; align:left}"
  txt(3) = "p.subtitle {font-size:20pt; color:black; text-align:center; font-weight:bold}"
  txt(4) = "p.btext {font-size:10pt; text-align:left}"
  txt(5) = "td {text-align:center; font-size:12pt}"
  txt(6) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
  txt(7) = "</STYLE></HEAD><BODY BGCOLOR="+quote+"WhiteSmoke"+quote+">"
  txt(8) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"
  txt(9) = "<TABLE width='990' border=0 cellpadding=6>"
  if (isfilepresent2(OUTDIR+"metrics.txt")) then
     txt(10) = "<tr><td width='300'><p class=btext><a href="+quote+"methodology.html"+quote+">Methodology</a> | <a href="+quote+"metrics.html"+quote+">Metrics Table</a><br>"+subtxt+"<br>"+namelist_txt+"<br>"
  else
     txt(10) = "<tr><td width='300'><p class=btext><a href="+quote+"methodology.html"+quote+">Methodology</a><br>"+subtxt+"<br>"+namelist_txt+"<br>"
  end if
  txt(11) = "Created: "+systemfunc("date")+"<br>CVDP Version "+VERSION+"</p></td>"
 
  txt(12) = "<td><p class=subtitle>"+WEBTITLE+"</p></td></tr></table><hr width=990 align=left>"
  delete([/subtxt,namelist_txt/])

  txt(13) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(14) = "<CAPTION>Means</CAPTION>"
  txt(15) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"sst.mean.ann.png","Annual")+"</td></tr>"    
  txt(16) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"tas.mean.ann.png","Annual")+"</td></tr>"    
  txt(17) = "<tr><td><b>PSL</b></td><td>"+table_link_setup(OUTDIR,"psl.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"psl.mean.ann.png","Annual")+"</td></tr>"    
  txt(18) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.mean.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"pr.mean.ann.png","Annual")+"</td></tr>"    
  txt(19) = "</table>"  
  
  txt(23) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(24) = "<CAPTION>Standard Deviations</CAPTION>"
  txt(25) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"sst.stddev.ann.png","Annual")+"</td></tr>"    
  txt(26) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"tas.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"tas.stddev.ann.png","Annual")+"</td></tr>"    
  txt(27) = "<tr><td><b>PSL</b></td><td>"+table_link_setup(OUTDIR,"psl.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"psl.stddev.ann.png","Annual")+"</td></tr>"    
  txt(28) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.stddev.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.mam.png","MAM")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.son.png","SON")+"</td><td>"+table_link_setup(OUTDIR,"pr.stddev.ann.png","Annual")+"</td></tr>"    
  txt(29) = "</table>"  
  
  txt(36) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(37) = "<CAPTION>Coupled Modes of Variability</CAPTION>"  
  txt(38) = "<tr><td><b>AMO</b></td>"  
  txt(39) = "<td>"+table_link_setup(OUTDIR,"amo.png","Pattern")+"</td>"
  txt(40) = "<td>"+table_link_setup(OUTDIR,"amo.timeseries.png","Timeseries")+"</td>"
  txt(41) = "<td>"+table_link_setup(OUTDIR,"amo.powspec.png","Power Spectra")+"</td></tr>"  
  txt(42) = "<tr><td><b>PDO</b></td>"  
  txt(43) = "<td>"+table_link_setup(OUTDIR,"pdo.png","Pattern")+"</td>"
  txt(44) = "<td>"+table_link_setup(OUTDIR,"pdo.timeseries.png","Timeseries")+"</td>"
  txt(45) = "<td>"+table_link_setup(OUTDIR,"pdo.powspec.png","Power Spectra")+"</td></tr>" 

  txt(46) = "<tr><td rowspan=5><b>ENSO</b></td>"
  txt(47) = "<td rowspan=3>Spatial Composites</td><td>"+table_link_setup(OUTDIR,"nino34.spatialcomp.jja0.png","JJA<sup>0</sup>")+"</td><td>"+table_link_setup(OUTDIR,"nino34.spatialcomp.son0.png","SON<sup>0</sup>")+"</td></tr>"
  txt(48) = "<tr><td>"+table_link_setup(OUTDIR,"nino34.spatialcomp.djf1.png","DJF<sup>+1</sup>")+"</td><td>"+table_link_setup(OUTDIR,"nino34.spatialcomp.mam1.png","MAM<sup>+1</sup>")+"</td></tr>"
  txt(49) = "<tr><td>"+table_link_setup(OUTDIR,"nino34.hov.elnino.png","El Ni&ntilde;o Hovm&ouml;ller")+"</td><td>"+table_link_setup(OUTDIR,"nino34.hov.lanina.png","La Ni&ntilde;a Hovm&ouml;ller")+"</td></tr>"
  txt(50) = "<tr><td rowspan=2>Ni&ntilde;o3.4</td><td>"+table_link_setup(OUTDIR,"nino34.timeseries.png","Timeseries")+"</td><td>"+table_link_setup(OUTDIR,"nino34.powspec.png","Power Spectra")+"</td></tr>"
  txt(51) = "<tr><td>"+table_link_setup(OUTDIR,"nino34.monstddev.png","Monthly  Std. Dev.")+"</td><td>"+table_link_setup(OUTDIR,"nino34.runstddev.png","Running Std. Dev.")+"</td></tr></table>"


  txt(58) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6><CAPTION>Atmospheric Modes of Variability</CAPTION>"  
  txt(60) = "<tr><td rowspan=3><b>NAM</b></td>"
  txt(61) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"nam.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nam.mam.png","MAM")+"</td>"
  txt(62) = "<td>"+table_link_setup(OUTDIR,"nam.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nam.son.png","SON")+"</td>"
  txt(63) = "<td>"+table_link_setup(OUTDIR,"nam.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nam.mon.png","Monthly")+"</td></tr>"
  txt(64) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.mam.png","MAM")+"</td>"
  txt(65) = "<td>"+table_link_setup(OUTDIR,"nam.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.son.png","SON")+"</td>"
  txt(66) = "<td>"+table_link_setup(OUTDIR,"nam.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nam.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(67) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"nam.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nam.sstreg.mam.png","MAM")+"</td>"
  txt(68) = "<td>"+table_link_setup(OUTDIR,"nam.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nam.sstreg.son.png","SON")+"</td>"
  txt(69) = "<td>"+table_link_setup(OUTDIR,"nam.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nam.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(70) = "<tr><td rowspan=3><b>SAM</b></td>"
  txt(71) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"sam.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sam.mam.png","MAM")+"</td>"
  txt(72) = "<td>"+table_link_setup(OUTDIR,"sam.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sam.son.png","SON")+"</td>"
  txt(73) = "<td>"+table_link_setup(OUTDIR,"sam.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"sam.mon.png","Monthly")+"</td></tr>"
  txt(74) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.mam.png","MAM")+"</td>"
  txt(75) = "<td>"+table_link_setup(OUTDIR,"sam.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.son.png","SON")+"</td>"
  txt(76) = "<td>"+table_link_setup(OUTDIR,"sam.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"sam.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(77) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"sam.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"sam.sstreg.mam.png","MAM")+"</td>"
  txt(78) = "<td>"+table_link_setup(OUTDIR,"sam.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"sam.sstreg.son.png","SON")+"</td>"
  txt(79) = "<td>"+table_link_setup(OUTDIR,"sam.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"sam.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(80) = "<tr><td rowspan=3><b>NAO</b></td>"
  txt(81) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"nao.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nao.mam.png","MAM")+"</td>"
  txt(82) = "<td>"+table_link_setup(OUTDIR,"nao.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nao.son.png","SON")+"</td>"
  txt(83) = "<td>"+table_link_setup(OUTDIR,"nao.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nao.mon.png","Monthly")+"</td></tr>"
  txt(84) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.mam.png","MAM")+"</td>"
  txt(85) = "<td>"+table_link_setup(OUTDIR,"nao.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.son.png","SON")+"</td>"
  txt(86) = "<td>"+table_link_setup(OUTDIR,"nao.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nao.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(87) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"nao.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"nao.sstreg.mam.png","MAM")+"</td>"
  txt(88) = "<td>"+table_link_setup(OUTDIR,"nao.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"nao.sstreg.son.png","SON")+"</td>"
  txt(89) = "<td>"+table_link_setup(OUTDIR,"nao.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"nao.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(90) = "<tr><td rowspan=3><b>PNA</b></td>"
  txt(91) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"pna.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pna.mam.png","MAM")+"</td>"
  txt(92) = "<td>"+table_link_setup(OUTDIR,"pna.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pna.son.png","SON")+"</td>"
  txt(93) = "<td>"+table_link_setup(OUTDIR,"pna.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"pna.mon.png","Monthly")+"</td></tr>"
  txt(94) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.mam.png","MAM")+"</td>"
  txt(95) = "<td>"+table_link_setup(OUTDIR,"pna.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.son.png","SON")+"</td>"
  txt(96) = "<td>"+table_link_setup(OUTDIR,"pna.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"pna.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(97) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"pna.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"pna.sstreg.mam.png","MAM")+"</td>"
  txt(98) = "<td>"+table_link_setup(OUTDIR,"pna.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"pna.sstreg.son.png","SON")+"</td>"
  txt(99) = "<td>"+table_link_setup(OUTDIR,"pna.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"pna.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(100) = "<tr><td rowspan=3><b>NPO</b></td>"
  txt(101) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"npo.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"npo.mam.png","MAM")+"</td>"
  txt(102) = "<td>"+table_link_setup(OUTDIR,"npo.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"npo.son.png","SON")+"</td>"
  txt(103) = "<td>"+table_link_setup(OUTDIR,"npo.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"npo.mon.png","Monthly")+"</td></tr>"
  txt(104) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.mam.png","MAM")+"</td>"
  txt(105) = "<td>"+table_link_setup(OUTDIR,"npo.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.son.png","SON")+"</td>"
  txt(106) = "<td>"+table_link_setup(OUTDIR,"npo.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"npo.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(107) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"npo.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"npo.sstreg.mam.png","MAM")+"</td>"
  txt(108) = "<td>"+table_link_setup(OUTDIR,"npo.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"npo.sstreg.son.png","SON")+"</td>"
  txt(109) = "<td>"+table_link_setup(OUTDIR,"npo.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"npo.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(110) = "<tr><td rowspan=3><b>PSA1</b></td>"
  txt(111) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"psa1.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa1.mam.png","MAM")+"</td>"
  txt(112) = "<td>"+table_link_setup(OUTDIR,"psa1.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa1.son.png","SON")+"</td>"
  txt(113) = "<td>"+table_link_setup(OUTDIR,"psa1.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa1.mon.png","Monthly")+"</td></tr>"
  txt(114) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.mam.png","MAM")+"</td>"
  txt(115) = "<td>"+table_link_setup(OUTDIR,"psa1.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.son.png","SON")+"</td>"
  txt(116) = "<td>"+table_link_setup(OUTDIR,"psa1.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa1.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(117) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"psa1.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa1.sstreg.mam.png","MAM")+"</td>"
  txt(118) = "<td>"+table_link_setup(OUTDIR,"psa1.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa1.sstreg.son.png","SON")+"</td>"
  txt(119) = "<td>"+table_link_setup(OUTDIR,"psa1.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa1.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(120) = "<tr><td rowspan=3><b>PSA2</b></td>"
  txt(121) = "<td>Patterns</td><td>"+table_link_setup(OUTDIR,"psa2.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.mam.png","MAM")+"</td>"
  txt(122) = "<td>"+table_link_setup(OUTDIR,"psa2.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.son.png","SON")+"</td>"
  txt(123) = "<td>"+table_link_setup(OUTDIR,"psa2.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.mon.png","Monthly")+"</td></tr>"
  txt(124) = "<tr><td>Timeseries</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.mam.png","MAM")+"</td>"
  txt(125) = "<td>"+table_link_setup(OUTDIR,"psa2.timeseries.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.son.png","SON")+"</td>"
  txt(126) = "<td>"+table_link_setup(OUTDIR,"psa2.timeseries.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(127) = "<tr><td>SST Regressions</td><td>"+table_link_setup(OUTDIR,"psa2.sstreg.djf.png","DJF")+"</td><td>"+table_link_setup(OUTDIR,"psa2.sstreg.mam.png","MAM")+"</td>"
  txt(128) = "<td>"+table_link_setup(OUTDIR,"psa2.sstreg.jja.png","JJA")+"</td><td>"+table_link_setup(OUTDIR,"psa2.sstreg.son.png","SON")+"</td>"
  txt(129) = "<td>"+table_link_setup(OUTDIR,"psa2.sstreg.ann.png","Annual")+"</td><td>"+table_link_setup(OUTDIR,"psa2.sstreg.mon.png","Monthly")+"</td></tr>"
  txt(130) = "</tr></table>"

  txt(132) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(133) = "<CAPTION>AMOC</CAPTION>"
  txt(134) = "<tr><td>"+table_link_setup(OUTDIR,"amoc.mean.ann.png","Means")+"</td>"
  txt(135) = "<td>"+table_link_setup(OUTDIR,"amoc.stddev.ann.png","Standard Deviations")+"</td>"
  txt(136) = "<td>"+table_link_setup(OUTDIR,"amoc.ann.png","Patterns")+"</td>"
  txt(137) = "<td>"+table_link_setup(OUTDIR,"amoc.timeseries.ann.png","Timeseries")+"</td></tr>"
  txt(138) = "<tr><td>"+table_link_setup(OUTDIR,"amoc.powspec.ann.png","Spectra")+"</td>"
  txt(139) = "<td>"+table_link_setup(OUTDIR,"amoc.tasreg.ann.png","TAS Regressions")+"</td>"
  txt(140) = "<td>"+table_link_setup(OUTDIR,"amoc.sstreg.ann.png","SST Regressions")+"</td>"
  txt(141) = "<td>"+table_link_setup(OUTDIR,"amoc_amo.leadlag.ann.png","AMO/AMOC Lag Correlations")+"</td></tr></table>" 

  
  txt(162) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(163) = "<CAPTION>Global Trend Maps</CAPTION>"
  txt(164) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.trends.djf.png","DJF")+"</td>"
  txt(165) = "<td>"+table_link_setup(OUTDIR,"sst.trends.mam.png","MAM")+"</td>"
  txt(166) = "<td>"+table_link_setup(OUTDIR,"sst.trends.jja.png","JJA")+"</td>"
  txt(167) = "<td>"+table_link_setup(OUTDIR,"sst.trends.son.png","SON")+"</td>"
  txt(168) = "<td>"+table_link_setup(OUTDIR,"sst.trends.ann.png","Annual")+"</td>"
  txt(169) = "<td>"+table_link_setup(OUTDIR,"sst.trends.mon.png","Monthly")+"</td></tr>"
  txt(170) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"air2m.trends.djf.png","DJF")+"</td>"
  txt(171) = "<td>"+table_link_setup(OUTDIR,"air2m.trends.mam.png","MAM")+"</td>"
  txt(172) = "<td>"+table_link_setup(OUTDIR,"air2m.trends.jja.png","JJA")+"</td>"
  txt(173) = "<td>"+table_link_setup(OUTDIR,"air2m.trends.son.png","SON")+"</td>"
  txt(174) = "<td>"+table_link_setup(OUTDIR,"air2m.trends.ann.png","Annual")+"</td>"
  txt(175) = "<td>"+table_link_setup(OUTDIR,"air2m.trends.mon.png","Monthly")+"</td></tr>"
  txt(176) = "<tr><td><b>PSL</b></td><td>"+table_link_setup(OUTDIR,"psl.trends.djf.png","DJF")+"</td>"
  txt(177) = "<td>"+table_link_setup(OUTDIR,"psl.trends.mam.png","MAM")+"</td>"
  txt(178) = "<td>"+table_link_setup(OUTDIR,"psl.trends.jja.png","JJA")+"</td>"
  txt(179) = "<td>"+table_link_setup(OUTDIR,"psl.trends.son.png","SON")+"</td>"
  txt(180) = "<td>"+table_link_setup(OUTDIR,"psl.trends.ann.png","Annual")+"</td>"
  txt(181) = "<td>"+table_link_setup(OUTDIR,"psl.trends.mon.png","Monthly")+"</td></tr>"
  txt(182) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.trends.djf.png","DJF")+"</td>"
  txt(183) = "<td>"+table_link_setup(OUTDIR,"pr.trends.mam.png","MAM")+"</td>"
  txt(184) = "<td>"+table_link_setup(OUTDIR,"pr.trends.jja.png","JJA")+"</td>"
  txt(185) = "<td>"+table_link_setup(OUTDIR,"pr.trends.son.png","SON")+"</td>"
  txt(186) = "<td>"+table_link_setup(OUTDIR,"pr.trends.ann.png","Annual")+"</td>"
  txt(187) = "<td>"+table_link_setup(OUTDIR,"pr.trends.mon.png","Monthly")+"</td></tr>" 
  txt(188) = "<tr><td><b>SND</b></td><td>"+table_link_setup(OUTDIR,"snd.trends.djf.png","DJF")+"</td>"
  txt(189) = "<td>"+table_link_setup(OUTDIR,"snd.trends.mam.png","MAM")+"</td>"
  txt(190) = "<td>"+table_link_setup(OUTDIR,"snd.trends.jja.png","JJA")+"</td>"
  txt(191) = "<td>"+table_link_setup(OUTDIR,"snd.trends.son.png","SON")+"</td>"
  txt(192) = "<td>"+table_link_setup(OUTDIR,"snd.trends.ann.png","Annual")+"</td>"
  txt(193) = "<td>"+table_link_setup(OUTDIR,"snd.trends.mon.png","Monthly")+"</td></tr></table>"
 
  txt(194) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(195) = "<CAPTION>Global Timeseries</CAPTION>"
  txt(196) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.timeseries.djf.png","DJF")+"</td>"
  txt(197) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.mam.png","MAM")+"</td>"
  txt(198) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.jja.png","JJA")+"</td>"
  txt(199) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.son.png","SON")+"</td>"
  txt(200) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.ann.png","Annual")+"</td>"
  txt(201) = "<td>"+table_link_setup(OUTDIR,"sst.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(202) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"air2m.timeseries.djf.png","DJF")+"</td>"
  txt(203) = "<td>"+table_link_setup(OUTDIR,"air2m.timeseries.mam.png","MAM")+"</td>"
  txt(204) = "<td>"+table_link_setup(OUTDIR,"air2m.timeseries.jja.png","JJA")+"</td>"
  txt(205) = "<td>"+table_link_setup(OUTDIR,"air2m.timeseries.son.png","SON")+"</td>"
  txt(206) = "<td>"+table_link_setup(OUTDIR,"air2m.timeseries.ann.png","Annual")+"</td>"
  txt(207) = "<td>"+table_link_setup(OUTDIR,"air2m.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(208) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.timeseries.djf.png","DJF")+"</td>"
  txt(209) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.mam.png","MAM")+"</td>"
  txt(210) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.jja.png","JJA")+"</td>"
  txt(211) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.son.png","SON")+"</td>"
  txt(212) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.ann.png","Annual")+"</td>"
  txt(213) = "<td>"+table_link_setup(OUTDIR,"pr.timeseries.mon.png","Monthly")+"</td></tr>"
  txt(214) = "</table>"

  txt(215) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(216) = "<CAPTION>Global Timeseries Running Trends (Monthly)</CAPTION>"
  txt(217) = "<tr><td><b>SST</b></td><td>"+table_link_setup(OUTDIR,"sst.8yr_runtrend.mon.png","8yr")+"</td>"
  txt(218) = "<td>"+table_link_setup(OUTDIR,"sst.10yr_runtrend.mon.png","10yr")+"</td>"
  txt(219) = "<td>"+table_link_setup(OUTDIR,"sst.12yr_runtrend.mon.png","12yr")+"</td>"
  txt(220) = "<td>"+table_link_setup(OUTDIR,"sst.14yr_runtrend.mon.png","14yr")+"</td>"
  txt(221) = "<td>"+table_link_setup(OUTDIR,"sst.16yr_runtrend.mon.png","16yr")+"</td></tr>"
  txt(222) = "<tr><td><b>TAS</b></td><td>"+table_link_setup(OUTDIR,"air2m.8yr_runtrend.mon.png","8yr")+"</td>"
  txt(223) = "<td>"+table_link_setup(OUTDIR,"air2m.10yr_runtrend.mon.png","10yr")+"</td>"
  txt(224) = "<td>"+table_link_setup(OUTDIR,"air2m.12yr_runtrend.mon.png","12yr")+"</td>"
  txt(225) = "<td>"+table_link_setup(OUTDIR,"air2m.14yr_runtrend.mon.png","14yr")+"</td>"
  txt(226) = "<td>"+table_link_setup(OUTDIR,"air2m.16yr_runtrend.mon.png","16yr")+"</td></tr>"
  txt(227) = "</table>"

 
  txt(240) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(241) = "<CAPTION>Climatological Zonal Averages</CAPTION>"
  txt(242) = "<tr><td><b>PR</b></td><td>"+table_link_setup(OUTDIR,"pr.za.djf.png","DJF")+"</td>"
  txt(243) = "<td>"+table_link_setup(OUTDIR,"pr.za.mam.png","MAM")+"</td>"
  txt(244) = "<td>"+table_link_setup(OUTDIR,"pr.za.jja.png","JJA")+"</td>"
  txt(245) = "<td>"+table_link_setup(OUTDIR,"pr.za.son.png","SON")+"</td>"
  txt(246) = "<td>"+table_link_setup(OUTDIR,"pr.za.ann.png","Annual")+"</td></tr>"
  txt(247) = "</table>"
  
  txt(249) = "<br><TABLE width='600' border=1 frame=void rules=all cellpadding=6>"
  txt(250) = "<CAPTION>Additional Indices</CAPTION>"
  txt(251) = "<tr><td>"+table_link_setup(OUTDIR,"npi.timeseries.ndjfm.png","North Pacific Index")+"</td>"
  txt(252) = "<td>"+table_link_setup(OUTDIR,"tna.timeseries.png","Tropical North Atlantic SST")+"</td>"
  txt(253) = "<td>"+table_link_setup(OUTDIR,"tsa.timeseries.png","Tropical South Atlantic SST")+"</td>"
  txt(254) = "<td>"+table_link_setup(OUTDIR,"tio.timeseries.png","Tropical Indian Ocean SST")+"</td></tr>"
  txt(255) = "<tr><td>"+table_link_setup(OUTDIR,"nino12.timeseries.png","ni&ntilde;o1+2 Timeseries")+"</td>"
  txt(256) = "<td>"+table_link_setup(OUTDIR,"nino3.timeseries.png","ni&ntilde;o3 Timeseries")+"</td>"
  txt(257) = "<td>"+table_link_setup(OUTDIR,"nino4.timeseries.png","ni&ntilde;o4 Timeseries")+"</td>"
  txt(258) = "<td>"+table_link_setup(OUTDIR,"iod.timeseries.png","Indian Ocean Dipole")+"</td></tr>"
  txt(259) = "</table>"
  
;  if (OT.eq."True") then
;     txt(265) = "<p class=btext>Output timeseries written to "+OUTDIR+"</p>"
;  else
     txt(265) = "<br>"
;  end if  
  txt(266) = "</body></html>"
  
  tt = ind(.not.ismissing(txt))
  txt2 = txt(:tt(dimsizes(tt)-1)+1)
  txt2 = where(ismissing(txt2),"",txt2)
  asciiwrite(OUTDIR+"index.html",txt2)
  delete([/txt,tt,txt2,quote,names,na,nsim_noobs,OT/])
;----------------------------------------------------------------------------
;--   Create calculation description webpage
;----------------------------------------------------------------------------
  txt = new(200,"string")
  quote = str_get_dq()
  
  txt(0) = "<HTML><HEAD>"
  txt(1) = "<TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
  txt(2) = "img.title {width:990; height:105; align:left}"
  txt(3) = "p.subtitle1 {font-size:12pt; color:gray}"
  txt(4) = "p.subtitle {font-size:16pt; color:black; font-weight:bold; text-align:left}"
  txt(5) = "p.btext {font-size:9pt}"
  txt(6) = "td {text-align:left; font-size:12pt}"
  txt(7) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
  txt(8) = "</STYLE></HEAD><BODY BGCOLOR="+quote+"WhiteSmoke"+quote+">"
  txt(9) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"
  txt(10)= "<p class=subtitle>"
  txt(11) = "Methodology</p><p><a href="+quote+"index.html"+quote+">Back to Diagnostics Plots</a></p>"  
  txt(12) = "<p><ul><li>TS is surface ("+quote+"skin"+quote+") temperature and is used in lieu of sea surface temperatures (SSTs), TAS is 2m air temperature and is equivalent to TREFHT.<li>PR is total precipitation and is equivalent to CESM variables PRECC+PRECL, PSL is sea level pressure, and SNOWDP is snow depth."
  txt(13) = "<li>The annual cycle is removed prior to every calculation by subtracting the long-term monthly means. Exceptions: The annual cycle is not removed for mean and standard deviation spatial maps.<li> Area-averages are always based on cosine of latitude weighting.<li>The following calculations used detrended data: standard deviations, nino3.4 spectra, and ENSO spatial composites."
  txt(14) = "<li>For visual clarity, the Y-axis <i>may</i> differ amongst individual panels on a particular plot.<li>For more information on observational datasets and climate indices, see the <a href="+quote+"http://climatedataguide.ucar.edu/"+quote+" target="+quote+"_blank"+quote+">Climate Data Guide</a>.</ul></p>"

  txt(43) = "<TABLE width='800' border=1 frame=void rules=all cellpadding=6>"
  txt(44) = "<tr><td><b>AMO</b> (Atlantic Multidecadal Oscillation)</td>"  
  txt(45) = "<td>Monthly index timeseries defined as North Atlantic (0:60&deg;N, 80&deg;W:0&deg;E) SST snomalies minus global (60&deg;S:60&deg;N) SST anomalies. Pattern created by regressing global SST anomalies onto index timeseries and smoothing with a 9-point spatial filter. Low pass-filtered timeseries (black curve) is based on a a 61-month running mean. Based on Trenberth, K. E., and D. J. Shea, 2006: Atlantic hurricanes and natural variability in 2005, <i>Geophys. Res. Lett.</i>, 33, L12704, doi:10.1029/2006GL026894.</td></tr>"
  txt(46) = "<tr><td><b>PDO</b> (Pacific Decadal Oscillation)</td>"  
  txt(47) = "<td>Monthly index timeseries defined as the leading principal component (PC) of North Pacific (20:70&deg;N, 110&deg;E:100&deg;W) area-weighted SST* anomalies, where SST* denotes that the global mean SST has been removed at each timestep. Pattern created by regressing global SST anomalies onto normalized PC timeseries. Low pass-filtered timeseries (black curve) is based on a a 61-month running mean. See Deser, C., M. A. Alexander, S. -P. Xie, and A. S. Phillips, 2010: Sea surface temperature variability: patterns and mechanisms. <i>Ann. Rev. Mar. Sci.</i>, 2010.2, 115-143, doi:10.1146/annurev-marine-120408-151453. Also see Mantua, N. J., S. R. Hare, Y. Zhang, J. M. Wallace, and R. Francis, 1997: A Pacific interdecadal climate oscillation with impacts on salmon production. <i>Bull. Amer. Met. Soc.</i>, 1069-1079.</td></tr>"
  txt(48) = "<tr><td><b>ENSO Spatial Composites</b></td>"
  txt(49) = "<td>The normalized December nino3.4 timeseries is used to composite all years greater than 1 standard deviation (El Ni&ntilde;o) and all years less that -1 standard deviation (La Ni&ntilde;a). The December nino3.4 timeseries is based on the December values of the monthly nino3.4 time series smoothed with a 3-point binomial filter. See Deser, C., A. S. Phillips, R. A. Tomas, Y. Okumura, M. A. Alexander, A. Capotondi, J. D. Scott, Y. -O. Kwon, and M. Ohba, 2012: ENSO and Pacific Decadal Variability in Community Climate System Model Version 4. <i>J. Climate</i>, 25, 2622-2651, doi: 10.1175/JCLI-D-11-00301.1. The red/blue shading on the Nino3.4 time series denotes positive/negative departures from the best-fit linear trend line.</td></tr>"   
  txt(50) = "<tr><td><b>NAM</b> (Northern Annular Mode)</td>"  
  txt(51) = "<td>Seasonal/annual PSL averages are formed, square root of the cosine of the latitude weighting is applied, and then the leading EOF and associated principal component (PC) timeseries are computed over 20:90&deg;N, 0:360&deg;E. Pattern created by regressing global PSL anomalies onto normalized PC timeseries. Based on Hurrell, J. W., and C. Deser, 2009: North Atlantic climate variability: The role of the North Atlantic Oscillation. <i>J. Mar. Syst.</i>, 78, 28-41, doi:10.1016/j.jmarsys.2008.11.026. Also see Thompson, D. W. J., and J. M. Wallace, 2000: Annular modes in the extratropical circulation. Part I: Month-to-month variability. <i>J. Climate</i>, 13, 1000-1016.</td></tr>"  
  txt(52) = "<tr><td><b>NAO</b> (North Atlantic Oscillation)</td>"  
  txt(53) = "<td>Seasonal/annual PSL averages are formed, square root of the cosine of latitude weighting is applied, and then the leading EOF and associated principal component (PC) timeseries are computed over 20:80&deg;N, 90&deg;W:40&deg;E. Pattern created by regressing global PSL anomalies onto normalized PC timeseries. Based on Hurrell, J. W. and C. Deser, 2009: North Atlantic climate variability: The role of the North Atlantic Oscillation. <i>J. Mar. Syst.</i>, 78, 28-41, doi:10.1016/j.jmarsys.2008.11.026.</td></tr>"  
  txt(54) = "<tr><td><b>SAM/PSA1/PSA2</b> (Southern Annular Mode, Pacific South American Patterns 1/2)</td>"  
  txt(55) = "<td>Seasonal/annual PSL averages are formed, square root of the cosine of latitude weighting is applied, and then the 1st (SAM), 2nd (PSA1) and 3rd (PSA2) EOFs and associated principal component (PC) timeseries are computed over 20:90&deg;S, 0:360&deg;E. Patterns created by regressing global PSL anomalies onto normalized PC timeseries. SAM calculation based on Thompson, D. W. J. and J.M. Wallace, 2000: Annular modes in the extratropical circulation. Part I: Month-to-month variability. <i>J. Climate</i>, 13, 1000-1016.</td></tr>"  
  txt(56) = "<tr><td><b>PNA/NPO</b> (Pacific North American Pattern, North Pacific Oscillation)</td>" 
  txt(57) = "<td>Seasonal/annual PSL averages are formed, the square root of the cosine of the latitude weighting is applied, and then the 1st (PNA) and 2nd (NPO) EOFs and associated principal component (PC) timeseries are computed over 20:85&deg;N, 120&deg;E:120&deg;W. Patterns created by regressing global PSL anomalies onto normalized PC timeseries.</td></tr>" 
  txt(58) = "<tr><td><b>Atmospheric Mode SST Regressions</b></td>"
  txt(59) = "<td>SST anomalies are regressed upon the normalized atmospheric mode timeseries.</td></tr>" 
  txt(60) = "<tr><td><b>Global Average Running Trends</b></td>"
  txt(61) = "<td>N-year running trends are computed by calculating the linear trend over the N-year interval beginning at each successive timestep.  For instance, for a global timeseries that runs from 1970-2012, the 8yr running trend value for January 1970 is the linear trend during January 1970 - December 1977, and the value for January 2005 is the linear trend during January 2005 - December 2012.</td></tr>"  
  txt(62) = "<tr><td><b>Metrics Table</b></td>"
  txt(63) = "<td>Area-weighted pattern correlations and rms differences are calculated between observations and each model simulation (regridded to match the observational grid) for 11 climate metrics. The Total Score column shows the average of the 11 pattern correlations (Z-transformed) and rms differences. Domains Used: Means, standard deviations, AMO, and PDO: Global. Hovmollers: Entire longitude/temporal range shown. NAM: Entire domain shown (20:90&deg;N). SAM: Entire domain shown (20:90&deg;S).</td></tr>" 
  txt(64) = "<tr><td><b>Power Spectra</b></td>"
  txt(65) = "<td>The best-fit first-order Markov red noise spectrum (red curve) and its 95% (blue curve) and 99% (green curve) confidence bounds are shown on each panel. Top X-axis shows the period (in years), and the bottom X-axis shows the frequency (cycles/mo). If calculated, the observational spectrum is overlaid in gray on each model spectrum. The spectra are displayed in variance-conserving form.</td></tr>"
  txt(66) = "<tr><td><b>Climatological Zonal Averages</b></td>"
  txt(67) = "<td>Climatological means are zonally averaged over the globe.</td></tr>"
  txt(68) = "<tr><td><b>NPI</b> (North Pacific Index)</td>"
  txt(69) = "<td>Winter (December-March) average PSL anomalies area-averaged over 30&deg;:65&deg;N, 160&deg;E:140&deg;W. Based on Trenberth, K. E. and J. W. Hurrell, 1994: Decadal atmosphere-ocean variations in the Pacific, <i>Climate Dynamics</i>, 9, 303-319.</td></tr>"  
  txt(70) = "<tr><td><b>Additional Indices</b></td>"
  txt(71) = "<td>Area-averaged anomalies as specified in the plot title. Red/blue shading denotes positive/negative departures from the best-fit linear trend line.</td></tr>"  
  txt(70) = "<tr><td><b>AMOC</b></td>"
  txt(71) = "<td>MOC annual averages are formed, weighted by the cosine of the latitude and by the depth of the vertical layer, and then the data is masked by setting all those areas to missing where the variance is less than 1.e-6. The leading EOF and associated principal component (PC) timeseries are then computed over the Atlantic basin from 33&deg;S:90&deg;N. The AMOC patterns are created by regressing Atlantic MOC anomalies onto normalized PC timeseries. The SST/TAS patterns are created by regressing global TAS/SST anomalies onto the normalized PC timeseries. A 15pt low-pass filter is applied to the AMO and AMOC PC timeseries prior to the lead/lag correlation calculations, and at least 90 years of data are required. Unlike the cited reference paper the data is never detrended. For CCSM4 and CESM1 the MOC variable is read in, the components are summed over, and the Atlantic Ocean + Mediterranean Sea + Labrador Sea + GIN Sea + Arctic Ocean + Hudson Bay transport region is selected. For CCSM2 and CCSM3 the same transport region is selected but only the Eulerian Mean component is used as that is all that is available. For CMIP5 data the atlantic_arctic_ocean basin is used, while for CMIP3 data the atlantic_ocean geo_region is used. See Danabasoglu, G., S. G. Yeager, Y. -O. Kwon, J. J. Tribbia, A. S. Phillips, and J. W. Hurrell, 2012. Variability of the Atlantic Meridional Overturning Circulation in CCSM4. <i>J. Climate</i>, 25, 5153-5172, doi: 10.1175/JCLI-D-11-00463.1.</td></tr>" 

  txt(75) = "</table>"
  txt(76) = ""
  
  txt(100) = "<br>Created "+systemfunc("date")
  txt(101) = "<p class=btext>CVDP Version "+VERSION+"</p>"
  txt(102) = "</body></html>"
  
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
  
     txt(0) = "<HTML><HEAD>"
     txt(1) = "<TITLE>Climate Variability Diagnostics Package</TITLE><STYLE>"
     txt(2) = "img.title {width:990; height:105; align:left}"
     txt(3) = "p.subtitle1 {font-size:15pt; color:gray}"
     txt(4) = "p.subtitle {font-size:20pt; color:black; text-align:center; font-weight:bold}"
     txt(5) = "p.btext {font-size:12pt; text-align:left}"
     txt(6) = "p.btext2 {font-size:9pt; text-align:left}"
     txt(7) = "td {text-align:center; font-size:12pt}"
     txt(8) = "CAPTION {text-align:left; font-size:14pt; font-weight:bold}"
     txt(9) = "</STYLE></HEAD><BODY BGCOLOR="+quote+"WhiteSmoke"+quote+">"
     txt(10) = "<img class=title src="+quote+"cas-cvdp.png"+quote+" alt="+quote+"CVDP logo"+quote+">"
     txt(11)= "<TABLE width='990' border=0 cellpadding=6>"
     txt(12) = "<tr><td width='300'><p class=btext><a href="+quote+"index.html"+quote+">Back to Diagnostics Plots</a></p></td>"
     txt(13) = "<td><p class=subtitle>"+WEBTITLE+"</p></td></tr></table><hr width=990 align=left><p class=subtitle1>Metrics Table</p>" 
     
     z = asciiread(OUTDIR+"metrics.txt",(/-1/),"string")
     nlines = dimsizes(z)
     txt(14) = "<PRE>"
     do gg = 0,nlines-1
        txt(15+gg) = z(gg)
     end do
     txt(15+nlines) = "</PRE>"
     txt(15+nlines+1) = "<p class=btext><a href="+quote+"namelist_obs"+quote+">Observations Used</a></p>"
     txt(15+nlines+2) = "Created "+systemfunc("date")
     txt(15+nlines+3) = "<p class=btext2>CVDP Version "+VERSION+"</p></body></html>"
  
     tt = ind(ismissing(txt))
     txt2 = txt(:tt(0)-1)
     asciiwrite(OUTDIR+"metrics.html",txt2)
     delete([/txt,txt2,tt,OUTDIR/])
  end if 
  print("Finished: webpage.ncl")  
end
