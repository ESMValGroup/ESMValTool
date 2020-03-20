;********************************************
; Lisp code for an NCL major mode
;********************************************
; Updated: Thu Mar  2 17:39:33 MST 2017
;
; Modified version for ESMValTool v2.0
; Mattia Righi, DLR (Germany)
;
; Revision 0.34
;  - Updated to include new functions, resources, and 
;    keywords added in NCL 6.0.0
; 
; Revision 0.33
; Changes to 0.32 by T. Corti, ETH Zurich and David Brown, 
; Changes between 0.2 and 0.3 by C. Schreck and A. Srock, University at Albany
; Changes between 0.1 and 0.2 Provided by Heiko Klein of Norway

; August 19 2003 Sylvia Murphy 
; National Center for Atmospheric Research
; Does text highlighting for NCL reserved words, built-in functions, 
; gsn* functions, contributed and shea-util functions, text, and comments.
; Does automatic indenting between begin and end statments, and within
; do loops and if statements.
; 
; Emacs has a lot more power that these functions. I do not use that
; functionality, so i did not spend any more time trying to add abbreviations,
; special keymaps etc.
;
; Updates in version 0.32
; Added Comment Handling (M-;). 
;  - Insert a comment at the end of the current line
;  - Alternatively comment/uncomment selected region
; Use syntactic fontification for comments and strings
; Correct fontification of strings containing a semicolon (;)
; Added highlightning for resources using font-lock-constant-face
; All documented functions are now highlighted (modification by D. Brown)
;
; Updates in version 0.3:
; Added more keywords (full list from NCL documentation)
; Changed color mapping (font-lock) settings:
;   - removed usage of font-lock-reference-face
;   - added usage of font-lock-builtin-face
;   - NCL built-in functions now use font-lock-builtin-face
;   - contributed and shea_util functions now use font-lock-function-face
;   - added boolean and value test keywords
;   - added keywords for beginning and ending arrays: (/ and /)
;   - all keywords now use font-lock-keyword-face
;   - explicitly fontifies strings with font-lock-string-face
; Changed syntax type of underscore to "word" instead of punctuation
; Updated "How to Use" instructions for ease of inclusion with Xemacs
;
; KNOWN PROBLEMS in version 0.32:
; 1) Comment Handling does not work in xemacs
; 2) Comments may not fontify on file open in xemacs 
;
; KNOWN PROBLEMS in version 0.3:
; 1) Comments with embedded strings don't initially fontify properly, but
;    do change if line modified somehow
; 2) Strings containing a semicolon (;) do not fontify properly
;
; KNOWN PROBLEMS THAT VERSION 0.2 fixed
; 1) Works with xemacs 21.*
; 2) Works with emacs 20.*

; KNOWN PROBLEMS in Version 0.1
; 1) Only partially works with emacs version 20.3.2
;    a) highlights only comments and text, and only after tabs
;    b) indentation appears to work
; 2) Does not work with xemacs
; 3) Not all NCL built-in functions are highlighted. I listed MY favorite 
;    ones.
; 4) Have not demonstrated how to change the indentation value in .emacs
; 5) The ncl-in-comment function does not work. Its calls are commented out.
;
;********************************************
; HOW TO USE
;********************************************
; 1) place this file somewhere on your local system e.g. ~your_home/bin

; 2) in your .emacs or .xemacs/custom.el file, add and properly modify //
; the following (without the comments):
  ;(setq auto-mode-alist (cons '("\.ncl$" . ncl-mode) auto-mode-alist))
  ;(autoload 'ncl-mode "LOCATION/ncl.el") 
  ;(add-hook 'ncl-mode-hook
  ;      (lambda ()  
  ;       )
  ;   )

; 3) setup display colors for font-lock.  You may also want to set default
; foreground and background colors.  Colors can be Xwindows names or #rrggbb.
; These should also go somewhere in your .emacs or .xemacs/custom.el file.
;     ; highlight comments
;         (set-face-foreground font-lock-comment-face "FireBrick")
;     ; highlight strings
;         (set-face-foreground font-lock-string-face "SlateGray")
;     ; highlight keywords, array descriptors, and tests                             
;         (set-face-foreground font-lock-keyword-face "Purple")
;     ; highlight built-in functions
;         (set-face-foreground font-lock-builtin-face "Blue")
;     ; highlight gsn* functions
;         (set-face-foreground font-lock-variable-name-face "SteelBlue")
;     ; highlight shea_util and contributed functions
;         (set-face-foreground font-lock-function-name-face  "CadetBlue")
;     ; highlight ESMValTool interface and lib functions
;         (set-face-foreground font-lock-type-face  "OrangeRed1") 
;     ; highlight resources
;         (set-face-foreground font-lock-constant-face  "ForestGreen")


;********************************************************************
(defvar ncl-mode-hook nil
  "*List of functions to call when entering ncl mode.")

(defvar ncl-font-lock-keywords
   '(


    ;; comments. the period (.) means a ; and any character after it except a 
    ;; newline while the asterisk (*) means repeated all occurrences.

    ; this is only for XEmacs!
    ("\\(;.*\\)" 1 font-lock-comment-face)
	
 
    ;; strings.  .*? means the shortest possible group of characters within
    ;; the quotes (only on one line)
    ;("\\(\".*?\"\\)" 1 font-lock-string-face )

    ;; NCL keywords
    ("\\<\\(begin\\|break\\|byte\\|character\\|continue\\|create\\|defaultapp\\|do\\|double\\|else\\|elseif\\|end\\|enumeric\\|external\\|file\\|float\\|function\\|getvalues\\|graphic\\|group\\|if\\|integer\\|int64\\|list\\|load\\|local\\|logical\\|long\\|new\\|_Missing\\|Missing\\|new\\|noparent\\|numeric\\|procedure\\|quit\\|QUIT\\|Quit\\|record\\|return\\|setvalues\\|short\\|snumeric\\|stop\\|string\\|then\\|ubyte\\|uint\\|uint64\\|ulong\\|ushort\\|while\\|\\)\\>" 1 font-lock-keyword-face)

    ;; Array definition chars and tests - couldn't get to work in list above...
    ("\\(\(\/\\)" 1 font-lock-keyword-face )
    ("\\(\/\)\\)" 1 font-lock-keyword-face )
    ("\\(\\\\\\\)" 1 font-lock-keyword-face )
    ("\\(->\\)" 1 font-lock-keyword-face )
    ("\\(\\.eq\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.ne\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.lt\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.le\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.gt\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.ge\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.and\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.or\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.not\\.\\)" 1 font-lock-keyword-face )
    ("\\(\\.xor\\.\\)" 1 font-lock-keyword-face )

    ;; ncl built-in functions
    ("\\<\\(abs\\|acos\\|addfile\\|addfiles\\|all\\|angmom_atm\\|any\\|area_conserve_remap\\|area_hi2lores\\|area_poly_sphere\\|asciiread\\|asciiwrite\\|asin\\|atan\\|atan2\\|attsetvalues\\|avg\\|betainc\\|bin_avg\\|bin_sum\\|bw_bandpass_filter\\|cancor\\|cbinread\\|cbinwrite\\|cd_calendar\\|cd_inv_calendar\\|cdfbin_p\\|cdfbin_pr\\|cdfbin_s\\|cdfbin_xn\\|cdfchi_p\\|cdfchi_x\\|cdfgam_p\\|cdfgam_x\\|cdfnor_p\\|cdfnor_x\\|cdft_p\\|cdft_t\\|ceemdan\\|ceil\\|center_finite_diff\\|center_finite_diff_n\\|cfftb\\|cfftf\\|cfftf_frq_reorder\\|charactertodouble\\|charactertofloat\\|charactertointeger\\|charactertolong\\|charactertoshort\\|charactertostring\\|chartodouble\\|chartofloat\\|chartoint\\|chartointeger\\|chartolong\\|chartoshort\\|chartostring\\|chiinv\\|cla_sq\\|clear\\|color_index_to_rgba\\|conform\\|conform_dims\\|cos\\|cosh\\|count_unique_values\\|count_unique_values_n\\|covcorm\\|covcorm_xy\\|craybinnumrec\\|craybinrecread\\|create_graphic\\|csa1\\|csa1d\\|csa1s\\|csa1x\\|csa1xd\\|csa1xs\\|csa2\\|csa2d\\|csa2l\\|csa2ld\\|csa2ls\\|csa2lx\\|csa2lxd\\|csa2lxs\\|csa2s\\|csa2x\\|csa2xd\\|csa2xs\\|csa3\\|csa3d\\|csa3l\\|csa3ld\\|csa3ls\\|csa3lx\\|csa3lxd\\|csa3lxs\\|csa3s\\|csa3x\\|csa3xd\\|csa3xs\\|csc2s\\|csgetp\\|css2c\\|cssetp\\|cssgrid\\|csstri\\|csvoro\\|cumsum\\|cz2ccm\\|datatondc\\|day_of_week\\|day_of_year\\|days_in_month\\|default_fillvalue\\|delete\\|depth_to_pres\\|destroy\\|determinant\\|dewtemp_trh\\|dgeevx_lapack\\|dim_acumrun_n\\|dim_avg\\|dim_avg_n\\|dim_avg_wgt\\|dim_avg_wgt_n\\|dim_cumsum\\|dim_cumsum_n\\|dim_gamfit_n\\|dim_gbits\\|dim_max\\|dim_max_n\\|dim_median\\|dim_median_n\\|dim_min\\|dim_min_n\\|dim_num\\|dim_num_n\\|dim_numrun_n\\|dim_pqsort\\|dim_pqsort_n\\|dim_product\\|dim_product_n\\|dim_rmsd\\|dim_rmsd_n\\|dim_rmvmean\\|dim_rmvmean_n\\|dim_rmvmed\\|dim_rmvmed_n\\|dim_spi_n\\|dim_standardize\\|dim_standardize_n\\|dim_stat4\\|dim_stat4_n\\|dim_stddev\\|dim_stddev_n\\|dim_sum\\|dim_sum_n\\|dim_sum_wgt\\|dim_sum_wgt_n\\|dim_variance\\|dim_variance_n\\|dimsizes\\|doubletobyte\\|doubletochar\\|doubletocharacter\\|doubletofloat\\|doubletoint\\|doubletointeger\\|doubletolong\\|doubletoshort\\|dpres_hybrid_ccm\\|dpres_plevel\\|draw\\|draw_color_palette\\|dsgetp\\|dsgrid2\\|dsgrid2d\\|dsgrid2s\\|dsgrid3\\|dsgrid3d\\|dsgrid3s\\|dspnt2\\|dspnt2d\\|dspnt2s\\|dspnt3\\|dspnt3d\\|dspnt3s\\|dssetp\\|dtrend\\|dtrend_msg\\|dtrend_msg_n\\|dtrend_n\\|dtrend_quadratic\\|dtrend_quadratic_msg_n\\|dv2uvF\\|dv2uvf\\|dv2uvG\\|\\)\\>" 1 font-lock-builtin-face)

    ("\\<\\(dv2uvg\\|dz_height\\|echo_off\\|echo_on\\|eemd\\|eof2data\\|eof2data_n\\|eof_varimax\\|eofcor\\|eofcor_pcmsg\\|eofcor_ts\\|eofcov\\|eofcov_pcmsg\\|eofcov_ts\\|eofunc\\|eofunc_n\\|eofunc_ts\\|eofunc_ts_n\\|eofunc_varimax\\|equiv_sample_size\\|erf\\|erfc\\|esacr\\|esacv\\|esccr\\|esccv\\|escorc\\|escorc_n\\|escovc\\|exit\\|exp\\|exp_tapersh\\|exp_tapersh_wgts\\|exp_tapershC\\|extval_mlegev\\|ezfftb\\|ezfftb_n\\|ezfftf\\|ezfftf_n\\|f2fosh\\|f2foshv\\|f2fsh\\|f2fshv\\|f2gsh\\|f2gshv\\|fabs\\|fbindirread\\|fbindirwrite\\|fbinnumrec\\|fbinread\\|fbinrecread\\|fbinrecwrite\\|fbinwrite\\|fft2db\\|fft2df\\|fftshift\\|fileattdef\\|filechunkdimdef\\|filedimdef\\|fileexists\\|filegrpdef\\|filevarattdef\\|filevarchunkdef\\|filevarcompressleveldef\\|filevardef\\|filevardimsizes\\|filwgts_lancos\\|filwgts_lanczos\\|filwgts_normal\\|floattobyte\\|floattochar\\|floattocharacter\\|floattoint\\|floattointeger\\|floattolong\\|floattoshort\\|floor\\|fluxEddy\\|fo2fsh\\|fo2fshv\\|fourier_info\\|frame\\|fspan\\|ftcurv\\|ftcurvd\\|ftcurvi\\|ftcurvp\\|ftcurvpi\\|ftcurvps\\|ftcurvs\\|ftest\\|ftgetp\\|ftkurv\\|ftkurvd\\|ftkurvp\\|ftkurvpd\\|ftsetp\\|ftsurf\\|g2fsh\\|g2fshv\\|g2gsh\\|g2gshv\\|gamma\\|gammainc\\|gaus\\|gaus_lobat\\|gaus_lobat_wgt\\|gc_aangle\\|gc_clkwise\\|gc_dangle\\|gc_inout\\|gc_latlon\\|gc_onarc\\|gc_pnt2gc\\|gc_qarea\\|gc_tarea\\|generate_2d_array\\|get_color_rgba\\|get_cpu_time\\|get_isolines\\|get_ncl_version\\|get_script_name\\|get_script_prefix_name\\|get_sphere_radius\\|get_unique_values\\|getbitsone\\|getenv\\|getfiledimsizes\\|getfilegrpnames\\|getfilepath\\|getfilevaratts\\|getfilevarchunkdimsizes\\|getfilevardims\\|getfilevardimsizes\\|getfilevarnames\\|getfilevartypes\\|getvaratts\\|getvardims\\|gradsf\\|gradsg\\|greg2jul\\|grid2triple\\|hlsrgb\\|hsvrgb\\|hydro\\|hyi2hyo\\|idsfft\\|igradsf\\|igradsF\\|igradsG\\|igradsg\\|ilapsf\\|ilapsF\\|ilapsg\\|ilapsG\\|ilapvf\\|ilapvg\\|ind\\|ind_resolve\\|int2p\\|int2p_n\\|integertobyte\\|integertochar\\|integertocharacter\\|integertoshort\\|inttobyte\\|inttochar\\|inttoshort\\|inverse_matrix\\|is_string_numeric\\|isatt\\|isbigendian\\|isbyte\\|ischar\\|iscoord\\|isdefined\\|isdim\\|isdimnamed\\|isdouble\\|isenumeric\\|isfile\\|isfilepresent\\|isfilevar\\|isfilevaratt\\|isfilevarcoord\\|isfilevardim\\|isfloat\\|isfunc\\|isgraphic\\|isint\\|isint64\\|isinteger\\|isleapyear\\|islogical\\|islong\\|ismissing\\|isnan_ieee\\|isnumeric\\|ispan\\|isproc\\|\\)\\>" 1 font-lock-builtin-face)

    ("\\<\\(isscalar\\|isshort\\|issnumeric\\|isstring\\|isubyte\\|isuint\\|isuint64\\|isulong\\|isunlimited\\|isunsigned\\|isushort\\|isvar\\|jul2greg\\|kmeans_as136\\|kolsm2_n\\|kron_product\\|lapsF\\|lapsf\\|lapsG\\|lapsg\\|lapvf\\|lapvg\\|latlon2utm\\|lclvl\\|lderuvf\\|lderuvg\\|linint1\\|linint1_n\\|linint2\\|linint2_points\\|linmsg\\|linmsg_n\\|linrood_latwgt\\|linrood_wgt\\|list_files\\|list_filevars\\|list_hlus\\|list_procfuncs\\|list_vars\\|ListAppend\\|ListCount\\|ListGetType\\|ListIndex\\|ListIndexFromName\\|ListPop\\|ListPush\\|ListSetType\\|loadscript\\|local_max\\|local_min\\|log\\|log10\\|longtobyte\\|longtochar\\|longtocharacter\\|longtoint\\|longtointeger\\|longtoshort\\|lspoly\\|lspoly_n\\|mask\\|max\\|maxind\\|min\\|minind\\|mixed_layer_depth\\|mixhum_ptd\\|mixhum_ptrh\\|mjo_cross_coh2pha\\|mjo_cross_segment\\|moc_globe_atl\\|monthday\\|namedcolor2rgb\\|namedcolor2rgba\\|natgrid\\|natgridd\\|natgrids\\|ncargpath\\|ncargversion\\|ndctodata\\|ndtooned\\|new\\|NewList\\|ngezlogo\\|nggcog\\|nggetp\\|nglogo\\|ngsetp\\|NhlAddAnnotation\\|NhlAddData\\|NhlAddOverlay\\|NhlAddPrimitive\\|NhlAppGetDefaultParentId\\|NhlChangeWorkstation\\|NhlClassName\\|NhlClearWorkstation\\|NhlDataPolygon\\|NhlDataPolyline\\|NhlDataPolymarker\\|NhlDataToNDC\\|NhlDestroy\\|NhlDraw\\|NhlFrame\\|NhlFreeColor\\|NhlGetBB\\|NhlGetClassResources\\|NhlGetErrorObjectId\\|NhlGetNamedColorIndex\\|NhlGetParentId\\|NhlGetParentWorkstation\\|NhlGetWorkspaceObjectId\\|NhlIsAllocatedColor\\|NhlIsApp\\|NhlIsDataComm\\|NhlIsDataItem\\|NhlIsDataSpec\\|NhlIsTransform\\|NhlIsView\\|NhlIsWorkstation\\|NhlName\\|NhlNDCPolygon\\|NhlNDCPolyline\\|NhlNDCPolymarker\\|NhlNDCToData\\|NhlNewColor\\|NhlNewDashPattern\\|NhlNewMarker\\|NhlPalGetDefined\\|NhlRemoveAnnotation\\|NhlRemoveData\\|NhlRemoveOverlay\\|NhlRemovePrimitive\\|NhlSetColor\\|NhlSetDashPattern\\|NhlSetMarker\\|NhlUpdateData\\|NhlUpdateWorkstation\\|nice_mnmxintvl\\|nngetaspectd\\|nngetaspects\\|nngetp\\|nngetsloped\\|nngetslopes\\|nngetwts\\|nngetwtsd\\|nnpnt\\|nnpntd\\|nnpntend\\|nnpntendd\\|nnpntinit\\|nnpntinitd\\|nnpntinits\\|nnpnts\\|nnsetp\\|num\\|obj_anal_ic\\|omega_ccm\\|onedtond\\|overlay\\|paleo_outline\\|pdfxy_bin\\|poisson_grid_fill\\|pop_remap\\|potmp_insitu_ocn\\|prcwater_dp\\|pres2hybrid\\|pres_hybrid_ccm\\|pres_hybrid_jra55\\|pres_sigma\\|print\\|print_table\\|printFileVarSummary\\|printVarSummary\\|product\\|pslec\\|pslhor\\|pslhyp\\|qsort\\|rand\\|random_chi\\|random_gamma\\|random_normal\\|random_setallseed\\|random_uniform\\|rcm2points\\|rcm2rgrid\\|rdsstoi\\|read_colormap_file\\|reg_multlin\\|regCoef\\|regcoef\\|regCoef_n\\|regline\\|relhum\\|relhum_ice\\|relhum_water\\|replace_ieeenan\\|reshape\\|reshape_ind\\|rgbhls\\|\\)\\>" 1 font-lock-builtin-face)

    ("\\<\\(rgbhsv\\|rgbyiq\\|rgrid2rcm\\|rhomb_trunc\\|rhomb_trunC\\|rip_cape_2d\\|rip_cape_3d\\|round\\|rtest\\|runave\\|runave_n\\|set_default_fillvalue\\|set_sphere_radius\\|setfileoption\\|sfvp2uvf\\|sfvp2uvg\\|shaeC\\|shaec\\|shagC\\|shagc\\|shgetnp\\|shgetp\\|shgrid\\|shorttobyte\\|shorttochar\\|shorttocharacter\\|show_ascii\\|shseC\\|shsec\\|shsetp\\|shsgC\\|shsgc\\|shsgc_R42\\|sigma2hybrid\\|simpeq\\|simpne\\|sin\\|sindex_yrmo\\|sinh\\|sizeof\\|sleep\\|smth9\\|snindex_yrmo\\|solve_linsys\\|span_color_indexes\\|span_color_rgba\\|span_named_colors\\|sparse_matrix_mult\\|spcorr\\|spcorr_n\\|specx_anal\\|specxy_anal\\|sprintf\\|sprinti\\|sqrt\\|sqsort\\|srand\\|stat2\\|stat4\\|stat_medrng\\|stat_trim\\|status_exit\\|stdatmus_p2tdz\\|stdatmus_z2tdp\\|stddev\\|str_capital\\|str_concat\\|str_fields_count\\|str_get_cols\\|str_get_dq\\|str_get_field\\|str_get_nl\\|str_get_sq\\|str_get_tab\\|str_index_of_substr\\|str_insert\\|str_is_blank\\|str_join\\|str_left_strip\\|str_lower\\|str_match\\|str_match_ic\\|str_match_ic_regex\\|str_match_ind\\|str_match_ind_ic\\|str_match_ind_ic_regex\\|str_match_ind_regex\\|str_match_regex\\|str_right_strip\\|str_split\\|str_split_by_length\\|str_split_csv\\|str_squeeze\\|str_strip\\|str_sub_str\\|str_switch\\|str_upper\\|stringtochar\\|stringtocharacter\\|stringtodouble\\|stringtofloat\\|stringtoint\\|stringtointeger\\|stringtolong\\|stringtoshort\\|strlen\\|student_t\\|sum\\|svd_lapack\\|svdcov\\|svdcov_sv\\|svdstd\\|svdstd_sv\\|system\\|systemfunc\\|tan\\|tanh\\|taper\\|taper_n\\|tdclrs\\|tdctri\\|tdcudp\\|tdcurv\\|tddtri\\|tdez2d\\|tdez3d\\|tdgetp\\|tdgrds\\|tdgrid\\|tdgtrs\\|tdinit\\|tditri\\|tdlbla\\|tdlblp\\|tdlbls\\|tdline\\|tdlndp\\|tdlnpa\\|tdlpdp\\|tdmtri\\|tdotri\\|tdpara\\|tdplch\\|tdprpa\\|tdprpi\\|tdprpt\\|tdsetp\\|tdsort\\|tdstri\\|tdstrs\\|tdttri\\|thornthwaite\\|tobyte\\|tochar\\|todouble\\|tofloat\\|toint\\|toint64\\|tointeger\\|tolong\\|toshort\\|tosigned\\|tostring\\|tostring_with_format\\|totype\\|toubyte\\|touint\\|touint64\\|toulong\\|tounsigned\\|toushort\\|trend_manken\\|tri_trunC\\|tri_trunc\\|triple2grid\\|triple2grid2d\\|trop_wmo\\|ttest\\|typeof\\|undef\\|unique_string\\|update\\|ushorttoint\\|ut_calendar\\|ut_calendar_fix\\|ut_inv_calendar\\|ut_inv_calendar_fix\\|utm2latlon\\|uv2dv_cfd\\|uv2dvf\\|uv2dvF\\|uv2dvg\\|uv2dvG\\|uv2sfvpF\\|uv2sfvpf\\|uv2sfvpG\\|uv2sfvpg\\|uv2vr_cfd\\|uv2vrdvF\\|uv2vrdvf\\|\\)\\>" 1 font-lock-builtin-face)

    ("\\<\\(uv2vrdvG\\|uv2vrdvg\\|uv2vrF\\|uv2vrf\\|uv2vrG\\|uv2vrg\\|v5d_close\\|v5d_create\\|v5d_setLowLev\\|v5d_setUnits\\|v5d_write\\|v5d_write_var\\|variance\\|vhaeC\\|vhaec\\|vhagC\\|vhagc\\|vhsec\\|vhseC\\|vhsgc\\|vhsgC\\|vibeta\\|vinth2p\\|vinth2p_ecmwf\\|vinth2p_ecmwf_nodes\\|vinth2p_nodes\\|vintp2p_ecmwf\\|vr2uvF\\|vr2uvf\\|vr2uvG\\|vr2uvg\\|vrdv2uvf\\|vrdv2uvF\\|vrdv2uvg\\|vrdv2uvG\\|wavelet\\|wavelet_default\\|weibull\\|wetbulb\\|wgt_area_smooth\\|wgt_areaave\\|wgt_areaave2\\|wgt_arearmse\\|wgt_arearmse2\\|wgt_areasum2\\|wgt_runave\\|wgt_runave_n\\|wgt_vert_avg_beta\\|wgt_volave\\|wgt_volave_ccm\\|wgt_volrmse\\|wgt_volrmse_ccm\\|where\\|wk_smooth121\\|wmbarb\\|wmbarbmap\\|wmdrft\\|wmgetp\\|wmlabs\\|wmsetp\\|wmstnm\\|wmvect\\|wmvectmap\\|wmvlbl\\|wrf_avo\\|wrf_cape_2d\\|wrf_cape_3d\\|wrf_dbz\\|wrf_eth\\|wrf_helicity\\|wrf_ij_to_ll\\|wrf_interp_1d\\|wrf_interp_2d_xy\\|wrf_interp_3d_z\\|wrf_latlon_to_ij\\|wrf_ll_to_ij\\|wrf_omega\\|wrf_pvo\\|wrf_rh\\|wrf_slp\\|wrf_smooth_2d\\|wrf_td\\|wrf_tk\\|wrf_updraft_helicity\\|wrf_uvmet\\|wrf_virtual_temp\\|wrf_wetbulb\\|wrf_wps_close_int\\|wrf_wps_open_int\\|wrf_wps_rddata_int\\|wrf_wps_rdhead_int\\|wrf_wps_read_int\\|wrf_wps_write_int\\|write_matrix\\|write_table\\|yiqrgb\\|z2geouv\\|zonal_mpsi\\|\\)\\>" 1 font-lock-builtin-face)

    ;; bootstrap functions
    ("\\<\\(bootstrap_correl\\|bootstrap_diff\\|bootstrap_estimate\\|bootstrap_regcoef\\|bootstrap_stat\\|\\)\\>" 1 font-lock-function-name-face)

    ;; contributed functions
    ("\\<\\(addfiles_GetVar\\|advect_variable\\|albedo_ccm\\|area_conserve_remap_Wrap\\|area_hi2lores_Wrap\\|array_append_record\\|assignFillValue\\|brunt_vaisala_atm\\|byte2flt\\|byte2flt_hdf\\|calcDayAnomTLL\\|calcMonAnomLLLT\\|calcMonAnomLLT\\|calcMonAnomTLL\\|calcMonAnomTLLL\\|calculate_daily_values\\|calculate_monthly_values\\|calculate_segment_values\\|cd_convert\\|changeCase\\|changeCaseChar\\|clmDayTLL\\|clmDayTLLL\\|clmMon2clmDay\\|clmMonLLLT\\|clmMonLLT\\|clmMonTLL\\|clmMonTLLL\\|closest_val\\|cohsq_c2p\\|cohsq_p2c\\|copy_VarAtts\\|copy_VarCoords\\|copy_VarCoords_1\\|copy_VarCoords_2\\|copy_VarMeta\\|copyatt\\|coriolis_param\\|crossp3\\|cshstringtolist\\|cssgrid_Wrap\\|dble2flt\\|decimalPlaces\\|delete_VarAtts\\|demod_cmplx\\|dim_avg_n_Wrap\\|dim_avg_wgt_n_Wrap\\|dim_avg_wgt_Wrap\\|dim_avg_Wrap\\|dim_cumsum_n_Wrap\\|dim_cumsum_Wrap\\|dim_max_n_Wrap\\|dim_maxind\\|dim_min_n_Wrap\\|dim_minind\\|dim_rmsd_n_Wrap\\|dim_rmsd_Wrap\\|dim_rmvmean_n_Wrap\\|dim_rmvmean_Wrap\\|dim_rmvmed_n_Wrap\\|dim_rmvmed_Wrap\\|dim_standardize_n_Wrap\\|dim_standardize_Wrap\\|dim_stddev_n_Wrap\\|dim_stddev_Wrap\\|dim_sum_n_Wrap\\|dim_sum_wgt_n_Wrap\\|dim_sum_wgt_Wrap\\|dim_sum_Wrap\\|dim_variance_n_Wrap\\|dim_variance_Wrap\\|dpres_plevel_Wrap\\|dtrend_leftdim\\|dv2uvF_Wrap\\|dv2uvG_Wrap\\|eady_growth_rate\\|eofcor_Wrap\\|eofcov_Wrap\\|eofunc_n_Wrap\\|eofunc_north\\|eofunc_ts_n_Wrap\\|eofunc_ts_Wrap\\|eofunc_varimax_reorder\\|eofunc_varimax_Wrap\\|eofunc_Wrap\\|epflux\\|epsZero\\|extract_globalatts_hdf5\\|f2fosh_Wrap\\|f2foshv_Wrap\\|f2fsh_Wrap\\|f2fshv_Wrap\\|f2gsh_Wrap\\|f2gshv_Wrap\\|fbindirSwap\\|fbinseqSwap1\\|fbinseqSwap2\\|fire_index_haines\\|flt2dble\\|flt2string\\|fo2fsh_Wrap\\|fo2fshv_Wrap\\|g2fsh_Wrap\\|g2fshv_Wrap\\|g2gsh_Wrap\\|g2gshv_Wrap\\|generate_sample_indices\\|generate_unique_indices\\|genNormalDist\\|get1Dindex_Collapse\\|get1Dindex_Exclude\\|get_d2r\\|get_file_suffix\\|get_pi\\|get_r2d\\|GetFillColor\\|GetFillColorIndex\\|getFillValue\\|getind_latlon2d\\|getVarDimNames\\|getVarFillValue\\|grad_latlon_cfd\\|grib_stime2itime\\|hyi2hyo_Wrap\\|ilapsF_Wrap\\|ilapsG_Wrap\\|ind_nearest_coord\\|indStrSubset\\|int2dble\\|int2flt\\|int2p_n_Wrap\\|int2p_Wrap\\|isMonotonic\\|isStrSubset\\|latent_heat_water\\|latGau\\|latGauWgt\\|latGlobeF\\|latGlobeFo\\|latRegWgt\\|linint1_n_Wrap\\|linint1_Wrap\\|linint2_points_Wrap\\|linint2_Wrap\\|local_max_1d\\|local_min_1d\\|lonFlip\\|lonGlobeF\\|lonGlobeFo\\|lonPivot\\|merge_levels_sfc\\|mod\\|month_to_annual\\|month_to_annual_weighted\\|month_to_season\\|month_to_season12\\|month_to_seasonN\\|monthly_total_to_daily_mean\\|nameDim\\|natgrid_Wrap\\|NewCosWeight\\|niceLatLon2D\\|NormCosWgtGlobe\\|numAsciiCol\\|numAsciiRow\\|numeric2int\\|obj_anal_ic_deprecated\\|obj_anal_ic_Wrap\\|omega_ccm_driver\\|omega_to_w\\|oneDtostring\\|pack_values\\|parse_globalatts_hdf5\\|pattern_cor\\|pdfx\\|pdfxy\\|pdfxy_conform\\|pot_temp\\|pot_temp_equiv\\|pot_vort_hybrid\\|pot_vort_isobaric\\|pres2hybrid_Wrap\\|print_clock\\|printMinMax\\|quadroots\\|rcm2points_Wrap\\|rcm2rgrid_Wrap\\|readAsciiHead\\|readAsciiTable\\|reg_multlin_stats\\|region_ind\\|regline_stats\\|relhum_ttd\\|replaceSingleChar\\|RGBtoCmap\\|rgrid2rcm_Wrap\\|rho_mwjf\\|rigrad_bruntv_atm\\|rm_single_dims\\|rmAnnCycle1D\\|\\)\\>" 1 font-lock-function-name-face)

    ("\\<\\(rmInsufData\\|rmMonAnnCycLLLT\\|rmMonAnnCycLLT\\|rmMonAnnCycTLL\\|runave_n_Wrap\\|runave_Wrap\\|satvpr_water_bolton\\|satvpr_water_stipanuk\\|short2flt\\|short2flt_hdf\\|shsgc_R42_Wrap\\|sign_f90\\|sign_matlab\\|smth9_Wrap\\|smthClmDayTLL\\|smthClmDayTLLL\\|SqrtCosWeight\\|stat_dispersion\\|static_stability\\|stdMonLLLT\\|stdMonLLT\\|stdMonTLL\\|stdMonTLLL\\|symMinMaxPlt\\|table_attach_columns\\|table_attach_rows\\|time_reassign\\|time_reassign_cv2var\\|time_to_newtime\\|time_to_newtime_fix\\|transpose\\|triple2grid_Wrap\\|ut_convert\\|ut_convert_fix\\|uv2dvF_Wrap\\|uv2dvG_Wrap\\|uv2vrF_Wrap\\|uv2vrG_Wrap\\|vapor_pres_rh\\|venn2_difference\\|venn2_intersection\\|venn2_union\\|vr2uvF_Wrap\\|vr2uvG_Wrap\\|w_to_omega\\|wallClockElapseTime\\|wave_number_spc\\|wetbulb_stull\\|wgt_areaave_Wrap\\|wgt_runave_leftdim\\|wgt_runave_n_Wrap\\|wgt_runave_Wrap\\|wgt_vertical_n\\|wind_component\\|wind_direction\\|wind_speed\\|wind_stats\\|yyyyddd_to_yyyymmdd\\|yyyymm_time\\|yyyymm_to_yyyyfrac\\|yyyymmdd_time\\|yyyymmdd_to_yyyyddd\\|yyyymmdd_to_yyyyfrac\\|yyyymmddhh_time\\|yyyymmddhh_to_yyyyfrac\\|zonal_mpsi_Wrap\\|zonalAve\\|\\)\\>" 1 font-lock-function-name-face)

    ;; crop functions
    ("\\<\\(actvpr_mnmx_fao56\\|actvpr_rhmean_fao56\\|crop_water_need\\|daylight_fao56\\|netlw_fao56\\|netrad_fao56\\|netsw_fao56\\|prsatm_tz_fao56\\|prsatm_z_fao56\\|psychro_fao56\\|radext_fao56\\|radsol2_fao56\\|radsol3_hargreaves_fao56\\|radsol_clrsky_fao56\\|radsol_fao56\\|refevt_hargreaves_fao56\\|refevt_penman_fao56\\|refevt_turc\\|refevt_turc_rh\\|rhum_fao56\\|satvpr_mean_fao56\\|satvpr_slope_fao56\\|satvpr_tdew_fao56\\|satvpr_temp_fao56\\|soil_heatflux_month_fao56\\|tdew_actvpr_fao56\\|u2_fao56\\|\\)\\>" 1 font-lock-function-name-face)

    ;; diagnostics functions
    ("\\<\\(band_pass_area_time\\|band_pass_area_time_plot\\|band_pass_hovmueller\\|band_pass_hovmueller_plot\\|band_pass_latlon_time\\|band_pass_latlon_time_plot\\|decomposeSymAsym\\|mjo_cross\\|mjo_cross_plot\\|mjo_phase_background\\|mjo_space_time_cross\\|mjo_spectra\\|mjo_spectra_season\\|mjo_wavenum_freq_season\\|mjo_wavenum_freq_season_plot\\|mjo_xcor_lag_ovly\\|mjo_xcor_lag_ovly_panel\\|mjo_xcor_lag_season\\|resolveWavesHayashi\\|wkSpaceTime\\|wkSpaceTime_cam\\|\\)\\>" 1 font-lock-function-name-face)

    ;; ESMF regridding
    ("\\<\\(curvilinear_to_SCRIP\\|ESMF_regrid\\|ESMF_regrid_gen_weights\\|ESMF_regrid_with_weights\\|latlon_to_SCRIP\\|rectilinear_to_SCRIP\\|unstructured_to_ESMF\\|\\)\\>" 1 font-lock-function-name-face)

    ;; ext_val functions
    ("\\<\\(extval_frechet\\|extval_gev\\|extval_gumbel\\|extval_mlegam\\|extval_pareto\\|extval_recurrence_table\\|extval_return_period\\|extval_return_prob\\|extval_weibull\\|\\)\\>" 1 font-lock-function-name-face)

    ;; heat_stress functions
    ("\\<\\(heat_apptemp\\|heat_discoi\\|heat_discoi_stull\\|heat_esidx_moran\\|heat_humidex\\|heat_index_nws\\|heat_swamp_cooleff\\|heat_thic_thip\\|heat_wbgt_inout\\|heat_wbgt_simplified\\|\\)\\>" 1 font-lock-function-name-face)

    ;; pop_remap functions
    ("\\<\\(PopLatLon\\|PopLatLonV\\|\\)\\>" 1 font-lock-function-name-face)

    ;; shea_util functions
    ("\\<\\(add90LatX\\|add90LatY\\|boxplot\\|ColorNegDashZeroPosContour\\|ColorShadeLeGeContour\\|drawNDCGrid\\|infoTimeStamp\\|landsea_mask\\|msgValOutline\\|pie_chart\\|plt_pdfxy\\|setColorContourClear\\|ShadeCOI\\|ShadeGeLeContour\\|ShadeGtContour\\|ShadeLtContour\\|ShadeLtGtContour\\|simple_legend\\|specx_ci\\|\\)\\>" 1 font-lock-function-name-face)

    ;; skewt functions
    ("\\<\\(skewT_BackGround\\|skewT_PlotData\\|\\)\\>" 1 font-lock-function-name-face)

    ;; user_contributed functions
    ("\\<\\(box_percentile_plot\\|calendar_decode2\\|calendar_decode2_fix\\|cd_inv_string\\|cd_string\\|kf_filter\\|run_cor\\|time_axis_labels\\|ut_string\\|ut_string_fix\\|\\)\\>" 1 font-lock-function-name-face)

    ;; wrf_arw functions
    ("\\<\\(wrf_contour\\|wrf_map\\|wrf_map_overlay\\|wrf_map_overlays\\|wrf_map_resources\\|wrf_map_zoom\\|wrf_overlay\\|wrf_overlays\\|wrf_user_getvar\\|wrf_user_ij_to_ll\\|wrf_user_intrp2d\\|wrf_user_intrp3d\\|wrf_user_latlon_to_ij\\|wrf_user_list_times\\|wrf_user_ll_to_ij\\|wrf_user_unstagger\\|wrf_user_vert_interp\\|wrf_vector\\|\\)\\>" 1 font-lock-function-name-face)

    ;; wrf_contributed functions
    ("\\<\\(wrf_mapres_c\\|wrf_times_c\\|\\)\\>" 1 font-lock-function-name-face)

    ;; wind_rose functions
    ("\\<\\(WindRoseBasic\\|WindRoseColor\\|WindRoseThickLine\\|\\)\\>" 1 font-lock-function-name-face)

    ;; gsn csm plot templates and special gsn functions
    ("\\<\\(gsn_add_annotation\\|gsn_add_polygon\\|gsn_add_polyline\\|gsn_add_polymarker\\|gsn_add_shapefile_polygons\\|gsn_add_shapefile_polylines\\|gsn_add_shapefile_polymarkers\\|gsn_add_text\\|gsn_attach_plots\\|gsn_blank_plot\\|gsn_contour\\|gsn_contour_map\\|gsn_contour_shade\\|gsn_coordinates\\|gsn_create_labelbar\\|gsn_create_legend\\|gsn_create_text\\|gsn_csm_attach_zonal_means\\|gsn_csm_blank_plot\\|gsn_csm_contour\\|gsn_csm_contour_map\\|gsn_csm_contour_map_ce\\|gsn_csm_contour_map_overlay\\|gsn_csm_contour_map_polar\\|gsn_csm_hov\\|gsn_csm_lat_time\\|gsn_csm_map\\|gsn_csm_map_ce\\|gsn_csm_map_polar\\|gsn_csm_pres_hgt\\|gsn_csm_pres_hgt_streamline\\|gsn_csm_pres_hgt_vector\\|gsn_csm_streamline\\|gsn_csm_streamline_contour_map\\|gsn_csm_streamline_contour_map_ce\\|gsn_csm_streamline_contour_map_polar\\|gsn_csm_streamline_map\\|gsn_csm_streamline_map_ce\\|gsn_csm_streamline_map_polar\\|gsn_csm_streamline_scalar\\|gsn_csm_streamline_scalar_map\\|gsn_csm_streamline_scalar_map_ce\\|gsn_csm_streamline_scalar_map_polar\\|gsn_csm_time_lat\\|gsn_csm_vector\\|gsn_csm_vector_map\\|gsn_csm_vector_map_ce\\|gsn_csm_vector_map_polar\\|gsn_csm_vector_scalar\\|gsn_csm_vector_scalar_map\\|gsn_csm_vector_scalar_map_ce\\|gsn_csm_vector_scalar_map_polar\\|gsn_csm_x2y\\|gsn_csm_x2y2\\|gsn_csm_xy\\|gsn_csm_xy2\\|gsn_csm_xy3\\|gsn_csm_y\\|gsn_define_colormap\\|gsn_draw_colormap\\|gsn_draw_named_colors\\|gsn_histogram\\|gsn_labelbar_ndc\\|gsn_legend_ndc\\|gsn_map\\|gsn_merge_colormaps\\|gsn_open_wks\\|gsn_panel\\|gsn_polygon\\|gsn_polygon_ndc\\|gsn_polyline\\|gsn_polyline_ndc\\|gsn_polymarker\\|gsn_polymarker_ndc\\|gsn_retrieve_colormap\\|gsn_reverse_colormap\\|gsn_streamline\\|gsn_streamline_map\\|gsn_streamline_scalar\\|gsn_streamline_scalar_map\\|gsn_table\\|gsn_text\\|gsn_text_ndc\\|gsn_vector\\|gsn_vector_map\\|gsn_vector_scalar\\|gsn_vector_scalar_map\\|gsn_xy\\|gsn_y\\|hsv2rgb\\|maximize_output\\|reset_device_coordinates\\|\\)\\>" 1 font-lock-variable-name-face)

    ;; ncl resources (the list is split in several lines to avoid too long regular expressions)
    ("\\<\\(amDataXF\\|amDataYF\\|amJust\\|amOn\\|amOrthogonalPosF\\|amParallelPosF\\|amResizeNotify\\|amSide\\|amTrackData\\|amViewId\\|amZone\\|appDefaultParent\\|appFileSuffix\\|appResources\\|appSysDir\\|appUsrDir\\|caCopyArrays\\|caXArray\\|caXCast\\|caXMaxV\\|caXMinV\\|caXMissingV\\|caYArray\\|caYCast\\|caYMaxV\\|caYMinV\\|caYMissingV\\|cnCellFillEdgeColor\\|cnCellFillMissingValEdgeColor\\|cnConpackParams\\|cnConstFEnableFill\\|cnConstFLabelAngleF\\|cnConstFLabelBackgroundColor\\|cnConstFLabelConstantSpacingF\\|cnConstFLabelFont\\|cnConstFLabelFontAspectF\\|cnConstFLabelFontColor\\|cnConstFLabelFontHeightF\\|cnConstFLabelFontQuality\\|cnConstFLabelFontThicknessF\\|cnConstFLabelFormat\\|cnConstFLabelFuncCode\\|cnConstFLabelJust\\|cnConstFLabelOn\\|cnConstFLabelOrthogonalPosF\\|cnConstFLabelParallelPosF\\|cnConstFLabelPerimColor\\|cnConstFLabelPerimOn\\|cnConstFLabelPerimSpaceF\\|cnConstFLabelPerimThicknessF\\|cnConstFLabelSide\\|cnConstFLabelString\\|cnConstFLabelTextDirection\\|cnConstFLabelZone\\|cnConstFUseInfoLabelRes\\|cnExplicitLabelBarLabelsOn\\|cnExplicitLegendLabelsOn\\|cnExplicitLineLabelsOn\\|cnFillBackgroundColor\\|cnFillColor\\|cnFillColors\\|cnFillDotSizeF\\|cnFillDrawOrder\\|cnFillMode\\|cnFillOn\\|cnFillOpacityF\\|cnFillPalette\\|cnFillPattern\\|cnFillPatterns\\|cnFillScaleF\\|cnFillScales\\|cnFixFillBleed\\|cnGridBoundFillColor\\|cnGridBoundFillPattern\\|cnGridBoundFillScaleF\\|cnGridBoundPerimColor\\|cnGridBoundPerimDashPattern\\|cnGridBoundPerimOn\\|cnGridBoundPerimThicknessF\\|cnHighLabelAngleF\\|cnHighLabelBackgroundColor\\|cnHighLabelConstantSpacingF\\|cnHighLabelCount\\|cnHighLabelFont\\|cnHighLabelFontAspectF\\|cnHighLabelFontColor\\|cnHighLabelFontHeightF\\|cnHighLabelFontQuality\\|cnHighLabelFontThicknessF\\|cnHighLabelFormat\\|cnHighLabelFuncCode\\|cnHighLabelPerimColor\\|cnHighLabelPerimOn\\|cnHighLabelPerimSpaceF\\|cnHighLabelPerimThicknessF\\|cnHighLabelString\\|cnHighLabelsOn\\|cnHighLowLabelOverlapMode\\|cnHighUseLineLabelRes\\|cnInfoLabelAngleF\\|cnInfoLabelBackgroundColor\\|cnInfoLabelConstantSpacingF\\|cnInfoLabelFont\\|cnInfoLabelFontAspectF\\|cnInfoLabelFontColor\\|cnInfoLabelFontHeightF\\|cnInfoLabelFontQuality\\|cnInfoLabelFontThicknessF\\|cnInfoLabelFormat\\|cnInfoLabelFuncCode\\|cnInfoLabelJust\\|cnInfoLabelOn\\|cnInfoLabelOrthogonalPosF\\|cnInfoLabelParallelPosF\\|cnInfoLabelPerimColor\\|cnInfoLabelPerimOn\\|cnInfoLabelPerimSpaceF\\|cnInfoLabelPerimThicknessF\\|cnInfoLabelSide\\|cnInfoLabelString\\|cnInfoLabelTextDirection\\|cnInfoLabelZone\\|cnLabelBarEndLabelsOn\\|cnLabelBarEndStyle\\|cnLabelDrawOrder\\|cnLabelMasking\\|cnLabelScaleFactorF\\|cnLabelScaleValueF\\|cnLabelScalingMode\\|cnLegendLevelFlags\\|cnLevelCount\\|cnLevelFlag\\|cnLevelFlags\\|cnLevelSelectionMode\\|cnLevelSpacingF\\|cnLevels\\|cnLineColor\\|cnLineColors\\|cnLineDashPattern\\|cnLineDashPatterns\\|cnLineDashSegLenF\\|cnLineDrawOrder\\|cnLineLabelAngleF\\|cnLineLabelBackgroundColor\\|cnLineLabelConstantSpacingF\\|cnLineLabelCount\\|cnLineLabelDensityF\\|cnLineLabelFont\\|cnLineLabelFontAspectF\\|cnLineLabelFontColor\\|cnLineLabelFontColors\\|cnLineLabelFontHeightF\\|cnLineLabelFontQuality\\|cnLineLabelFontThicknessF\\|cnLineLabelFormat\\|cnLineLabelFuncCode\\|cnLineLabelInterval\\|cnLineLabelPerimColor\\|cnLineLabelPerimOn\\|cnLineLabelPerimSpaceF\\|cnLineLabelPerimThicknessF\\|cnLineLabelPlacementMode\\|cnLineLabelStrings\\|cnLineLabelsOn\\|cnLinePalette\\|cnLineThicknessF\\|cnLineThicknesses\\|cnLinesOn\\|cnLowLabelAngleF\\|cnLowLabelBackgroundColor\\|cnLowLabelConstantSpacingF\\|cnLowLabelCount\\|cnLowLabelFont\\|cnLowLabelFontAspectF\\|cnLowLabelFontColor\\|cnLowLabelFontHeightF\\|cnLowLabelFontQuality\\|cnLowLabelFontThicknessF\\|cnLowLabelFormat\\|cnLowLabelFuncCode\\|cnLowLabelPerimColor\\|cnLowLabelPerimOn\\|cnLowLabelPerimSpaceF\\|cnLowLabelPerimThicknessF\\|cnLowLabelString\\|cnLowLabelsOn\\|cnLowUseHighLabelRes\\|cnMaxDataValueFormat\\|cnMaxLevelCount\\|cnMaxLevelValF\\|cnMaxPointDistanceF\\|cnMinLevelValF\\|cnMissingValFillColor\\|cnMissingValFillPattern\\|cnMissingValFillScaleF\\|cnMissingValPerimColor\\|cnMissingValPerimDashPattern\\|cnMissingValPerimGridBoundOn\\|cnMissingValPerimOn\\|cnMissingValPerimThicknessF\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(cnMonoFillColor\\|cnMonoFillPattern\\|cnMonoFillScale\\|cnMonoLevelFlag\\|cnMonoLineColor\\|cnMonoLineDashPattern\\|cnMonoLineLabelFontColor\\|cnMonoLineThickness\\|cnNoDataLabelOn\\|cnNoDataLabelString\\|cnOutOfRangeFillColor\\|cnOutOfRangeFillPattern\\|cnOutOfRangeFillScaleF\\|cnOutOfRangePerimColor\\|cnOutOfRangePerimDashPattern\\|cnOutOfRangePerimOn\\|cnOutOfRangePerimThicknessF\\|cnRasterCellSizeF\\|cnRasterMinCellSizeF\\|cnRasterModeOn\\|cnRasterSampleFactorF\\|cnRasterSmoothingOn\\|cnScalarFieldData\\|cnSmoothingDistanceF\\|cnSmoothingOn\\|cnSmoothingTensionF\\|cnSpanFillPalette\\|cnSpanLinePalette\\|ctCopyTables\\|ctXElementSize\\|ctXMaxV\\|ctXMinV\\|ctXMissingV\\|ctXTable\\|ctXTableLengths\\|ctXTableType\\|ctYElementSize\\|ctYMaxV\\|ctYMinV\\|ctYMissingV\\|ctYTable\\|ctYTableLengths\\|ctYTableType\\|dcDelayCompute\\|errBuffer\\|errFileName\\|errFilePtr\\|errLevel\\|errPrint\\|errUnitNumber\\|gsClipOn\\|gsColors\\|gsEdgeColor\\|gsEdgeDashPattern\\|gsEdgeDashSegLenF\\|gsEdgeThicknessF\\|gsEdgesOn\\|gsFillBackgroundColor\\|gsFillColor\\|gsFillDotSizeF\\|gsFillIndex\\|gsFillLineThicknessF\\|gsFillOpacityF\\|gsFillScaleF\\|gsFont\\|gsFontAspectF\\|gsFontColor\\|gsFontHeightF\\|gsFontOpacityF\\|gsFontQuality\\|gsFontThicknessF\\|gsLineColor\\|gsLineDashPattern\\|gsLineDashSegLenF\\|gsLineLabelConstantSpacingF\\|gsLineLabelFont\\|gsLineLabelFontAspectF\\|gsLineLabelFontColor\\|gsLineLabelFontHeightF\\|gsLineLabelFontQuality\\|gsLineLabelFontThicknessF\\|gsLineLabelFuncCode\\|gsLineLabelString\\|gsLineOpacityF\\|gsLineThicknessF\\|gsMarkerColor\\|gsMarkerIndex\\|gsMarkerOpacityF\\|gsMarkerSizeF\\|gsMarkerThicknessF\\|gsSegments\\|gsTextAngleF\\|gsTextConstantSpacingF\\|gsTextDirection\\|gsTextFuncCode\\|gsTextJustification\\|gsnAboveYRefLineBarColors\\|gsnAboveYRefLineBarFillScales\\|gsnAboveYRefLineBarPatterns\\|gsnAboveYRefLineColor\\|gsnAddCyclic\\|gsnAttachBorderOn\\|gsnAttachPlotsXAxis\\|gsnBelowYRefLineBarColors\\|gsnBelowYRefLineBarFillScales\\|gsnBelowYRefLineBarPatterns\\|gsnBelowYRefLineColor\\|gsnBoxMargin\\|gsnCenterString\\|gsnCenterStringFontColor\\|gsnCenterStringFontHeightF\\|gsnCenterStringFuncCode\\|gsnCenterStringOrthogonalPosF\\|gsnCenterStringParallelPosF\\|gsnContourLineThicknessesScale\\|gsnContourNegLineDashPattern\\|gsnContourPosLineDashPattern\\|gsnContourZeroLineThicknessF\\|gsnDebugWriteFileName\\|gsnDraw\\|gsnFrame\\|gsnHistogramBarColors\\|gsnHistogramBarWidthPercent\\|gsnHistogramBinIntervals\\|gsnHistogramBinMissing\\|gsnHistogramBinWidth\\|gsnHistogramClassIntervals\\|gsnHistogramCompare\\|gsnHistogramComputePercentages\\|gsnHistogramComputePercentagesNoMissing\\|gsnHistogramDiscreteBinValues\\|gsnHistogramDiscreteClassValues\\|gsnHistogramHorizontal\\|gsnHistogramMinMaxBinsOn\\|gsnHistogramNumberOfBins\\|gsnHistogramPercentSign\\|gsnHistogramSelectNiceIntervals\\|gsnLeftString\\|gsnLeftStringFontColor\\|gsnLeftStringFontHeightF\\|gsnLeftStringFuncCode\\|gsnLeftStringOrthogonalPosF\\|gsnLeftStringParallelPosF\\|gsnLeftXRefLineBarColors\\|gsnLeftXRefLineBarFillScales\\|gsnLeftXRefLineBarPatterns\\|gsnLeftXRefLineColor\\|gsnMajorLatSpacing\\|gsnMajorLonSpacing\\|gsnMaskLambertConformal\\|gsnMaskLambertConformalOutlineOn\\|gsnMaximize\\|gsnMinorLatSpacing\\|gsnMinorLonSpacing\\|gsnPanelBottom\\|gsnPanelCenter\\|gsnPanelDebug\\|gsnPanelFigureStrings\\|gsnPanelFigureStringsBackgroundFillColor\\|gsnPanelFigureStringsFontHeightF\\|gsnPanelFigureStringsJust\\|gsnPanelFigureStringsPerimOn\\|gsnPanelLabelBar\\|gsnPanelLeft\\|gsnPanelMainFont\\|gsnPanelMainFontColor\\|gsnPanelMainFontHeightF\\|gsnPanelMainPosXF\\|gsnPanelMainPosYF\\|gsnPanelMainString\\|gsnPanelRight\\|gsnPanelRowSpec\\|gsnPanelScalePlotIndex\\|gsnPanelTop\\|gsnPanelXF\\|gsnPanelXWhiteSpacePercent\\|gsnPanelYF\\|gsnPanelYWhiteSpacePercent\\|gsnPaperHeight\\|gsnPaperMargin\\|gsnPaperOrientation\\|gsnPaperWidth\\|gsnPolar\\|gsnPolarLabelDistance\\|gsnPolarLabelFont\\|gsnPolarLabelFontHeightF\\|gsnPolarLabelSpacing\\|gsnPolarTime\\|gsnPolarUT\\|gsnRightString\\|gsnRightStringFontColor\\|gsnRightStringFontHeightF\\|gsnRightStringFuncCode\\|gsnRightStringOrthogonalPosF\\|gsnRightStringParallelPosF\\|gsnRightXRefLineBarColors\\|gsnRightXRefLineBarFillScales\\|gsnRightXRefLineBarPatterns\\|gsnRightXRefLineColor\\|gsnScalarContour\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(gsnScale\\|gsnShape\\|gsnSpreadColorEnd\\|gsnSpreadColorStart\\|gsnSpreadColors\\|gsnStringFont\\|gsnStringFontColor\\|gsnStringFontHeightF\\|gsnStringFuncCode\\|gsnTickMarksOn\\|gsnXAxisIrregular2Linear\\|gsnXAxisIrregular2Log\\|gsnXRefLine\\|gsnXRefLineColor\\|gsnXRefLineColors\\|gsnXRefLineDashPattern\\|gsnXRefLineDashPatterns\\|gsnXRefLineThicknessF\\|gsnXRefLineThicknesses\\|gsnXYAboveFillColors\\|gsnXYBarChart\\|gsnXYBarChartBarWidth\\|gsnXYBarChartColors\\|gsnXYBarChartColors2\\|gsnXYBarChartFillDotSizeF\\|gsnXYBarChartFillLineThicknessF\\|gsnXYBarChartFillOpacityF\\|gsnXYBarChartFillScaleF\\|gsnXYBarChartOutlineOnly\\|gsnXYBarChartOutlineThicknessF\\|gsnXYBarChartPatterns\\|gsnXYBarChartPatterns2\\|gsnXYBelowFillColors\\|gsnXYFillColors\\|gsnXYFillOpacities\\|gsnXYLeftFillColors\\|gsnXYRightFillColors\\|gsnYAxisIrregular2Linear\\|gsnYAxisIrregular2Log\\|gsnYRefLine\\|gsnYRefLineColor\\|gsnYRefLineColors\\|gsnYRefLineDashPattern\\|gsnYRefLineDashPatterns\\|gsnYRefLineThicknessF\\|gsnYRefLineThicknesses\\|gsnZonalMean\\|gsnZonalMeanXMaxF\\|gsnZonalMeanXMinF\\|gsnZonalMeanYRefLine\\|lbAutoManage\\|lbBottomMarginF\\|lbBoxCount\\|lbBoxEndCapStyle\\|lbBoxFractions\\|lbBoxLineColor\\|lbBoxLineDashPattern\\|lbBoxLineDashSegLenF\\|lbBoxLineThicknessF\\|lbBoxLinesOn\\|lbBoxMajorExtentF\\|lbBoxMinorExtentF\\|lbBoxSeparatorLinesOn\\|lbBoxSizing\\|lbFillBackground\\|lbFillColor\\|lbFillColors\\|lbFillDotSizeF\\|lbFillLineThicknessF\\|lbFillOpacityF\\|lbFillPattern\\|lbFillPatterns\\|lbFillScaleF\\|lbFillScales\\|lbJustification\\|lbLabelAlignment\\|lbLabelAngleF\\|lbLabelAutoStride\\|lbLabelBarOn\\|lbLabelConstantSpacingF\\|lbLabelDirection\\|lbLabelFont\\|lbLabelFontAspectF\\|lbLabelFontColor\\|lbLabelFontHeightF\\|lbLabelFontQuality\\|lbLabelFontThicknessF\\|lbLabelFuncCode\\|lbLabelJust\\|lbLabelOffsetF\\|lbLabelPosition\\|lbLabelStride\\|lbLabelStrings\\|lbLabelsOn\\|lbLeftMarginF\\|lbMaxLabelLenF\\|lbMinLabelSpacingF\\|lbMonoFillColor\\|lbMonoFillPattern\\|lbMonoFillScale\\|lbOrientation\\|lbOverrideFillOpacity\\|lbPerimColor\\|lbPerimDashPattern\\|lbPerimDashSegLenF\\|lbPerimFill\\|lbPerimFillColor\\|lbPerimOn\\|lbPerimThicknessF\\|lbRasterFillOn\\|lbRightMarginF\\|lbTitleAngleF\\|lbTitleConstantSpacingF\\|lbTitleDirection\\|lbTitleExtentF\\|lbTitleFont\\|lbTitleFontAspectF\\|lbTitleFontColor\\|lbTitleFontHeightF\\|lbTitleFontQuality\\|lbTitleFontThicknessF\\|lbTitleFuncCode\\|lbTitleJust\\|lbTitleOffsetF\\|lbTitleOn\\|lbTitlePosition\\|lbTitleString\\|lbTopMarginF\\|lgAutoManage\\|lgBottomMarginF\\|lgBoxBackground\\|lgBoxLineColor\\|lgBoxLineDashPattern\\|lgBoxLineDashSegLenF\\|lgBoxLineThicknessF\\|lgBoxLinesOn\\|lgBoxMajorExtentF\\|lgBoxMinorExtentF\\|lgDashIndex\\|lgDashIndexes\\|lgItemCount\\|lgItemOrder\\|lgItemPlacement\\|lgItemPositions\\|lgItemType\\|lgItemTypes\\|lgJustification\\|lgLabelAlignment\\|lgLabelAngleF\\|lgLabelAutoStride\\|lgLabelConstantSpacingF\\|lgLabelDirection\\|lgLabelFont\\|lgLabelFontAspectF\\|lgLabelFontColor\\|lgLabelFontHeightF\\|lgLabelFontQuality\\|lgLabelFontThicknessF\\|lgLabelFuncCode\\|lgLabelJust\\|lgLabelOffsetF\\|lgLabelPosition\\|lgLabelStride\\|lgLabelStrings\\|lgLabelsOn\\|lgLeftMarginF\\|lgLegendOn\\|lgLineColor\\|lgLineColors\\|lgLineDashSegLenF\\|lgLineDashSegLens\\|lgLineLabelConstantSpacingF\\|lgLineLabelFont\\|lgLineLabelFontAspectF\\|lgLineLabelFontColor\\|lgLineLabelFontColors\\|lgLineLabelFontHeightF\\|lgLineLabelFontHeights\\|lgLineLabelFontQuality\\|lgLineLabelFontThicknessF\\|lgLineLabelFuncCode\\|lgLineLabelStrings\\|lgLineLabelsOn\\|lgLineThicknessF\\|lgLineThicknesses\\|lgMarkerColor\\|lgMarkerColors\\|lgMarkerIndex\\|lgMarkerIndexes\\|lgMarkerSizeF\\|lgMarkerSizes\\|lgMarkerThicknessF\\|lgMarkerThicknesses\\|lgMonoDashIndex\\|lgMonoItemType\\|lgMonoLineColor\\|lgMonoLineDashSegLen\\|lgMonoLineLabelFontColor\\|lgMonoLineLabelFontHeight\\|lgMonoLineThickness\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(lgMonoMarkerColor\\|lgMonoMarkerIndex\\|lgMonoMarkerSize\\|lgMonoMarkerThickness\\|lgOrientation\\|lgPerimColor\\|lgPerimDashPattern\\|lgPerimDashSegLenF\\|lgPerimFill\\|lgPerimFillColor\\|lgPerimOn\\|lgPerimThicknessF\\|lgRightMarginF\\|lgTitleAngleF\\|lgTitleConstantSpacingF\\|lgTitleDirection\\|lgTitleExtentF\\|lgTitleFont\\|lgTitleFontAspectF\\|lgTitleFontColor\\|lgTitleFontHeightF\\|lgTitleFontQuality\\|lgTitleFontThicknessF\\|lgTitleFuncCode\\|lgTitleJust\\|lgTitleOffsetF\\|lgTitleOn\\|lgTitlePosition\\|lgTitleString\\|lgTopMarginF\\|mpAreaGroupCount\\|mpAreaMaskingOn\\|mpAreaNames\\|mpAreaTypes\\|mpBottomAngleF\\|mpBottomMapPosF\\|mpBottomNDCF\\|mpBottomNPCF\\|mpBottomPointLatF\\|mpBottomPointLonF\\|mpBottomWindowF\\|mpCenterLatF\\|mpCenterLonF\\|mpCenterRotF\\|mpCountyLineColor\\|mpCountyLineDashPattern\\|mpCountyLineDashSegLenF\\|mpCountyLineThicknessF\\|mpDataBaseVersion\\|mpDataResolution\\|mpDataSetName\\|mpDefaultFillColor\\|mpDefaultFillPattern\\|mpDefaultFillScaleF\\|mpDynamicAreaGroups\\|mpEllipticalBoundary\\|mpFillAreaSpecifiers\\|mpFillBoundarySets\\|mpFillColor\\|mpFillColors\\|mpFillDotSizeF\\|mpFillDrawOrder\\|mpFillOn\\|mpFillPatternBackground\\|mpFillPattern\\|mpFillPatterns\\|mpFillScaleF\\|mpFillScales\\|mpFixedAreaGroups\\|mpGeophysicalLineColor\\|mpGeophysicalLineDashPattern\\|mpGeophysicalLineDashSegLenF\\|mpGeophysicalLineThicknessF\\|mpGreatCircleLinesOn\\|mpGridAndLimbDrawOrder\\|mpGridAndLimbOn\\|mpGridLatSpacingF\\|mpGridLineColor\\|mpGridLineDashPattern\\|mpGridLineDashSegLenF\\|mpGridLineThicknessF\\|mpGridLonSpacingF\\|mpGridMaskMode\\|mpGridMaxLatF\\|mpGridPolarLonSpacingF\\|mpGridSpacingF\\|mpInlandWaterFillColor\\|mpInlandWaterFillPattern\\|mpInlandWaterFillScaleF\\|mpLabelDrawOrder\\|mpLabelFontColor\\|mpLabelFontHeightF\\|mpLabelsOn\\|mpLambertMeridianF\\|mpLambertParallel1F\\|mpLambertParallel2F\\|mpLandFillColor\\|mpLandFillPattern\\|mpLandFillScaleF\\|mpLeftAngleF\\|mpLeftCornerLatF\\|mpLeftCornerLonF\\|mpLeftMapPosF\\|mpLeftNDCF\\|mpLeftNPCF\\|mpLeftPointLatF\\|mpLeftPointLonF\\|mpLeftWindowF\\|mpLimbLineColor\\|mpLimbLineDashPattern\\|mpLimbLineDashSegLenF\\|mpLimbLineThicknessF\\|mpLimitMode\\|mpMaskAreaSpecifiers\\|mpMaskOutlineSpecifiers\\|mpMaxLatF\\|mpMaxLonF\\|mpMinLatF\\|mpMinLonF\\|mpMonoFillColor\\|mpMonoFillPattern\\|mpMonoFillScale\\|mpNationalLineColor\\|mpNationalLineDashPattern\\|mpNationalLineDashSegLenF\\|mpNationalLineThicknessF\\|mpOceanFillColor\\|mpOceanFillPattern\\|mpOceanFillScaleF\\|mpOutlineBoundarySets\\|mpOutlineDrawOrder\\|mpOutlineMaskingOn\\|mpOutlineOn\\|mpOutlineSpecifiers\\|mpPerimDrawOrder\\|mpPerimLineColor\\|mpPerimLineDashPattern\\|mpPerimLineDashSegLenF\\|mpPerimLineThicknessF\\|mpPerimOn\\|mpPolyMode\\|mpProjection\\|mpProvincialLineColor\\|mpProvincialLineDashPattern\\|mpProvincialLineDashSegLenF\\|mpProvincialLineThicknessF\\|mpRelativeCenterLat\\|mpRelativeCenterLon\\|mpRightAngleF\\|mpRightCornerLatF\\|mpRightCornerLonF\\|mpRightMapPosF\\|mpRightNDCF\\|mpRightNPCF\\|mpRightPointLatF\\|mpRightPointLonF\\|mpRightWindowF\\|mpSatelliteAngle1F\\|mpSatelliteAngle2F\\|mpSatelliteDistF\\|mpShapeMode\\|mpSpecifiedFillColors\\|mpSpecifiedFillDirectIndexing\\|mpSpecifiedFillPatterns\\|mpSpecifiedFillPriority\\|mpSpecifiedFillScales\\|mpTopAngleF\\|mpTopMapPosF\\|mpTopNDCF\\|mpTopNPCF\\|mpTopPointLatF\\|mpTopPointLonF\\|mpTopWindowF\\|mpUSStateLineColor\\|mpUSStateLineDashPattern\\|mpUSStateLineDashSegLenF\\|mpUSStateLineThicknessF\\|pmAnnoManagers\\|pmAnnoViews\\|pmLabelBarDisplayMode\\|pmLabelBarHeightF\\|pmLabelBarKeepAspect\\|pmLabelBarOrthogonalPosF\\|pmLabelBarParallelPosF\\|pmLabelBarSide\\|pmLabelBarWidthF\\|pmLabelBarZone\\|pmLegendDisplayMode\\|pmLegendHeightF\\|pmLegendKeepAspect\\|pmLegendOrthogonalPosF\\|pmLegendParallelPosF\\|pmLegendSide\\|pmLegendWidthF\\|pmLegendZone\\|pmOverlaySequenceIds\\|pmTickMarkDisplayMode\\|pmTickMarkZone\\|pmTitleDisplayMode\\|pmTitleZone\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(prGraphicStyle\\|prPolyType\\|prXArray\\|prYArray\\|sfCopyData\\|sfCopyData\\|sfDataArray\\|sfDataArray\\|sfDataMaxV\\|sfDataMaxV\\|sfDataMinV\\|sfDataMinV\\|sfElementNodes\\|sfExchangeDimensions\\|sfFirstNodeIndex\\|sfMissingValueV\\|sfMissingValueV\\|sfXArray\\|sfXArray\\|sfXCActualEndF\\|sfXCActualEndF\\|sfXCActualStartF\\|sfXCActualStartF\\|sfXCEndIndex\\|sfXCEndSubsetV\\|sfXCEndV\\|sfXCStartIndex\\|sfXCStartSubsetV\\|sfXCStartV\\|sfXCStride\\|sfXCellBounds\\|sfYArray\\|sfYArray\\|sfYCActualEndF\\|sfYCActualEndF\\|sfYCActualStartF\\|sfYCActualStartF\\|sfYCEndIndex\\|sfYCEndSubsetV\\|sfYCEndV\\|sfYCStartIndex\\|sfYCStartSubsetV\\|sfYCStartV\\|sfYCStride\\|sfYCellBounds\\|stArrowLengthF\\|stArrowStride\\|stCrossoverCheckCount\\|stExplicitLabelBarLabelsOn\\|stLabelBarEndLabelsOn\\|stLabelFormat\\|stLengthCheckCount\\|stLevelColors\\|stLevelCount\\|stLevelPalette\\|stLevelSelectionMode\\|stLevelSpacingF\\|stLevels\\|stLineColor\\|stLineOpacityF\\|stLineStartStride\\|stLineThicknessF\\|stMapDirection\\|stMaxLevelCount\\|stMaxLevelValF\\|stMinArrowSpacingF\\|stMinDistanceF\\|stMinLevelValF\\|stMinLineSpacingF\\|stMinStepFactorF\\|stMonoLineColor\\|stNoDataLabelOn\\|stNoDataLabelString\\|stScalarFieldData\\|stScalarMissingValColor\\|stSpanLevelPalette\\|stStepSizeF\\|stStreamlineDrawOrder\\|stUseScalarArray\\|stVectorFieldData\\|stZeroFLabelAngleF\\|stZeroFLabelBackgroundColor\\|stZeroFLabelConstantSpacingF\\|stZeroFLabelFont\\|stZeroFLabelFontAspectF\\|stZeroFLabelFontColor\\|stZeroFLabelFontHeightF\\|stZeroFLabelFontQuality\\|stZeroFLabelFontThicknessF\\|stZeroFLabelFuncCode\\|stZeroFLabelJust\\|stZeroFLabelOn\\|stZeroFLabelOrthogonalPosF\\|stZeroFLabelParallelPosF\\|stZeroFLabelPerimColor\\|stZeroFLabelPerimOn\\|stZeroFLabelPerimSpaceF\\|stZeroFLabelPerimThicknessF\\|stZeroFLabelSide\\|stZeroFLabelString\\|stZeroFLabelTextDirection\\|stZeroFLabelZone\\|tfDoNDCOverlay\\|tfPlotManagerOn\\|tfPolyDrawList\\|tfPolyDrawOrder\\|tiDeltaF\\|tiMainAngleF\\|tiMainConstantSpacingF\\|tiMainDirection\\|tiMainFont\\|tiMainFontAspectF\\|tiMainFontColor\\|tiMainFontHeightF\\|tiMainFontQuality\\|tiMainFontThicknessF\\|tiMainFuncCode\\|tiMainJust\\|tiMainOffsetXF\\|tiMainOffsetYF\\|tiMainOn\\|tiMainPosition\\|tiMainSide\\|tiMainString\\|tiUseMainAttributes\\|tiXAxisAngleF\\|tiXAxisConstantSpacingF\\|tiXAxisDirection\\|tiXAxisFont\\|tiXAxisFontAspectF\\|tiXAxisFontColor\\|tiXAxisFontHeightF\\|tiXAxisFontQuality\\|tiXAxisFontThicknessF\\|tiXAxisFuncCode\\|tiXAxisJust\\|tiXAxisOffsetXF\\|tiXAxisOffsetYF\\|tiXAxisOn\\|tiXAxisPosition\\|tiXAxisSide\\|tiXAxisString\\|tiYAxisAngleF\\|tiYAxisConstantSpacingF\\|tiYAxisDirection\\|tiYAxisFont\\|tiYAxisFontAspectF\\|tiYAxisFontColor\\|tiYAxisFontHeightF\\|tiYAxisFontQuality\\|tiYAxisFontThicknessF\\|tiYAxisFuncCode\\|tiYAxisJust\\|tiYAxisOffsetXF\\|tiYAxisOffsetYF\\|tiYAxisOn\\|tiYAxisPosition\\|tiYAxisSide\\|tiYAxisString\\|tmBorderLineColor\\|tmBorderThicknessF\\|tmEqualizeXYSizes\\|tmLabelAutoStride\\|tmSciNoteCutoff\\|tmXBAutoPrecision\\|tmXBBorderOn\\|tmXBDataLeftF\\|tmXBDataRightF\\|tmXBFormat\\|tmXBIrrTensionF\\|tmXBIrregularPoints\\|tmXBLabelAngleF\\|tmXBLabelConstantSpacingF\\|tmXBLabelDeltaF\\|tmXBLabelDirection\\|tmXBLabelFont\\|tmXBLabelFontAspectF\\|tmXBLabelFontColor\\|tmXBLabelFontHeightF\\|tmXBLabelFontQuality\\|tmXBLabelFontThicknessF\\|tmXBLabelFuncCode\\|tmXBLabelJust\\|tmXBLabelStride\\|tmXBLabels\\|tmXBLabelsOn\\|tmXBMajorLengthF\\|tmXBMajorLineColor\\|tmXBMajorOutwardLengthF\\|tmXBMajorThicknessF\\|tmXBMaxLabelLenF\\|tmXBMaxTicks\\|tmXBMinLabelSpacingF\\|tmXBMinorLengthF\\|tmXBMinorLineColor\\|tmXBMinorOn\\|tmXBMinorOutwardLengthF\\|tmXBMinorPerMajor\\|tmXBMinorThicknessF\\|tmXBMinorValues\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(tmXBMode\\|tmXBOn\\|tmXBPrecision\\|tmXBStyle\\|tmXBTickEndF\\|tmXBTickSpacingF\\|tmXBTickStartF\\|tmXBValues\\|tmXMajorGrid\\|tmXMajorGridLineColor\\|tmXMajorGridLineDashPattern\\|tmXMajorGridThicknessF\\|tmXMinorGrid\\|tmXMinorGridLineColor\\|tmXMinorGridLineDashPattern\\|tmXMinorGridThicknessF\\|tmXTAutoPrecision\\|tmXTBorderOn\\|tmXTDataLeftF\\|tmXTDataRightF\\|tmXTFormat\\|tmXTIrrTensionF\\|tmXTIrregularPoints\\|tmXTLabelAngleF\\|tmXTLabelConstantSpacingF\\|tmXTLabelDeltaF\\|tmXTLabelDirection\\|tmXTLabelFont\\|tmXTLabelFontAspectF\\|tmXTLabelFontColor\\|tmXTLabelFontHeightF\\|tmXTLabelFontQuality\\|tmXTLabelFontThicknessF\\|tmXTLabelFuncCode\\|tmXTLabelJust\\|tmXTLabelStride\\|tmXTLabels\\|tmXTLabelsOn\\|tmXTMajorLengthF\\|tmXTMajorLineColor\\|tmXTMajorOutwardLengthF\\|tmXTMajorThicknessF\\|tmXTMaxLabelLenF\\|tmXTMaxTicks\\|tmXTMinLabelSpacingF\\|tmXTMinorLengthF\\|tmXTMinorLineColor\\|tmXTMinorOn\\|tmXTMinorOutwardLengthF\\|tmXTMinorPerMajor\\|tmXTMinorThicknessF\\|tmXTMinorValues\\|tmXTMode\\|tmXTOn\\|tmXTPrecision\\|tmXTStyle\\|tmXTTickEndF\\|tmXTTickSpacingF\\|tmXTTickStartF\\|tmXTValues\\|tmXUseBottom\\|tmYLAutoPrecision\\|tmYLBorderOn\\|tmYLDataBottomF\\|tmYLDataTopF\\|tmYLFormat\\|tmYLIrrTensionF\\|tmYLIrregularPoints\\|tmYLLabelAngleF\\|tmYLLabelConstantSpacingF\\|tmYLLabelDeltaF\\|tmYLLabelDirection\\|tmYLLabelFont\\|tmYLLabelFontAspectF\\|tmYLLabelFontColor\\|tmYLLabelFontHeightF\\|tmYLLabelFontQuality\\|tmYLLabelFontThicknessF\\|tmYLLabelFuncCode\\|tmYLLabelJust\\|tmYLLabelStride\\|tmYLLabels\\|tmYLLabelsOn\\|tmYLMajorLengthF\\|tmYLMajorLineColor\\|tmYLMajorOutwardLengthF\\|tmYLMajorThicknessF\\|tmYLMaxLabelLenF\\|tmYLMaxTicks\\|tmYLMinLabelSpacingF\\|tmYLMinorLengthF\\|tmYLMinorLineColor\\|tmYLMinorOn\\|tmYLMinorOutwardLengthF\\|tmYLMinorPerMajor\\|tmYLMinorThicknessF\\|tmYLMinorValues\\|tmYLMode\\|tmYLOn\\|tmYLPrecision\\|tmYLStyle\\|tmYLTickEndF\\|tmYLTickSpacingF\\|tmYLTickStartF\\|tmYLValues\\|tmYMajorGrid\\|tmYMajorGridLineColor\\|tmYMajorGridLineDashPattern\\|tmYMajorGridThicknessF\\|tmYMinorGrid\\|tmYMinorGridLineColor\\|tmYMinorGridLineDashPattern\\|tmYMinorGridThicknessF\\|tmYRAutoPrecision\\|tmYRBorderOn\\|tmYRDataBottomF\\|tmYRDataTopF\\|tmYRFormat\\|tmYRIrrTensionF\\|tmYRIrregularPoints\\|tmYRLabelAngleF\\|tmYRLabelConstantSpacingF\\|tmYRLabelDeltaF\\|tmYRLabelDirection\\|tmYRLabelFont\\|tmYRLabelFontAspectF\\|tmYRLabelFontColor\\|tmYRLabelFontHeightF\\|tmYRLabelFontQuality\\|tmYRLabelFontThicknessF\\|tmYRLabelFuncCode\\|tmYRLabelJust\\|tmYRLabelStride\\|tmYRLabels\\|tmYRLabelsOn\\|tmYRMajorLengthF\\|tmYRMajorLineColor\\|tmYRMajorOutwardLengthF\\|tmYRMajorThicknessF\\|tmYRMaxLabelLenF\\|tmYRMaxTicks\\|tmYRMinLabelSpacingF\\|tmYRMinorLengthF\\|tmYRMinorLineColor\\|tmYRMinorOn\\|tmYRMinorOutwardLengthF\\|tmYRMinorPerMajor\\|tmYRMinorThicknessF\\|tmYRMinorValues\\|tmYRMode\\|tmYROn\\|tmYRPrecision\\|tmYRStyle\\|tmYRTickEndF\\|tmYRTickSpacingF\\|tmYRTickStartF\\|tmYRValues\\|tmYUseLeft\\|trGridType\\|trLineInterpolationOn\\|trXAxisType\\|trXCoordPoints\\|trXInterPoints\\|trXLog\\|trXMaxF\\|trXMinF\\|trXReverse\\|trXSamples\\|trXTensionF\\|trYAxisType\\|trYCoordPoints\\|trYInterPoints\\|trYLog\\|trYMaxF\\|trYMinF\\|trYReverse\\|trYSamples\\|trYTensionF\\|txAngleF\\|txBackgroundFillColor\\|txConstantSpacingF\\|txDirection\\|txFont\\|txFontAspectF\\|txFontColor\\|txFontHeightF\\|txFontOpacityF\\|txFontQuality\\|txFontThicknessF\\|txFuncCode\\|txJust\\|txPerimColor\\|txPerimDashLengthF\\|txPerimDashPattern\\|txPerimOn\\|txPerimSpaceF\\|txPerimThicknessF\\|txPosXF\\|txPosYF\\|txString\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(vcExplicitLabelBarLabelsOn\\|vcFillArrowEdgeColor\\|vcFillArrowEdgeThicknessF\\|vcFillArrowFillColor\\|vcFillArrowHeadInteriorXF\\|vcFillArrowHeadMinFracXF\\|vcFillArrowHeadMinFracYF\\|vcFillArrowHeadXF\\|vcFillArrowHeadYF\\|vcFillArrowMinFracWidthF\\|vcFillArrowWidthF\\|vcFillArrowsOn\\|vcFillOverEdge\\|vcGlyphOpacityF\\|vcGlyphStyle\\|vcLabelBarEndLabelsOn\\|vcLabelFontColor\\|vcLabelFontHeightF\\|vcLabelsOn\\|vcLabelsUseVectorColor\\|vcLevelColors\\|vcLevelCount\\|vcLevelPalette\\|vcLevelSelectionMode\\|vcLevelSpacingF\\|vcLevels\\|vcLineArrowColor\\|vcLineArrowHeadMaxSizeF\\|vcLineArrowHeadMinSizeF\\|vcLineArrowThicknessF\\|vcMagnitudeFormat\\|vcMagnitudeScaleFactorF\\|vcMagnitudeScaleValueF\\|vcMagnitudeScalingMode\\|vcMapDirection\\|vcMaxLevelCount\\|vcMaxLevelValF\\|vcMaxMagnitudeF\\|vcMinAnnoAngleF\\|vcMinAnnoArrowAngleF\\|vcMinAnnoArrowEdgeColor\\|vcMinAnnoArrowFillColor\\|vcMinAnnoArrowLineColor\\|vcMinAnnoArrowMinOffsetF\\|vcMinAnnoArrowSpaceF\\|vcMinAnnoArrowUseVecColor\\|vcMinAnnoBackgroundColor\\|vcMinAnnoConstantSpacingF\\|vcMinAnnoExplicitMagnitudeF\\|vcMinAnnoFont\\|vcMinAnnoFontAspectF\\|vcMinAnnoFontColor\\|vcMinAnnoFontHeightF\\|vcMinAnnoFontQuality\\|vcMinAnnoFontThicknessF\\|vcMinAnnoFuncCode\\|vcMinAnnoJust\\|vcMinAnnoOn\\|vcMinAnnoOrientation\\|vcMinAnnoOrthogonalPosF\\|vcMinAnnoParallelPosF\\|vcMinAnnoPerimColor\\|vcMinAnnoPerimOn\\|vcMinAnnoPerimSpaceF\\|vcMinAnnoPerimThicknessF\\|vcMinAnnoSide\\|vcMinAnnoString1\\|vcMinAnnoString1On\\|vcMinAnnoString2\\|vcMinAnnoString2On\\|vcMinAnnoTextDirection\\|vcMinAnnoZone\\|vcMinDistanceF\\|vcMinFracLengthF\\|vcMinLevelValF\\|vcMinMagnitudeF\\|vcMonoFillArrowEdgeColor\\|vcMonoFillArrowFillColor\\|vcMonoLineArrowColor\\|vcMonoWindBarbColor\\|vcNoDataLabelOn\\|vcNoDataLabelString\\|vcPositionMode\\|vcRefAnnoAngleF\\|vcRefAnnoArrowAngleF\\|vcRefAnnoArrowEdgeColor\\|vcRefAnnoArrowFillColor\\|vcRefAnnoArrowLineColor\\|vcRefAnnoArrowMinOffsetF\\|vcRefAnnoArrowSpaceF\\|vcRefAnnoArrowUseVecColor\\|vcRefAnnoBackgroundColor\\|vcRefAnnoConstantSpacingF\\|vcRefAnnoExplicitMagnitudeF\\|vcRefAnnoFont\\|vcRefAnnoFontAspectF\\|vcRefAnnoFontColor\\|vcRefAnnoFontHeightF\\|vcRefAnnoFontQuality\\|vcRefAnnoFontThicknessF\\|vcRefAnnoFuncCode\\|vcRefAnnoJust\\|vcRefAnnoOn\\|vcRefAnnoOrientation\\|vcRefAnnoOrthogonalPosF\\|vcRefAnnoParallelPosF\\|vcRefAnnoPerimColor\\|vcRefAnnoPerimOn\\|vcRefAnnoPerimSpaceF\\|vcRefAnnoPerimThicknessF\\|vcRefAnnoSide\\|vcRefAnnoString1\\|vcRefAnnoString1On\\|vcRefAnnoString2\\|vcRefAnnoString2On\\|vcRefAnnoTextDirection\\|vcRefAnnoZone\\|vcRefLengthF\\|vcRefMagnitudeF\\|vcScalarFieldData\\|vcScalarMissingValColor\\|vcScalarValueFormat\\|vcScalarValueScaleFactorF\\|vcScalarValueScaleValueF\\|vcScalarValueScalingMode\\|vcSpanLevelPalette\\|vcUseRefAnnoRes\\|vcUseScalarArray\\|vcVectorDrawOrder\\|vcVectorFieldData\\|vcWindBarbCalmCircleSizeF\\|vcWindBarbColor\\|vcWindBarbLineThicknessF\\|vcWindBarbScaleFactorF\\|vcWindBarbTickAngleF\\|vcWindBarbTickLengthF\\|vcWindBarbTickSpacingF\\|vcZeroFLabelAngleF\\|vcZeroFLabelBackgroundColor\\|vcZeroFLabelConstantSpacingF\\|vcZeroFLabelFont\\|vcZeroFLabelFontAspectF\\|vcZeroFLabelFontColor\\|vcZeroFLabelFontHeightF\\|vcZeroFLabelFontQuality\\|vcZeroFLabelFontThicknessF\\|vcZeroFLabelFuncCode\\|vcZeroFLabelJust\\|vcZeroFLabelOn\\|vcZeroFLabelOrthogonalPosF\\|vcZeroFLabelParallelPosF\\|vcZeroFLabelPerimColor\\|vcZeroFLabelPerimOn\\|vcZeroFLabelPerimSpaceF\\|vcZeroFLabelPerimThicknessF\\|vcZeroFLabelSide\\|vcZeroFLabelString\\|vcZeroFLabelTextDirection\\|vcZeroFLabelZone\\|vfCopyData\\|vfDataArray\\|vfExchangeDimensions\\|vfExchangeUVData\\|vfMagMaxV\\|vfMagMinV\\|vfMissingUValueV\\|vfMissingVValueV\\|vfPolarData\\|vfSingleMissingValue\\|vfUDataArray\\|vfUMaxV\\|vfUMinV\\|vfVDataArray\\|vfVMaxV\\|vfVMinV\\|vfXArray\\|vfXCActualEndF\\|vfXCActualStartF\\|vfXCEndIndex\\|vfXCEndSubsetV\\|vfXCEndV\\|vfXCStartIndex\\|vfXCStartSubsetV\\|vfXCStartV\\|vfXCStride\\|vfYArray\\|vfYCActualEndF\\|vfYCActualStartF\\|vfYCEndIndex\\|vfYCEndSubsetV\\|vfYCEndV\\|vfYCStartIndex\\|vfYCStartSubsetV\\|vfYCStartV\\|vfYCStride\\|vpAnnoManagerId\\|vpClipOn\\|vpHeightF\\|vpKeepAspect\\|vpOn\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(vpUseSegments\\|vpWidthF\\|vpXF\\|vpYF\\|wkAntiAlias\\|wkBackgroundColor\\|wkBackgroundOpacityF\\|wkColorMapLen\\|wkColorMap\\|wkColorModel\\|wkColorModel\\|wkDashTableLength\\|wkDefGraphicStyleId\\|wkDeviceLowerX\\|wkDeviceLowerX\\|wkDeviceLowerX\\|wkDeviceLowerY\\|wkDeviceLowerY\\|wkDeviceLowerY\\|wkDeviceUpperX\\|wkDeviceUpperX\\|wkDeviceUpperX\\|wkDeviceUpperY\\|wkDeviceUpperY\\|wkDeviceUpperY\\|wkFileName\\|wkFileName\\|wkFillTableLength\\|wkForegroundColor\\|wkFormat\\|wkFormat\\|wkFullBackground\\|wkFullBackground\\|wkGksWorkId\\|wkHeight\\|wkMarkerTableLength\\|wkMetaName\\|wkOrientation\\|wkOrientation\\|wkOrientation\\|wkPDFFileName\\|wkPDFFormat\\|wkPDFResolution\\|wkPSFileName\\|wkPSFormat\\|wkPSResolution\\|wkPaperHeightF\\|wkPaperHeightF\\|wkPaperHeightF\\|wkPaperSize\\|wkPaperSize\\|wkPaperSize\\|wkPaperWidthF\\|wkPaperWidthF\\|wkPaperWidthF\\|wkTopLevelViews\\|wkViews\\|wkVisualType\\|wkVisualType\\|wkWidth\\|wsCurrentSize\\|wsMaximumSize\\|wsThresholdSize\\|xyComputeXMax\\|xyComputeXMin\\|xyComputeYMax\\|xyComputeYMin\\|xyCoordData\\|xyCoordDataSpec\\|xyCurveDrawOrder\\|xyDashPattern\\|xyDashPatterns\\|xyExplicitLabels\\|xyExplicitLegendLabels\\|xyLabelMode\\|xyLineColor\\|xyLineColors\\|xyLineDashSegLenF\\|xyLineLabelConstantSpacingF\\|xyLineLabelFont\\|xyLineLabelFontAspectF\\|xyLineLabelFontColor\\|xyLineLabelFontColors\\|xyLineLabelFontHeightF\\|xyLineLabelFontQuality\\|xyLineLabelFontThicknessF\\|xyLineLabelFuncCode\\|xyLineOpacities\\|xyLineOpacityF\\|xyLineThicknessF\\|xyLineThicknesses\\|xyMarkLineMode\\|xyMarkLineModes\\|xyMarker\\|xyMarkerColor\\|xyMarkerColors\\|xyMarkerOpacities\\|xyMarkerOpacityF\\|xyMarkerSizeF\\|xyMarkerSizes\\|xyMarkerThicknessF\\|xyMarkerThicknesses\\|xyMarkers\\|xyMonoDashPattern\\|xyMonoLineColor\\|xyMonoLineLabelFontColor\\|xyMonoLineThickness\\|xyMonoMarkLineMode\\|xyMonoMarker\\|xyMonoMarkerColor\\|xyMonoMarkerSize\\|xyMonoMarkerThickness\\|xyXIrrTensionF\\|xyXIrregularPoints\\|xyXStyle\\|xyYIrrTensionF\\|xyYIrregularPoints\\|xyYStyle\\|\\)\\>" 1 font-lock-constant-face)

    ;; ESMValTool graphic resources (same as above, with res_ prepended)
    ("\\<\\(res_amDataXF\\|res_amDataYF\\|res_amJust\\|res_amOn\\|res_amOrthogonalPosF\\|res_amParallelPosF\\|res_amResizeNotify\\|res_amSide\\|res_amTrackData\\|res_amViewId\\|res_amZone\\|res_appDefaultParent\\|res_appFileSuffix\\|res_appResources\\|res_appSysDir\\|res_appUsrDir\\|res_caCopyArrays\\|res_caXArray\\|res_caXCast\\|res_caXMaxV\\|res_caXMinV\\|res_caXMissingV\\|res_caYArray\\|res_caYCast\\|res_caYMaxV\\|res_caYMinV\\|res_caYMissingV\\|res_cnCellFillEdgeColor\\|res_cnCellFillMissingValEdgeColor\\|res_cnConpackParams\\|res_cnConstFEnableFill\\|res_cnConstFLabelAngleF\\|res_cnConstFLabelBackgroundColor\\|res_cnConstFLabelConstantSpacingF\\|res_cnConstFLabelFont\\|res_cnConstFLabelFontAspectF\\|res_cnConstFLabelFontColor\\|res_cnConstFLabelFontHeightF\\|res_cnConstFLabelFontQuality\\|res_cnConstFLabelFontThicknessF\\|res_cnConstFLabelFormat\\|res_cnConstFLabelFuncCode\\|res_cnConstFLabelJust\\|res_cnConstFLabelOn\\|res_cnConstFLabelOrthogonalPosF\\|res_cnConstFLabelParallelPosF\\|res_cnConstFLabelPerimColor\\|res_cnConstFLabelPerimOn\\|res_cnConstFLabelPerimSpaceF\\|res_cnConstFLabelPerimThicknessF\\|res_cnConstFLabelSide\\|res_cnConstFLabelString\\|res_cnConstFLabelTextDirection\\|res_cnConstFLabelZone\\|res_cnConstFUseInfoLabelRes\\|res_cnExplicitLabelBarLabelsOn\\|res_cnExplicitLegendLabelsOn\\|res_cnExplicitLineLabelsOn\\|res_cnFillBackgroundColor\\|res_cnFillColor\\|res_cnFillColors\\|res_cnFillDotSizeF\\|res_cnFillDrawOrder\\|res_cnFillMode\\|res_cnFillOn\\|res_cnFillOpacityF\\|res_cnFillPalette\\|res_cnFillPattern\\|res_cnFillPatterns\\|res_cnFillScaleF\\|res_cnFillScales\\|res_cnFixFillBleed\\|res_cnGridBoundFillColor\\|res_cnGridBoundFillPattern\\|res_cnGridBoundFillScaleF\\|res_cnGridBoundPerimColor\\|res_cnGridBoundPerimDashPattern\\|res_cnGridBoundPerimOn\\|res_cnGridBoundPerimThicknessF\\|res_cnHighLabelAngleF\\|res_cnHighLabelBackgroundColor\\|res_cnHighLabelConstantSpacingF\\|res_cnHighLabelCount\\|res_cnHighLabelFont\\|res_cnHighLabelFontAspectF\\|res_cnHighLabelFontColor\\|res_cnHighLabelFontHeightF\\|res_cnHighLabelFontQuality\\|res_cnHighLabelFontThicknessF\\|res_cnHighLabelFormat\\|res_cnHighLabelFuncCode\\|res_cnHighLabelPerimColor\\|res_cnHighLabelPerimOn\\|res_cnHighLabelPerimSpaceF\\|res_cnHighLabelPerimThicknessF\\|res_cnHighLabelString\\|res_cnHighLabelsOn\\|res_cnHighLowLabelOverlapMode\\|res_cnHighUseLineLabelRes\\|res_cnInfoLabelAngleF\\|res_cnInfoLabelBackgroundColor\\|res_cnInfoLabelConstantSpacingF\\|res_cnInfoLabelFont\\|res_cnInfoLabelFontAspectF\\|res_cnInfoLabelFontColor\\|res_cnInfoLabelFontHeightF\\|res_cnInfoLabelFontQuality\\|res_cnInfoLabelFontThicknessF\\|res_cnInfoLabelFormat\\|res_cnInfoLabelFuncCode\\|res_cnInfoLabelJust\\|res_cnInfoLabelOn\\|res_cnInfoLabelOrthogonalPosF\\|res_cnInfoLabelParallelPosF\\|res_cnInfoLabelPerimColor\\|res_cnInfoLabelPerimOn\\|res_cnInfoLabelPerimSpaceF\\|res_cnInfoLabelPerimThicknessF\\|res_cnInfoLabelSide\\|res_cnInfoLabelString\\|res_cnInfoLabelTextDirection\\|res_cnInfoLabelZone\\|res_cnLabelBarEndLabelsOn\\|res_cnLabelBarEndStyle\\|res_cnLabelDrawOrder\\|res_cnLabelMasking\\|res_cnLabelScaleFactorF\\|res_cnLabelScaleValueF\\|res_cnLabelScalingMode\\|res_cnLegendLevelFlags\\|res_cnLevelCount\\|res_cnLevelFlag\\|res_cnLevelFlags\\|res_cnLevelSelectionMode\\|res_cnLevelSpacingF\\|res_cnLevels\\|res_cnLineColor\\|res_cnLineColors\\|res_cnLineDashPattern\\|res_cnLineDashPatterns\\|res_cnLineDashSegLenF\\|res_cnLineDrawOrder\\|res_cnLineLabelAngleF\\|res_cnLineLabelBackgroundColor\\|res_cnLineLabelConstantSpacingF\\|res_cnLineLabelCount\\|res_cnLineLabelDensityF\\|res_cnLineLabelFont\\|res_cnLineLabelFontAspectF\\|res_cnLineLabelFontColor\\|res_cnLineLabelFontColors\\|res_cnLineLabelFontHeightF\\|res_cnLineLabelFontQuality\\|res_cnLineLabelFontThicknessF\\|res_cnLineLabelFormat\\|res_cnLineLabelFuncCode\\|res_cnLineLabelInterval\\|res_cnLineLabelPerimColor\\|res_cnLineLabelPerimOn\\|res_cnLineLabelPerimSpaceF\\|res_cnLineLabelPerimThicknessF\\|res_cnLineLabelPlacementMode\\|res_cnLineLabelStrings\\|res_cnLineLabelsOn\\|res_cnLinePalette\\|res_cnLineThicknessF\\|res_cnLineThicknesses\\|res_cnLinesOn\\|res_cnLowLabelAngleF\\|res_cnLowLabelBackgroundColor\\|res_cnLowLabelConstantSpacingF\\|res_cnLowLabelCount\\|res_cnLowLabelFont\\|res_cnLowLabelFontAspectF\\|res_cnLowLabelFontColor\\|res_cnLowLabelFontHeightF\\|res_cnLowLabelFontQuality\\|res_cnLowLabelFontThicknessF\\|res_cnLowLabelFormat\\|res_cnLowLabelFuncCode\\|res_cnLowLabelPerimColor\\|res_cnLowLabelPerimOn\\|res_cnLowLabelPerimSpaceF\\|res_cnLowLabelPerimThicknessF\\|res_cnLowLabelString\\|res_cnLowLabelsOn\\|res_cnLowUseHighLabelRes\\|res_cnMaxDataValueFormat\\|res_cnMaxLevelCount\\|res_cnMaxLevelValF\\|res_cnMaxPointDistanceF\\|res_cnMinLevelValF\\|res_cnMissingValFillColor\\|res_cnMissingValFillPattern\\|res_cnMissingValFillScaleF\\|res_cnMissingValPerimColor\\|res_cnMissingValPerimDashPattern\\|res_cnMissingValPerimGridBoundOn\\|res_cnMissingValPerimOn\\|res_cnMissingValPerimThicknessF\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_cnMonoFillColor\\|res_cnMonoFillPattern\\|res_cnMonoFillScale\\|res_cnMonoLevelFlag\\|res_cnMonoLineColor\\|res_cnMonoLineDashPattern\\|res_cnMonoLineLabelFontColor\\|res_cnMonoLineThickness\\|res_cnNoDataLabelOn\\|res_cnNoDataLabelString\\|res_cnOutOfRangeFillColor\\|res_cnOutOfRangeFillPattern\\|res_cnOutOfRangeFillScaleF\\|res_cnOutOfRangePerimColor\\|res_cnOutOfRangePerimDashPattern\\|res_cnOutOfRangePerimOn\\|res_cnOutOfRangePerimThicknessF\\|res_cnRasterCellSizeF\\|res_cnRasterMinCellSizeF\\|res_cnRasterModeOn\\|res_cnRasterSampleFactorF\\|res_cnRasterSmoothingOn\\|res_cnScalarFieldData\\|res_cnSmoothingDistanceF\\|res_cnSmoothingOn\\|res_cnSmoothingTensionF\\|res_cnSpanFillPalette\\|res_cnSpanLinePalette\\|res_ctCopyTables\\|res_ctXElementSize\\|res_ctXMaxV\\|res_ctXMinV\\|res_ctXMissingV\\|res_ctXTable\\|res_ctXTableLengths\\|res_ctXTableType\\|res_ctYElementSize\\|res_ctYMaxV\\|res_ctYMinV\\|res_ctYMissingV\\|res_ctYTable\\|res_ctYTableLengths\\|res_ctYTableType\\|res_dcDelayCompute\\|res_errBuffer\\|res_errFileName\\|res_errFilePtr\\|res_errLevel\\|res_errPrint\\|res_errUnitNumber\\|res_gsClipOn\\|res_gsColors\\|res_gsEdgeColor\\|res_gsEdgeDashPattern\\|res_gsEdgeDashSegLenF\\|res_gsEdgeThicknessF\\|res_gsEdgesOn\\|res_gsFillBackgroundColor\\|res_gsFillColor\\|res_gsFillDotSizeF\\|res_gsFillIndex\\|res_gsFillLineThicknessF\\|res_gsFillOpacityF\\|res_gsFillScaleF\\|res_gsFont\\|res_gsFontAspectF\\|res_gsFontColor\\|res_gsFontHeightF\\|res_gsFontOpacityF\\|res_gsFontQuality\\|res_gsFontThicknessF\\|res_gsLineColor\\|res_gsLineDashPattern\\|res_gsLineDashSegLenF\\|res_gsLineLabelConstantSpacingF\\|res_gsLineLabelFont\\|res_gsLineLabelFontAspectF\\|res_gsLineLabelFontColor\\|res_gsLineLabelFontHeightF\\|res_gsLineLabelFontQuality\\|res_gsLineLabelFontThicknessF\\|res_gsLineLabelFuncCode\\|res_gsLineLabelString\\|res_gsLineOpacityF\\|res_gsLineThicknessF\\|res_gsMarkerColor\\|res_gsMarkerIndex\\|res_gsMarkerOpacityF\\|res_gsMarkerSizeF\\|res_gsMarkerThicknessF\\|res_gsSegments\\|res_gsTextAngleF\\|res_gsTextConstantSpacingF\\|res_gsTextDirection\\|res_gsTextFuncCode\\|res_gsTextJustification\\|res_gsnAboveYRefLineBarColors\\|res_gsnAboveYRefLineBarFillScales\\|res_gsnAboveYRefLineBarPatterns\\|res_gsnAboveYRefLineColor\\|res_gsnAddCyclic\\|res_gsnAttachBorderOn\\|res_gsnAttachPlotsXAxis\\|res_gsnBelowYRefLineBarColors\\|res_gsnBelowYRefLineBarFillScales\\|res_gsnBelowYRefLineBarPatterns\\|res_gsnBelowYRefLineColor\\|res_gsnBoxMargin\\|res_gsnCenterString\\|res_gsnCenterStringFontColor\\|res_gsnCenterStringFontHeightF\\|res_gsnCenterStringFuncCode\\|res_gsnCenterStringOrthogonalPosF\\|res_gsnCenterStringParallelPosF\\|res_gsnContourLineThicknessesScale\\|res_gsnContourNegLineDashPattern\\|res_gsnContourPosLineDashPattern\\|res_gsnContourZeroLineThicknessF\\|res_gsnDebugWriteFileName\\|res_gsnDraw\\|res_gsnFrame\\|res_gsnHistogramBarColors\\|res_gsnHistogramBarWidthPercent\\|res_gsnHistogramBinIntervals\\|res_gsnHistogramBinMissing\\|res_gsnHistogramBinWidth\\|res_gsnHistogramClassIntervals\\|res_gsnHistogramCompare\\|res_gsnHistogramComputePercentages\\|res_gsnHistogramComputePercentagesNoMissing\\|res_gsnHistogramDiscreteBinValues\\|res_gsnHistogramDiscreteClassValues\\|res_gsnHistogramHorizontal\\|res_gsnHistogramMinMaxBinsOn\\|res_gsnHistogramNumberOfBins\\|res_gsnHistogramPercentSign\\|res_gsnHistogramSelectNiceIntervals\\|res_gsnLeftString\\|res_gsnLeftStringFontColor\\|res_gsnLeftStringFontHeightF\\|res_gsnLeftStringFuncCode\\|res_gsnLeftStringOrthogonalPosF\\|res_gsnLeftStringParallelPosF\\|res_gsnLeftXRefLineBarColors\\|res_gsnLeftXRefLineBarFillScales\\|res_gsnLeftXRefLineBarPatterns\\|res_gsnLeftXRefLineColor\\|res_gsnMajorLatSpacing\\|res_gsnMajorLonSpacing\\|res_gsnMaskLambertConformal\\|res_gsnMaskLambertConformalOutlineOn\\|res_gsnMaximize\\|res_gsnMinorLatSpacing\\|res_gsnMinorLonSpacing\\|res_gsnPanelBottom\\|res_gsnPanelCenter\\|res_gsnPanelDebug\\|res_gsnPanelFigureStrings\\|res_gsnPanelFigureStringsBackgroundFillColor\\|res_gsnPanelFigureStringsFontHeightF\\|res_gsnPanelFigureStringsJust\\|res_gsnPanelFigureStringsPerimOn\\|res_gsnPanelLabelBar\\|res_gsnPanelLeft\\|res_gsnPanelMainFont\\|res_gsnPanelMainFontColor\\|res_gsnPanelMainFontHeightF\\|res_gsnPanelMainPosXF\\|res_gsnPanelMainPosYF\\|res_gsnPanelMainString\\|res_gsnPanelRight\\|res_gsnPanelRowSpec\\|res_gsnPanelScalePlotIndex\\|res_gsnPanelTop\\|res_gsnPanelXF\\|res_gsnPanelXWhiteSpacePercent\\|res_gsnPanelYF\\|res_gsnPanelYWhiteSpacePercent\\|res_gsnPaperHeight\\|res_gsnPaperMargin\\|res_gsnPaperOrientation\\|res_gsnPaperWidth\\|res_gsnPolar\\|res_gsnPolarLabelDistance\\|res_gsnPolarLabelFont\\|res_gsnPolarLabelFontHeightF\\|res_gsnPolarLabelSpacing\\|res_gsnPolarTime\\|res_gsnPolarUT\\|res_gsnRightString\\|res_gsnRightStringFontColor\\|res_gsnRightStringFontHeightF\\|res_gsnRightStringFuncCode\\|res_gsnRightStringOrthogonalPosF\\|res_gsnRightStringParallelPosF\\|res_gsnRightXRefLineBarColors\\|res_gsnRightXRefLineBarFillScales\\|res_gsnRightXRefLineBarPatterns\\|res_gsnRightXRefLineColor\\|res_gsnScalarContour\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_gsnScale\\|res_gsnShape\\|res_gsnSpreadColorEnd\\|res_gsnSpreadColorStart\\|res_gsnSpreadColors\\|res_gsnStringFont\\|res_gsnStringFontColor\\|res_gsnStringFontHeightF\\|res_gsnStringFuncCode\\|res_gsnTickMarksOn\\|res_gsnXAxisIrregular2Linear\\|res_gsnXAxisIrregular2Log\\|res_gsnXRefLine\\|res_gsnXRefLineColor\\|res_gsnXRefLineColors\\|res_gsnXRefLineDashPattern\\|res_gsnXRefLineDashPatterns\\|res_gsnXRefLineThicknessF\\|res_gsnXRefLineThicknesses\\|res_gsnXYAboveFillColors\\|res_gsnXYBarChart\\|res_gsnXYBarChartBarWidth\\|res_gsnXYBarChartColors\\|res_gsnXYBarChartColors2\\|res_gsnXYBarChartFillDotSizeF\\|res_gsnXYBarChartFillLineThicknessF\\|res_gsnXYBarChartFillOpacityF\\|res_gsnXYBarChartFillScaleF\\|res_gsnXYBarChartOutlineOnly\\|res_gsnXYBarChartOutlineThicknessF\\|res_gsnXYBarChartPatterns\\|res_gsnXYBarChartPatterns2\\|res_gsnXYBelowFillColors\\|res_gsnXYFillColors\\|res_gsnXYFillOpacities\\|res_gsnXYLeftFillColors\\|res_gsnXYRightFillColors\\|res_gsnYAxisIrregular2Linear\\|res_gsnYAxisIrregular2Log\\|res_gsnYRefLine\\|res_gsnYRefLineColor\\|res_gsnYRefLineColors\\|res_gsnYRefLineDashPattern\\|res_gsnYRefLineDashPatterns\\|res_gsnYRefLineThicknessF\\|res_gsnYRefLineThicknesses\\|res_gsnZonalMean\\|res_gsnZonalMeanXMaxF\\|res_gsnZonalMeanXMinF\\|res_gsnZonalMeanYRefLine\\|res_lbAutoManage\\|res_lbBottomMarginF\\|res_lbBoxCount\\|res_lbBoxEndCapStyle\\|res_lbBoxFractions\\|res_lbBoxLineColor\\|res_lbBoxLineDashPattern\\|res_lbBoxLineDashSegLenF\\|res_lbBoxLineThicknessF\\|res_lbBoxLinesOn\\|res_lbBoxMajorExtentF\\|res_lbBoxMinorExtentF\\|res_lbBoxSeparatorLinesOn\\|res_lbBoxSizing\\|res_lbFillBackground\\|res_lbFillColor\\|res_lbFillColors\\|res_lbFillDotSizeF\\|res_lbFillLineThicknessF\\|res_lbFillOpacityF\\|res_lbFillPattern\\|res_lbFillPatterns\\|res_lbFillScaleF\\|res_lbFillScales\\|res_lbJustification\\|res_lbLabelAlignment\\|res_lbLabelAngleF\\|res_lbLabelAutoStride\\|res_lbLabelBarOn\\|res_lbLabelConstantSpacingF\\|res_lbLabelDirection\\|res_lbLabelFont\\|res_lbLabelFontAspectF\\|res_lbLabelFontColor\\|res_lbLabelFontHeightF\\|res_lbLabelFontQuality\\|res_lbLabelFontThicknessF\\|res_lbLabelFuncCode\\|res_lbLabelJust\\|res_lbLabelOffsetF\\|res_lbLabelPosition\\|res_lbLabelStride\\|res_lbLabelStrings\\|res_lbLabelsOn\\|res_lbLeftMarginF\\|res_lbMaxLabelLenF\\|res_lbMinLabelSpacingF\\|res_lbMonoFillColor\\|res_lbMonoFillPattern\\|res_lbMonoFillScale\\|res_lbOrientation\\|res_lbOverrideFillOpacity\\|res_lbPerimColor\\|res_lbPerimDashPattern\\|res_lbPerimDashSegLenF\\|res_lbPerimFill\\|res_lbPerimFillColor\\|res_lbPerimOn\\|res_lbPerimThicknessF\\|res_lbRasterFillOn\\|res_lbRightMarginF\\|res_lbTitleAngleF\\|res_lbTitleConstantSpacingF\\|res_lbTitleDirection\\|res_lbTitleExtentF\\|res_lbTitleFont\\|res_lbTitleFontAspectF\\|res_lbTitleFontColor\\|res_lbTitleFontHeightF\\|res_lbTitleFontQuality\\|res_lbTitleFontThicknessF\\|res_lbTitleFuncCode\\|res_lbTitleJust\\|res_lbTitleOffsetF\\|res_lbTitleOn\\|res_lbTitlePosition\\|res_lbTitleString\\|res_lbTopMarginF\\|res_lgAutoManage\\|res_lgBottomMarginF\\|res_lgBoxBackground\\|res_lgBoxLineColor\\|res_lgBoxLineDashPattern\\|res_lgBoxLineDashSegLenF\\|res_lgBoxLineThicknessF\\|res_lgBoxLinesOn\\|res_lgBoxMajorExtentF\\|res_lgBoxMinorExtentF\\|res_lgDashIndex\\|res_lgDashIndexes\\|res_lgItemCount\\|res_lgItemOrder\\|res_lgItemPlacement\\|res_lgItemPositions\\|res_lgItemType\\|res_lgItemTypes\\|res_lgJustification\\|res_lgLabelAlignment\\|res_lgLabelAngleF\\|res_lgLabelAutoStride\\|res_lgLabelConstantSpacingF\\|res_lgLabelDirection\\|res_lgLabelFont\\|res_lgLabelFontAspectF\\|res_lgLabelFontColor\\|res_lgLabelFontHeightF\\|res_lgLabelFontQuality\\|res_lgLabelFontThicknessF\\|res_lgLabelFuncCode\\|res_lgLabelJust\\|res_lgLabelOffsetF\\|res_lgLabelPosition\\|res_lgLabelStride\\|res_lgLabelStrings\\|res_lgLabelsOn\\|res_lgLeftMarginF\\|res_lgLegendOn\\|res_lgLineColor\\|res_lgLineColors\\|res_lgLineDashSegLenF\\|res_lgLineDashSegLens\\|res_lgLineLabelConstantSpacingF\\|res_lgLineLabelFont\\|res_lgLineLabelFontAspectF\\|res_lgLineLabelFontColor\\|res_lgLineLabelFontColors\\|res_lgLineLabelFontHeightF\\|res_lgLineLabelFontHeights\\|res_lgLineLabelFontQuality\\|res_lgLineLabelFontThicknessF\\|res_lgLineLabelFuncCode\\|res_lgLineLabelStrings\\|res_lgLineLabelsOn\\|res_lgLineThicknessF\\|res_lgLineThicknesses\\|res_lgMarkerColor\\|res_lgMarkerColors\\|res_lgMarkerIndex\\|res_lgMarkerIndexes\\|res_lgMarkerSizeF\\|res_lgMarkerSizes\\|res_lgMarkerThicknessF\\|res_lgMarkerThicknesses\\|res_lgMonoDashIndex\\|res_lgMonoItemType\\|res_lgMonoLineColor\\|res_lgMonoLineDashSegLen\\|res_lgMonoLineLabelFontColor\\|res_lgMonoLineLabelFontHeight\\|res_lgMonoLineThickness\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_lgMonoMarkerColor\\|res_lgMonoMarkerIndex\\|res_lgMonoMarkerSize\\|res_lgMonoMarkerThickness\\|res_lgOrientation\\|res_lgPerimColor\\|res_lgPerimDashPattern\\|res_lgPerimDashSegLenF\\|res_lgPerimFill\\|res_lgPerimFillColor\\|res_lgPerimOn\\|res_lgPerimThicknessF\\|res_lgRightMarginF\\|res_lgTitleAngleF\\|res_lgTitleConstantSpacingF\\|res_lgTitleDirection\\|res_lgTitleExtentF\\|res_lgTitleFont\\|res_lgTitleFontAspectF\\|res_lgTitleFontColor\\|res_lgTitleFontHeightF\\|res_lgTitleFontQuality\\|res_lgTitleFontThicknessF\\|res_lgTitleFuncCode\\|res_lgTitleJust\\|res_lgTitleOffsetF\\|res_lgTitleOn\\|res_lgTitlePosition\\|res_lgTitleString\\|res_lgTopMarginF\\|res_mpAreaGroupCount\\|res_mpAreaMaskingOn\\|res_mpAreaNames\\|res_mpAreaTypes\\|res_mpBottomAngleF\\|res_mpBottomMapPosF\\|res_mpBottomNDCF\\|res_mpBottomNPCF\\|res_mpBottomPointLatF\\|res_mpBottomPointLonF\\|res_mpBottomWindowF\\|res_mpCenterLatF\\|res_mpCenterLonF\\|res_mpCenterRotF\\|res_mpCountyLineColor\\|res_mpCountyLineDashPattern\\|res_mpCountyLineDashSegLenF\\|res_mpCountyLineThicknessF\\|res_mpDataBaseVersion\\|res_mpDataResolution\\|res_mpDataSetName\\|res_mpDefaultFillColor\\|res_mpDefaultFillPattern\\|res_mpDefaultFillScaleF\\|res_mpDynamicAreaGroups\\|res_mpEllipticalBoundary\\|res_mpFillAreaSpecifiers\\|res_mpFillBoundarySets\\|res_mpFillColor\\|res_mpFillColors\\|res_mpFillDotSizeF\\|res_mpFillDrawOrder\\|res_mpFillOn\\|res_mpFillPatternBackground\\|res_mpFillPattern\\|res_mpFillPatterns\\|res_mpFillScaleF\\|res_mpFillScales\\|res_mpFixedAreaGroups\\|res_mpGeophysicalLineColor\\|res_mpGeophysicalLineDashPattern\\|res_mpGeophysicalLineDashSegLenF\\|res_mpGeophysicalLineThicknessF\\|res_mpGreatCircleLinesOn\\|res_mpGridAndLimbDrawOrder\\|res_mpGridAndLimbOn\\|res_mpGridLatSpacingF\\|res_mpGridLineColor\\|res_mpGridLineDashPattern\\|res_mpGridLineDashSegLenF\\|res_mpGridLineThicknessF\\|res_mpGridLonSpacingF\\|res_mpGridMaskMode\\|res_mpGridMaxLatF\\|res_mpGridPolarLonSpacingF\\|res_mpGridSpacingF\\|res_mpInlandWaterFillColor\\|res_mpInlandWaterFillPattern\\|res_mpInlandWaterFillScaleF\\|res_mpLabelDrawOrder\\|res_mpLabelFontColor\\|res_mpLabelFontHeightF\\|res_mpLabelsOn\\|res_mpLambertMeridianF\\|res_mpLambertParallel1F\\|res_mpLambertParallel2F\\|res_mpLandFillColor\\|res_mpLandFillPattern\\|res_mpLandFillScaleF\\|res_mpLeftAngleF\\|res_mpLeftCornerLatF\\|res_mpLeftCornerLonF\\|res_mpLeftMapPosF\\|res_mpLeftNDCF\\|res_mpLeftNPCF\\|res_mpLeftPointLatF\\|res_mpLeftPointLonF\\|res_mpLeftWindowF\\|res_mpLimbLineColor\\|res_mpLimbLineDashPattern\\|res_mpLimbLineDashSegLenF\\|res_mpLimbLineThicknessF\\|res_mpLimitMode\\|res_mpMaskAreaSpecifiers\\|res_mpMaskOutlineSpecifiers\\|res_mpMaxLatF\\|res_mpMaxLonF\\|res_mpMinLatF\\|res_mpMinLonF\\|res_mpMonoFillColor\\|res_mpMonoFillPattern\\|res_mpMonoFillScale\\|res_mpNationalLineColor\\|res_mpNationalLineDashPattern\\|res_mpNationalLineDashSegLenF\\|res_mpNationalLineThicknessF\\|res_mpOceanFillColor\\|res_mpOceanFillPattern\\|res_mpOceanFillScaleF\\|res_mpOutlineBoundarySets\\|res_mpOutlineDrawOrder\\|res_mpOutlineMaskingOn\\|res_mpOutlineOn\\|res_mpOutlineSpecifiers\\|res_mpPerimDrawOrder\\|res_mpPerimLineColor\\|res_mpPerimLineDashPattern\\|res_mpPerimLineDashSegLenF\\|res_mpPerimLineThicknessF\\|res_mpPerimOn\\|res_mpPolyMode\\|res_mpProjection\\|res_mpProvincialLineColor\\|res_mpProvincialLineDashPattern\\|res_mpProvincialLineDashSegLenF\\|res_mpProvincialLineThicknessF\\|res_mpRelativeCenterLat\\|res_mpRelativeCenterLon\\|res_mpRightAngleF\\|res_mpRightCornerLatF\\|res_mpRightCornerLonF\\|res_mpRightMapPosF\\|res_mpRightNDCF\\|res_mpRightNPCF\\|res_mpRightPointLatF\\|res_mpRightPointLonF\\|res_mpRightWindowF\\|res_mpSatelliteAngle1F\\|res_mpSatelliteAngle2F\\|res_mpSatelliteDistF\\|res_mpShapeMode\\|res_mpSpecifiedFillColors\\|res_mpSpecifiedFillDirectIndexing\\|res_mpSpecifiedFillPatterns\\|res_mpSpecifiedFillPriority\\|res_mpSpecifiedFillScales\\|res_mpTopAngleF\\|res_mpTopMapPosF\\|res_mpTopNDCF\\|res_mpTopNPCF\\|res_mpTopPointLatF\\|res_mpTopPointLonF\\|res_mpTopWindowF\\|res_mpUSStateLineColor\\|res_mpUSStateLineDashPattern\\|res_mpUSStateLineDashSegLenF\\|res_mpUSStateLineThicknessF\\|res_pmAnnoManagers\\|res_pmAnnoViews\\|res_pmLabelBarDisplayMode\\|res_pmLabelBarHeightF\\|res_pmLabelBarKeepAspect\\|res_pmLabelBarOrthogonalPosF\\|res_pmLabelBarParallelPosF\\|res_pmLabelBarSide\\|res_pmLabelBarWidthF\\|res_pmLabelBarZone\\|res_pmLegendDisplayMode\\|res_pmLegendHeightF\\|res_pmLegendKeepAspect\\|res_pmLegendOrthogonalPosF\\|res_pmLegendParallelPosF\\|res_pmLegendSide\\|res_pmLegendWidthF\\|res_pmLegendZone\\|res_pmOverlaySequenceIds\\|res_pmTickMarkDisplayMode\\|res_pmTickMarkZone\\|res_pmTitleDisplayMode\\|res_pmTitleZone\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_prGraphicStyle\\|res_prPolyType\\|res_prXArray\\|res_prYArray\\|res_sfCopyData\\|res_sfCopyData\\|res_sfDataArray\\|res_sfDataArray\\|res_sfDataMaxV\\|res_sfDataMaxV\\|res_sfDataMinV\\|res_sfDataMinV\\|res_sfElementNodes\\|res_sfExchangeDimensions\\|res_sfFirstNodeIndex\\|res_sfMissingValueV\\|res_sfMissingValueV\\|res_sfXArray\\|res_sfXArray\\|res_sfXCActualEndF\\|res_sfXCActualEndF\\|res_sfXCActualStartF\\|res_sfXCActualStartF\\|res_sfXCEndIndex\\|res_sfXCEndSubsetV\\|res_sfXCEndV\\|res_sfXCStartIndex\\|res_sfXCStartSubsetV\\|res_sfXCStartV\\|res_sfXCStride\\|res_sfXCellBounds\\|res_sfYArray\\|res_sfYArray\\|res_sfYCActualEndF\\|res_sfYCActualEndF\\|res_sfYCActualStartF\\|res_sfYCActualStartF\\|res_sfYCEndIndex\\|res_sfYCEndSubsetV\\|res_sfYCEndV\\|res_sfYCStartIndex\\|res_sfYCStartSubsetV\\|res_sfYCStartV\\|res_sfYCStride\\|res_sfYCellBounds\\|res_stArrowLengthF\\|res_stArrowStride\\|res_stCrossoverCheckCount\\|res_stExplicitLabelBarLabelsOn\\|res_stLabelBarEndLabelsOn\\|res_stLabelFormat\\|res_stLengthCheckCount\\|res_stLevelColors\\|res_stLevelCount\\|res_stLevelPalette\\|res_stLevelSelectionMode\\|res_stLevelSpacingF\\|res_stLevels\\|res_stLineColor\\|res_stLineOpacityF\\|res_stLineStartStride\\|res_stLineThicknessF\\|res_stMapDirection\\|res_stMaxLevelCount\\|res_stMaxLevelValF\\|res_stMinArrowSpacingF\\|res_stMinDistanceF\\|res_stMinLevelValF\\|res_stMinLineSpacingF\\|res_stMinStepFactorF\\|res_stMonoLineColor\\|res_stNoDataLabelOn\\|res_stNoDataLabelString\\|res_stScalarFieldData\\|res_stScalarMissingValColor\\|res_stSpanLevelPalette\\|res_stStepSizeF\\|res_stStreamlineDrawOrder\\|res_stUseScalarArray\\|res_stVectorFieldData\\|res_stZeroFLabelAngleF\\|res_stZeroFLabelBackgroundColor\\|res_stZeroFLabelConstantSpacingF\\|res_stZeroFLabelFont\\|res_stZeroFLabelFontAspectF\\|res_stZeroFLabelFontColor\\|res_stZeroFLabelFontHeightF\\|res_stZeroFLabelFontQuality\\|res_stZeroFLabelFontThicknessF\\|res_stZeroFLabelFuncCode\\|res_stZeroFLabelJust\\|res_stZeroFLabelOn\\|res_stZeroFLabelOrthogonalPosF\\|res_stZeroFLabelParallelPosF\\|res_stZeroFLabelPerimColor\\|res_stZeroFLabelPerimOn\\|res_stZeroFLabelPerimSpaceF\\|res_stZeroFLabelPerimThicknessF\\|res_stZeroFLabelSide\\|res_stZeroFLabelString\\|res_stZeroFLabelTextDirection\\|res_stZeroFLabelZone\\|res_tfDoNDCOverlay\\|res_tfPlotManagerOn\\|res_tfPolyDrawList\\|res_tfPolyDrawOrder\\|res_tiDeltaF\\|res_tiMainAngleF\\|res_tiMainConstantSpacingF\\|res_tiMainDirection\\|res_tiMainFont\\|res_tiMainFontAspectF\\|res_tiMainFontColor\\|res_tiMainFontHeightF\\|res_tiMainFontQuality\\|res_tiMainFontThicknessF\\|res_tiMainFuncCode\\|res_tiMainJust\\|res_tiMainOffsetXF\\|res_tiMainOffsetYF\\|res_tiMainOn\\|res_tiMainPosition\\|res_tiMainSide\\|res_tiMainString\\|res_tiUseMainAttributes\\|res_tiXAxisAngleF\\|res_tiXAxisConstantSpacingF\\|res_tiXAxisDirection\\|res_tiXAxisFont\\|res_tiXAxisFontAspectF\\|res_tiXAxisFontColor\\|res_tiXAxisFontHeightF\\|res_tiXAxisFontQuality\\|res_tiXAxisFontThicknessF\\|res_tiXAxisFuncCode\\|res_tiXAxisJust\\|res_tiXAxisOffsetXF\\|res_tiXAxisOffsetYF\\|res_tiXAxisOn\\|res_tiXAxisPosition\\|res_tiXAxisSide\\|res_tiXAxisString\\|res_tiYAxisAngleF\\|res_tiYAxisConstantSpacingF\\|res_tiYAxisDirection\\|res_tiYAxisFont\\|res_tiYAxisFontAspectF\\|res_tiYAxisFontColor\\|res_tiYAxisFontHeightF\\|res_tiYAxisFontQuality\\|res_tiYAxisFontThicknessF\\|res_tiYAxisFuncCode\\|res_tiYAxisJust\\|res_tiYAxisOffsetXF\\|res_tiYAxisOffsetYF\\|res_tiYAxisOn\\|res_tiYAxisPosition\\|res_tiYAxisSide\\|res_tiYAxisString\\|res_tmBorderLineColor\\|res_tmBorderThicknessF\\|res_tmEqualizeXYSizes\\|res_tmLabelAutoStride\\|res_tmSciNoteCutoff\\|res_tmXBAutoPrecision\\|res_tmXBBorderOn\\|res_tmXBDataLeftF\\|res_tmXBDataRightF\\|res_tmXBFormat\\|res_tmXBIrrTensionF\\|res_tmXBIrregularPoints\\|res_tmXBLabelAngleF\\|res_tmXBLabelConstantSpacingF\\|res_tmXBLabelDeltaF\\|res_tmXBLabelDirection\\|res_tmXBLabelFont\\|res_tmXBLabelFontAspectF\\|res_tmXBLabelFontColor\\|res_tmXBLabelFontHeightF\\|res_tmXBLabelFontQuality\\|res_tmXBLabelFontThicknessF\\|res_tmXBLabelFuncCode\\|res_tmXBLabelJust\\|res_tmXBLabelStride\\|res_tmXBLabels\\|res_tmXBLabelsOn\\|res_tmXBMajorLengthF\\|res_tmXBMajorLineColor\\|res_tmXBMajorOutwardLengthF\\|res_tmXBMajorThicknessF\\|res_tmXBMaxLabelLenF\\|res_tmXBMaxTicks\\|res_tmXBMinLabelSpacingF\\|res_tmXBMinorLengthF\\|res_tmXBMinorLineColor\\|res_tmXBMinorOn\\|res_tmXBMinorOutwardLengthF\\|res_tmXBMinorPerMajor\\|res_tmXBMinorThicknessF\\|res_tmXBMinorValues\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_tmXBMode\\|res_tmXBOn\\|res_tmXBPrecision\\|res_tmXBStyle\\|res_tmXBTickEndF\\|res_tmXBTickSpacingF\\|res_tmXBTickStartF\\|res_tmXBValues\\|res_tmXMajorGrid\\|res_tmXMajorGridLineColor\\|res_tmXMajorGridLineDashPattern\\|res_tmXMajorGridThicknessF\\|res_tmXMinorGrid\\|res_tmXMinorGridLineColor\\|res_tmXMinorGridLineDashPattern\\|res_tmXMinorGridThicknessF\\|res_tmXTAutoPrecision\\|res_tmXTBorderOn\\|res_tmXTDataLeftF\\|res_tmXTDataRightF\\|res_tmXTFormat\\|res_tmXTIrrTensionF\\|res_tmXTIrregularPoints\\|res_tmXTLabelAngleF\\|res_tmXTLabelConstantSpacingF\\|res_tmXTLabelDeltaF\\|res_tmXTLabelDirection\\|res_tmXTLabelFont\\|res_tmXTLabelFontAspectF\\|res_tmXTLabelFontColor\\|res_tmXTLabelFontHeightF\\|res_tmXTLabelFontQuality\\|res_tmXTLabelFontThicknessF\\|res_tmXTLabelFuncCode\\|res_tmXTLabelJust\\|res_tmXTLabelStride\\|res_tmXTLabels\\|res_tmXTLabelsOn\\|res_tmXTMajorLengthF\\|res_tmXTMajorLineColor\\|res_tmXTMajorOutwardLengthF\\|res_tmXTMajorThicknessF\\|res_tmXTMaxLabelLenF\\|res_tmXTMaxTicks\\|res_tmXTMinLabelSpacingF\\|res_tmXTMinorLengthF\\|res_tmXTMinorLineColor\\|res_tmXTMinorOn\\|res_tmXTMinorOutwardLengthF\\|res_tmXTMinorPerMajor\\|res_tmXTMinorThicknessF\\|res_tmXTMinorValues\\|res_tmXTMode\\|res_tmXTOn\\|res_tmXTPrecision\\|res_tmXTStyle\\|res_tmXTTickEndF\\|res_tmXTTickSpacingF\\|res_tmXTTickStartF\\|res_tmXTValues\\|res_tmXUseBottom\\|res_tmYLAutoPrecision\\|res_tmYLBorderOn\\|res_tmYLDataBottomF\\|res_tmYLDataTopF\\|res_tmYLFormat\\|res_tmYLIrrTensionF\\|res_tmYLIrregularPoints\\|res_tmYLLabelAngleF\\|res_tmYLLabelConstantSpacingF\\|res_tmYLLabelDeltaF\\|res_tmYLLabelDirection\\|res_tmYLLabelFont\\|res_tmYLLabelFontAspectF\\|res_tmYLLabelFontColor\\|res_tmYLLabelFontHeightF\\|res_tmYLLabelFontQuality\\|res_tmYLLabelFontThicknessF\\|res_tmYLLabelFuncCode\\|res_tmYLLabelJust\\|res_tmYLLabelStride\\|res_tmYLLabels\\|res_tmYLLabelsOn\\|res_tmYLMajorLengthF\\|res_tmYLMajorLineColor\\|res_tmYLMajorOutwardLengthF\\|res_tmYLMajorThicknessF\\|res_tmYLMaxLabelLenF\\|res_tmYLMaxTicks\\|res_tmYLMinLabelSpacingF\\|res_tmYLMinorLengthF\\|res_tmYLMinorLineColor\\|res_tmYLMinorOn\\|res_tmYLMinorOutwardLengthF\\|res_tmYLMinorPerMajor\\|res_tmYLMinorThicknessF\\|res_tmYLMinorValues\\|res_tmYLMode\\|res_tmYLOn\\|res_tmYLPrecision\\|res_tmYLStyle\\|res_tmYLTickEndF\\|res_tmYLTickSpacingF\\|res_tmYLTickStartF\\|res_tmYLValues\\|res_tmYMajorGrid\\|res_tmYMajorGridLineColor\\|res_tmYMajorGridLineDashPattern\\|res_tmYMajorGridThicknessF\\|res_tmYMinorGrid\\|res_tmYMinorGridLineColor\\|res_tmYMinorGridLineDashPattern\\|res_tmYMinorGridThicknessF\\|res_tmYRAutoPrecision\\|res_tmYRBorderOn\\|res_tmYRDataBottomF\\|res_tmYRDataTopF\\|res_tmYRFormat\\|res_tmYRIrrTensionF\\|res_tmYRIrregularPoints\\|res_tmYRLabelAngleF\\|res_tmYRLabelConstantSpacingF\\|res_tmYRLabelDeltaF\\|res_tmYRLabelDirection\\|res_tmYRLabelFont\\|res_tmYRLabelFontAspectF\\|res_tmYRLabelFontColor\\|res_tmYRLabelFontHeightF\\|res_tmYRLabelFontQuality\\|res_tmYRLabelFontThicknessF\\|res_tmYRLabelFuncCode\\|res_tmYRLabelJust\\|res_tmYRLabelStride\\|res_tmYRLabels\\|res_tmYRLabelsOn\\|res_tmYRMajorLengthF\\|res_tmYRMajorLineColor\\|res_tmYRMajorOutwardLengthF\\|res_tmYRMajorThicknessF\\|res_tmYRMaxLabelLenF\\|res_tmYRMaxTicks\\|res_tmYRMinLabelSpacingF\\|res_tmYRMinorLengthF\\|res_tmYRMinorLineColor\\|res_tmYRMinorOn\\|res_tmYRMinorOutwardLengthF\\|res_tmYRMinorPerMajor\\|res_tmYRMinorThicknessF\\|res_tmYRMinorValues\\|res_tmYRMode\\|res_tmYROn\\|res_tmYRPrecision\\|res_tmYRStyle\\|res_tmYRTickEndF\\|res_tmYRTickSpacingF\\|res_tmYRTickStartF\\|res_tmYRValues\\|res_tmYUseLeft\\|res_trGridType\\|res_trLineInterpolationOn\\|res_trXAxisType\\|res_trXCoordPoints\\|res_trXInterPoints\\|res_trXLog\\|res_trXMaxF\\|res_trXMinF\\|res_trXReverse\\|res_trXSamples\\|res_trXTensionF\\|res_trYAxisType\\|res_trYCoordPoints\\|res_trYInterPoints\\|res_trYLog\\|res_trYMaxF\\|res_trYMinF\\|res_trYReverse\\|res_trYSamples\\|res_trYTensionF\\|res_txAngleF\\|res_txBackgroundFillColor\\|res_txConstantSpacingF\\|res_txDirection\\|res_txFont\\|res_txFontAspectF\\|res_txFontColor\\|res_txFontHeightF\\|res_txFontOpacityF\\|res_txFontQuality\\|res_txFontThicknessF\\|res_txFuncCode\\|res_txJust\\|res_txPerimColor\\|res_txPerimDashLengthF\\|res_txPerimDashPattern\\|res_txPerimOn\\|res_txPerimSpaceF\\|res_txPerimThicknessF\\|res_txPosXF\\|res_txPosYF\\|res_txString\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_vcExplicitLabelBarLabelsOn\\|res_vcFillArrowEdgeColor\\|res_vcFillArrowEdgeThicknessF\\|res_vcFillArrowFillColor\\|res_vcFillArrowHeadInteriorXF\\|res_vcFillArrowHeadMinFracXF\\|res_vcFillArrowHeadMinFracYF\\|res_vcFillArrowHeadXF\\|res_vcFillArrowHeadYF\\|res_vcFillArrowMinFracWidthF\\|res_vcFillArrowWidthF\\|res_vcFillArrowsOn\\|res_vcFillOverEdge\\|res_vcGlyphOpacityF\\|res_vcGlyphStyle\\|res_vcLabelBarEndLabelsOn\\|res_vcLabelFontColor\\|res_vcLabelFontHeightF\\|res_vcLabelsOn\\|res_vcLabelsUseVectorColor\\|res_vcLevelColors\\|res_vcLevelCount\\|res_vcLevelPalette\\|res_vcLevelSelectionMode\\|res_vcLevelSpacingF\\|res_vcLevels\\|res_vcLineArrowColor\\|res_vcLineArrowHeadMaxSizeF\\|res_vcLineArrowHeadMinSizeF\\|res_vcLineArrowThicknessF\\|res_vcMagnitudeFormat\\|res_vcMagnitudeScaleFactorF\\|res_vcMagnitudeScaleValueF\\|res_vcMagnitudeScalingMode\\|res_vcMapDirection\\|res_vcMaxLevelCount\\|res_vcMaxLevelValF\\|res_vcMaxMagnitudeF\\|res_vcMinAnnoAngleF\\|res_vcMinAnnoArrowAngleF\\|res_vcMinAnnoArrowEdgeColor\\|res_vcMinAnnoArrowFillColor\\|res_vcMinAnnoArrowLineColor\\|res_vcMinAnnoArrowMinOffsetF\\|res_vcMinAnnoArrowSpaceF\\|res_vcMinAnnoArrowUseVecColor\\|res_vcMinAnnoBackgroundColor\\|res_vcMinAnnoConstantSpacingF\\|res_vcMinAnnoExplicitMagnitudeF\\|res_vcMinAnnoFont\\|res_vcMinAnnoFontAspectF\\|res_vcMinAnnoFontColor\\|res_vcMinAnnoFontHeightF\\|res_vcMinAnnoFontQuality\\|res_vcMinAnnoFontThicknessF\\|res_vcMinAnnoFuncCode\\|res_vcMinAnnoJust\\|res_vcMinAnnoOn\\|res_vcMinAnnoOrientation\\|res_vcMinAnnoOrthogonalPosF\\|res_vcMinAnnoParallelPosF\\|res_vcMinAnnoPerimColor\\|res_vcMinAnnoPerimOn\\|res_vcMinAnnoPerimSpaceF\\|res_vcMinAnnoPerimThicknessF\\|res_vcMinAnnoSide\\|res_vcMinAnnoString1\\|res_vcMinAnnoString1On\\|res_vcMinAnnoString2\\|res_vcMinAnnoString2On\\|res_vcMinAnnoTextDirection\\|res_vcMinAnnoZone\\|res_vcMinDistanceF\\|res_vcMinFracLengthF\\|res_vcMinLevelValF\\|res_vcMinMagnitudeF\\|res_vcMonoFillArrowEdgeColor\\|res_vcMonoFillArrowFillColor\\|res_vcMonoLineArrowColor\\|res_vcMonoWindBarbColor\\|res_vcNoDataLabelOn\\|res_vcNoDataLabelString\\|res_vcPositionMode\\|res_vcRefAnnoAngleF\\|res_vcRefAnnoArrowAngleF\\|res_vcRefAnnoArrowEdgeColor\\|res_vcRefAnnoArrowFillColor\\|res_vcRefAnnoArrowLineColor\\|res_vcRefAnnoArrowMinOffsetF\\|res_vcRefAnnoArrowSpaceF\\|res_vcRefAnnoArrowUseVecColor\\|res_vcRefAnnoBackgroundColor\\|res_vcRefAnnoConstantSpacingF\\|res_vcRefAnnoExplicitMagnitudeF\\|res_vcRefAnnoFont\\|res_vcRefAnnoFontAspectF\\|res_vcRefAnnoFontColor\\|res_vcRefAnnoFontHeightF\\|res_vcRefAnnoFontQuality\\|res_vcRefAnnoFontThicknessF\\|res_vcRefAnnoFuncCode\\|res_vcRefAnnoJust\\|res_vcRefAnnoOn\\|res_vcRefAnnoOrientation\\|res_vcRefAnnoOrthogonalPosF\\|res_vcRefAnnoParallelPosF\\|res_vcRefAnnoPerimColor\\|res_vcRefAnnoPerimOn\\|res_vcRefAnnoPerimSpaceF\\|res_vcRefAnnoPerimThicknessF\\|res_vcRefAnnoSide\\|res_vcRefAnnoString1\\|res_vcRefAnnoString1On\\|res_vcRefAnnoString2\\|res_vcRefAnnoString2On\\|res_vcRefAnnoTextDirection\\|res_vcRefAnnoZone\\|res_vcRefLengthF\\|res_vcRefMagnitudeF\\|res_vcScalarFieldData\\|res_vcScalarMissingValColor\\|res_vcScalarValueFormat\\|res_vcScalarValueScaleFactorF\\|res_vcScalarValueScaleValueF\\|res_vcScalarValueScalingMode\\|res_vcSpanLevelPalette\\|res_vcUseRefAnnoRes\\|res_vcUseScalarArray\\|res_vcVectorDrawOrder\\|res_vcVectorFieldData\\|res_vcWindBarbCalmCircleSizeF\\|res_vcWindBarbColor\\|res_vcWindBarbLineThicknessF\\|res_vcWindBarbScaleFactorF\\|res_vcWindBarbTickAngleF\\|res_vcWindBarbTickLengthF\\|res_vcWindBarbTickSpacingF\\|res_vcZeroFLabelAngleF\\|res_vcZeroFLabelBackgroundColor\\|res_vcZeroFLabelConstantSpacingF\\|res_vcZeroFLabelFont\\|res_vcZeroFLabelFontAspectF\\|res_vcZeroFLabelFontColor\\|res_vcZeroFLabelFontHeightF\\|res_vcZeroFLabelFontQuality\\|res_vcZeroFLabelFontThicknessF\\|res_vcZeroFLabelFuncCode\\|res_vcZeroFLabelJust\\|res_vcZeroFLabelOn\\|res_vcZeroFLabelOrthogonalPosF\\|res_vcZeroFLabelParallelPosF\\|res_vcZeroFLabelPerimColor\\|res_vcZeroFLabelPerimOn\\|res_vcZeroFLabelPerimSpaceF\\|res_vcZeroFLabelPerimThicknessF\\|res_vcZeroFLabelSide\\|res_vcZeroFLabelString\\|res_vcZeroFLabelTextDirection\\|res_vcZeroFLabelZone\\|res_vfCopyData\\|res_vfDataArray\\|res_vfExchangeDimensions\\|res_vfExchangeUVData\\|res_vfMagMaxV\\|res_vfMagMinV\\|res_vfMissingUValueV\\|res_vfMissingVValueV\\|res_vfPolarData\\|res_vfSingleMissingValue\\|res_vfUDataArray\\|res_vfUMaxV\\|res_vfUMinV\\|res_vfVDataArray\\|res_vfVMaxV\\|res_vfVMinV\\|res_vfXArray\\|res_vfXCActualEndF\\|res_vfXCActualStartF\\|res_vfXCEndIndex\\|res_vfXCEndSubsetV\\|res_vfXCEndV\\|res_vfXCStartIndex\\|res_vfXCStartSubsetV\\|res_vfXCStartV\\|res_vfXCStride\\|res_vfYArray\\|res_vfYCActualEndF\\|res_vfYCActualStartF\\|res_vfYCEndIndex\\|res_vfYCEndSubsetV\\|res_vfYCEndV\\|res_vfYCStartIndex\\|res_vfYCStartSubsetV\\|res_vfYCStartV\\|res_vfYCStride\\|res_vpAnnoManagerId\\|res_vpClipOn\\|res_vpHeightF\\|res_vpKeepAspect\\|res_vpOn\\|\\)\\>" 1 font-lock-constant-face)

    ("\\<\\(res_vpUseSegments\\|res_vpWidthF\\|res_vpXF\\|res_vpYF\\|res_wkAntiAlias\\|res_wkBackgroundColor\\|res_wkBackgroundOpacityF\\|res_wkColorMapLen\\|res_wkColorMap\\|res_wkColorModel\\|res_wkColorModel\\|res_wkDashTableLength\\|res_wkDefGraphicStyleId\\|res_wkDeviceLowerX\\|res_wkDeviceLowerX\\|res_wkDeviceLowerX\\|res_wkDeviceLowerY\\|res_wkDeviceLowerY\\|res_wkDeviceLowerY\\|res_wkDeviceUpperX\\|res_wkDeviceUpperX\\|res_wkDeviceUpperX\\|res_wkDeviceUpperY\\|res_wkDeviceUpperY\\|res_wkDeviceUpperY\\|res_wkFileName\\|res_wkFileName\\|res_wkFillTableLength\\|res_wkForegroundColor\\|res_wkFormat\\|res_wkFormat\\|res_wkFullBackground\\|res_wkFullBackground\\|res_wkGksWorkId\\|res_wkHeight\\|res_wkMarkerTableLength\\|res_wkMetaName\\|res_wkOrientation\\|res_wkOrientation\\|res_wkOrientation\\|res_wkPDFFileName\\|res_wkPDFFormat\\|res_wkPDFResolution\\|res_wkPSFileName\\|res_wkPSFormat\\|res_wkPSResolution\\|res_wkPaperHeightF\\|res_wkPaperHeightF\\|res_wkPaperHeightF\\|res_wkPaperSize\\|res_wkPaperSize\\|res_wkPaperSize\\|res_wkPaperWidthF\\|res_wkPaperWidthF\\|res_wkPaperWidthF\\|res_wkTopLevelViews\\|res_wkViews\\|res_wkVisualType\\|res_wkVisualType\\|res_wkWidth\\|res_wsCurrentSize\\|res_wsMaximumSize\\|res_wsThresholdSize\\|res_xyComputeXMax\\|res_xyComputeXMin\\|res_xyComputeYMax\\|res_xyComputeYMin\\|res_xyCoordData\\|res_xyCoordDataSpec\\|res_xyCurveDrawOrder\\|res_xyDashPattern\\|res_xyDashPatterns\\|res_xyExplicitLabels\\|res_xyExplicitLegendLabels\\|res_xyLabelMode\\|res_xyLineColor\\|res_xyLineColors\\|res_xyLineDashSegLenF\\|res_xyLineLabelConstantSpacingF\\|res_xyLineLabelFont\\|res_xyLineLabelFontAspectF\\|res_xyLineLabelFontColor\\|res_xyLineLabelFontColors\\|res_xyLineLabelFontHeightF\\|res_xyLineLabelFontQuality\\|res_xyLineLabelFontThicknessF\\|res_xyLineLabelFuncCode\\|res_xyLineOpacities\\|res_xyLineOpacityF\\|res_xyLineThicknessF\\|res_xyLineThicknesses\\|res_xyMarkLineMode\\|res_xyMarkLineModes\\|res_xyMarker\\|res_xyMarkerColor\\|res_xyMarkerColors\\|res_xyMarkerOpacities\\|res_xyMarkerOpacityF\\|res_xyMarkerSizeF\\|res_xyMarkerSizes\\|res_xyMarkerThicknessF\\|res_xyMarkerThicknesses\\|res_xyMarkers\\|res_xyMonoDashPattern\\|res_xyMonoLineColor\\|res_xyMonoLineLabelFontColor\\|res_xyMonoLineThickness\\|res_xyMonoMarkLineMode\\|res_xyMonoMarker\\|res_xyMonoMarkerColor\\|res_xyMonoMarkerSize\\|res_xyMonoMarkerThickness\\|res_xyXIrrTensionF\\|res_xyXIrregularPoints\\|res_xyXStyle\\|res_xyYIrrTensionF\\|res_xyYIrregularPoints\\|res_xyYStyle\\|\\)\\>" 1 font-lock-constant-face)

    ;; ESMValTool interface, shared and utility scripts
    ("\\<\\(read_data\\|select_metadata_by_atts\\|select_metadata_by_name\\|metadata_att_as_array\\|bname\\|basename\\|att2var\\|att2var_default\\|get_ncdf_name\\|get_ncdf_dir\\|ncdf_read\\|ncdf_define\\|ncdf_write\\|ncdf_att\\|copy_CoordNames_n\\|extend_var_at\\|remove_index\\|set_default_att\\|empty_str\\|log_info\\|log_debug\\|enter_msg\\|leave_msg\\|error_msg\\|tool_stop\\|exit_if_missing_atts\\|log_provenance\\|taylor_plot\\|contour_map\\|contour_map_polar\\|contour_map_ce\\|add_markers_to_map\\|get_title_suffix\\|remove_attrs\\|plot_two_by_one\\|plot_three_by_one_diff\\|two_by_one\\|three_by_one_diff\\|plot_three_by_one_vector\\|three_by_one_vector\\|plot_multipanel\\|multipanel\\|plot_multipanel_vector\\|multipanel_vector\\|seasonal_plot\\|xy_plot_wrapper\\|ts_line_wrapper\\|pr_u850_mean_plot\\|mjo_xcor_lag_plot\\|mjo_pr_ua_vari_plot\\|mjo_unvari_eof_plot\\|get_title_suffix\\|remove_attrs\\|plot_two_by_one\\|plot_three_by_one_diff\\|two_by_one\\|three_by_one_diff\\|plot_three_by_one_vector\\|three_by_one_vector\\|plot_multipanel\\|multipanel\\|plot_multipanel_vector\\|multipanel_vector\\|seasonal_plot\\|xy_plot_wrapper\\|ts_line_wrapper\\|xy_line_overlap\\|plot_precip_domain\\|precip_domain\\|month_sel\\|lat_names\\|add_line\\|add_scatt\\|add_legend\\|calcRegCoeffs\\|genZonalMeans\\|calcMeanAnnCycleMonthly\\|calcMeanAnnCycleAnnual\\|rmMeanAnnCycle\\|apfiltersmooth\\|smoothAnomalies\\|clmMon2clmDayn\\|scatterplot\\|scatterplot3D\\|scatterplot_markers\\|zonalmean_profile\\|contourplot\\|portrait_plot\\|circle_plot\\|profile_plev\\|aerosol_profile\\|aerosol_sizedist\\|xy_line\\|xy_line_anom\\|timeseries_station\\|cycle_plot\\|errorbar_plot\\|create_legend_lines\\|output_type\\|copy_VarAtt_sel\\|panelling\\|get_plot_dir\\|get_outfile_name\\|get_wks\\|add_markers\\|add_num_markers\\|add_errorbar\\|horizontal_whiskers\\|add_prediction_error\\|mjo_wave_freq_plot\\|addHorVertLinesCross_extended\\|mjo_cross_spectra_plot\\|mjo_ceof_plot\\|mjo_life_cycle_plot\\|vector_scalar_map_polar\\|project_style\\|place_debuginfo\\|place_description\\|gsnColorRange\\|format_units\\|set_log_ticks\\|sort_alphabetically\\|legend_lines\\|legend_markers\\|roi\\|extract_area\\|gridcell_area\\|map_area\\|area_operations\\|select_region\\|make_latlon2D\\|cdo_remapdis\\|guestimate_average_grid_area\\|get_lower_limits\\|get_upper_limits\\|is_regional\\|esmf_conserve_wrapper\\|rect2rect_interp\\|plev_lat_interp\\|get_dataset_minus_ref\\|esmf_conserve_wrapper_time\\|regrid_3D_to_rectilinear_grid\\|get_start_year\\|get_end_year\\|convert_units\\|UNIQ\\|union\\|set_inclusive_OR\\|intersection\\|is_array_subset\\|relative_complement\\|set_symmetric_difference\\|dim_stddev_wgt_Wrap\\|time_operations\\|calc_season_index\\|extract_season\\|month_to_season_extended\\|coswgt_areaave\\|coswgt_arearmse\\|coswgt_pattern_cor\\|interannual_variability\\|calculate_metric\\|normalize_metric\\|distrib_stats\\|lognormal_dist\\|add_labelbar\\|create_empty_array\\|data_read_in\\|data_read_in_ocean_MOC\\|data_read_in_ice\\|y_axis_check\\|check_custom_climo\\|isfilepresent2\\|table_link_setup\\|set_varAtts\\|create_timec\\|format_time\\|format_plev\\|format_lev\\|format_lat\\|format_lon\\|format_coords\\|read_cmor\\|format_variable\\|guess_bounds_time\\|guess_bounds_lev\\|guess_bounds_lat\\|guess_bounds_lon\\|guess_coord_bounds\\|set_global_atts\\|write_nc\\|write_nc_profile\\|set_size_array\\|process_EBAS_data\\|\\)\\>" 1 font-lock-type-face)

    ) 
  "words used in ncl-mode highlighting"
  )

(put 'ncl-mode 'font-lock-defaults 'ncl-font-lock-keywords)
;;************************************************
;; some variables used in the creation of ncl-mode
;;************************************************
(defvar ncl-mode-map nil
  "Keymap used in NCL mode.")
(defvar ncl-startup-message t
  "*Non-nil displays a startup message when `ncl-mode' is first called.")
(defconst ncl-mode-version "0.32")
;;************************************************
;; syntax table
;;************************************************
;; characters are preceeded by a ?.
;; "." indicates punctuation
;; "_" indicates a symbol
;; "\"" indicates a string (must escape the double-quote)
;; "<" indicates a comment
;; "w" indicates a word character

(defvar ncl-mode-syntax-table nil
  "Syntax table in use in `ncl-mode' buffers.")
(if ncl-mode-syntax-table ()
  (setq ncl-mode-syntax-table (make-syntax-table))
  (modify-syntax-entry ?\;  "<"  ncl-mode-syntax-table)
  (modify-syntax-entry ?+   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?-   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?*   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?/   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?^   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?#   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?=   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?%   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?<   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?>   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?\'  "\"" ncl-mode-syntax-table)
  (modify-syntax-entry ?\"  "\"" ncl-mode-syntax-table)
  (modify-syntax-entry ?\\  "." ncl-mode-syntax-table)
  (modify-syntax-entry ?_   "w"  ncl-mode-syntax-table)
  (modify-syntax-entry ?{   "\(}"  ncl-mode-syntax-table)
  (modify-syntax-entry ?}   "\){"  ncl-mode-syntax-table)
  (modify-syntax-entry ?$   "_"  ncl-mode-syntax-table)
  (modify-syntax-entry ?.   "."  ncl-mode-syntax-table)
  (modify-syntax-entry ?\n  ">"  ncl-mode-syntax-table)
  (modify-syntax-entry ?\f  ">"  ncl-mode-syntax-table))

(defvar ncl-find-symbol-syntax-table nil
  "Syntax table that treats symbol characters as word characters.")

(if ncl-find-symbol-syntax-table ()
  (setq ncl-find-symbol-syntax-table
        (copy-syntax-table ncl-mode-syntax-table))
  )
;;****************************************************************************
;; keymap
;;****************************************************************************
(defvar ncl-mode-map nil
  "Keymap used in NCL mode.")
(if ncl-mode-map ()
  (setq ncl-mode-map (make-sparse-keymap))
  (define-key ncl-mode-map "\t"       'ncl-indent-line))
;;****************************************************************************
;; indenting variables
;;****************************************************************************
(defvar ncl-main-block-indent 2
  "*Extra indentation for the main block of code. That is the block between 
the begin statement and the end statement.")

(defvar ncl-main-block-end -2
  "*The offset that places the end statement back on the left margin. This is
the negative of `ncl-main-block-indent`")

(defvar ncl-block-indent 2
  "*Extra indentation for do loops.")

(defvar ncl-block-end -2
  "*The offset that places the `end do` statement back to it origin.")

(defconst ncl-comment-line-start-skip "^[ \t]*;"
  "Regexp to match the start of a full-line comment. That is the 
_beginning_ of a line containing a comment delmiter `\;' preceded 
only by whitespace.")

;; defconst are constants that never change 
;; the \` matches only those at the beginning of the buffer and no other
(defconst ncl-begin "\\<\\(begin\\)\\>\\|\\`" 
  "Regular expression to find the begin statement.")

;; the \' matches only those at the end of the buffer and no other
(defconst ncl-end "\\<\\(^end$\\)\\>\\|\\'" 
  "Regular expression to find the line that indicates the end of a
script.")

(defconst ncl-begin-do "^[ /t]*do" 
  "Regular expression to find the beginning of a do loop.")

(defconst ncl-else "^[ /t]*else" 
  "Regular expression to find an else statment.")

(defconst ncl-begin-if "^[ /t]*if" 
  "Regular expression to find the beginning of a if statment.")

(defconst ncl-enddo "end[ ]do"
  "Regular expression to find the end of a do loop")

(defconst ncl-endif "end[ ]if"
  "Regular expression to find the end of a if statement")

(defconst ncl-identifier "[a-zA-Z][a-zA-Z0-9$_]+"
  "Regular expression matching an NCL identifier.")

(defconst ncl-label (concat ncl-identifier ":")
  "Regular expression matching NCL labels.")

(defvar ncl-no-change-comment ";;"
  "*The indentation of a comment that starts with this regular
expression will not be changed. Note that the indentation of a comment
at the beginning of a line is never changed.")

(defvar ncl-begin-line-comment nil
  "*A comment anchored at the beginning of line.
A comment matching this regular expression will not have its
indentation changed.  If nil the default is \"^\;\", i.e., any line
beginning with a \"\;\".  Expressions for comments at the beginning of
the line should begin with \"^\".")

(defvar ncl-code-comment ";;[^;]"
  "*A comment that starts with this regular expression on a line by
itself is indented as if it is a part of NCL code.  As a result if
the comment is not preceded by whitespace it is unchanged.")
;;****************************************************************************
;; indenting functions
;;****************************************************************************
(defun ncl-beginning-of-statement ()
  "Move to beginning of the current statement. Skips back past statement 
continuations. Point is placed at the beginning of the line whether or not 
this is an actual statement."
  (if (save-excursion (forward-line -1) (ncl-is-continuation-line))
      (ncl-previous-statement)
    (beginning-of-line)))

(defun ncl-end-of-statement ()
  "Moves point to the end of the current NCL statement. If not in a statement 
just moves to end of line. Returns position."
  (interactive)
  (while (and (ncl-is-continuation-line)
              (= (forward-line 1) 0)))
  (end-of-line) (point))

(defun ncl-previous-statement ()
  "Moves point to beginning of the previous statement. Returns t if the 
current line before moving is the beginning of the first non-comment 
statement in the file, and nil otherwise."
  (interactive)
  (let (first-statement)
    (if (not (= (forward-line -1) 0))
        ;; first line in file
        t
      ;; skip blank lines, label lines, include lines and line comments
      (while (and 
              ;; The current statement is the first statement until we
              ;; reach another statement.
              (setq first-statement
                    (or 
                     (looking-at ncl-comment-line-start-skip)
                     (looking-at "[ \t]*$")
                     (looking-at (concat "[ \t]*" ncl-label "[ \t]*$"))
                     (looking-at "^@")))
              (= (forward-line -1) 0)))
      ;; skip continuation lines
      (while (and 
              (save-excursion
                (forward-line -1)
                (ncl-is-continuation-line))
              (= (forward-line -1) 0)))
      first-statement)))

(defun ncl-is-continuation-line ()
  "Tests if current line is continuation line."
  (save-excursion
    (ncl-look-at "\\<\\$")))

(defun ncl-look-at (regexp &optional cont beg)
  "Searches current line from current point for the regular expression REGEXP.
If optional argument CONT is non-nil, searches to the end of the current 
statement. If optional arg BEG is non-nil, search starts from the beginning 
of the current statement. Ignores matches that end in a comment or inside a 
string expression. Returns point if successful, nil otherwise.  This function 
produces unexpected results if REGEXP contains quotes or a comment delimiter. 
The search is case insensitive.  If successful leaves point after the match, 
otherwise, does not move point."
  (let ((here (point))
        (old-syntax-table (syntax-table))
        (case-fold-search t)
        eos
        found)
    (set-syntax-table ncl-find-symbol-syntax-table)
    (setq eos
          (if cont
              (save-excursion (ncl-end-of-statement) (point))
            (save-excursion (end-of-line) (point))))
    (if beg (ncl-beginning-of-statement))
    (while (and (setq found (re-search-forward regexp eos t))
                (ncl-quoted)))
    (set-syntax-table old-syntax-table)
    (if (not found) (goto-char here))
    found))
 
(defun ncl-in-quote ()
  "Returns location of the opening quote if point is in a NCL string constant,
nil otherwise. Ignores comment delimiters on the current line. Properly 
handles nested quotation marks and octal constants - a double quote followed 
by an octal digit."
;;; Treat an octal inside an apostrophe to be a normal string. Treat a
;;; double quote followed by an octal digit to be an octal constant
;;; rather than a string. Therefore, there is no terminating double
;;; quote.
  (save-excursion
    ;; Because single and double quotes can quote each other we must
    ;; search for the string start from the beginning of line.
    (let* ((start (point))
           (eol (progn (end-of-line) (point)))
           (bq (progn (beginning-of-line) (point)))
           (endq (point))
           (data (match-data))
           delim
           found)
          (while  (< endq start)
            ;; Find string start
            ;; Don't find an octal constant beginning with a double quote
            (if (re-search-forward "\"[^0-7]\\|'\\|\"$" eol 'lim)
                ;; Find the string end. In NCL, two consecutive delimiters 
		;; after the start of a string act as an escape for the 
                ;; delimiter in the string. Two consecutive delimiters alone 
		;; (i.e., not after the start of a string) is the the 
		;; null string.
                (progn 
                  ;; Move to position after quote
                  (goto-char (1+ (match-beginning 0)))
                  (setq bq (1- (point)))
                  ;; Get the string delimiter
                  (setq delim (char-to-string (preceding-char)))
                  ;; Check for null string
                  (if (looking-at delim)
                      (progn (setq endq (point)) (forward-char 1))
                    ;; Look for next unpaired delimiter
                    (setq found (search-forward delim eol 'lim))
                    (while (looking-at delim)
                      (forward-char 1)
                      (setq found (search-forward delim eol 'lim)))
                    (if found
                        (setq endq (- (point) 1))
                      (setq endq (point)))
                    ))
              (progn (setq bq (point)) (setq endq (point)))))
          (store-match-data data)
      ;; return string beginning position or nil
      (if (> start bq) bq))))


(defun ncl-quoted ()
  "Returns t if point is in a comment or quoted string. nil otherwise."
;  (or (ncl-in-comment) (ncl-in-quote)))
  (or (ncl-in-quote)))

(defun ncl-in-comment ()
  "Returns t if point is inside a comment, nil otherwise."
  (save-excursion
    (let ((here (point)))
      (and (ncl-goto-comment) (> here (point))))))

(defun ncl-goto-comment ()
  "Move to start of comment delimiter on current line. Moves to end of line if
there is no comment delimiter. Ignores comment delimiters in strings. Returns 
point if comment found and nil otherwise."
  (let ((eos (progn (end-of-line) (point)))
        (data (match-data))
        found)
    ;; Look for first comment delimiter not in a string
    (beginning-of-line)
    (setq found (search-forward comment-start eos 'lim))
    (while (and found (ncl-in-quote))
      (setq found (search-forward comment-start eos 'lim)))
    (store-match-data data)
    (and found (not (ncl-in-quote))
         (progn
           (backward-char 1)
           (point)))))

(defun ncl-current-statement-indent ()
  "Return indentation of the current statement. If in a statement, moves to 
beginning of statement before finding indent."
  (ncl-beginning-of-statement)
  (ncl-current-indent))

(defun ncl-current-indent ()
  "Return the column of the indentation of the current line.  Skips any 
whitespace. Returns 0 if the end-of-line follows the whitespace."
  (save-excursion
    (beginning-of-line)
    (skip-chars-forward " \t")
    ;; if we are at the end of blank line return 0
    (cond ((eolp) 0)
          ((current-column)))))

(defun ncl-calculate-indent ()
  "Return appropriate indentation for current line as NCL code."
  (save-excursion
    (beginning-of-line)
    (cond 
     ;; if line is "begin" do nothing and exit
     ((ncl-look-at ncl-begin) 0)
     ;; calculate indent based on previous and current statements
     (t (let ((the-indent
	      ;; calculate indent based on previous statement
	      (save-excursion
		(cond
		 ;; retreive the previous statement
		 ( (ncl-previous-statement) 0)

		 ;; indent if previous statment is begin
		 ((ncl-look-at ncl-begin t)
		  (+ (ncl-current-statement-indent) ncl-main-block-indent))
		 
		 ;; indent if previous statment is do 
		 ((ncl-look-at ncl-begin-do t)
		  (+ (ncl-current-statement-indent) ncl-block-indent))

		 ;; indent if previous statment is if 
		 ((ncl-look-at ncl-begin-if t)
		  (+ (ncl-current-statement-indent) ncl-block-indent))

		 ;; indent if previous statment is else
		 ((ncl-look-at ncl-else t)
		  (+ (ncl-current-statement-indent) ncl-block-indent))

		 ((ncl-current-statement-indent))))))
	  ;; adjust the indentation based on the current statement
	  (cond
	   ;; do loop
	   ((ncl-look-at ncl-enddo t)
	    (+ the-indent ncl-block-end))
	   ;; if statement
	   ((ncl-look-at ncl-endif t)
	    (+ the-indent ncl-block-end))
	   ;; else statement
	   ((ncl-look-at ncl-else t)
	    (+ the-indent ncl-block-end))

	   ;; End block
	   ((ncl-look-at ncl-end t)
	    (+ the-indent ncl-main-block-end)) ;; end gets negative indent
	   (the-indent))

	  )))))

(defun ncl-indent-to (col &optional min)
  "Indent from point with spaces until column COL. Inserts space before 
markers at point."
  (if (not min) (setq min 0))
  (insert-before-markers
   (make-string (max min (- col (current-column))) ? )))

(defun ncl-indent-left-margin (col)
  "Indent the current line to column COL. Indents such that first 
non-whitespace character is at column COL. Inserts spaces before markers at 
point."
  (save-excursion
    (beginning-of-line)
    (delete-horizontal-space)
    (ncl-indent-to col)))

(defun ncl-comment-hook ()
  "Compute indent for the beginning of the NCL comment delimiter."
  (if (or (looking-at ncl-no-change-comment)
          (if ncl-begin-line-comment
              (looking-at ncl-begin-line-comment)
              (looking-at "^\;")))
      (current-column)
    (if (looking-at ncl-code-comment)
        (if (save-excursion (skip-chars-backward " \t") (bolp))
            ;; On line by itself, indent as code
            (let ((tem (ncl-calculate-indent)))
              (if (listp tem) (car tem) tem))
          ;; after code - do not change
          (current-column))
      (skip-chars-backward " \t")
      (max (if (bolp) 0 (1+ (current-column)))
           comment-column))))

(defun ncl-indent-line ()
  "Indents current NCL line as code or as a comment."
  (interactive)
  ;; Move point out of left margin.
  (if (save-excursion
        (skip-chars-backward " \t")
        (bolp))
      (skip-chars-forward " \t"))
  (let ((mloc (point-marker)))
    (save-excursion
      (beginning-of-line)
      (if (looking-at ncl-comment-line-start-skip)
          ;; Indentation for a line comment
          (progn
            (skip-chars-forward " \t")
            (ncl-indent-left-margin (ncl-comment-hook)))
        ;; Indent for code line
        (beginning-of-line)
        (if (or
             ;; a label line
             (looking-at (concat "^" ncl-label "[ \t]*$"))
             ;; a batch command
             (looking-at "^[ \t]*@"))
            ;; leave flush left
            nil
          ;; indent the line
          (ncl-indent-left-margin (ncl-calculate-indent)))
        ;; Adjust parallel comment
;        (end-of-line) 
;        (if (ncl-in-comment)
;            (indent-for-comment))
	))
    (goto-char mloc)
    ;; Get rid of marker
    (set-marker mloc nil)
    ))

;;****************************************************************************
;; the command to comment/uncomment text
;;****************************************************************************
(defun ncl-comment-dwim (arg)
"Comment or uncomment current line or region in a smart way.
For detail, see `comment-dwim'."
   (interactive "*P")
   (require 'newcomment)
   (let ((deactivate-mark nil) (comment-start ";") (comment-end ""))
     (comment-dwim arg)))

;;****************************************************************************
;; define ncl mode
;;****************************************************************************
(defun ncl-mode ()
  "Major mode for editing NCL .ncl files"
  (interactive)
  (kill-all-local-variables)
  (setq major-mode 'ncl-mode)
  (setq mode-name "NCL")

  ;; modify the keymap
  (define-key ncl-mode-map [remap comment-dwim] 'ncl-comment-dwim)

  (if ncl-startup-message
      (message "Emacs NCL mode version %s." ncl-mode-version)
    ) 
;**************************
;; indentation
;**************************
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'ncl-indent-line)
  (use-local-map ncl-mode-map)
;**************************
;; these ensure syntax hightlighting
;**************************
;; font-lock setup for various emacs: XEmacs, Emacs 19.29+, Emacs <19.29.
;; taken from html-helper-mode, adapted to ncl
  (cond	((string-match "XEmacs\\|Lucid" (emacs-version)) ; XEmacs/Lucid
	 (put major-mode 'font-lock-keywords-case-fold-search t)
	 (put major-mode 'font-lock-syntactic-keywords t)
         (put major-mode 'font-lock-maximum-decoration 2)
	 )
        ;; not sure if this is correct
	;; XEmacs (19.13, at least) guesses the rest correctly.
	;; If any older XEmacs don't, then tell me.
	;;
	((string-lessp "19.28.89" emacs-version) ; Emacs 19.29 and later
	 (make-local-variable 'font-lock-defaults)
	 (setq font-lock-defaults '(ncl-font-lock-keywords nil t)))
	;;
	(t ; Emacs 19.28 and older
	 (make-local-variable 'font-lock-keywords-case-fold-search)
	 (make-local-variable 'font-lock-keywords)
	 ;;(make-local-variable 'font-lock-no-comments)
	 (setq font-lock-keywords-case-fold-search t)
	 (setq font-lock-keywords ncl-font-lock-keywords)
	 ;;(setq font-lock-no-comments t)
         ))

  (font-lock-mode 1)
;  (setq font-lock-maximum-decoration t)
;  (make-local-variable 'font-lock-defaults)
;  (setq font-lock-defaults 'ncl-keywords)
;  (make-local-variable 'comment-start)
;  (setq comment-start ";")
; turn this on if debuging this code
   (setq debug_on_error t)
  (set-syntax-table ncl-mode-syntax-table)
  (run-hooks 'ncl-mode-hook)
  )
;;************************************************************************
  (provide 'ncl)
