; Settings used by most EAMonsoon plot_types
plot_type_info = True
plot_type_info@season = "JJAS"
plot_type_info@ref_model = (/"TRMM-3B43-v7-0.25deg"/)

; Settings specific for the precip basic plot_type
plot_type_info@latrange_basic = (/-20, 40/)
plot_type_info@lonrange_basic = (/30, 150/)

plot_type_info@cn_levels_mean_basic = (/2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5/)
plot_type_info@cn_levels_mean_diff_basic = (/-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15/)
plot_type_info@cn_levels_stddev_basic = (/0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0/)
plot_type_info@cn_levels_stddev_diff_basic = (/-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7/)


; Settings specific for the precip seasonal plot_type
plot_type_info@latrange_seasonal = (/5, 30/)
plot_type_info@lonrange_seasonal = (/65, 95/)


; Settings shared by the precip global domain/intensity plot_type
plot_type_info@summer_season = "MJJAS"
plot_type_info@winter_season = "NDJFM"

; Settings specific for the precip global domain plot_type
plot_type_info@latrange_domain = (/-50, 50/)

plot_type_info@cn_levels_mean_domain = (/0.01, 1.01/)
plot_type_info@cn_levels_mean_diff_domain = (/-0.99, 0.01, 1.01/)

plot_type_info@colors_domain = (/"white", "blue"/)
plot_type_info@colors_diff_domain = (/"red", "white", "blue"/)

plot_type_info@high_intensity_cutoff_mmday = 2.0


; Settings specific for the precip intensity plot_type
plot_type_info@latrange_intensity = (/-50, 50/)
plot_type_info@cn_levels_mean_intensity = (/-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0/)
plot_type_info@cn_levels_mean_diff_intensity = (/-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0/)
