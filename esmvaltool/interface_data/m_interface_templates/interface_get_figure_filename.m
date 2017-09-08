function figure_name = interface_get_figure_filename(diag_script_base,...
                                                             variable,...
                                                           field_type,...
                                                             aux_info,...
                                                              idx_mod)
%                                                   return val [1] : string
%   Arguments:
%           @brief Construct a figure output file name
%           @param diags_script_base  -  The running diag script running withouth its suffix
%         @param variable  -  Current variable
%         @param field_type  -  Current field type
%         @param aux_info  -  User supplied info to put in figure filename
%         @param idx_mod  -  Current model number
    load('m_interface.mat');
    aux_sep = '_';  % Skip auxiliary separator if no aux_info
    if (length(aux_info) == 0)
        aux_sep = '';
    end

    if (idx_mod == -1)
        figure_name = [diag_script_base, '_', variable,...
                                         '_', field_type,...
                                         aux_sep, aux_info];
    else
        figure_name = [diag_script_base, '_', variable,...
                                         '_', field_type,...
                                         '_', aux_info, aux_sep,...
                                         figfiles_suffix{idx_mod}];
        figure_name = strrep(figure_name, '${MIP}', models_mip{idx_mod});
        figure_name = strrep(figure_name, '${EXP}', models_experiment{idx_mod});
    end
end
