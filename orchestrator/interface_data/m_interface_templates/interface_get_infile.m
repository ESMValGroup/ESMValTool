function curr_file = interface_get_infile(variable, field, idx)
%                        return val [1] : string
% Arguments:
%         @brief Reconstructs the current (idx) input filename
%         @param diag_script_base  -  The running the diag script running withouth its suffix
%         @param variable  -  Current variable
%         @param field_type  -  Current field type
    load('m_interface.mat');
    curr_file = strrep(fullpaths{idx}, '${FIELD}', field);
    curr_file = strrep(curr_file, '${VARIABLE}', variable);
    curr_file = strrep(curr_file, '${MIP}', models_mip{idx});
    curr_file = strrep(curr_file, '${EXP}', models_experiment{idx});
end
