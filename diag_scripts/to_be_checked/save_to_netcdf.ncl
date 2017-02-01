;
;  Save a derived field to a netcdf file
;
;  2008/11/26   C. Fischer   created
;  2009/01/27   CAF      changed output filename and climo directory structure
;  2009/06/24   CAF      added date range to output filename
;
load "interface_data/filesinfo"
load "interface_data/modelinfo"
begin
    orig_field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_dir = getenv("ESMValTool_plot_dir")
    wkdir = plot_dir + "/save_to_netcdf/"
    print("save_to_netcdf: Saving to netcdf file" + " " \
               + orig_field_number + " " + variable)
    data_new = True
    do i = 0, dimsizes(models@dir) - 1, 1
        filepath = interface_get_infiles(variable, orig_field_number, i)
        filename = systemfunc("basename " + filepath)
        out_filepath = wkdir + filename
        error = systemfunc("rm -f " + out_filepath)
        fout = addfile(out_filepath, "c")
        A = extract_data(i, data, -1, 0, 0)
        fout->$variable$ = A
        delete(A)
    end do
end
