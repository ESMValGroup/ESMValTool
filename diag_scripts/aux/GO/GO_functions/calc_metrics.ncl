load "./interface_data/ncl.interface"

undef("calc_metrics")
procedure calc_metrics(models, dim_MOD, refModel, iref, aYdata)
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering calc_metrics", verbosity, 4)
    ;; logic:  test 3 metrics.
    ;; set up arrays
    gm = new((/dim_MOD/), "float")
    gc = new((/dim_MOD/), "float")
    gv = new((/dim_MOD/), "float")
    gsum = new((/dim_MOD/), "float")
    gradetx = new((/dim_MOD + 1/), "string")

    ng = 3.

    ;; select 'obs' (as a model)
    do imod = 0, dim_MOD - 1
        info_output("Models are: " + models@name(imod) + refModel, verbosity, 4)
        if (models@name(imod) .eq. refModel) then
            iref = imod
        end if
    end do

    ;; loop for models (also calculate for 'obs' should be perfect)
    ref_avg = avg(aYdata(iref, :))

    ref_anom = aYdata(iref, :) - ref_avg

    do imod = 0, dim_MOD - 1
        ;; remove mean
        tmp = avg(aYdata(imod, :))
        ts = aYdata(imod, :) - tmp
        delete(tmp)

        ;; A: 'mean' following waugh & eyring
        ;; or select an acceptable 'delta'
        dd = ng

        gm(imod) = 1. - abs(avg(aYdata(imod, :)) - ref_avg) / stddev(ref_anom) / dd

        if (gm(imod) .lt. 0.) then
            gm(imod) = 0.
        end if

        ;; B: "correlation"
        ;; need to be able to match years!
        gg = ind(.not. ismissing(ref_anom))
        gc(imod) = (esccr(ts(gg), ref_anom(gg), 0) + 1.) / 2.

        ;; ; C: "variance"
        gv(imod) = 1. - abs((stddev(ts) - stddev(ref_anom))) / stddev(ref_anom) / ng
        if (gv(imod) .lt. 0.) then
            gv(imod) = 0.
        end if

        gsum(imod) = (gm(imod) + gc(imod) + gv(imod)) / 3.
        gradetx(imod) = sprintf("%5.2f", gm(imod)) + " " + sprintf("%5.2f", gc(imod)) \
                                                   + " " + sprintf("%5.2f", gv(imod)) + " " \
                                                   + sprintf("%5.2f", gsum(imod)) + " "\
                                                   + models@name(imod) + " " \
                                                   + models@experiment(imod)
    end do
    gradetx(dim_MOD) = " GM   " + " GC   " + " GV   " + " GSUM " + " Model" + "Experiment"
    info_output(">>>>>>>> Leaving calc_metrics", verbosity, 4)
end
