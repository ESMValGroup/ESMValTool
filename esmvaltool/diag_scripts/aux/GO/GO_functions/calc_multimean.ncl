load "./interface_data/ncl.interface"

undef("calc_multimean")
function calc_multimean(models, Xdata, Ydata,  dim_MOD, nMax, dimX)
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering calc_multimean", verbosity, 4)
    tmpn = models@name
    delete(models@name)
    models@name = new((/dim_MOD + 1/), "string")
    models@name(0) = "MEAN"
    models@name(1:dim_MOD) = tmpn
    delete(tmpn)

    tmpx = Xdata
    tmpy = Ydata

    ;; create bigger arrays (mean goes in index 0)  and reload
    Xdata1 = new((/dim_MOD + 1, nMax/), "float")
    Ydata1 = new((/dim_MOD + 1, nMax/), "float")
    Xdata1(1:dim_MOD, :) = (/tmpx(:, :)/)
    Ydata1(1:dim_MOD, :) = (/tmpy(:, :)/)

    ;; copy X array from model 0
    Xdata1(0, :) = tmpx(0, :)
    ;; average Y array for each time...
    do it = 0, dimX(0) - 1
        tmpx1d = ndtooned(tmpx)
        tmpy1d = ndtooned(tmpy)
        ;; Calculate average of all models at time it
        Ydata1(0, it) = avg(tmpy1d(ind(tmpx1d .eq. tmpx(0, it))))
    end do

    dim_MOD = dim_MOD + 1
    info_output(">>>>>>>> Leaving calc_multimean", verbosity, 4)
    return ([/Xdata1, Ydata1/])
end
