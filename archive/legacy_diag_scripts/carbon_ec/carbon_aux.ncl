; #############################################################################
; # EXTERNAL CARBON FUNCTIONS
; #############################################################################
;  Please consider using rf extending existing routines before adding new ones.
;  Check the header of each routine for documentation.
;
;  Contents:
;     function cumulate_carbon: cummulates the variable over time
;     function Cond_PDF: calculates the normal and conditional PDF
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/latlon.ncl"

; #############################################################################
undef("cumulate_carbon")
function cumulate_carbon(NBP[*][*][*]:numeric)

;  Arguments:
;     NBP            : carbon flux variable
;
;  Return value:
;     array containing cumulated carbon fluxes
;
;  Description:
;     Calculates cumulated carbon fluxes
;
;  Caveats
;
;  Modification history:
;     20191118-zechlau_sabrina: written.
;
local funcname, scriptname, verbosity, idx_new, NBP_tmp, NBP, y

begin

  funcname = "cumulate_carbon"
  scriptname = "diag_scripts/carbon_ec/carbon_aux.ncl"
  enter_msg(scriptname, funcname)

    ; Integrate
    idx_new = dimsizes(NBP)
    NBP_tmp = new(dimsizes(NBP), typeof(NBP))

    if (dimsizes(idx_new).eq.3) then
      NBP_tmp(0, :, :) = 0.0
      do y = 0, dimsizes(NBP(:, 0, 0)) - 2
        NBP_tmp(y + 1, :, :) = NBP_tmp(y, :, :) + NBP(y, :, :)
      end do

      copy_VarAtts(NBP, NBP_tmp)
      copy_VarCoords(NBP, NBP_tmp)

    elseif (dimsizes(idx_new).eq.1) then
      NBP_tmp(0) = 0.0
      do y = 0, dimsizes(NBP) - 2
        NBP_tmp(y + 1) = NBP_tmp(y:) + NBP(y)
      end do
      copy_VarAtts(NBP_tmp, NBP)
      copy_VarCoords(NBP_tmp, NBP)

    else
      error_msg("f", scriptname, funcname, \
                "rank " + dimsizes(idx_new) + " of NBP not supported")
    end if
    delete(y)

  return(NBP_tmp)
end

; #######################################################################
undef("Cond_PDF")
function Cond_PDF(xfit[*]:numeric,
                  yfit[*]:numeric,
                  Xval[*]:numeric,
                  Yval[*]:numeric,
                  MSE[*]:numeric,
                  obs_std[*]:numeric
                  )

; Arguments:
;  * xfit    : 1-D array of the X variable (from regressionline y=a*x+b)
;  * yfit    : 1-D array of the Y variable (from regressionline y=a*x+b)
;  * Xval    : 1-D array of the X variable (from models)
;  * Yval    : 1-D array of the Y variable (from models)
;  * MSE     : 1-D array of the mean squared error of regressionline
;  * obs_std : 1-D array of the conditional values
;
; Description:
;  * calculates the normal (Py_pr) and conditional distribution (Py) following
;    Wilks, D.S. - Statistical Methods in the Atmospheric Science
;
; Return value (PDF):
;    PDF@y     = X values of distribution
;    PDF@Py    = Y values of normal distribution
;    PDF@Py_pr = Y values of conditional distribution
;
; Modification history:
;    * 20140128-zechlau_sabrina: written.

local Pi, nfitx, mfity, miny, maxy, dy, y, sigma, \
      Px, dx, Pxy, Py_given_x, y, Py, Py_pr
begin

  ; Return variable as logical
  PDF = True

  Pi    = 3.142
  nfitx = dimsizes(xfit)

  miny  = min(Yval)-1.1*(max(Yval)-min(Yval))
  maxy  = max(Yval)+1.1*(max(Yval)-min(Yval))

  mfity = fspan(0, 500, 500)
  dy    = (maxy-miny)/dimsizes(mfity)
  y     = new(dimsizes(mfity), "float")
  y     = miny+dy*mfity

  ; calculate +/- 1 sima confidence limits of regressionline
  sigma = new(nfitx, "float")
  do n = 0, nfitx-1
    sigma(n) = sqrt(MSE) * sqrt(1.0 + \
                                (1.0 / dimsizes(Xval)) + \
                                (sum((xfit(n) - avg(Xval)) ^ 2) \
                                 / (dimsizes(Xval) * sum(Xval ^ 2))))
  end do

  if (obs_std(0).ne.0) then
    ; Calculate PDF for conditional values
    Px = new(nfitx, "float")
    Px!0 = "x"
    Px&x = xfit
    dx = xfit(1)-xfit(0)
    Px = [1 / (obs_std(1) * sqrt(2 * Pi))] * exp(-0.5 * (((xfit - obs_std(0)) \
                                                         / obs_std(1)) ^ 2))

    ; Calculate contours of probability in (x,y) space
    Pxy   = new((/nfitx, dimsizes(mfity)/), "float")
    Pxy!0 = "x"
    Pxy&x = xfit
    Pxy!1 = "y"
    Pxy&y = y
    Py    = new(dimsizes(mfity), "float")
    Py!0  = "y"
    Py&y  = y

    do m = 0, dimsizes(mfity)-1
      Py(m) = 0.0
      do n = 0, nfitx-1
        Py_given_x = (1/(sigma(n) * sqrt(2 * Pi))) * \
                     exp(-0.5 * (((y(m) - yfit(n)) / sigma(n)) ^ 2))
        Pxy(n, m) = Px(n) * Py_given_x
        Py(m) = Py(m) + Pxy(n, m)*dx
      end do
    end do

    ; Asign return values conditional PDF
    PDF@Py = Py
  end if

  ; Calculate mean and stdev of (equal model weight) prior
  mn_pr  = avg(Yval(:))
  std_pr = stddev(Yval(:))

  Py_pr = new(dimsizes(mfity), "float")
  Py_pr = (1 / (std_pr * sqrt(2 * Pi))) * exp(-0.5 * (((y - mn_pr) / \
                                                      std_pr) ^ 2))

  ; Assign return values for prior PDF
  PDF@y  = y
  PDF@Py_pr = Py_pr
  PDF@sigma = sigma

  return(PDF)

end
