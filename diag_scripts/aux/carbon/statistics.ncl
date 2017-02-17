;#######################################################################
;# GENERAL ROUTINES FOR STATISTICS                                     #
;#######################################################################
;# Please consider using of extending existing routines before adding  #
;# new ones. Check the header of each routine for documentation.       #
;#                                                                     #
;# Contents:                                                           #
;#                                                                     #
;#    * function Cond_PDF: calculates the normal and conditional PDF   #
;#    * function seas_cycle: analyzes the seasonal cycle of an periodic#
;#                           variable                                  #    
;#                                                                     #
;#######################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"

;#######################################################################
undef("Cond_PDF")
function Cond_PDF(xfit[*]:numeric,   \ ; x-values (from regressionline y=a*x+b)
                  yfit[*]:numeric,   \ ; y-values (from regressionline y=a*x+b)
                  Xval[*]:numeric,   \ ; x-values (from models)
                  Yval[*]:numeric,   \ ; y-values (from models)
                  MSE[*]:numeric,    \ ; mean sqaured error of regressionline
                  obs_std[*]:numeric \ ; conditional values
                  )
;;
;; Arguments:
;;  * xfit    : x-values (from regressionline y=a*x+b)
;;  * yfit    : y-values (from regressionline y=a*x+b)
;;  * Xval    : x-values (from models)
;;  * Yval    : y-values (from models)
;;  * MSE     : mean sqaured error of regressionline
;;  * obs_std : conditional values
;;
;; Description:
;;  * calculates the normal (Py_pr) and conditional dirstibution (Py) following
;;    Wilks, D.S. - Statistical Methods in the Atmospheric Sciences
;;
;; Return value (PDF):
;;    PDF@y     = X values of distribution
;;    PDF@Py    = Y values of normal distribution
;;    PDF@Py_pr = Y values of conditional distribution
;;
;; Arguments:
;;  xfit = 1-D array of the X variable (from regressionline y=a*x+b)  
;;  yfit = 1-D array of the Y variable (from regressionline y=a*x+b)  
;;  Xval = 1-D array of the Y variable                                
;;  Yval = 1-D array of the Y variable                                
;;  MSE  =  1-D array of the mean sqaured error of regressionline      
;;  obs_std = 1-D array of the conditional values                     
;;                                                                      
;; Modification history:
;;    * 20140128 written by S. Wenzel (sabrina.wenzel@dlr.de)
;;
local Pi,nfitx,mfity,miny,maxy,dy,y,sigma,Px,dx,Pxy,Py_given_x,y,Py,Py_pr
begin

  ;;Return variable as logical
  PDF = True

  Pi    = 3.142
  nfitx = dimsizes(xfit)
  
  miny  = min(Yval)-1.1*(max(Yval)-min(Yval))
  maxy  = max(Yval)+1.1*(max(Yval)-min(Yval))

  mfity = fspan(0,500,500)
  dy    = (maxy-miny)/dimsizes(mfity)
  y     = new(dimsizes(mfity),"float")
  y     = miny+dy*mfity

  ;;calculate +/- 1 sima confidence limits of regressionline 
  sigma = new(nfitx,"float")
  do n =0,nfitx-1
    sigma(n) = sqrt(MSE)*sqrt(1.0+(1.0/dimsizes(Xval))+(sum((xfit(n)-avg(Xval))^2)/(dimsizes(Xval)*sum(Xval^2))))
  end do

  if (obs_std(0).ne.0) then
    ;;Calculate PDF for conditional values
    Px = new(nfitx,"float");xfit ;
    Px!0 = "x"
    Px&x = xfit
    dx = xfit(1)-xfit(0)
    Px = [1/((obs_std(1))*sqrt(2*Pi))] * exp(-0.5*(((xfit-obs_std(0))/obs_std(1))^2))

    ;;Calculate contours of probability in (x,y) space
    Pxy   = new((/nfitx,dimsizes(mfity)/),"float")
    Pxy!0 = "x"
    Pxy&x = xfit
    Pxy!1 = "y"
    Pxy&y = y
    Py    = new(dimsizes(mfity),"float")
    Py!0  = "y"
    Py&y  = y

    do m = 0,dimsizes(mfity)-1
      Py(m) = 0.0
      do n = 0,nfitx-1 
        Py_given_x = (1/(sigma(n)*sqrt(2*Pi))) * exp(-0.5*(((y(m)-yfit(n))/sigma(n))^2))
        Pxy(n,m) = Px(n)*Py_given_x
        Py(m) = Py(m)+Pxy(n,m)*dx
      end do
    end do

    ;;Asign return values conditional PDF
    PDF@Py = Py
  end if

  ;;Calculate mean and stdev of (equal model weight) prior
  mn_pr  = avg(Yval(:))
  std_pr = stddev(Yval(:))

  Py_pr = new(dimsizes(mfity),"float");y
  Py_pr = (1/(std_pr*sqrt(2*Pi))) * exp(-0.5*(((y-mn_pr)/std_pr)^2))

  ;;Asign return values for prior PDF
  PDF@y  = y
  PDF@Py_pr = Py_pr
  PDF@sigma = sigma

  return(PDF)

end

;########################################################################
undef("seas_cycle")
function seas_cycle(\;Xann[*]:numeric, \ ; x-values annual
                   Xdata[*]:numeric, \ ; x-values monthly
                   Ydata[*]:numeric, \ ; y-values 
                     opt[1]:string)
;;
;; Author: Sabrina Wenzel (DLR, Germany)
;;
;; Description:
;;  * analyses the sasonal cycle of an de-trendet periodic timeseries eg CO2
;;    - calculates the relative seasonal amplitude of the timeseries
;;    - finds the zero crossing point (phase) in Days Of Year (DOY)
;;    - integrates over growing season (negative values of cycle)
;;
;; Return value (SEAS_CYC):
;;    SEAS_CYC@Xamp = X values of relative amplitude
;;    SEAS_CYC@Yamp = Y values of relative amplitude
;;    SEAS_CYC@Xday = X values of phase (zero crossing)
;;    SEAS_CYC@Yday = Y values of phase (zero crossing)
;;    SEAS_CYC@Yintg = Y values of integral for growing season
;;    SEAS_CYC@Yintd = Y values of integral for dormant season
;;
;; Arguments:
;;  Xdata = 1-D array of x-values monthly
;;  Ydata = 1-D array of y-values 
;;  opt   = possible options: 
;;         'amplitude' only calculates the relative amplitude of the cycle
;;         'integral'  only calculates the integral over growing season
;;         'phase'     only calculates the phase of the cycle
;;         'all'       calculates all
;;                                                                      
;; Modification history:
;;    * 20140218 written by S. Wenzel (sabrina.wenzel@dlr.de)
;;
local y, aa, n, ttt, ttx, txx, ymin, ymax, y2, aY2, aX2, aY1, aX1, \
      alpha, tmpxh, y1
begin

  ;;Return variable as logical
  SEAS_CYC = True

  ;;define arrays for output
  tmp  = ispan(toint(min(Xdata(:))),toint(max(Xdata(:))),1) * 1. 
  nyr  = dimsizes(tmp) 
  Xann = tmp
  delete(tmp)
  Xamp = new(nyr,float)
  Yamp = new(nyr,float)
  Xampi = new(nyr,float)
  Yampi = new(nyr,float)
  Xmax = new(nyr,float)
  Ymax = new(nyr,float)
  Xmin = new(nyr,float)
  Ymin = new(nyr,float)
  Xday = new(nyr,float)
  Yday = new(nyr,float)
  Yintg = new(nyr,float)
  Yintd = new(nyr,float) 
  Ygrow = new(nyr,float) 
  Xmon = tofloat(Ydata&time) ;units must be in 'days since JJJJ-MM-DD hh:mm:ss'
  Ydata&time = Xdata

  ;;analyze each year seperatly:
  do y = 0,nyr-1 ;0,dimsizes(aYdetr(imod,:))-1
    n = 0
    ttt = ind(toint(Ydata&time).eq.toint(Xann(y))) 
    aa  = Ydata(ttt) 
    n   = dimsizes(aa)

    if (num(ismissing(aa)).le.5) then ;.not.all(ismissing(aa)).or.
      if (max(aa).ne.min(aa)) then

       ;;growth rate:
        if (opt .eq. "growth")then    
          if (n(0) .eq. 12 .and. min(ttt) .ge. 12) then
            tmi = min(ttt)
            Ygrow(y-1) = avg(Ydata(imod,tmi-2:tmi+1))-avg(Ydata(imod,tmi-14:tmi-11))
          end if
        end if

        if (opt .eq. "ampint" .or. opt .eq. "all")then
            tmp = dim_cumsum(aa,0)/12
            delete(aa)
            aa = tmp
            delete(tmp)
          if (n(0) .eq. 12) then
            Yampi(y) = max(aa)-min(aa)  ;rel amplitude eq max + min of seas cycl
            ttx = maxind(aa)
            txx = y*12+ttx
            Xampi(y) = Xdata(txx)                 ;find maching x value
            delete(ttx)
            delete(txx)
          end if
        end if

        ;;relativ seasonal amplitude:
        if (opt .eq. "amplitude" .or. opt .eq. "all")then
          if (n(0) .eq. 12) then
            Yamp(y) = max(aa)-min(aa)  ;rel amplitude eq max + min of seas cycl
            ttx = maxind(aa)
            txx = y*12+ttx
            ;Xamp(y) = Xdata(txx)                 ;find maching x value
            delete(ttx)
            delete(txx)
          end if
        end if  ;opt 'amplitude'

        ;;relativ seasonal amplitude:
        if (opt .eq. "maximum" .or. opt .eq. "all")then
          if (n(0) .gt. 9.) then
            Ymax(y) = max(aa(0:9))               ;max of seas cycl
            ttx = maxind(aa(0:9))
            txx = y*12+ttx
            Xmax(y) = Xdata(txx)                 ;find maching x value
            delete(ttx)
            delete(txx)
          end if
        end if  ;opt 'maximum'

        ;;relativ seasonal amplitude:
        if (opt .eq. "minimum" .or. opt .eq. "all")then
          if (n(0) .gt. 9.) then
            Ymin(y) = min(aa(0:9))               ;max of seas cycl
            ttx = maxind(aa(0:9))
            txx = y*12+ttx
            Xmin(y) = Xdata(txx)                 ;find maching x value
            delete(ttx)
            delete(txx)
          end if
        end if  ;opt 'minimum'

        ;;linearly interpolate in values to get zero crossing point:
        if (opt .eq. "phase" .or. opt .eq. "all")then
          if (n(0).gt.9. .and. y.ne.0) then
            ymin = maxind(aa(0:9))               ;find min and max in
            ymax = minind(aa(0:9))               ;the seasonal cycle

            if (ymax.gt.ymin) then
              y2  = ind(aa(ymin:ymax).lt.0.0)    ;only consider values lt 0.0
              aY2 = aa(min(y2+ymin))             ;find the closest neg. value  to 0.0
              ttx = min(y2+ymin)
              txx = y*12+ttx
              aX2 = Xmon(txx)-Xmon(y*12)

              if (.not.ismissing(aa(ttx-1))) then
                aY1 = aa(ttx-1)                  ;at ttx-1 is closest pos. value to 0.0
                aX1 = Xmon(txx-1)-Xmon(y*12)
              
                if (aX2.ne.aX1) then
                  alpha = (aY2-aY1)/(aX2-aX1)    ;interpolate between ponts
                  tmpxh = aX1 - aY1/alpha        ;before and after zero crossing
                  Yday(y) = tmpxh
                  Xday(y) = Xdata(txx)
                  delete(tmpxh)
                  delete(alpha)
                end if
                delete(aY1)
                delete(aX1)
              end if

              delete(ttx)
              delete(txx)
              delete(aX2)
              delete(aY2)
              delete(ymax)
              delete(ymin)
              delete(y2)
            end if
          end if
        end if  ;opt 'phase'

      end if
    end if

    ;;integrate over growing saison:
    if (opt .eq. "integral" .or. opt .eq. "all")then
      if (n(0) .gt. 9.) then
        y1 = ind(aa(1:n-1).lt.0.0)                  ;only consider values lt 0.0 for integral
        if (.not.all(ismissing(y1))) then
          Yintg(y) = sum(aa);(y1))
        end if
        delete(y1)
        y1 = ind(aa(1:n-1).gt.0.0)                  ;only consider values gt 0.0 for integral
        if (.not.all(ismissing(y1))) then
          Yintd(y) = sum(aa(y1))
        end if
        delete(y1)
      end if
    end if
    delete(n)
    delete(aa) 
    delete(ttt)

  end do ;annual loop 

  SEAS_CYC@Xamp = Xamp
  SEAS_CYC@Yamp = Yamp
  SEAS_CYC@Xmax = Xmax
  SEAS_CYC@Ymax = Ymax
  SEAS_CYC@Xmin = Xmin
  SEAS_CYC@Ymin = Ymin
  SEAS_CYC@Xday = Xday
  SEAS_CYC@Yday = Yday
  SEAS_CYC@Yintg = Yintg
  SEAS_CYC@Yintd = Yintd
  SEAS_CYC@Ygrow = Ygrow
  SEAS_CYC@Yai   = Yampi

  return(SEAS_CYC)
end
