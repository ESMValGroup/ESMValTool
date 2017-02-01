;;#############################################################################
undef("notextoverlap")
procedure notextoverlap(wks[1]:graphic,
                        plot[1]:graphic)

;;
;; Arguments
;;    plot: plot identifier, e.g. created by gsn_*.
;;    xdata: array of x values, same size as ydata.
;;    ydata: array of y values, same size as xdata.
;;
;; Description
;;    Adjusts text position to be inside the box.
;;    Checks for overlapping text strings.
;;    If text strings overlap than chose another 'txJust' position
;;
;; Caveats
;;
;; References
;;    https://www.ncl.ucar.edu/Applications/Scripts/text_17.ncl
;;    https://www.ncl.ucar.edu/Applications/Scripts/text_10.ncl
;;
;; Modification history
;;    20150722-A_wenz_sa: written.
;;
local xndc, yndc, rgt, lft, top, bot, width, height, nstrs, removed_list, \
      rcount, ibox1, ibox2, vpx, vpy, vpw, vph, nstrs
begin

    funcname = "notextoverlap"
    scriptname = "plot_scripts/ncl/aux_plotting.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    getvalues plot
        "pmAnnoViews" : text_ids
    end getvalues

    nstrs  = dimsizes(text_ids)
    tx     = new(4, float)  ; X top left, top right, bottom right, bottom left
    ty     = new(4, float)  ; Y top left, top right, bottom right, bottom left
    rgt    = new(nstrs, float)  ; hold right position of text box.
    lft    = new(nstrs, float)  ; "    left     "     "   "    " 
    top    = new(nstrs, float)  ; "    top      "     "   "    " 
    bot    = new(nstrs, float)  ; "    bottom   "     "   "    " 
    width  = new(nstrs, float)  ; width of each text box
    height = new(nstrs, float)  ; height of each text box
    text   = new(nstrs, string) 

    do i = 0, dimsizes(text_ids) - 1
        if (NhlClassName(text_ids(i)).eq."textItemClass") then

            ;; Retrieve the text box location and size
            getvalues text_ids(i)
                "txString"  : text(i)
                "vpXF"      : vpx
                "vpYF"      : vpy
                "vpWidthF"  : width(i)
                "vpHeightF" : height(i)
            end getvalues

            print(i + "  " + text(i) + "  " + vpx + " "+ vpy + " "+ width(i) + " "+ height(i))

            ;; This is the box that encloses the text in NDC space
            xndc = (/vpx, vpx + width(i), vpx + width(i), vpx/)
            yndc = (/vpy, vpy, vpy - height(i), vpy - height(i)/)

            ;; Convert box fron NDC coordinates to plot coordinates
            NhlNDCToData(plot, xndc, yndc, tx, ty)
    
            print("---->" + tx(0) + "  "  + tx(1) + "  "  + tx(2) + "  "  + tx(3))
            print("---->" + ty(0) + "  "  + ty(1) + "  "  + ty(2) + "  "  + ty(3))

            ;; Adjust the text that it fits inside the box
            if (all(ismissing((/tx(1), tx(2), ty(1), ty(2)/))))
                setvalues text_ids(i)
                    "txJust" : "CenterRight"
                end setvalues
            end if
            delete([/xndc,yndc/])
            delete([/vpx,vpy/])
        end if
    end do
    NhlUpdateData(plot)
    draw(plot)
    frame(wks)

    ;; Convert the x/y center of each box to NDC coordinates
    xndc = new(nstrs, float)  ;  X center of box in NDC coords
    yndc = new(nstrs, float)  ;  Y center of box in NDC coords
    datatondc(plot, tofloat(xdata), tofloat(ydata), xndc, yndc)

    ;; Calculate the four corners of each text box in NDC coordinates
    top = yndc + height/2.
    bot = yndc - height/2.
    lft = xndc - width/2.
    rgt = xndc + width/2.
    
    removed_list = new(nstrs,integer)
    rcount = 0                         ; Number of canged text boxes.

;;---ibox1 is index of box that we're checking to see if it overlaps.
;;   ibox2 is the index of the box that we're checking box ibox1
;;   against.
    do ibox1 = 0,nstrs-1
        ibox2 = 0
        overlap_found = False 

        do while(ibox2.lt.nstrs.and..not.overlap_found)

;;---Check if any one of the corners represented by box "ibox1" is inside
;;   the box representd by "ibox2".

            if (ibox1.ne.ibox2.and.(.not.any(ibox2.eq.removed_list))) then
                if ((top(ibox1).le.top(ibox2).and.top(ibox1).ge.bot(ibox2).and. \
                    lft(ibox1).le.rgt(ibox2).and.lft(ibox1).ge.lft(ibox2)).or. \
                    (bot(ibox1).le.top(ibox2).and.bot(ibox1).ge.bot(ibox2).and. \
                    lft(ibox1).le.rgt(ibox2).and.lft(ibox1).ge.lft(ibox2)).or. \
                    (bot(ibox1).le.top(ibox2).and.bot(ibox1).ge.bot(ibox2).and. \
                    rgt(ibox1).le.rgt(ibox2).and.rgt(ibox1).ge.lft(ibox2)).or. \
                    (top(ibox1).le.top(ibox2).and.top(ibox1).ge.bot(ibox2).and. \
                    rgt(ibox1).le.rgt(ibox2).and.rgt(ibox1).ge.lft(ibox2))) then
                    
                    overlap_found = True         ; Mark that an overlap has been found.
                    
                    getvalues text_ids(ibox1)
                    "txJust" : text_just
                    end getvalues

                    if (any(text_just .eq. (/"TopLeft","CenterRight"/))) then
                        setvalues text_ids(ibox1)
                        "txJust" : "TopRight"
                        end setvalues
                    else
                        setvalues text_ids(ibox1)
                        "txJust" : "CenterRight"
                        end setvalues
                    end if

                    removed_list(rcount) = ibox1    ; Update removed_list and
                    rcount = rcount + 1             ; its counter.
                    delete(text_just)
                end if
            end if
            ibox2 = ibox2 + 1               ; Update the box counter
        end do
    end do

;;---Replace the bad strings
    NhlUpdateData(plot)
    
    leave_msg(scriptname, funcname, 4)

end
