;
; Requires: none
; Test case for Hisako
;
info = True
info@derived = False
info@long_name = "Surf Temp (radiative)"
info@units = "K"
info@pcp_cnlvl = (/210,220,230,240,250,255,260,265,270,275,280,285,290,295,300/)
info@pcp_dcnlvl = (/-15,-12,-9,-7,-5,-3,-1,0,1,3,5,7,9,12,15/) 
info@scp_cnlvl = (/210,220,230,240,250,260,270,275,280,285,290,295,300,305,310/)
info@scp_dcnlvl = (/-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,12/)

info@acp_cnlvl = (/210,220,230,240,250,260,270,275,280,285,290,295,300,305,310/)
info@acp_dcnlvl = (/-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,12/)

;info@sdp_cnlvl = (/210,220,230,240,250,260,270,275,280,285,290,295,300,305,310/)
info@sdp_cnlvl = (/-12,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,12/)
info@rgb_file = "amwg.rgb"
