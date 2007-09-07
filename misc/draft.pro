PRO draft
;+
; NAME:     draft 
;
;
;
; PURPOSE:  add a draft sign on your plot
;
;
;
; CATEGORY: plotting
;
;
;
; CALLING SEQUENCE:
;           draft
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;           2005 - AB - creation
;-

string = "DRAFT"
angle = ATAN((!P.CLIP[3]-!P.CLIP[1])*1.0/$
             (!P.CLIP[2]-!P.CLIP[0]))

charsize=SQRT((!P.CLIP[2]-!P.CLIP[0])^2+(!P.CLIP[3]-!P.CLIP[1])^2)/ $
  !D.X_CH_SIZE/(STRLEN(string)+12)

X_CENTER = (!P.CLIP[2]-!P.CLIP[0])/2 + $
  !P.CLIP[0]+!D.X_CH_SIZE*charsize*cos(angle)/2
Y_CENTER = (!P.CLIP[3]-!P.CLIP[1])/2 + $
  !P.CLIP[1]-!D.Y_CH_SIZE*charsize*sin(angle)/2

DEVICE, /HELVETICA, /BOLD, FONT_INDEX=4

xyouts,X_CENTER,Y_CENTER,$
 "!4"+string+"!N",$
  CHARSIZE=charsize,ALIGNMENT=0.5,/DEVICE,$
  ORIENTATION=angle/!pi*180,$
  CHARTHICK=charsize,color=196,FONT=4

END
