PRO copyright, STRING=string
;+
; NAME:     copyright
;
;
;
; PURPOSE:  add a copyright info on your plot
;
;
;
; CATEGORY: plotting
;
;
;
; CALLING SEQUENCE: 
;           copyright,[STRING=string]
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
;           string : the string to print (default: $USERNAME)
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
;        2005 - AB - creation
;-

IF NOT KEYWORD_SET(string) THEN $
  string = GETENV("USERNAME")

xyouts,!P.CLIP[2],!P.CLIP[3]+!D.Y_CH_SIZE/2,$
  "© ("+(STRSPLIT(systime()," ",/EXTRACT))[4]+") "+string , $
  CHARSIZE=0.8,ALIGNMENT=1,/DEVICE

END
