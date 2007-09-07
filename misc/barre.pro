PRO barre,X,Y,TEXT,LOG=log
;+
; NAME: 
;         barre
;
;
; PURPOSE:
;        draw a bar
;
;
; CATEGORY:
;        plot
;
;
; CALLING SEQUENCE:
;        barre,X,Y,TEXT
;
;
; INPUTS:
;       X : abscissa, two dimentionnal array [X_start, X_end]
;       Y : ordinate, float number/integer
;    TEXT : legend
;
; OPTIONAL INPUTS:
;     LOG : If the plot is in log-log format
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
;
;-


oplot,X,Y*[1,1],thick=!p.thick*2

IF KEYWORD_SET(log) THEN BEGIN
    XYOUTS,10.^(TOTAL(ALOG10(X)*[-1,1])/2+ALOG10(X[0])),Y*1.1,TEXT,$
      ALIGNMENT=0.5

RETURN
END


XYOUTS,TOTAL(X*[-1.,1.])/2+X[0],10.^(ALOG10(Y)*1.1),TEXT, $
  ALIGNMENT=0.5
RETURN

END
