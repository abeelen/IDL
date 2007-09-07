pro plot_histo, data, histo,bins,$
                NBINS=nbins,BINSIZE=binsize,$
                MAX=max,MIN=min,$
                OVERPLOT=overplot,NOPLOT=noplot,_EXTRA=extra,COLOR=color

;+
; NAME: plot_histo
;
;
;
; PURPOSE:
;         plot histogram of an array correctly 
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
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
;
;-

IF NOT KEYWORD_SET(nbins) THEN nbins = 20
IF NOT KEYWORD_SET(binsize) THEN binsize=(MAX(data,/NAN)-MIN(data,/NAN))/nbins

histo=[0, histogram(data,binsize=binsize, max=max,min=min,locations=X,/NAN), 0]
bins = [X[0], X, X[N_ELEMENTS(X)-1] ]+[-1,INDGEN(N_ELEMENTS(X))*0+1,3]*(X[1]-X[0])/2

IF KEYWORD_SET(NOPLOT) THEN $
  RETURN

IF NOT KEYWORD_SET(overplot) THEN $
PLOT,bins,histo,/XSTYLE,/NODATA,PSYM=10,_EXTRA=extra

OPLOT,bins,histo,PSYM=10,_EXTRA=extra,color=color

return
end
