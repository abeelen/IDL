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
  
  ON_ERROR, 2 
  
  IF NOT KEYWORD_SET(nbins)   THEN nbins = 20
  IF NOT KEYWORD_SET(max)     THEN max = MAX(data,/NAN)
  IF NOT KEYWORD_SET(min)     THEN min = MIN(data,/NAN)
  IF NOT KEYWORD_SET(binsize) THEN binsize=(max-min)/(nbins-1)
  
  histo=[0, histogram(data,binsize=binsize, max=max,min=min,locations=X,/NAN), 0]
  bins = [X[0], X, X[N_ELEMENTS(X)-1] ]+[-1,INDGEN(N_ELEMENTS(X))*0+1,3]*(X[1]-X[0])/2
  
  IF KEYWORD_SET(NOPLOT) THEN $
     RETURN
  
  IF NOT KEYWORD_SET(overplot) THEN $
     PLOT,bins,histo,/XSTYLE,/NODATA,PSYM=10,_EXTRA=extra
  
  OPLOT,bins,histo,PSYM=10,_EXTRA=extra,color=color
  
  return
end
