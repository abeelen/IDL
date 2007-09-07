PRO plot_sdss_spectrum,spectrum, linetable, TITLE=title,_EXTRA=_extra
;+
; NAME:
;           plot_sdss_spectrum
;
;
; PURPOSE:  plot an SDSS spectrum the "SDSS way"
;          
;
;
; CATEGORY: plotting
;
;
;
; CALLING SEQUENCE:
;           plot_sdss_spectrum,spectrum,[TITLE=]
;
;
; INPUTS:
;          specrum : the spectrum read by read_sdss_specrum
;
;
; OPTIONAL INPUTS:
;          linetable : description of the line detected in the spectrum 
;          title     : the title of the plot
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
;          read_sdss_spectrum
;          textoidl
;          SXPAR (ASTRON lib)
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;       2006 - AB - created
;-


IF NOT KEYWORD_SET(Title) THEN BEGIN 
    ra = SXPAR(spectrum.header,'RADEG')
    dec = SXPAR(spectrum.header,'DECDEG')
    MJD = SXPAR(spectrum.header,'MJD')
    Plate = SXPAR(spectrum.header,'PLATEID')
    Fiber = SXPAR(spectrum.header,'FIBERID')
    
    TITLE = "RA=" + STRING(ra,   FORMAT='(F10.5)')+$
      ", DEC="    + STRING(dec,  FORMAT='(F10.5)')+$
      ", MJD="    + STRING(MJD,  FORMAT='(I5)')+$
      ", Plate="  + STRING(Plate,FORMAT='(I4)')+$
      ", Fiber="  + STRING(Fiber,FORMAT='(I3)')
ENDIF

plot,spectrum.wavelength,smooth(spectrum.count,8,/NAN),$
  XTITLE='Wavelength [A]',$
  YTITLE=TEXTOIDL(' F_\lambda (10^{-17} erg cm^{-2} s^{-1} A^{-1}'),$
  TITLE=TITLE,$
  /XSTYLE,_EXTRA=_extra

oplot,spectrum.wavelength,smooth(spectrum.counterr,8,/NAN),color=4

FOR iLine = 0,N_ELEMENTS(linetable)-1 DO BEGIN
    IF linetable[iLine].wave NE -9999 AND $
      ($
        (!X.TYPE EQ 0 AND $
         MIN(!X.CRANGE) LT linetable[iLine].wave AND MAX(!X.CRANGE) GT linetable[iLine].wave) $
        OR $
        (!X.TYPE EQ 1 AND $
         10.^MIN(!X.CRANGE) LT linetable[iLine].wave AND 10.^MAX(!X.CRANGE) GT linetable[iLine].wave) $
      ) THEN BEGIN
        
        close = WHERE($
                       ABS(spectrum.wavelength-linetable[iLine].wave) EQ $
                       MIN(ABS(spectrum.wavelength-linetable[iLine].wave),/NAN)$
                     )
        IF !Y.TYPE THEN BEGIN
            maxY=10.^MAX(!Y.CRANGE)*0.85
        ENDIF ELSE BEGIN
            maxY=MAX(!Y.CRANGE)
        ENDELSE
        
        OPLOT,[1,1]*linetable[iLine].wave,$
          [ spectrum.count[close],maxY]+[2,-2],$
          color=5,linestyle=1

        XYOUTS,linetable[iLine].wave,$
          maxY-1,$
          linetable[iLine].name,$
          ALIGNMENT=0,$
          ORIENTATION=90,$
          color=5

    ENDIF
ENDFOR

END
