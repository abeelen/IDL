PRO fit_line, filename, result, perror
;+
; NAME:    
;         fit_line
;
;
;
; PURPOSE: 
;          fit a gaussian line on a spectrum obtained by the task UVFIT
;          in GILDAS
;
;
;
; CATEGORY: 
;          fit/plot
;
;
;
; CALLING SEQUENCE: 
;          fit_line, filename [ , result]
;
;
;
; INPUTS:  
;          filename : the name of the fits file containing the spectrum
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
;          result : the result of the fit (MPFIT format)
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
;        2004, A. Beelen created
;        2005, A. Beelen cleaning
;-

tab = mrdfits(filename,0,header,/SILENT)

velocity  = tab[*,3]
intensity = tab[*,11]
error     = tab[*,12]

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], parname:"", $
                     limits:[0.D,0]}, 4)
parinfo(*).parname  = ["int_flux", "position", "FWHM", "continuum"] 
parinfo(*).value    = [0, 0, 300,0]*1.D

fcnargs = {X:velocity, $
           Y:intensity, $
           ERR:error}


; line in emission
parinfo[0].limited = [1,0]
parinfo[0].limits  = [0,0]

; guess line position
parinfo[1].value = velocity[WHERE(intensity EQ MAX(intensity))]

result = mpfit("mp_gauss", parinfo=parinfo, $
               functargs=fcnargs, $
               MAXITER=2000,$
               STATUS=status,PERROR=perror,/QUIET)

print,"Gaussian fit Result : "
print,FORMAT='("  int flux  : ",F9.6," +- ",F9.6," Jy km/s")',result[0],perror[0]
print,FORMAT='("  base      : ",F9.6," +- ",F9.6," mJy")',result[3]*1.d3,perror[3]*1.d3
print,FORMAT='("  position  : ",F9.5," +- ",F9.5," km/s")',result[1],perror[1]
print,FORMAT='("  FWHM      : ",F9.3," +- ",F9.3," km/s")',result[2],perror[2]

; Plotting routines

ploterror,velocity-result[1],intensity*1.d3,error*1.d3,/nodata,$
  XTITLE=TEXTOIDL(' \Delta v (km/s)'),$
  YTITLE="Flux density (mJy)",$
  XSTYLE=9

oploterror,velocity-result[1],intensity*1.d3,error*1.d3,$
  psym=10,errcolor=1,color=1

oplot,!X.CRANGE,[0,0],LINESTYLE=2
oplot,[0,0],!Y.CRANGE,LINESTYLE=1

X = indgen(100,/FLOAT)/100*(max(!X.CRANGE)-min(!X.CRANGE))+min(!X.CRANGE)
temp_result = result
temp_result[1]=0
gauss_plus_cont,X,temp_result,F
OPLOT,X,F*1.d3,linestyle=2,color=5

END
