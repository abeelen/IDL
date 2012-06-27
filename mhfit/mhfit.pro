;+
; NAME: 
;      MHFIT
;
;
; PURPOSE:
;      Perform Metropolis Hasting Markov Chain to perfom fit
;
;
; CATEGORY:
;      Fitting
;
;
; CALLING SEQUENCE:
;      parms = MHFIT(MYFUNCT, start_parms, INCOVAR=incovar
;                    FUNCTARGS=fcnargs, MAXITER=maxiter, 
;                    QUIET=quiet, COVAR=covar, perror=perror, $
;                    NPRINT=nprint, ITERPROC=iterproc, $
;                    PARINFO=parinfo, SCALE=scale, $
;                    CHAINS=chains, ACCEPT=accept, lnL = lnL, $
;                    SAVE_STEP=save_step, RESTORE_STEP=restore_step
;
; INPUTS:
;   MYFUNCT - a string variable containing the name of the function to
;             be minimized.  The function should return the weighted
;             deviations between the model and the data, as described
;             above.
;
;             For EXTERNAL evaluation of functions, this parameter
;             should be set to a value of "(EXTERNAL)".
;
;   START_PARAMS - An array of starting values for each of the
;                  parameters of the model.  The number of parameters
;                  should be fewer than the number of measurements.
;                  Also, the parameters should have the same data type
;                  as the measurements (double is preferred).
;
;                  This parameter is optional if the PARINFO keyword
;                  is used (but see PARINFO).  The PARINFO keyword
;                  provides a mechanism to fix or constrain individual
;                  parameters.  If both START_PARAMS and PARINFO are
;                  passed, then the starting *value* is taken from
;                  START_PARAMS, but the *constraints* are taken from
;                  PARINFO.
;   INCOVAR      - An NxN Matric containing the input covariance
;                  Matrix for the sampling.
; 
;   FUNCTARGS - A structure which contains the parameters to be passed
;               to the user-supplied function specified by MYFUNCT via
;               the _EXTRA mechanism.  This is the way you can pass
;               additional data to your user-supplied function without
;               using common blocks.
;
;               Consider the following example:
;                if FUNCTARGS = { XVAL:[1.D,2,3], YVAL:[1.D,4,9],
;                                 ERRVAL:[1.D,1,1] }
;                then the user supplied function should be declared
;                like this:
;                FUNCTION MYFUNCT, P, XVAL=x, YVAL=y, ERRVAL=err
;
;               By default, no extra parameters are passed to the
;               user-supplied function, but your function should
;               accept *at least* one keyword parameter.  [ This is to
;               accomodate a limitation in IDL's _EXTRA
;               parameter-passing mechanism. ]
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;   MAXITER - The maximum number of iterations to perform.  If the
;             number is exceeded, then the STATUS value is set to 5
;             and MPFIT returns.
;             Default: 10000 iterations
;
;   QUIET    - set this keyword when no textual output should be printed
;              by MHFIT
;
;
;   PARINFO  - See MPFIT documention, accepted parameters :
;     .VALUE - the starting parameter value (but see the START_PARAMS
;              parameter for more information).
;  
;     .FIXED - a boolean value, whether the parameter is to be held
;              fixed or not.  Fixed parameters are not varied by
;              MPFIT, but are passed on to MYFUNCT for evaluation.
;  
;     .LIMITED - a two-element boolean array.  If the first/second
;                element is set, then the parameter is bounded on the
;                lower/upper side.  A parameter can be bounded on both
;                sides.  Both LIMITED and LIMITS must be given
;                together.
;  
;     .LIMITS - a two-element float or double array.  Gives the
;               parameter limits on the lower and upper sides,
;               respectively.  Zero, one or two of these values can be
;               set, depending on the values of LIMITED.  Both LIMITED
;               and LIMITS must be given together.
;  
;     .PARNAME - a string, giving the name of the parameter.  The
;                fitting code of MPFIT does not use this tag in any
;                way.  However, the default ITERPROC will print the
;                parameter name if available.
;
;     .TIED - a string expression which "ties" the parameter to other
;             free or fixed parameters.  Any expression involving
;             constants and the parameter array P are permitted.
;             Example: if parameter 2 is always to be twice parameter
;             1 then use the following: parinfo(2).tied = '2 * P(1)'.
;             Since they are totally constrained, tied parameters are
;             considered to be fixed; no errors are computed for them.
;             [ NOTE: the PARNAME can't be used in expressions. ]
;  
;
;     .MPPRINT - if set to 1, then the default ITERPROC will print the
;                parameter value.  If set to 0, the parameter value
;                will not be printed.  This tag can be used to
;                selectively print only a few parameter values out of
;                many.  Default: 1 (all parameters printed)
;
;     .MPFORMAT - IDL format string to print the parameter within
;                 ITERPROC.  Default: '(G20.6)' An empty string will
;                 also use the default.
;
;
;
; OUTPUTS:
;  
;
;
; OPTIONAL OUTPUTS:
;
;   COVAR - the covariance matrix for the set of parameters returned
;           by MHFIT. 
;
;   PERROR - The formal 1-sigma errors in each parameter, computed
;            from the covariance matrix.  
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
;   MHFIT need to have access to Markwardt IDL Library and the IDL Astronomy User's Library
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;   2012-06 Beelen. A Creation
;-
PRO mhfit_dummy
  ;; Enclose in a procedure so these are not defined in the main level
  FORWARD_FUNCTION myfunct_fft, mhfit_step, mhfit, mpfit_call, mpfit_enorm, mpfit
END


FUNCTION MYFUNCT_FFT, p,dp,  k=K, lnPk=lnPk, Err=err
;; Helper function to fit the fft of the chains

  P_0   = p[0]
  k_s   = p[1]
  alpha = p[2]
  

  w = (k_s/k)^alpha

  model = P_0*w/(w+1)
  
  IF N_PARAMS() GT 1 THEN BEGIN
     DP = DBLARR(N_ELEMENTS(k), 3)
     DP[*,0] = w/(w+1)
     DP[*,1] = P_0 * alpha / k_s * w / (w+1)^2
     DP[*,2] = P_0 * w * ALOG(k_s/k) / (w+1)^2
  ENDIF

  return, (lnPk-ALOG(model))
END

PRO MHFIT_ACHAIN, chains, parinfo
;; Analyse chains for convergence (FFT fits)

s    = size(chains)
npar = s[1]
N    = s[2]

k     = 2*!pi*(FINDGEN(N/2+1))/N

if (N lt 2000 ) then jmax=99     ;this can be varied. Should be about 10 x j* 
if (N ge 2000 ) then jmax=999    ;this can be varied. Should be about 10 x j* 
if (N ge 10000) then jmax=N/2.5  ;this can be varied. Should be about 10 x j* 

fft_chains =  (ABS(FFT(STANDARDIZE(chains,/DOUBLE), -1, DIMENSION=2)))^2*N
fft_chains = fft_chains[*,0:N/2+1]


;; Finding non fixed/tied parameters
mpfit_parinfo, parinfo, tagnames, 'TIED',   ptied, default='', n=npar
mpfit_parinfo, parinfo, tagnames, 'FIXED', pfixed, default=0,  n=npar
goodChains = WHERE(pfixed EQ 0 AND STRTRIM(ptied,2) EQ '', nGoodChains)

erase
multiplot,[ROUND(SQRT(nGoodChains*1.)),CEIL(nGoodChains*1./ROUND(SQRT(nGoodChains*1.)))], mxtitle="k"
PRINT, "parname","j*","r","conv", FORMAT='(A20,1X,A16,A22,4X,A4)'
PRINT, REPLICATE("-",67), FORMAT='(67A1)'

FOR iPar = 0, nGoodChains-1 DO BEGIN
   plot_oo, k[1:*], fft_chains[goodChains[iPar],1:*],/XSTYLE,YTITLE="P(k) for "+parinfo[goodChains[iPar]].parname

   
   lparinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 3)
   lparinfo[0].limited = [1,0]
   lparinfo[0].limits  = [MIN(fft_chains[goodChains[iPar],1:*]),0]
   
   lparinfo[1].limited = [1,1]
   lparinfo[1].limits  = MINMAX(k[1:*])

   lparinfo[2].limited = [1,1]
   lparinfo[2].limits  = [0,5]
   
   ;; first guesses
   lparinfo.value = [exp(MEAN(ALOG(fft_chains[goodChains[iPar],1:10]))),200./N*2*!pi, 2]
   
   ;; two pass fit to avoid small scale artefacts (Dunkley et al. 2004)
;;   fcnargs = {k: k[1:jmax],lnPk:ALOG(REFORM(fft_chains[goodChains[iPar],1:jmax])), ERR:k[1:jmax]*0+1}
   fcnargs = {k: k[1:jmax],lnPk:ALOG(REFORM(fft_chains[goodChains[iPar],1:jmax])), ERR:k[1:jmax]*0+1}
   params = MPFIT('MYFUNCT_FFT',PARINFO=lparinfo,FUNCTARGS=fcnargs,AUTODERIVATIVE=0,PERROR=perror,/QUIET)
   lparinfo.value = params
   oplot, k, params[0]*((params[1]/k)^params[2])/((params[1]/k)^params[2]+1),color=2
;;    oplot, [1,1]*jmax*2*!pi/N,10.^!Y.CRANGE, linestyle=2, color=2

   ;; select data based on the new jmax = 10 j*
   jmax = 10*params[1]*N/(2*!pi) < (N/2+1)-1
   fcnargs = {k: k[1:jmax],lnPk:ALOG(REFORM(fft_chains[goodChains[iPar],1:jmax])), ERR:k[1:jmax]*0+1}
   params = MPFIT('MYFUNCT_FFT',PARINFO=lparinfo,FUNCTARGS=fcnargs,AUTODERIVATIVE=0,PERROR=perror,/QUIET)
   oplot, k, params[0]*((params[1]/k)^params[2])/((params[1]/k)^params[2]+1),color=3
   oplot, [1,1]*jmax*2*!pi/N,10.^!Y.CRANGE, linestyle=2, color=3
   
   js = [params[1], perror[1]]*N/(2*!pi)
   r  = [params[0], perror[0]]/N
   print, parinfo[goodChains[iPar]].parname, js, r, (js[0] GT 20)?'*':'-', (r[0] LT 0.01)?'*':'-', $
          FORMAT='(A20,1X,F7.1,"+-",F4.2,3X, G10.2,"+-",G10.2,2X,": ",2A)'

   multiplot
ENDFOR

multiplot,/default

END

PRO MHFIT_PELIPSE, center, covar, nsigma=nsigma, color=color
;; Helper function to plot an elipse

  IF NOT KEYWORD_SET(nsigma) THEN nsigma = [1]


  ;; Set up to draw ellipsis later
  nPoint= 100
  theta    = 2D0*!PI*dindgen(nPoint)/(nPoint-1)
  x_circle = cos(theta)
  y_circle = sin(theta)
  
;; Find eigen values and vector of the sub covariance matrix 
  eigenvalues = sqrt(eigenql(covar,eigenvectors=eigenvector,/double))
  
  FOR iSigma = 0, N_ELEMENTS(nsigma)-1 DO BEGIN 
     x_ellipse = center[0]+$
                 x_circle*nsigma[iSigma]*eigenvalues[0]*eigenvector[0,0]+ $
                 y_circle*nsigma[iSigma]*eigenvalues[1]*eigenvector[0,1]
     y_ellipse = center[1]+$
                 x_circle*nsigma[iSigma]*eigenvalues[0]*eigenvector[1,0]+ $
                 y_circle*nsigma[iSigma]*eigenvalues[1]*eigenvector[1,1]
     
     
     oplot,x_ellipse, y_ellipse, color=color
  ENDFOR


END

PRO MHFIT_PCHAIN, chains, parinfo, nsigma=nsigma, nbins = nbins, clouds=clouds, ct=ct,ROBUST=robust, FIT=fit, SWIDTH=swidth, PCOVAR=Pcovar, AUTO_LIMITS=auto_limits, ACCEPT=accept, LNL=LnL, ORIG=orig, INCOVAR=incovar
;; Ploting markov chains.... 

  IF NOT KEYWORD_SET(nsigma) THEN nsigma = [1,2,3]
  IF NOT KEYWORD_SET(nbins)  THEN nbins  = 200
  IF NOT KEYWORD_SET(ct)     THEN ct  = 0
  IF NOT KEYWORD_SET(swidth) THEN swidth = 1

  npar = N_ELEMENTS(chains[*,0])


  ;; Setting up color table (forcing color 0 to white)
  TVLCT, r_orig, g_orig, b_orig,/get
  LOADCT, ct,/SILENT
  TVLCT, r_ct, g_ct, b_ct,/get
  IF !D.NAME EQ 'PS' THEN Begin
     r_ct[0] = 255b
     g_ct[0] = 255b
     b_ct[0] = 255b
  ENDIF
  TVLCT, r_orig, g_orig, b_orig

  ;; Finding non fixed/tied parameters
  mpfit_parinfo, parinfo, tagnames, 'TIED',   ptied, default='', n=npar
  mpfit_parinfo, parinfo, tagnames, 'FIXED', pfixed, default=0,  n=npar
  goodChains = WHERE(pfixed EQ 0 AND STRTRIM(ptied,2) EQ '', nGoodChains)

  IF nGoodChains EQ 0 THEN $
     MESSAGE, 'EE - All parameters are either fixed or tied'
  
  parname = 'P('+strtrim(LINDGEN(npar),2)+')'
  IF  N_ELEMENTS(parinfo) GT 0 THEN BEGIN
     parinfo_tags = tag_names(parinfo)
     wh = where(parinfo_tags EQ 'PARNAME', ct)
     IF ct EQ 1 THEN BEGIN
        wh = where(parinfo.parname NE '', ct)
        IF ct GT 0 THEN $
           parname(wh) = strmid(parinfo(wh).parname,0,25)
     ENDIF
  ENDIF



  MHFIT_BSTAT, chains[goodChains,*], params, perror, covar, ROBUST=robust, FIT=FIT
  
  nSigma_plot = 6

  erase
  MULTIPLOT,[nGoodChains,nGoodChains],/ROWMAJOR
  FOR I=0, nGoodChains-1 DO BEGIN 
     
     ;; by defaults nSigma_plot sigma limits around the mean ...
     XRANGE = params[I] + [-1,1]*perror[I]*nSigma_plot

     ;; ... or use the limits from the parameters
     IF NOT KEYWORD_SET(auto_limits) THEN BEGIN
        mpfit_parinfo, parinfo, tagnames, 'LIMITED', limited, status=st1
        mpfit_parinfo, parinfo, tagnames, 'LIMITS',  limits,  status=st2
        
        IF st1 EQ 1 AND st2 EQ 1 THEN BEGIN
           parlim = WHERE(limited EQ 1, nLim)
           XRANGE[parlim] = limits[parlim]
        ENDIF
     ENDIF
     
     FOR J=0, nGoodChains-1 DO BEGIN 
        ;; by defaults nSigma_plot sigma limits around the mean ...
        YRANGE = params[J] + [-1,1]*perror[J]*nSigma_plot

        ;; ... or use the limits from the parameters
        IF NOT KEYWORD_SET(auto_limits) THEN BEGIN
           mpfit_parinfo, parinfo, tagnames, 'LIMITED', limited, status=st1
           mpfit_parinfo, parinfo, tagnames, 'LIMITS',  limits,  status=st2
           
           IF st1 EQ 1 AND st2 EQ 1 THEN BEGIN
              parlim = WHERE(limited EQ 1, nLim)
              YRANGE[parlim] = limits[parlim]
           ENDIF
        ENDIF

        IF I GT J THEN BEGIN
           ;; Upper triangle -> SKIP
           multiplot
           CONTINUE
        ENDIF

        IF I EQ 0 THEN $
           YTITLE=parname[goodChains[J]] $
        ELSE $
           YTITLE=""
   
        IF J EQ nGoodChains-1 THEN $
           XTITLE=parname[goodChains[I]] $
        ELSE $
           XTITLE=""
                         

        IF I EQ J THEN BEGIN
           YRANGE=[0,1]
           plot,[0],[0],/NODATA, $
                XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE
           ;; Draw the histogram...
           yy = histogram(chains[goodChains[I],*], NBINS=nbins, LOCATIONS=xx, MIN=MIN(XRANGE), MAX=MAX(XRANGE))
           oplot, xx, yy*1.D/MAX(yy),psym=10
           IF KEYWORD_SET(pcovar) THEN BEGIN 
              ;; overplot the corresponding gaussian
              xx = DINDGEN(nBins*10)/(nBins*10-1)*(MAX(XRANGE)-MIN(XRANGE))+MIN(XRANGE)
              yy = exp(-(xx-params[I])^2/(2*perror[I]^2))
              oplot,xx,yy, color=5
           ENDIF
           IF KEYWORD_SET(orig) THEN BEGIN
              ;; overplot starting value
              oplot, [1,1]*parinfo[goodChains[I]].value,[0,1],linestyle=1, color=6
              IF parinfo[goodChains[I]].limited[0] EQ 1 THEN $
                 oplot, [1,1]*parinfo[goodChains[I]].limits[0], !Y.CRANGE, linestyle=2,color=6
              IF parinfo[goodChains[I]].limited[1] EQ 1 THEN $
                 oplot, [1,1]*parinfo[goodChains[I]].limits[1], !Y.CRANGE, linestyle=2,color=6

              IF KEYWORD_SET(INCOVAR) THEN BEGIN
                 xx = DINDGEN(nBins*10)/(nBins*10-1)*(MAX(XRANGE)-MIN(XRANGE))+MIN(XRANGE)
                 yy = exp(-(xx-parinfo[goodChains[I]].value)^2/(2*incovar[goodChains[I], goodChains[I]]) )
                 oplot, xx, yy, color=6
              ENDIF
           ENDIF
           multiplot
        ENDIF
        IF I LT J THEN BEGIN
           plot,[0],[0],/NODATA, $
                XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE

           ;; plot the points 
           IF KEYWORD_SET(clouds) THEN BEGIN 
              oplot, chains[goodChains[I],*],chains[goodChains[J],*],psym=3
           ENDIF ELSE BEGIN
              width  = (MAX(XRANGE)-MIN(XRANGE))/nbins
              height = (MAX(YRANGE)-MIN(YRANGE))/nbins
              
              bitmap = HIST_2D(chains[goodChains[I],*],chains[goodChains[J],*], $
                               BIN1=width, BIN2=height, $
                               MIN1=MIN(XRANGE)+width/2,  MAX1=MAX(XRANGE)-width/2,$
                               MIN2=MIN(YRANGE)+height/2, MAX2=MAX(YRANGE)-height/2)
              
              IF swidth GT 1 THEN $
                 bitmap = SMOOTH(bitmap, swidth, /EDGE_TRUNCATE)
              bitmap[WHERE(bitmap EQ 0B)] = 0
              TVLCT, r_ct, g_ct, b_ct
              IMDISP, bitmap, /USEPOS, POSITION=!p.position, BOTTOM=0
              TVLCT, r_orig, g_orig, b_orig
              ;; REDRAW the plot
              plot,[0],[0],/NODATA, $
                   XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                   YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE
              
              ;; contour, bitmap, NLEVELS=3,XRANGE=XRANGE, /XSTYLE,YRANGE=YRANGE,/YSTYLE,POSITION=!p.position
              
           ENDELSE
           
           IF KEYWORD_SET(pcovar) THEN BEGIN
              
              ;; overplot an ellipse
              IF covar[I,J] NE 0 AND $
                 covar[J,I] NE 0 THEN BEGIN
                 
                 lcovar = [ [ covar[I,I], covar[J,I]],$
                            [ covar[I,J], covar[J,J]] ]
              
                 
                 MHFIT_PELIPSE, [params[I],params[J]], lcovar, nsigma=nsigma, color=5
              ENDIF
           ENDIF

           IF KEYWORD_SET(orig) THEN BEGIN
              ;; overplot starting value
              oplot, [parinfo[goodChains[I]].value],[parinfo[goodChains[J]].value],psym=2, color=6
              IF parinfo[goodChains[I]].limited[0] EQ 1 THEN $
                 oplot, [1,1]*parinfo[goodChains[I]].limits[0], !Y.CRANGE, linestyle=2,color=6
              IF parinfo[goodChains[I]].limited[1] EQ 1 THEN $
                 oplot, [1,1]*parinfo[goodChains[I]].limits[1], !Y.CRANGE, linestyle=2,color=6
              IF parinfo[goodChains[J]].limited[0] EQ 1 THEN $
                 oplot, !X.CRANGE, [1,1]*parinfo[goodChains[J]].limits[0], linestyle=2,color=6
              IF parinfo[goodChains[J]].limited[1] EQ 1 THEN $
                 oplot, !X.CRANGE, [1,1]*parinfo[goodChains[J]].limits[1], linestyle=2,color=6

              IF KEYWORD_SET(INCOVAR) THEN BEGIN
                 
                 IF incovar[goodChains[I],goodChains[J]] NE 0 AND $
                    incovar[goodChains[J],goodChains[I]] NE 0 THEN BEGIN
                    lcovar = [ [ incovar[goodChains[I],goodChains[I]], incovar[goodChains[J],goodChains[I]]],$
                               [ incovar[goodChains[I],goodChains[J]], incovar[goodChains[J],goodChains[J]]] ]
                    MHFIT_PELIPSE, [parinfo[goodChains[I]].value,parinfo[goodChains[J]].value], lcovar, nsigma=nsigma, color=6
                 ENDIF
                 
              ENDIF
           ENDIF

           multiplot
        ENDIF
        
     ENDFOR
  ENDFOR
  multiplot, /default

  IF KEYWORD_SET(accept) OR KEYWORD_SET(LnL) THEN BEGIN
   
     !p.charsize /= 2 
     multiplot,[8,8]
     FOR I=0, 5 DO $
        multiplot
     IF KEYWORD_SET(accept) THEN BEGIN
        IF NOT KEYWORD_SET(LnL) THEN $
           multiplot, /DOXAXIS, /DOYAXIS $
        ELSE $
           multiplot, /DOYAXIS
                       
        plot, FINDGEN(N_ELEMENTS(accept))/accept,/XSTYLE,YTITLE="Accept rate"
        FOR I=0, 6 DO multiplot
   ENDIF
     IF KEYWORD_SET(LnL) THEN BEGIN 
        multiplot, /DOXAXIS, /DOYAXIS
        plot, LnL,/XSTYLE,YTITLE="Log L"
     ENDIF
     multiplot,/default
     !p.charsize *= 2

ENDIF

END

PRO MHFIT_PCOV, params, covar, parinfo, nsigma=nsigma
;; Ploting covariance matrix

  IF NOT KEYWORD_SET(nsigma) THEN nsigma = [1,2,3]

  ;; To draw ellipsis later
  nPoint= 100
  theta    = 2D0*!PI*dindgen(nPoint)/(nPoint-1)
  x_circle = cos(theta)
  y_circle = sin(theta)

  npar = N_ELEMENTS(params)

  parname = 'P('+strtrim(LINDGEN(npar),2)+')'
  IF  N_ELEMENTS(parinfo) GT 0 THEN BEGIN
     parinfo_tags = tag_names(parinfo)
     wh = where(parinfo_tags EQ 'PARNAME', ct)
     IF ct EQ 1 THEN BEGIN
        wh = where(parinfo.parname NE '', ct)
        IF ct GT 0 THEN $
           parname(wh) = strmid(parinfo(wh).parname,0,25)
     ENDIF
  ENDIF


  erase
  MULTIPLOT,[npar,npar],/ROWMAJOR
  FOR I=0, npar-1 DO BEGIN 
     ;; 4 sigma
     XRANGE = params[I] + [-1,1]*4*sqrt(covar[I,I])
     FOR J=0, npar-1 DO BEGIN 
        YRANGE = params[J] + [-1,1]*4*sqrt(covar[J,J])
        IF I GT J THEN BEGIN
           ;; SKIP
           multiplot
           CONTINUE
        ENDIF

        IF I EQ 0 THEN $
           YTITLE=parname[J] $
        ELSE $
           YTITLE=""
   
        IF J EQ npar-1 THEN $
           XTITLE=parname[I] $
        ELSE $
           XTITLE=""
                         

        IF I EQ J THEN BEGIN
           YRANGE=[0,1]
           plot,[0],[0],/NODATA, $
                XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE

           ;; Draw a gaussian
           xx = DINDGEN(nPoint)/(nPoint-1)*(MAX(XRANGE)-MIN(XRANGE))+MIN(XRANGE)
           yy = exp(-(xx-params[I])^2/(2*covar[I,I]))
           oplot,xx,yy
           multiplot
           
        ENDIF
        IF I LT J THEN BEGIN
           plot,[0],[0],/NODATA, $
                XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE

           ;; Draw an ellipse
           lcovar = [ [ covar[I,I], covar[J,I]],$
                      [ covar[I,J], covar[J,J]] ]
;; Find eigen values and vector of the sub covariance matrix 
           eigenvalues = sqrt(eigenql(lcovar,eigenvectors=eigenvector,/double))
           
           FOR iSigma = 0, N_ELEMENTS(nsigma)-1 DO BEGIN 
              x_ellipse = params[0]+$
                          x_circle*nsigma[iSigma]*eigenvalues[0]*eigenvector[0,0]+ $
                          y_circle*nsigma[iSigma]*eigenvalues[1]*eigenvector[0,1]
              y_ellipse = params[1]+$
                          x_circle*nsigma[iSigma]*eigenvalues[0]*eigenvector[1,0]+ $
                          y_circle*nsigma[iSigma]*eigenvalues[1]*eigenvector[1,1]
              
              
              oplot,x_ellipse, y_ellipse
           ENDFOR
           multiplot
        ENDIF
        
     ENDFOR
  ENDFOR
  
multiplot, /default
END

PRO MHFIT_PCENT, chains, percentils, PERCENT=percent,QUIET=quiet, PARINFO=parinfo
;; compute percentil statistics

  IF NOT KEYWORD_SET(percent) THEN percent = [erfc(1./sqrt(2))/2, 0.5, 1-erfc(1./sqrt(2))/2]
 
  s = size(chains)
  nChain   = s[1]
  nIter    = s[2]
  nPercent = N_ELEMENTS(percent)

  percentils = DBLARR(nChain, nPercent)

  FOR iChain=0, nChain-1 DO BEGIN
     lchain = chains[iChain,*]
     lchain = lchain[SORT(lchain)]
     percentils[iChain,*] = lchain[(nIter+1)*percent]
  ENDFOR

  IF NOT KEYWORD_SET(quiet) AND N_ELEMENTS(percent) EQ 3 THEN BEGIN
     FOR iChain=0, nChain-1 DO BEGIN
        mean = percentils[iChain,1]
        err  = ABS(percentils[iChain,[2,0]]-mean)

        parname = 'P('+strtrim(iChain,2)+')'
        IF KEYWORD_SET(parinfo) THEN BEGIN
           parinfo_tags = tag_names(parinfo)
           wh = where(parinfo_tags EQ 'PARNAME', ct)
           if ct EQ 1 THEN BEGIN
              parname = strmid(parinfo[iChain].parname,0,25)
           endif
        ENDIF
        PRINT,  parname, mean, err, FORMAT='(A25," = ",G10.3, " + ", G10.3, " - ", G10.3)'
        
     ENDFOR
     
  ENDIF


END 

PRO MHFIT_BSTAT, chains, params, perror, covar, ROBUST=robust, FIT=FIT
;; Basic statistic on the chains (mean/variance/covariance)

  s       = size(chains)
  npar    = s[1]
  maxiter = s[2]
  
  params = DBLARR(npar)
  covar  = DBLARR(npar,npar)
  perror = DBLARR(npar)
  
  IF NOT KEYWORD_SET(robust) AND NOT KEYWORD_SET(fit) THEN BEGIN
     FOR I=0, npar-1 DO $
        params[I] = TOTAL(chains[I,*])/maxiter
     
     FOR I=0, npar-1 DO BEGIN
        FOR J=I, npar-1 DO BEGIN
                                ; Adapted from
                                ; Numerically-stable "two-pass" formula, which offers less
                                ; round-off error. Page 613, Numerical Recipes in C.
           resid_I = chains[I,*]-params[I]
           resid_J = chains[J,*]-params[J]
           covar[I,J] = ( TOTAL( resid_I*resid_J) - $
                          TOTAL(resid_I)*TOTAL(resid_J)/maxiter ) /(maxiter-1)
           covar[J,I] = covar[I,J]
        ENDFOR
     ENDFOR
     
     ;; Compute errors in parameters
     i = lindgen(npar)
     wh = where(covar(i,i) GE 0, ct)
     if ct GT 0 THEN $
        perror(wh) = sqrt(covar(wh, wh))
  ENDIF
  IF KEYWORD_SET(fit) THEN BEGIN 
     FOR I=0, npar-1 DO BEGIN
        yy = histogram(chains[I,*], NBINS=200,LOCATIONS=xx, MIN=MIN(chains[I,*]), MAX=MAX(chains[I,*]))
        dummy = gaussfit(xx, yy, results, nterms=3)
        params[I]  = results[1]
        perror[I]  = results[2]
        covar[I,I] = results[2]^2
     ENDFOR
  ENDIF

  IF KEYWORD_SET(robust) THEN BEGIN
     FOR I=0, npar-1 DO BEGIN
        mean = BIWEIGHT_MEAN(chains[I,*],sigma)
        params[I] = mean
        perror[I] = sigma
        covar[I,I] = sigma^2
     ENDFOR
     
  ENDIF
  

END

FUNCTION MHFIT_STEP, xall, covar, qulim, ulim, qllim, llim, ifree=ifree, NOCHOLESKY=NOCHOLESKY

  ;; Drawing values from the multivariate normal distribution
  ;; described by the covariance matrix covar
  ;; See http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution
  
  ;; Cholesky decomposition lcovar = LL^T if needed
  ;; Otherwise covar is supposed to be L already (avoid multiple decomposition)
  IF NOT KEYWORD_SET(noCholesky) THEN BEGIN
     lower_L = covar
     LA_CHOLDC, lower_L, /DOUBLE
     ;; Zero out the lower triangular portion to get L^T
     FOR I=0, N_ELEMENTS(covar[*,0])-2 DO lower_L[I,I+1:*] = 0
     covar = lower_L
  ENDIF

  if n_elements(ifree)  EQ 0 THEN ifree = lindgen(n_elements(xall))

  xnew = xall
  npar = N_ELEMENTS(xall)
  nfree = N_ELEMENTS(ifree)
  badStep  = 1

  buff_vec = DBLARR(npar)
  WHILE(badStep GT 0) DO BEGIN
     
     buff_free       = RANDOMN(seed, nfree)
     buff_vec[ifree] = buff_free
     step_vec        = covar#buff_vec

     ;; Reverse the sign of the step if we are up against the
     ;; parameters limit
     wh = WHERE( (qulim AND (xall[ifree] GT ulim-step_vec[ifree])) OR $
                 (qllim AND (xall[ifree] LT llim-step_vec[ifree])), badStep)
     IF badStep GT 0 THEN step_vec[wh] = -step_vec[wh]
     
     ;; Redraw the step if still outside bound
     wh = WHERE((qulim AND (xall[ifree] GT ulim-step_vec[ifree]) OR $
                 qllim AND (xall[ifree] LT llim-step_vec[ifree])), badStep)
  ENDWHILE

  RETURN, xall+step_vec
  

END

FUNCTION MHFIT, fcn, xall, INCOVAR=incovar, FUNCTARGS=fcnargs, $
                MAXITER=maxiter, COVAR=covar, perror=perror, $
                nprint=nprint, iterproc=iterproc, $
                PARINFO=parinfo, quiet=quiet, $
                QUERY=query, SCALE=scale, $
                CHAINS=chains, accept=accept, lnL = lnL, $
                SAVE_STEP=save_step, RESTORE_STEP=restore_step
;; Main function

  IF keyword_set(query) THEN return, 1

  IF n_elements(iterproc) EQ 0 THEN iterproc = 'MPFIT_DEFITER'
  IF n_elements(maxiter)  EQ 0 THEN maxiter  = 10000L
  IF n_elements(nprint)   EQ 0 THEN nprint   = maxiter/5

  IF KEYWORD_SET(restore_step) THEN RESTORE, restore_step

  IF n_params() EQ 0 OR NOT KEYWORD_SET(incovar) THEN BEGIN
      message, "USAGE: PARMS = MPFIT('MYFUNCT', START_PARAMS, INCOVAR=COVAR ... )", /info
      return, !values.d_nan
   ENDIF

  status = 0L
  errmsg = ''
  
  ;; Detect MPFIT and crash IF it was not found
  catch, catcherror
  IF catcherror NE 0 THEN BEGIN
     MPFIT_NOTFOUND:
     catch, /cancel
     message, 'ERROR: the required function MPFIT must be in your IDL path', /info
     return, !values.d_nan
  ENDIF
  IF mpfit(/query) NE 1 THEN goto, MPFIT_NOTFOUND
  catch, /cancel
  
  IF n_params() EQ 0 THEN BEGIN
     message, "USAGE: PARMS = MHFIT('MYFUNCT', X, Y, ERR, "+ $
              "START_PARAMS, ... )", /info
     return, !values.d_nan
  ENDIF
  
  common mpfit_config, mpconfig
  mpconfig = {fastnorm: keyword_set(fastnorm), proc: 0, nfev: 0L, damp: 0}
  common mpfit_machar, machvals
  
  
  
  ;; Parse FCN function name - be sure it is a scalar string
  sz = size(fcn)
  IF sz[0] NE 0 THEN BEGIN
     FCN_NAME:
     errmsg = 'ERROR: MYFUNCT must be a scalar string'
     goto, TERMINATE
  ENDIF
  IF sz[sz[0]+1] NE 7 THEN goto, FCN_NAME
  
  isext = 0
  catch_msg = 'parsing input parameters'
  ;; Parameters can either be stored in parinfo, or x.  Parinfo takes
  ;; precedence IF it exists.
  IF n_elements(xall) EQ 0 AND n_elements(parinfo) EQ 0 THEN BEGIN
     errmsg = 'ERROR: must pass parameters in P or PARINFO'
     goto, TERMINATE
  ENDIF
  
  ;; Be sure that PARINFO is of the right type
  IF n_elements(parinfo) GT 0 THEN BEGIN
     parinfo_size = size(parinfo)
     IF parinfo_size[parinfo_size[0]+1] NE 8 THEN BEGIN
        errmsg = 'ERROR: PARINFO must be a structure.'
        goto, TERMINATE
     ENDIF
     IF n_elements(xall) GT 0 AND n_elements(xall) NE n_elements(parinfo) $
     THEN BEGIN
        errmsg = 'ERROR: number of elements in PARINFO and P must agree'
        goto, TERMINATE
     ENDIF
  ENDIF
  
  ;; IF the parameters were not specIFied at the command line, THEN
  ;; extract them from PARINFO
  IF n_elements(xall) EQ 0 THEN BEGIN
     mpfit_parinfo, parinfo, tagnames, 'VALUE', xall, status=status
     IF status EQ 0 THEN BEGIN
        errmsg = 'ERROR: either P or PARINFO(*).VALUE must be supplied.'
        goto, TERMINATE
     ENDIF
     
     sz = size(xall)
     ;; Convert to double IF parameters are not float or double
     IF sz(sz(0)+1) NE 4 AND sz(sz(0)+1) NE 5 THEN $
        xall = double(xall)
  ENDIF
  npar = n_elements(xall)
  
  ;; TIED parameters?
  mpfit_parinfo, parinfo, tagnames, 'TIED', ptied, default='', n=npar
  ptied = strtrim(ptied, 2)
  wh = where(ptied NE '', qanytied) 
  qanytied = qanytied GT 0
  mpconfig = create_struct(mpconfig, 'QANYTIED', qanytied, 'PTIED', ptied)
  
  ;; FIXED parameters ?
  mpfit_parinfo, parinfo, tagnames, 'FIXED', pfixed, default=0, n=npar
  pfixed = pfixed EQ 1
  pfixed = pfixed OR (ptied NE '') ;; Tied parameters are also effectively fixed
  
  
  ;; Finish up the free parameters
  iFree = where(pfixed NE 1, nfree)
  IF nfree EQ 0 THEN BEGIN
     errmsg = 'ERROR: no free parameters'
     goto, TERMINATE
  ENDIF

  ;; scale factor for optimap step-size (see Dunkley et al 2004)
  IF NOT KEYWORD_SET(SCALE)    THEN scale = 2.4/sqrt(nfree)

  ;; scale is for the standard deviation, so square for the
  ;; covariances...
  lcovar = scale^2*incovar

  ;; Cholesky decomposition lcovar = LL^T
  ;; Done once for all future iteration as the a priori multivariate
  ;; normal distribution does not change with iteration
  lower_L = lcovar
  LA_CHOLDC, lower_L,/DOUBLE
  ;; Zero out the lower triangular portion to get L^T
  FOR I=0, N_ELEMENTS(lower_L[*,0])-2 DO lower_L[I,I+1:*] = 0
  

  ;; Compose only VARYING parameters
  xnew = xall      ;; xnew is the set of parameters to be returned
  x = xnew[iFree]  ;; x is the set of free parameters
  
  ;; LIMITED parameters ?
  mpfit_parinfo, parinfo, tagnames, 'LIMITED', limited, status=st1
  mpfit_parinfo, parinfo, tagnames, 'LIMITS',  limits,  status=st2
  IF st1 EQ 1 AND st2 EQ 1 THEN BEGIN
     
     ;; Error checking on limits in parinfo
     wh = where((limited[0,*] AND xall LT limits[0,*]) OR $
                (limited[1,*] AND xall GT limits[1,*]), ct)
     IF ct GT 0 THEN BEGIN
        errmsg = 'ERROR: parameters are not within PARINFO limits'
        goto, TERMINATE
     ENDIF
     wh = where(limited[0,*] AND limited[1,*] AND $
                limits[0,*] GE limits[1,*] AND $
                pfixed EQ 0, ct)
     IF ct GT 0 THEN BEGIN
        errmsg = 'ERROR: PARINFO parameter limits are not consistent'
        goto, TERMINATE
     ENDIF
     
     
     ;; Transfer structure values to local variables
     qulim = limited[1, iFree]
     ulim  = limits [1, iFree]
     qllim = limited[0, iFree]
     llim  = limits [0, iFree]
     
     wh = where(qulim OR qllim, ct)
     IF ct GT 0 THEN qanylim = 1 else qanylim = 0
     
  ENDIF else BEGIN
     
     ;; Fill in local variables with dummy values
     qulim = lonarr(nfree)
     ulim  = x * 0.
     qllim = qulim
     llim  = x * 0.
     qanylim = 0
     
  endelse
  
  ;; Initialize the number of parameters pegged at a hard limit value
  wh = where((qulim AND (x EQ ulim)) OR (qllim AND (x EQ llim)), npegged)
  
  n = n_elements(x)

  
  common mpfit_error, mperr
  
  mperr = 0
  catch_msg = 'calling '+fcn
  fvec = mpfit_call(fcn, xnew, _EXTRA=fcnargs)
  IFlag = mperr
  IF IFlag LT 0 THEN BEGIN
     errmsg = 'ERROR: first call to "'+fcn+'" failed'
     goto, TERMINATE
  ENDIF
  
  catch_msg = 'calling MPFIT_SETMACHAR'
  sz = size(fvec[0])
  isdouble = (sz[sz[0]+1] EQ 5)
  
  mpfit_setmachar, double=isdouble
  MACHEP0 = machvals.machep
  DWARF   = machvals.minnum
  
  szx = size(x)
  ;; The parameters and the squared deviations should have the same
  ;; type.  Otherwise the MACHAR-based evaluation will fail.
  catch_msg = 'checking parameter data'
  tp = szx[szx[0]+1]
  IF tp NE 4 AND tp NE 5 THEN BEGIN
     IF NOT keyword_set(quiet) THEN BEGIN
        message, 'WARNING: input parameters must be at least FLOAT', /info
        message, '         (converting parameters to FLOAT)', /info
     ENDIF
     x = float(x)
     xnew = float(x)
     szx = size(x)
  ENDIF
  IF isdouble AND tp NE 5 THEN BEGIN
     IF NOT keyword_set(quiet) THEN BEGIN
        message, 'WARNING: data is DOUBLE but parameters are FLOAT', /info
        message, '         (converting parameters to DOUBLE)', /info
     ENDIF
     x = double(x)
     xnew = double(xnew)
  ENDIF
  
  m = n_elements(fvec)
  IF (m LT n) THEN BEGIN
     errmsg = 'ERROR: number of parameters must not exceed data'
     goto, TERMINATE
  ENDIF
  
  iAccept = 1L  ;; The accepted iteration index
  iLoop   = 1L  ;; The actual iteration index

  ;; Start of the MCMC....
     
  ;; IF requested, call fcn to enable printing of iterates
  IF qanytied THEN mpfit_tie, xnew, ptied
  dof = (n_elements(fvec) - nfree) > 1L

  ;; chains contains parameters for all accepted iterations
  ;; lnL  containts the log likelihood for all accepted iterations
  ;; accept containts the actual iteration index, use to compute
  ;; acceptance ratio
  chains = DBLARR(npar, maxiter)
  lnL    = DBLARR(maxiter)
  accept = DBLARR(maxiter)

  ;; Likelihood defined as \propto exp (-chi2/2)
  ;; So log Likelihood : 
  prev = -1.D * (mpfit_enorm(fvec)^2) / 2

  chains[*,iAccept-1] = xnew
  accept[iAccept-1]   = iLoop
  lnL[iAccept-1]      = prev

  startTime = systime(/julian)

  IF KEYWORD_SET(restore_step) THEN BEGIN 
     RESTORE, restore_step
     good = WHERE(accept NE 0, nGood)
     iLoop   = good[nGood-1]
     iAccept = accept[iLoop]
     xnew    = chains[*,nGood-1]
  ENDIF

  WHILE(iAccept LT maxiter) DO BEGIN

     acceptRate = iAccept*1./iLoop


     IF nprint GT 0 AND iterproc NE '' THEN BEGIN
        catch_msg = 'calling '+iterproc
        IFlag = 0L
        IF iLoop MOD nprint EQ 0 THEN BEGIN

;; DEBUG
           IF KEYWORD_SET(save_step) THEN  $
              save, FILENAME=date_conv( startTime, 'FITS')+'_chains.dat', parinfo, incovar, fcnargs, chains, accept, lnL, startTime, maxiter, /COMPRESS
;; DEBUG
           

           mperr = 0
           xnew0 = xnew
           
           call_procedure, iterproc, fcn, xnew, iAccept, -2*prev, $
                           FUNCTARGS=fcnargs, parinfo=parinfo, quiet=quiet, $
                           dof=dof, _EXTRA=iterargs
           IF NOT KEYWORD_SET(QUIET) THEN $
              call_procedure, 'mpfit_defprint',  acceptRate , $
                              FORMAT='("acceptance ratio = " , G15.8)', _EXTRA=iterargs
           IFlag = mperr
           
           ;; Check for user termination
           IF IFlag LT 0 THEN BEGIN  
              errmsg = 'WARNING: premature termination by "'+iterproc+'"'
              goto, TERMINATE
           ENDIF
           
           ;; IF parameters were changed (grrr..) THEN re-tie
           IF max(abs(xnew0-xnew)) GT 0 THEN BEGIN
              IF qanytied THEN mpfit_tie, xnew, ptied
              x = xnew(iFree)
           ENDIF
           
        ENDIF
     ENDIF

     ;; xnew = mhfit_step(xall, lcovar, qulim, ulim, qllim, llim, iFree=iFree)
     xnew = mhfit_step(xall, lower_L, qulim, ulim, qllim, llim, iFree=iFree,/NOCHOLESKY)
     
     ;; Compute new likelihood
     mperr = 0
     catch_msg = 'calling '+fcn
     fvec = mpfit_call(fcn, xnew, _EXTRA=fcnargs)
     IFlag = mperr
     IF IFlag LT 0 THEN BEGIN
        errmsg = 'ERROR: first call to "'+fcn+'" failed'
        goto, TERMINATE
     ENDIF


     cur = -1.D * (mpfit_enorm(fvec)^2) / 2

;; L_cur / L_prev GT random

     IF ( (ALOG(RANDOMU(seed)) + prev - cur) LT 0 ) THEN BEGIN

        prev  = cur
        xall  = xnew
        iAccept += 1
        chains[*,iAccept-1] = xnew
        accept[iAccept-1]   = iLoop
        lnL[iAccept-1]      = prev

     ENDIF

     iLoop+=1
  ENDWHILE

  ;; Basic Chain statistic
  MHFIT_BSTAT, chains, params, perror, covar
  
  RETURN, params

TERMINATE:

END
