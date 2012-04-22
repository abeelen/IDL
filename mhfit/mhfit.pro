FORWARD_FUNCTION myfunct_fft, mhfit_step, mhfit, mpfit_call, mpfit_enorm, mpfit

FUNCTION MYFUNCT_FFT, p, k=K, Pk=pk, Err=err
  P_0   = p[0]
  k_s   = p[1]
  alpha = p[2]
  

  model = P_0*((k_s/k)^alpha)/((k_s/k)^alpha+1)
  return, (Pk-model)
END

PRO MHFIT_ACHAIN, chains, parinfo

s    = size(chains)
npar = s[1]
N    = s[2]

k     = 2*!pi*(FINDGEN(N/2+1))/N

fft_chains =  (ABS(FFT(chains, -1, DIMENSION=2)))^2*N
;; keep only the positive part
fft_chains = fft_chains[*,0:N/2+1]

erase
multiplot,[1,npar], mxtitle="k"
FOR iPar = 0, npar-1 DO BEGIN
   plot_oo, k[1:*], fft_chains[iPar,1:*],/XSTYLE,YTITLE="P(k) for "+parinfo[iPar].parname
   fcnargs = {k: k[1:*],Pk:fft_chains[iPar,1:*] , ERR:k[1:*]*0+1}
  
   parms = MPFIT('MYFUNCT_FFT',[fft_chains[iPar,1],1d-1, 1.5],FUNCTARGS=fcnargs,/QUIET)
   oplot, k, parms[0]*((parms[1]/k)^parms[2])/((parms[1]/k)^parms[2]+1),color=2
   multiplot
   print, parinfo[iPar].parname, parms[1]*N/(2*!pi)
ENDFOR

multiplot,/default

END

PRO MHFIT_PCHAIN, chains, parinfo, nsigma=nsigma, nbins = nbins, noclouds=noclouds, ct=ct,NOCOVAR=nocovar, SWIDTH=swidth
;; Ploting markov chains.... 

  IF NOT KEYWORD_SET(nsigma) THEN nsigma = [1,2,3]
  IF NOT KEYWORD_SET(nbins)  THEN nbins  = 200
  IF NOT KEYWORD_SET(ct)     THEN ct  = 0
  IF NOT KEYWORD_SET(swidth) THEN swidth = 1

  TVLCT, r_orig, g_orig, b_orig,/get
  LOADCT, ct
  TVLCT, r_ct, g_ct, b_ct,/get
  IF !D.NAME EQ 'PS' THEN Begin
     r_ct[0] = 255b
     g_ct[0] = 255b
     b_ct[0] = 255b
  ENDIF
     
  TVLCT, r_orig, g_orig, b_orig

  ;; To draw ellipsis later
  nPoint= 100
  theta    = 2D0*!PI*dindgen(nPoint)/(nPoint-1)
  x_circle = cos(theta)
  y_circle = sin(theta)


  MHFIT_BSTAT, chains, params, perror, covar

  npar = N_ELEMENTS(params)

  erase
  MULTIPLOT,[npar,npar],/ROWMAJOR
  FOR I=0, npar-1 DO BEGIN 
     ;; 4 sigma
     XRANGE = params[I] + [-1,1]*4*sqrt(covar[I,I])
     parlim = WHERE(parinfo[I].limited EQ 1, nLim)
     IF nLim NE 0 THEN $
        XRANGE[parlim] = parinfo[I].limits[parlim]
     
     FOR J=0, npar-1 DO BEGIN 
        YRANGE = params[J] + [-1,1]*4*sqrt(covar[J,J])
        parlim = WHERE(parinfo[J].limited EQ 1, nLim)
        IF nLim NE 0 THEN $
           YRANGE[parlim] = parinfo[J].limits[parlim]


        IF I GT J THEN BEGIN
           ;; SKIP
           multiplot
           CONTINUE
        ENDIF

        IF I EQ 0 THEN $
           YTITLE=parinfo[J].parname $
        ELSE $
           YTITLE=""
   
        IF J EQ npar-1 THEN $
           XTITLE=parinfo[I].parname $
        ELSE $
           XTITLE=""
                         

        IF I EQ J THEN BEGIN
           YRANGE=[0,1]
           plot,[0],[0],/NODATA, $
                XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE
           ;; Draw the histogram...
           plot_histo, chains[I,*], nbins=nbins, yy, xx,/NOPLOT
           oplot, xx, yy*1.D/MAX(yy),psym=10

           IF NOT KEYWORD_SET(nocovar) THEN BEGIN
              ;; overplot the corresponding gaussian
              xx = DINDGEN(nPoint)/(nPoint-1)*(MAX(XRANGE)-MIN(XRANGE))+MIN(XRANGE)
              yy = exp(-(xx-params[I])^2/(2*covar[I,I]))
              oplot,xx,yy, color=5
           ENDIF
           multiplot
        ENDIF
        IF I LT J THEN BEGIN
           plot,[0],[0],/NODATA, $
                XRANGE=XRANGE,/XSTYLE, XTITLE=XTITLE, $
                YRANGE=YRANGE,/YSTYLE, YTITLE=YTITLE

           ;; plot the points 
           IF NOT KEYWORD_SET(noclouds) THEN BEGIN 
              oplot, chains[I,*],chains[J,*],psym=3
           ENDIF ELSE BEGIN
              width  = (MAX(XRANGE)-MIN(XRANGE))/nbins
              height = (MAX(YRANGE)-MIN(YRANGE))/nbins
              
              bitmap = HIST_2D(chains[I,*],chains[J,*], $
                               BIN1=width, BIN2=height, $
                               MIN1=MIN(XRANGE), MAX1=MAX(XRANGE),$
                               MIN2=MIN(YRANGE), MAX2=MAX(YRANGE))
              
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

           IF NOT KEYWORD_SET(nocovar) THEN BEGIN
              
              ;; overplot an ellipse
              
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
                 
                 
                 oplot,x_ellipse, y_ellipse, color=5
              ENDFOR
           ENDIF
           multiplot
        ENDIF
        
     ENDFOR
  ENDFOR
  
multiplot, /default
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
           YTITLE=parinfo[J].parname $
        ELSE $
           YTITLE=""
   
        IF J EQ npar-1 THEN $
           XTITLE=parinfo[I].parname $
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

PRO MHFIT_BSTAT, chains, params, perror, covar
;; Basic statistic on the chains

  s       = size(chains)
  npar    = s[1]
  maxiter = s[2]
  
  params = DBLARR(npar)
  covar  = DBLARR(npar,npar)
  
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
  perror = replicate(abs(covar(0))*0., npar)
  wh = where(covar(i,i) GE 0, ct)
  if ct GT 0 then $
     perror(wh) = sqrt(covar(wh, wh))
  
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

  if n_elements(ifree)  EQ 0 then ifree = lindgen(n_elements(xall))

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
                  maxiter=maxiter, COVAR=covar, perror=perror, $
                  nprint=nprint, iterproc=iterproc, $
                  acceptRate=acceptRate, PARINFO=parinfo, quiet=quiet, $
                  QUERY=query, SCALE=scale, CHAINS=chains

  if keyword_set(query) then return, 1

  if n_elements(iterproc) EQ 0 then iterproc = 'MPFIT_DEFITER'
  if n_elements(maxiter)  EQ 0 then maxiter  = 10000L
  if n_elements(nprint)   EQ 0 then nprint   = maxiter/10

  if n_params() EQ 0 OR NOT KEYWORD_SET(incovar) then begin
      message, "USAGE: PARMS = MPFIT('MYFUNCT', START_PARAMS, INCOVAR=COVAR ... )", /info
      return, !values.d_nan
   endif

  status = 0L
  errmsg = ''
  
  ;; Detect MPFIT and crash if it was not found
  catch, catcherror
  if catcherror NE 0 then begin
     MPFIT_NOTFOUND:
     catch, /cancel
     message, 'ERROR: the required function MPFIT must be in your IDL path', /info
     return, !values.d_nan
  endif
  if mpfit(/query) NE 1 then goto, MPFIT_NOTFOUND
  catch, /cancel
  
  if n_params() EQ 0 then begin
     message, "USAGE: PARMS = MHFIT('MYFUNCT', X, Y, ERR, "+ $
              "START_PARAMS, ... )", /info
     return, !values.d_nan
  endif
  
  common mpfit_config, mpconfig
  mpconfig = {fastnorm: keyword_set(fastnorm), proc: 0, nfev: 0L, damp: 0}
  common mpfit_machar, machvals
  
  
  
  ;; Parse FCN function name - be sure it is a scalar string
  sz = size(fcn)
  if sz[0] NE 0 then begin
     FCN_NAME:
     errmsg = 'ERROR: MYFUNCT must be a scalar string'
     goto, TERMINATE
  endif
  if sz[sz[0]+1] NE 7 then goto, FCN_NAME
  
  isext = 0
  catch_msg = 'parsing input parameters'
  ;; Parameters can either be stored in parinfo, or x.  Parinfo takes
  ;; precedence if it exists.
  if n_elements(xall) EQ 0 AND n_elements(parinfo) EQ 0 then begin
     errmsg = 'ERROR: must pass parameters in P or PARINFO'
     goto, TERMINATE
  endif
  
  ;; Be sure that PARINFO is of the right type
  if n_elements(parinfo) GT 0 then begin
     parinfo_size = size(parinfo)
     if parinfo_size[parinfo_size[0]+1] NE 8 then begin
        errmsg = 'ERROR: PARINFO must be a structure.'
        goto, TERMINATE
     endif
     if n_elements(xall) GT 0 AND n_elements(xall) NE n_elements(parinfo) $
     then begin
        errmsg = 'ERROR: number of elements in PARINFO and P must agree'
        goto, TERMINATE
     endif
  endif
  
  ;; If the parameters were not specified at the command line, then
  ;; extract them from PARINFO
  if n_elements(xall) EQ 0 then begin
     mpfit_parinfo, parinfo, tagnames, 'VALUE', xall, status=status
     if status EQ 0 then begin
        errmsg = 'ERROR: either P or PARINFO(*).VALUE must be supplied.'
        goto, TERMINATE
     endif
     
     sz = size(xall)
     ;; Convert to double if parameters are not float or double
     if sz(sz(0)+1) NE 4 AND sz(sz(0)+1) NE 5 then $
        xall = double(xall)
  endif
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
  ifree = where(pfixed NE 1, nfree)
  if nfree EQ 0 then begin
     errmsg = 'ERROR: no free parameters'
     goto, TERMINATE
  endif
  
  IF NOT KEYWORD_SET(SCALE) THEN scale = 2.4/sqrt(nfree)
  lcovar = scale*incovar

  ;; Cholesky decomposition lcovar = LL^T
  ;; Done once for all future iteration as the a priori multivariate
  ;; normal distribution does not change with iteration
  lower_L = lcovar
  LA_CHOLDC, lower_L,/DOUBLE
  ;; Zero out the lower triangular portion to get L^T
  FOR I=0, N_ELEMENTS(lower_L[*,0])-2 DO lower_L[I,I+1:*] = 0
  

  ;; Compose only VARYING parameters
  xnew = xall      ;; xnew is the set of parameters to be returned
  x = xnew[ifree]  ;; x is the set of free parameters
  
  ;; LIMITED parameters ?
  mpfit_parinfo, parinfo, tagnames, 'LIMITED', limited, status=st1
  mpfit_parinfo, parinfo, tagnames, 'LIMITS',  limits,  status=st2
  if st1 EQ 1 AND st2 EQ 1 then begin
     
     ;; Error checking on limits in parinfo
     wh = where((limited[0,*] AND xall LT limits[0,*]) OR $
                (limited[1,*] AND xall GT limits[1,*]), ct)
     if ct GT 0 then begin
        errmsg = 'ERROR: parameters are not within PARINFO limits'
        goto, TERMINATE
     endif
     wh = where(limited[0,*] AND limited[1,*] AND $
                limits[0,*] GE limits[1,*] AND $
                pfixed EQ 0, ct)
     if ct GT 0 then begin
        errmsg = 'ERROR: PARINFO parameter limits are not consistent'
        goto, TERMINATE
     endif
     
     
     ;; Transfer structure values to local variables
     qulim = limited[1, ifree]
     ulim  = limits [1, ifree]
     qllim = limited[0, ifree]
     llim  = limits [0, ifree]
     
     wh = where(qulim OR qllim, ct)
     if ct GT 0 then qanylim = 1 else qanylim = 0
     
  endif else begin
     
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
  iflag = mperr
  if iflag LT 0 then begin
     errmsg = 'ERROR: first call to "'+fcn+'" failed'
     goto, TERMINATE
  endif
  
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
  if tp NE 4 AND tp NE 5 then begin
     if NOT keyword_set(quiet) then begin
        message, 'WARNING: input parameters must be at least FLOAT', /info
        message, '         (converting parameters to FLOAT)', /info
     endif
     x = float(x)
     xnew = float(x)
     szx = size(x)
  endif
  if isdouble AND tp NE 5 then begin
     if NOT keyword_set(quiet) then begin
        message, 'WARNING: data is DOUBLE but parameters are FLOAT', /info
        message, '         (converting parameters to DOUBLE)', /info
     endif
     x = double(x)
     xnew = double(xnew)
  endif
  
  m = n_elements(fvec)
  if (m LT n) then begin
     errmsg = 'ERROR: number of parameters must not exceed data'
     goto, TERMINATE
  endif
  
  ;; Likelihood defined as \propto exp (-chi2/2)
  ;; So log Likelihood : 
  lnL = -1.D * (mpfit_enorm(fvec)^2) / 2
  iter = 1L
  iLoop = 1L

  ;; Start of the MCMC....
     
  ;; If requested, call fcn to enable printing of iterates
  if qanytied then mpfit_tie, xnew, ptied
  dof = (n_elements(fvec) - nfree) > 1L

  prev = lnL
  chains = DBLARR(npar, maxiter)
  chains[*,iter-1] = xnew
  
  WHILE(iter LT maxiter) DO BEGIN

     acceptRate = iter*1./iLoop

     if nprint GT 0 AND iterproc NE '' then begin
        catch_msg = 'calling '+iterproc
        iflag = 0L
        if iLoop MOD nprint EQ 0 then begin
           mperr = 0
           xnew0 = xnew
           
           call_procedure, iterproc, fcn, xnew, iter, -2*lnL, $
                           FUNCTARGS=fcnargs, parinfo=parinfo, quiet=quiet, $
                           dof=dof, _EXTRA=iterargs
           IF NOT KEYWORD_SET(QUIET) THEN $
              call_procedure, 'mpfit_defprint',  acceptRate , $
                              FORMAT='("acceptance ratio = " , G15.8)', _EXTRA=iterargs
           iflag = mperr
           
           ;; Check for user termination
           if iflag LT 0 then begin  
              errmsg = 'WARNING: premature termination by "'+iterproc+'"'
              goto, TERMINATE
           endif
           
           ;; If parameters were changed (grrr..) then re-tie
           if max(abs(xnew0-xnew)) GT 0 then begin
              if qanytied then mpfit_tie, xnew, ptied
              x = xnew(ifree)
           endif
           
        endif
     endif

     ;; xnew = mhfit_step(xall, lcovar, qulim, ulim, qllim, llim, ifree=ifree)
     xnew = mhfit_step(xall, lower_L, qulim, ulim, qllim, llim, ifree=ifree,/NOCHOLESKY)
     
     ;; Compute new likelihood
     mperr = 0
     catch_msg = 'calling '+fcn
     fvec = mpfit_call(fcn, xnew, _EXTRA=fcnargs)
     iflag = mperr
     if iflag LT 0 then begin
        errmsg = 'ERROR: first call to "'+fcn+'" failed'
        goto, TERMINATE
     endif

     lnL = -1.D * (mpfit_enorm(fvec)^2) / 2
     cur = lnL

;; L_cur / L_prev GT random

     IF ( (ALOG(RANDOMU(seed)) + prev - cur) LT 0 ) THEN BEGIN

        prev  = cur
        xall  = xnew
        iter += 1
        chains[*,iter-1] = xnew

     ENDIF

     iLoop+=1
  ENDWHILE

  ;; Basic Chain statistic
  MHFIT_BSTAT, chains, params, perror, covar
  
  RETURN, params

TERMINATE:

END
