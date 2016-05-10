FUNCTION MYFUNCT_2, p, X=x, Y=y, ERR=err
                                ; Parameter values are passed in "p"
  model = P[0]+P[1]*x
  return, (y-model)/err
END

FUNCTION MYFUNCT_3, p, X=x, Y=y, ERR=err
                                ; Parameter values are passed in "p"
  model = P[0]+P[1]*x+P[2]*x^2
  return, (y-model)/err
END

FUNCTION MYPRIOR_2, p, _EXTRA=_extra
  RETURN, 0
END

FUNCTION MYPRIOR_3, p, _EXTRA=_extra

  IF P[1] LT 1 OR P[1] GT 5 THEN $
     RETURN, !VALUES.F_NAN
  
  IF ABS(P[2]) GT 1d-2 THEN $
     RETURN, !VALUES.F_NAN

  mu_p0 = 0.98
  sigma_p0 = 0.1
  prior_p0 = 1./(sqrt(2*!pi)*sigma_p0)*EXP(-(P[0]-mu_p0)^2/(2*sigma_p0^2))

  RETURN, prior_p0
  ;; prior_p2 = 0; -0.0001
  ;; sigma_p2 = 0.01
  ;; prior_p2 = 1./(sqrt(2*!pi)*sigma_p2)*EXP(-(P[2]-prior_p2)^2/(2*sigma_p2^2))

  ;; RETURN, (prior_p0+prior_p2)/2
  
END

FUNCTION GEN_PRIOR_3, N
  
  mu_p0 = 0.98
  sigma_p0 = 0.1
  
  p0 = RANDOMN(seed, N)*sigma_p0+mu_p0
  
  min_p1 = 1
  max_p1 = 5
  
  p1 = RANDOMU(seed,N)*(max_p1-min_p1)+min_p1
  
  min_p2 = -1d-2
  max_p2 = +1d-2
  
  p2 = RANDOMU(seed,N)*(max_p2-min_p2)+min_p2
  
  RETURN, TRANSPOSE([[p0],[p1],[p2]])
  
END


FUNCTION GEN_PRIOR_2, N
  
  mu_p0 = 0.98
  sigma_p0 = 0.1
  
  p0 = RANDOMN(seed, N)*sigma_p0+mu_p0
  
  min_p1 = 1
  max_p1 = 5
  
  p1 = RANDOMU(seed,N)*(max_p1-min_p1)+min_p1
  
  RETURN, TRANSPOSE([[p0],[p1]])
  
END


RESOLVE_ROUTINE,'mhfit',/IS_FUNCTION,/COMPILE_FULL_FILE

n = 50
x = FINDGEN(n)
s = FLTARR(n)+5
noise = RANDOMN(2, n)*s

y = -0.001*x^2 + 2*x+1  

p0 = [1.,1]
fa = {X:x, Y:Y+noise, ERR:s}
incovar_2 = [[1.,0], $
           [0.,1.]]


params_2 = mhfit('myfunct_2',p0, incovar=incovar_2, functargs = fa, $
               perror=perror_2, covar=outcovar_2,/QUIET, maxiter=1000)

params_2b = mhfit('myfunct_2',params_2, incovar=outcovar_2, functargs = fa, $
                  perror=perror_2b, chains=chains_2b, covar=outcovar_2b,/QUIET, maxiter=2000, LnL=LnL_2b)

parinfo_2  = REPLICATE({parname: '', value:0.D, fixed:0, limited:[0,0], limits:[0.D,0],tied:''},2)
parinfo_2.value = params_2b
parinfo_2[0].parname = "origin"
parinfo_2[1].parname = "slope"
parinfo_2[0].fixed = 0

window,2
mhfit_pchain, chains_2b, parinfo_2, /auto,/pcovar,/oresult, LnL=LnL_2b

;; parinfo_2[0].fixed = 1
;; params_2c = mhfit('myfunct_2', parinfo=parinfo_2, incovar=outcovar_2b, functargs = fa, $
;;                   perror=perror_2b, chains=chains_2c, covar=outcovar_2c,/QUIET, maxiter=2000)

;; window,1
;; mhfit_pchain, chains_2c, parinfo_2, /auto,/pcovar


N_particules = 1000
logZ_2 = NESTED_sampling( 'myfunct_2', GEN_PRIOR_2(N_particules), functargs=fa, incovar=outcovar_2b, maxiter=50, information=information_2, posterior_sample=posterior_sample_2)



;;; 3 parameters

parinfo_3  = REPLICATE({parname: '', value:0.D, fixed:0, limited:[0,0], limits:[0.D,0],tied:''},3)
parinfo_3[[0,1]].value = params_2b
parinfo_3[0].parname = "origin"
parinfo_3[1].parname = "slope"
parinfo_3[2].parname = "square"

incovar_3 = [[1.,0.   ,0.],$
             [0.,1.d-3,0.],$
             [0.,0.   ,1.d-4]]

params_3 = mhfit('myfunct_3',parinfo=parinfo_3, incovar=incovar_3, functargs = fa, $
                 perror=perror_3, chains=chains_3, covar=outcovar_3,/QUIET, maxiter=2000)

params_3b = mhfit('myfunct_3',parinfo=parinfo_3, incovar=outcovar_3, functargs = fa, $
                 perror=perror_3b, chains=chains_3b, covar=outcovar_3b,/QUIET, maxiter=2000, LnL=Lnl_3b)


window,3
mhfit_pchain, chains_3, parinfo_3, /auto,/pcovar,/oresult, LnL=LnL_3b


;; params_3_p = mhfit('myfunct_3', parinfo=parinfo_3, incovar=outcovar_3, functargs = fa, $
;;                   perror=perror_3_p, chains=chains_3_p, covar=outcovar_3_p,/QUIET, maxiter=2000, LnL=LnL_3_p, logLikelihoodPrior_fcn='myprior_3')

;; mhfit_pchain, chains_3_p, parinfo_3, /auto,/pcovar,/oresult, LnL=LnL_3_p

logZ_3 = NESTED_sampling( 'myfunct_3', GEN_PRIOR_3(N_particules), functargs=fa, incovar=outcovar_3b, maxiter=50, information=information_3, posterior_sample=posterior_sample_3)


window,0
ploterror, x,y+noise,s,psym=2
oplot, x,y,color=5
oplot, x, params_2[0]+x*params_2[1],color=2
oplot, x, params_2b[0]+x*params_2b[1],color=3
oplot, x, params_3b[0]+x*params_3b[1]+x^2*params_3b[2], color=4
;; oplot, x, params_2c[0]+x*params_2c[1],color=4

PRINT, "Order 1"
FOR I=0, 1 DO BEGIN 
   print, params_2b[I], "+-", perror_2b[I]
ENDFOR
PRINT, "CHI2R : ", TOTAL(MYFUNCT_2( params_2b, X=fa.x, Y=fa.y, ERR=fa.err)^2)/(N-2)
PRINT, "logZ : ", logZ_2, " +- ", SQRT(information_2/N_particules)

PRINT, "Order 2"
FOR I=0, 2 DO BEGIN 
   print, params_3b[I], "+-", perror_3b[I]
ENDFOR
PRINT, "CHI2R : ", TOTAL(MYFUNCT_3( params_3b, X=fa.x, Y=fa.y, ERR=fa.err)^2)/(N-3)
PRINT, "logZ : ", logZ_3, " +- ", SQRT(information_3/N_particules)


;; PRINT, "Order 2 w prior"
;; FOR I=0, 2 DO BEGIN 
;;    print, params_3_p[I], "+-", perror_3_p[I]
;; ENDFOR
;; PRINT, "CHI2R : ", TOTAL(MYFUNCT_3( params_3_p, X=fa.x, Y=fa.y, ERR=fa.err)^2)/(N-3)



END
