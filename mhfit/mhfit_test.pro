FUNCTION MYFUNCT, p, X=x, Y=y, ERR=err
                                ; Parameter values are passed in "p"
  model = P[0]+P[1]*x
  return, (y-model)/err
END


n = 10
x = INDGEN(n)
y = 2*x+1                       ;+RANDOMN(seed,n)
s = FLTARR(n)+1

p0 = [1.,1]
fa = {X:x, Y:Y, ERR:s}
incovar = [[1.,0], $
           [0.,1.]]

result = mhfit('myfunct',p0, incovar=incovar, functargs = fa, $
               perror=perror, covar=outcovar,/QUIET)
print, result
print, perror

result2 = mhfit('myfunct',result, incovar=outcovar, functargs = fa, $
                chains=chains, covar=outcovar2,/QUIET)

mhfit_pchain, chains


parinfo  = REPLICATE({parname: '', value:0.D, fixed:0, limited:[0,0], limits:[0.D,0],tied:''},2)
parinfo.value = result2
parinfo[0].parname = "origin"
parinfo[1].parname = "slope"
parinfo[0].fixed = 1

result3 = mhfit('myfunct',parinfo=parinfo, incovar=outcovar2, functargs = fa, $
                chains=chains, covar=outcovar3,/QUIET)

mhfit_pchain, chains, parinfo, /auto,/pcovar


END
