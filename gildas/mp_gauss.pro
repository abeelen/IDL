PRO GAUSS_plus_CONT,X,p,F

; p[0] = INT
; p[1] = position offset
; p[2] = FWHM
; p[3] = Continuum

sigma = P[2]/(2*SQRT(2*ALOG(2)))
Z = (X-P[1])/sigma
EZ = EXP(-Z^2/2.)

F = P[0]*EZ/(sigma*SQRT(2*!pi))+P[3]

END

FUNCTION MP_GAUSS, p ,X = x,Y = y ,ERR = err

gauss_plus_cont,x, p, model

return, (y-model)/err
END



