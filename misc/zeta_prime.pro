FUNCTION zeta_prime, s
;+
; NAME:
;       ZETA_PRIME
;
;
; PURPOSE:
;       compute the derivative of the zeta function
;
;
; CATEGORY:
;       math
;
;
; CALLING SEQUENCE:
;       result = zeta_prime(s)
;
;
; INPUTS:
;       s 
;
;
; OPTIONAL INPUTS:
;       None.
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
;          A. Beelen, 2002,        IAS
;-

Rs = DOUBLE(s)
IF (Rs LE 1) THEN BEGIN
PRINT, "Error : Real part of argument must be greater than 1" 
RETURN,0
ENDIF

;o = 1.d-12
; On majore totalement le nombre de terme pour avoir une certaine precision
;number_of_term = 5*CEIL((1./o)^(1./s))

number_of_term = 65535U
serie_k = (DINDGEN(number_of_term)+2)
returned = DINDGEN(N_ELEMENTS(s))
IF N_ELEMENTS(returned) GT 1 THEN BEGIN
    FOR I=0,N_ELEMENTS(returned)-1 DO BEGIN
        returned[I] = -TOTAL(ALOG(serie_k)/((serie_k)^(s[I])))
    ENDFOR 
ENDIF ELSE BEGIN
    returned = -TOTAL(ALOG(serie_k)/((serie_k)^(s)))
ENDELSE



RETURN, returned
END
