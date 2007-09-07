FUNCTION zeta, s
;+
; NAME:
;       ZETA
;
;
; PURPOSE:
;       compute the zeta function
;
;
; CATEGORY:
;       math
;
;
; CALLING SEQUENCE:
;       result = zeta(s)
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
IF (WHERE(Rs LE 1))[0] NE -1 THEN BEGIN
PRINT, "Error : Real part of argument must be greater than 1" 
RETURN,0
ENDIF

;o = 1.d-12

; On majore totalement le nombre de terme pour avoir une certaine precision
;number_of_term = 5*CEIL((1./o)^(1./s))

number_of_term = 65535U
serie_k = (DINDGEN(number_of_term)+1)
returned = DINDGEN(N_ELEMENTS(s))

IF N_ELEMENTS(returned) GT 1 THEN BEGIN
    FOR I=0,N_ELEMENTS(returned)-1 DO BEGIN
        returned[I] = TOTAL(serie_k^(-s[I]))
    ENDFOR 
ENDIF ELSE BEGIN
    returned = TOTAL(serie_k^(-s))
ENDELSE



RETURN,returned 
END
