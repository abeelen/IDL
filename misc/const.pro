;+
; NAME:
;      CONST
;
;
; PURPOSE:
;      Define some physical constant in SI
;
;
; CATEGORY:
;      Definition
;
;
; CALLING SEQUENCE:
;      @const
;
;
; INPUTS:
;      None.
;
;
; OPTIONAL INPUTS:
;      None.
;
;
; KEYWORD PARAMETERS:
;      None.
;
;
; OUTPUTS:
;      systems variables
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;      None.
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
;      A. Beelen, 2003            IAS
;-

DEFSYSV, '!const',1                                          ; is constant already present

DEFSYSV, '!c_lum',2.99792458D8                               ; Velocity of light (in vacuo) (m s^{-1})
DEFSYSV, '!G', 6.672D-11                                     ; Gravitational constant (m^3 kg^{-1} s^{-2})
DEFSYSV, '!AU',1.4959787D11                                  ; Astronomical Unit of Distance (m)
DEFSYSV, '!pc',3.0857D16                                     ; Parsec (=AU/sin 1") (m)
DEFSYSV, '!M_sol',1.9891D30                                  ; Solar Mass (kg)
DEFSYSV, '!L_sol',3.85D26                                    ; Solar Luminosity (W)
DEFSYSV, '!mb_sol',4.74                                      ; Absolute bolometric magnitude of the sun (Allen's Astro Quantities)
DEFSYSV, '!h_c',6.6260693D-34                                   ; Planck's Constant (J s = m^2 kg s^{-2})
DEFSYSV, '!k_b',1.3806505D-23                                   ; Boltzmann constant (J K^{-1} =  m^2 kg s^{-3} K^{-1})
DEFSYSV, '!sigma',5.6696D-8                                  ; Stefan Boltzmann constant (W m^{-2} K^{-4})

DEFSYSV, '!erg',1.D-7                                        ; erg unit (J)
DEFSYSV, '!jansky',1.D-26                                    ; janski unit (W m^{-2} Hz^{-1})
DEFSYSV, '!pc2m',3.0857678D16                                ; parsec to meter coefficient
DEFSYSV, '!arcsec2Deg',60.^(-2)                              ; arcsec to deg
DEFSYSV, '!arcsec2radian',60.^(-2)*180.^(-1)*!Dpi            ; arcsec to radian
DEFSYSV, '!arcsec2sr',60.^(-4)*180.^(-2)*!Dpi^(2)            ; arcsec2 to steradian
DEFSYSV, '!deg2sr',180.^(-2)*!Dpi^(2)                        ; deg2 to steradian
