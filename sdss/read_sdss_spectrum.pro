PRO read_sdss_spectrum, filename, spectrum, linetable
;+
; NAME:       read_sdss_spectrum
;
;
;
; PURPOSE:    read a sdss spectrum into a IDL structure used by other routine
;
;
;
; CATEGORY:   I/O
;
;
;
; CALLING SEQUENCE:
;             read_sdss_spectrum, filename, spectrum, [linetable]
;
;
; INPUTS:
;             filename : the input filename with the complete path
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
;             spectrum  : the IDL structure
;
;
;
; OPTIONAL OUTPUTS:
;             linetable : the list of line detected in the spectrum in
;                         IDL structure
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
;        MRDFITS, SXPAR (ASTRON Libs)
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;         2006 - AB - created
;-


fitsfile = MRDFITS(filename,0,header,/silent)

spectrum            = fitsfile[*,0]
spectrumNOcontinuum = fitsfile[*,1]
error               = fitsfile[*,2]
bitmask             = fitsfile[*,3]

COEFF0 = SXPAR(header,'COEFF0')
COEFF1 = SXPAR(header,'COEFF1')
Z      = SXPAR(header,'Z')
wavelength = 10.^(COEFF0+COEFF1*INDGEN(N_ELEMENTS(spectrum)))

spectrum = {$
             header:          header,$
             wavelength:      wavelength,$
             count:           spectrum,$
             counterr:        error,$
             countNObaseline: spectrumNOcontinuum, $
             redshift:        Z $

           }

; Line ID from http://www.sdss.org/dr2/dm/flatFiles/spSpec.html

restWave = [1033.82, 1215.67, 1240.81, 1305.53, 1335.31, $
            1399.8, 1549.48, 1640.4, 1665.85, 1857.4, $
            1908.734, 2326.0, 2439.5, 2799.117, 3346.79, $
            3426.85, 3727.092, 3729.875, 3798.976, 3836.47, $
            3889.0, 3934.777, 3969.588, 4072.3, 4102.89, $
            4305.61, 4341.68, 4364.436, 4862.68, 4960.295, $
            5008.240, 5176.7, 5895.6, 6302.046, 6365.536, $
            6549.86, 6564.61, 6585.27, 6707.89, 6718.29, $
            6732.67, 8500.36, 8544.44, 8664.52 ]

Name = ['OVI', 'Ly_a', 'NV', 'OI', 'CII', $
        'SiIV+OIV', 'CIV', 'HeII', 'OIII', 'AlIII', $
        'CIII', 'CII', 'NeIV', 'MgII', 'NeV', $
        'NeV', 'OII', 'OII', 'Hh', 'Hy', $
        'HeI', 'K', 'H', 'SII', 'Hd',$
        'G' , 'Hg', 'OIII', TEXTOIDL('H_\beta'), 'OIII', $
        'OIII', 'Mg', 'Na', 'OI', 'OI', $
        'NII', 'Ha', 'NII', 'Li', 'SII', $
        'SII', 'CaII', 'CaII', 'CaII']

; read the fitted line
fitstable = MRDFITS(filename,2,header,/silent)
linetable = REPLICATE({line,$
                       name:      " ",$
                       wave:      0.0,$
                       waveerr:   0.0,$
                       sigma:     0.0,$
                       sigmaerr:  0.0,$
                       height:    0.0,$
                       heighterr: 0.0,$
                       z:         0.0,$
                       zeer:      0.0 $
                      }$
                      ,N_ELEMENTS(fitstable))

FOR I=0,N_ELEMENTS(fitstable)-1 DO BEGIN
    linetable[I] = {name:      Name[WHERE(FITSTABLE[I].RESTWAVE EQ restWave)],$
                    wave:      fitstable[I].wave,$
                    waveerr:   fitstable[I].waveerr,$
                    sigma:     fitstable[I].sigma,$
                    sigmaerr:  fitstable[I].sigmaerr,$
                    height:    fitstable[I].height,$
                    heighterr: fitstable[I].heighterr,$
                    z:         fitstable[I].z,$
                    zeer:      fitstable[I].zerr $
                   }

ENDFOR


END
