FUNCTION filename_sdss, dr3struct
;+
; NAME:
;         filename_sdss
;
;
; PURPOSE:
;         construct the fits filename containing a spectra from one
;         structure returned by MRDFITS on the SDSS QSO Catalog
;         (contain the .plate, .MJD_SPEC, .fiberid elements)
;
;
; CATEGORY:
;         I/O
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;        dr3struct : a structure containing .plate .MJD_SPEC and
;                    .fiberid elements
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
;        the name of the fits file containing that spectra
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
;         2006 - AB - created
;-

return, STRING(dr3struct.plate,FORMAT='(I4.4)')+$
  "/1d/spSpec-"+STRING(dr3struct.MJD_SPEC,FORMAT='(I5.5)')+$
  "-"+STRING(dr3struct.plate,FORMAT='(I4.4)')+$
  "-"+STRING(dr3struct.fiberid,FORMAT='(I3.3)')+$
  ".fit"
END
