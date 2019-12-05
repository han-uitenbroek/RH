PRO readBRS
;+
; NAME:
;	READBRS
;
; PURPOSE:
;	This routine reads the record structure for the background
;       opacities. For each wavelength index this structure specifies
;       the record number of the backgound opacities (all spatial locations)
;       for that wavelength.
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       READBRS
;
; INPUTS:
;	brs.out (by default, hardwired).
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	atmos.backgr_rn:  Record numbers of background opacity (and
;                         emissivity and scattering coefficient) for
;                         each wavelength index. The opacity data itself
;                         is usually stored in background.dat
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       geometrycommon, atmoscommon
; SIDE EFFECTS:
; RESTRICTIONS:
;       The geometry of the atmosphere should be known (passed via COMMON block
;       geometrycommon), and the atmos structure should exist.
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Thu Jan  3 14:51:14 2008 --
;-

@geometry.common
@atmos.common

  openr, unit, "brs.out", /GET_LUN, /XDR

  atmosID = ''
  Nspace = 0L  &  Nspect = 0L
  readu, unit, atmosID, Nspace, Nspect

  hasline     = intarr(Nspect)
  ispolarized = intarr(Nspect)

  IF (atmos.moving OR tag_present(atmos, 'STOKES')) THEN $
   backgrrecno = lonarr(2*Nspect * geometry.Nrays) $
  ELSE $
   backgrrecno = lonarr(Nspect)

  readu, unit, hasline, ispolarized, backgrrecno

  backgrflags = $
   replicate({hasline: hasline[0], ispolarized: ispolarized[0]}, Nspect)
  backgrflags.hasline     = hasline
  backgrflags.ispolarized = ispolarized

  atmos = create_struct(atmos, 'backgrflags', backgrflags, $
                        'backgrrecno', backgrrecno)
  free_lun, unit
END
