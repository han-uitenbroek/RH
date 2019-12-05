FUNCTION read2datmos, fileName, BFILE=Bfile
;+
; NAME:
;	READ2DATMOS
;
; PURPOSE:
;	This routine reads 2-dimensional model atmospheres in RHSC2D format.
;       (see rhsc2d/readatmos.c for format definition).
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       atmos = READ2DATMOS(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read.
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	atmos:  Model atmosphere in atmos structure.
;                (see routine READATMOS.PRO).
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Sep 15 16:42:10 2004 --
;-

  Nx = 0L  &  Nz = 0L  &  NHydr = 0L

  openr, unit, fileName, /GET_LUN, /XDR

  readu, unit, Nx, Nz, NHydr
  point_lun, unit, 0

  atmos = {Nx: Nx,  Nz: Nz,  NHydr: NHydr, $
           boundary: lonarr(3),  dx: dblarr(Nx),  z: dblarr(Nz), $
           T: dblarr(Nx, Nz),  n_elec: dblarr(Nx, Nz),  $
           vturb: dblarr(Nx, Nz),  vx: dblarr(Nx, Nz),  $
           vz: dblarr(Nx, Nz), $
           nH: dblarr(Nx, Nz, NHydr)}

  readu, unit, atmos
  free_lun, unit

  IF (keyword_set(BFILE)) THEN BEGIN
    openr, unit2, Bfile, /GET_LUN, /XDR

    Bdata = dblarr(Nx, Nz, 3)
    readu, unit2, Bdata
    free_lun, unit2

    atmos = create_struct(atmos, 'B', Bdata[*, *, 0], $
                          'gamma', Bdata[*, *, 1], 'chi', Bdata[*, *, 2])
  ENDIF


  return, atmos
END
