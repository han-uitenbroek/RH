FUNCTION read3datmos, fileName, BFILE=Bfile
;+
; NAME:
;	READ3DATMOS
;
; PURPOSE:
;	This routine reads 3-dimensional model atmospheres in RHCC3D format.
;       (see rhcc2d/readatmos.c for format definition).
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       atmos = READ3DATMOS(fileName)
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
;   --- Last modified: Wed Sep 16 16:36:33 2009 --
;-

  Nx = 0L  &  Ny = 0L  &  Nz = 0L  &  NHydr = 0L

  openr, unit, fileName, /GET_LUN, /XDR

  readu, unit, Nx, Ny, Nz, NHydr
  point_lun, unit, 0

  atmos = {Nx: Nx,  Ny: Ny,  Nz: Nz,  NHydr: NHydr, $
           boundary: lonarr(2),  dx: 0.0D0,  dy: 0.0D0,  z: dblarr(Nz), $
           T: dblarr(Nx, Ny, Nz),  n_elec: dblarr(Nx, Ny, Nz),  $
           vturb: dblarr(Nx, Ny, Nz),  vx: dblarr(Nx, Ny, Nz),  $
           vy: dblarr(Nx, Ny, Nz),  vz: dblarr(Nx, Ny, Nz), $
           nH: dblarr(Nx, Ny, Nz, NHydr)}

  readu, unit, atmos
  free_lun, unit

  IF (keyword_set(BFILE)) THEN BEGIN
    openr, unit2, Bfile, /GET_LUN, /XDR

    Bdata = dblarr(Nx, Ny, Nz, 3)
    readu, unit2, Bdata
    free_lun, unit2

    atmos = create_struct(atmos, 'B', Bdata[*, *, *, 0], $
                          'gamma', Bdata[*, *, *, 1], $
                          'chi', Bdata[*, *, *, 2])
  ENDIF

  return, atmos
END
