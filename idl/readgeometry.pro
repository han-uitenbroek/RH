; -------- file: -------------------------- readgeometry.pro---------- ;

; -------- begin -------------------------- readGeometry.pro --------- ;

FUNCTION readGeometry, fileName

;+
; NAME:
;	READGEOMETRY
;
; PURPOSE:
;	This routine reads the geometry output file from the RH radiative
;       transfer programs (RHF1D, RHSC2D, RHCC3D, RHSPHERE). 
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       succes = READGEOMETRY(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read.
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	geometry:  (Passed via the geometrycommon COMMON block).
;                  Geometry of the model atmosphere. Format is defined
;                  in readgeometry.c.
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       geometrycommon, filescommon
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Fri Mar  5 14:32:43 1999 --
;-

@geometry.common
@files.common

  WHILE (NOT existFile(fileName, UNIT=geometryUnit, /XDR)) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new geometry file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.out', TITLE='Geometric data file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0
  ENDWHILE

  ;; --- Store name in common block Files_Common --     -------------- ;

  geometryFile = fileName

  type = 0L  &  Nrays = 0L
  readu, geometryUnit, type, Nrays

  ;; --- Structure definitions --                       -------------- ;

  CASE type OF
    0: BEGIN
      geometryType = "ONE_D_PLANE"
      Ndep = 0L
      readu, geometryUnit, Ndep
      geometry = { Nrays: Nrays,  Ndep: Ndep, $
                   xmu: dblarr(Nrays),  wmu: dblarr(Nrays), $
                   height: dblarr(Ndep),  cmass: dblarr(Ndep), $
                   tau500: dblarr(Ndep), vz: dblarr(Ndep) }
    END
    1: BEGIN
      geometryType = "TWO_D_PLANE"
      Nx = 0L  &  Nz = 0L
      readu, geometryUnit, Nx, Nz
      geometry = { Nrays: Nrays,  Nx: Nx,  Nz: Nz,  angleSet: 0L, $
                   xmu: dblarr(Nrays),  ymu: dblarr(Nrays), $
                   wmu: dblarr(Nrays), x: dblarr(Nx), z: dblarr(Nz), $
                   vx: dblarr(Nx, Nz), vz: dblarr(Nx, Nz) }
    END
    3: BEGIN
      geometryType = "THREE_D_PLANE"
      Nx = 0L  &  Ny = 0L  &  Nz = 0L
      readu, geometryUnit, Nx, Ny, Nz
      geometry = { Nrays: Nrays,  Nx: Nx,  Ny: Ny,  Nz: Nz,  angleSet: 0L, $
                   xmu: dblarr(Nrays),  ymu: dblarr(Nrays), $
                   wmu: dblarr(Nrays), $
                   dx: 0.0D+0,  dy: 0.0D+0,  z: dblarr(Nz), $
                   vx: dblarr(Nx, Ny, Nz), vy: dblarr(Nx, Ny, Nz), $
                   vz: dblarr(Nx, Ny, Nz)}
    END
    2: BEGIN
      geometryType = "SPHERICAL_SYMMETRIC"
      Nradius = 0L  &  Ncore = 0L
      readu, geometryUnit, Nradius, Ncore
      geometry = { Nrays: Nrays, Nradius: Nradius,  Ncore: Ncore, $
                   Radius: 0.0D+0,  xmu: dblarr(Nrays),  wmu: dblarr(Nrays), $
                   r: dblarr(Nradius),  cmass: dblarr(Nradius), $
                   tau500: dblarr(Nradius), vr: dblarr(Nradius) }
    END
    ELSE: BEGIN
      typeStr = string(type, FORMAT='(I3)')
      answer = dialog_message(/ERROR, "Not a supported geometry: " + typeStr)
      return, 0
    END
  ENDCASE

  ;; --- Read rest in one go --                         -------------- ;

  point_lun, geometryUnit, 0
  readu, geometryUnit, type, geometry
  free_lun, geometryUnit
  return, 1
END
; -------- end ---------------------------- readGeometry.pro --------- ;
