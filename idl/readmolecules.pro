; -------- begin -------------------------- readMolPops.pro ---------- ;

FUNCTION readMolPops, molecule, fileName

@geometry.common
@files.common

  IF (n_tags(molecule) LE 0) THEN return, 0

  CASE geometryType OF
    "ONE_D_PLANE":         nv = dblarr(geometry.Ndep, molecule.Nv)
    "SPHERICAL_SYMMETRIC": nv = dblarr(geometry.Nradius, molecule.Nv)

    "TWO_D_PLANE":         nv = dblarr(geometry.Nx, geometry.Nz, $
                                       molecule.Nv)
    "THREE_D_PLANE":       nv = dblarr(geometry.Nx, geometry.Ny, $
                                       geometry.Nz, molecule.Nv)
  ENDCASE

  f = file_search(fileName, COUNT=Nf)
  IF (Nf EQ 0) THEN return, 0

  openr, lun, fileName, /XDR, /GET_LUN

  atmosID = ''
  Nvibr = 0L  &  Nspace = 0L
  readu, lun, atmosID, Nvibr, Nspace, nv
  IF (Nvibr NE molecule.Nv) THEN BEGIN
    print, "Error: Nvibr != molecule.Nv in readMolPops for molecule ", $
     molecule.ID
    return, 0
  ENDIF

  molecule.nv_ptr = ptr_new(nv)
  readu, lun, nv
  molecule.nvstar_ptr = ptr_new(nv)

  free_lun, lun
  return, 1
END
; -------- end ---------------------------- readMolPops.pro ---------- ;

; -------- begin -------------------------- readmolecules.pro -------- ;

FUNCTION readmolecules, fileName

;+
; NAME:
;	READMOLECULES
;
; PURPOSE:
;	This routine reads molecular output data.
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       succes = READMOLECULES(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read (usually molecules.out).
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;       succes:     Returns 1 IF fileName exists AND was succesfully read,
;                   0 otherwise.
;	molecules:  Array of molecule structures (passed via COMMON block
;                   atmoscommon).
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       geometry.common, atmos.common
; SIDE EFFECTS:
; RESTRICTIONS:
;       The geometry of the atmosphere should be known (via geometrycommon).
;       See also READGEOMETRY.PRO.
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Sep 20 13:39:59 2006 --
;-

@geometry.common
@atmos.common

  WHILE (NOT existFile(fileName, UNIT=molUnit, /XDR)) DO BEGIN
   answer = dialog_message( /QUESTION, "Find molecular data file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.out', TITLE='Molecular data file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0

  ENDWHILE

  CASE geometryType OF
    "ONE_D_PLANE":         n = dblarr(geometry.Ndep)
    "SPHERICAL_SYMMETRIC": n = dblarr(geometry.Nradius)
    "TWO_D_PLANE":         n = dblarr(geometry.Nx, geometry.Nz)
    "THREE_D_PLANE":       n = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)
  ENDCASE


  Nmolecules = 0L
  readu, molUnit, Nmolecules

  molstruct = {ID: '', Nv: 0L, NJ: 0L, $
               Ediss: 0.0D+0, Tmin: 0.0D+0, Tmax: 0.0D+0, n_ptr: ptr_new(), $
               E_ptr: ptr_new(), nv_ptr: ptr_new(), nvstar_ptr: ptr_new()}

  molecules = replicate(molstruct, Nmolecules)
  FOR m=0, Nmolecules-1 DO BEGIN
    readstructure, molstruct, molUnit, N_MAX_TAGS=6
    molecules[m] = molstruct
    readu, molUnit, n
    molecules[m].n_ptr = ptr_new(n)

    IF (molecules[m].Nv GT 0  AND  molecules[m].NJ GT 0) THEN BEGIN
      E = dblarr(molecules[m].NJ, molecules[m].Nv)
      readu, molUnit, E
      molecules[m].E_ptr = ptr_new(E, /NO_COPY)

      popsFile = string(FORMAT='("pops_mol.", A, ".out")', molecules[m].ID)
      r = readMolPops(molstruct, popsFile)
      molecules[m].nv_ptr = molstruct.nv_ptr
      molecules[m].nvstar_ptr = molstruct.nvstar_ptr
    ENDIF
  ENDFOR

  nHmin = n
  readu, molUnit, nHmin

  free_lun, molUnit
  return, 1
END
; -------- end ---------------------------- readmolecules.pro -------- ;
