; -------- file: -------------------------- readatmos.pro ------------ ;

; -------- begin -------------------------- readAtmos.pro ------------ ;

FUNCTION readAtmos, fileName

;+
; NAME:
;	READDATMOS
;
; PURPOSE:
;	This routine reads atmosphere output files from the RH radiative
;       transfer programs (RHF1D, RHSC2D, RHCC3D, RHSPHERE).
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       succes = READDATMOS(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read.
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;       succes: IF the file exists AND is read succesfully 1 is returned,
;               0 otherwise.
;	atmos:  (Passed via the atmoscommon COMMON block) Model atmosphere
;               in atmos structure.
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       geometrycommon, atmoscommon, filescommon
;
; SIDE EFFECTS:
; RESTRICTIONS:
;       The file's geometry should be known via the geometry structure
;       in the geometrycommon COMMON block (see READGEOMETRY.PRO)
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Fri Dec  2 14:05:48 2005 --
;-

@geometry.common
@atmos.common
@files.common


  WHILE (NOT existFile( fileName, UNIT=atmosUnit, /XDR )) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new atmosphere file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.out', TITLE='Atmospheric data file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0

  ENDWHILE

  ;; --- Store name in common block Files_Common --     -------------- ;

  atmosFile = fileName
  NHydr = 0L  &  Nelem = 0L
  readu, atmosUnit, NHydr, Nelem

  ;; --- Structure definitions --                       -------------- ;

  elem = {id: string(FORMAT='(A2)', ''), weight: 0.0D+0,  abund: 0.0D+0}

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      Ndep = geometry.Ndep
      atmos = {NHydr: NHydr,  Nelem: Nelem, moving: 0L, $
               T: dblarr(Ndep),  n_elec: dblarr(Ndep), $
               vturb: dblarr(Ndep),  nH: dblarr(Ndep, NHydr),  ID: '', $
               elements: replicate(elem, Nelem)}
    END
    "TWO_D_PLANE": BEGIN
      Nx = geometry.Nx  &  Nz = geometry.Nz
      atmos = {NHydr: NHydr,  Nelem: Nelem, moving: 0L, $
               T: dblarr(Nx, Nz),  n_elec: dblarr(Nx, Nz), $
               vturb: dblarr(Nx, Nz),  nH: dblarr(Nx, Nz, NHydr),  ID: '', $
               elements: replicate(elem, Nelem)}
    END
    "THREE_D_PLANE": BEGIN
      Nx = geometry.Nx  &  Ny = geometry.Ny  &  Nz = geometry.Nz
      atmos = {NHydr: NHydr,  Nelem: Nelem, moving: 0L, $
               T: dblarr(Nx, Ny, Nz),  n_elec: dblarr(Nx, Ny, Nz), $
               vturb: dblarr(Nx, Ny, Nz),  nH: dblarr(Nx, Ny, Nz, NHydr), $
               ID: '',  elements: replicate(elem, Nelem)}
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      Nradius = geometry.Nradius
      atmos = {NHydr: NHydr,  Nelem: Nelem, moving: 0L, $
               T: dblarr(Nradius),  n_elec: dblarr(Nradius), $
               vturb: dblarr(Nradius),  nH: dblarr(Nradius, NHydr),  ID: '', $
               elements: replicate(elem, Nelem)}
    END
  ENDCASE
  
  ;; --- Read in one go --                              -------------- ;

  point_lun, atmosUnit, 0
  readu, atmosUnit, atmos

  ;; --- If magnetic fields were present in the 1-D case -- ---------- ;

  IF (geometryType NE "SPHERICAL_SYMMETRIC") THEN BEGIN
    Ntags = n_tags(atmos)
    on_ioerror, no_Stokes
    stokes = 0L
    readu, atmosUnit, stokes
    atmos = create_struct(atmos, 'stokes', stokes, 'B', atmos.T, $
                          'gamma_B', atmos.T, 'chi_B', atmos.T)

    ;; --- Read the field and the two position angles -- ------------- ;

    readstructure, atmos, atmosUnit, N_START=Ntags+1, N_MAX_TAGS=3
    on_ioerror, NULL
  ENDIF
  
no_stokes:
  free_lun, atmosUnit

  readMetals, metalFile
  IF (NOT readMolecules(moleculeFile)) THEN molecules = 0

  ;; --- Read background record structure --            -------------- ;

  readbrs

  return, 1
END
; -------- end ---------------------------- readAtmos.pro ------------ ;
