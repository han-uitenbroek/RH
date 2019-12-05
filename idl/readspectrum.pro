; -------- file: -------------------------- readspectrum.pro --------- ;

; -------- begin -------------------------- readSpectrum.pro --------- ;

FUNCTION readSpectrum, fileName

;+
; NAME:
;	READSPECTRUM
;
; PURPOSE:
;	This routine reads the output spectral data FOR the RH radiative
;       transfer routines RHF1D, RHSC2D, RHCC3D, RHSPHERE.
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       succes = READSPECTRUM(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read.
;       asrs.out:  Filename with active set record structure (Default,
;                  hardwired). this quantity is the ActiveSet analogon
;                  to the Background Record Structure (see READBRS.PRO).
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	spectrum:  (Communicated via spectrumcommon COMMON block).
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       geometrycommon, atmoscommon, spectrumcommon, filescommon
; SIDE EFFECTS:
; RESTRICTIONS:
;       The geometry of the atmosphere must be known (via geometrycommon)
;       and it must be known whether the atmosphere is moving or not
;       (via atmos in atmoscommon).
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Tue Jun 30 10:52:25 2009 --
;-

@input.common
@geometry.common
@atmos.common
@spectrum.common
@files.common


  WHILE (NOT existFile(fileName, UNIT=spectrumUnit, /XDR)) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new spectrum file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.out', $
                               TITLE='Spectroscopic data file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0

  ENDWHILE

  ;; --- Store name in common block Files_Common --     -------------- ;

  spectrumFile = fileName

  Nspect = 0L
  readu, spectrumUnit, Nspect

  ;; --- Structure definitions --                       -------------- ;

  CASE geometryType OF 
    "ONE_D_PLANE": $
     spectrum = {Nspect: Nspect,  lambda: dblarr(Nspect), $
                 I: dblarr(geometry.Nrays, Nspect), $
                 vacuum_to_air: 0L,  air_limit: 0.0D+0}
    "TWO_D_PLANE": $
     spectrum = {Nspect: Nspect,  lambda: dblarr(Nspect), $
                 I: dblarr(geometry.Nx, geometry.Nrays, Nspect), $
                 vacuum_to_air: 0L,  air_limit: 0.0D+0}
    "THREE_D_PLANE": $
     spectrum = {Nspect: Nspect,  lambda: dblarr(Nspect), $
                 I: dblarr(geometry.Nx, geometry.Ny, $
                           geometry.Nrays, Nspect), $
                 vacuum_to_air: 0L,  air_limit: 0.0D+0}
    "SPHERICAL_SYMMETRIC": $
     spectrum = {Nspect: Nspect,  lambda: dblarr(Nspect), $
                 I: dblarr(geometry.Nrays, Nspect), $
                 vacuum_to_air: 0L,  air_limit: 0.0D+0}
  ENDCASE

  ;; --- Read in one go --                              -------------- ;

  point_lun, spectrumUnit, 0
  readu, spectrumUnit, spectrum

  IF (tag_present(atmos, 'STOKES') OR inputData.backgr_pol) THEN BEGIN

    ;; --- Read Stokes profiles if magnetic field present -- --------- ;

    Ntags = n_tags(spectrum)
    spectrum = create_struct(spectrum, 'stokes_Q', spectrum.I, $
                             'stokes_U', spectrum.I, 'stokes_V', spectrum.I)
    readstructure, spectrum, spectrumUnit, N_START=Ntags, N_MAX_TAGS=3
  ENDIF

  free_lun, spectrumUnit

  ;; --- Permute indices so that spectral index runs fastest, and
  ;;     angular index runs slowest --                  -------------- ;

  CASE (geometryType) OF
    "ONE_D_PLANE": BEGIN
      spectrum = replace_tag(spectrum, 'I', transpose(spectrum.I))
      IF (tag_present(atmos, 'STOKES') OR inputData.backgr_pol) THEN BEGIN
        spectrum = $
         replace_tag(spectrum, 'STOKES_Q', transpose(spectrum.Stokes_Q))
        spectrum = $
         replace_tag(spectrum, 'STOKES_U', transpose(spectrum.Stokes_U))
        spectrum = $
         replace_tag(spectrum, 'STOKES_V', transpose(spectrum.Stokes_V))
      ENDIF
    END
    "TWO_D_PLANE": BEGIN
      spectrum = replace_tag(spectrum, $
                             'I', reform(transpose(spectrum.I, [2, 0, 1])))

     IF (tag_present(atmos, 'STOKES') OR inputData.backgr_pol) THEN BEGIN
        spectrum = replace_tag(spectrum, 'STOKES_Q', $
                               reform(transpose(spectrum.Stokes_Q, [2, 0, 1])))
        spectrum = replace_tag(spectrum, 'STOKES_U', $
                               reform(transpose(spectrum.Stokes_U, [2, 0, 1])))
        spectrum = replace_tag(spectrum, 'STOKES_V', $
                               reform(transpose(spectrum.Stokes_V, [2, 0, 1])))
      ENDIF
    END
    "THREE_D_PLANE": BEGIN
      spectrum = replace_tag(spectrum, 'I', reform(spectrum.I))
      IF (tag_present(atmos, 'STOKES') OR inputData.backgr_pol) THEN BEGIN
        spectrum = replace_tag(spectrum, 'STOKES_Q', reform(spectrum.Stokes_Q))
        spectrum = replace_tag(spectrum, 'STOKES_U', reform(spectrum.Stokes_U))
        spectrum = replace_tag(spectrum, 'STOKES_V', reform(spectrum.Stokes_V))
      ENDIF
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      spectrum = replace_tag(spectrum, 'I', transpose(spectrum.I))
     END
   ENDCASE

  ;; Check is needed when opacity.out is absent (possible when
  ;; filename for OPACITY_OUT is set to "none" in keyword.input

  files = file_search('asrs.out', COUNT=Nfile)
  IF (Nfile GT 0) THEN BEGIN
    openr, unit, 'asrs.out', /GET_LUN, /XDR
    IF (atmos.moving OR tag_present(atmos, 'STOKES') OR $
        inputData.PRD_angle_dep) THEN $
     as_rn = lonarr(spectrum.Nspect * geometry.Nrays) $
    ELSE $
     as_rn = lonarr(spectrum.Nspect)
    readu, unit, as_rn
    spectrum = create_struct(spectrum, 'as_rn', as_rn)
    free_lun, unit
  ENDIF
  
  return, 1
END
; -------- end ---------------------------- readSpectrum.pro --------- ;

