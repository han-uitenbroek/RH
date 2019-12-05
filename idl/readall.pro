;+
; NAME:
;	READALL
;
; PURPOSE:
;	This routine reads geometry, atmosphere, background, atomic data,
;       and spectral output data.
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       .r readall
;
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	geometry, atmos, spectrum, atom
;                (all communicated via COMMON blocks).
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       filescommon, geometrycommon, atmoscommon, spectrumcommon
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Fri Apr 17 16:24:53 2009 --
;-

  print, "-- Reading inputData ...."
  IF (readInput('input.out')) THEN $
   print, "     Success" $
  ELSE $
   print, "     Failed"
@input.common

  print, "-- Reading geometry ...."
  IF (readgeometry('geometry.out')) THEN $
   print, "     Success" $
  ELSE $
   print, "     Failed"
@geometry.common

  print, "-- Reading atmosphere ...."

@files.common

  metalFile  = "metals.out"  &  moleculeFile = "molecules.out"
  IF (readatmos('atmos.out')) THEN $
   print, "     Success" $
  ELSE $
   print, "     Failed"
@atmos.common

  print, "-- Reading spectrum ...."
  IF (readspectrum('spectrum.out')) THEN $
   print, "     Success" $
  ELSE $
   print, "     Failed"
  print, "-- Reading flux ...."
  IF (readflux('flux.out')) THEN $
   print, "     Success" $
  ELSE $
   print, "     Failed"
@spectrum.common

  print, "-- Reading atomic models ...."
  atom_files = file_search("atom.*.out", COUNT=Nfile)
  print, "Found ", Nfile, " atoms: "
  IF (Nfile GT 0) THEN BEGIN
    print, "  atom_files: ", atom_files
    print, " Read atoms with: atom = readatom(atom_files[n])"
  ENDIF

END
