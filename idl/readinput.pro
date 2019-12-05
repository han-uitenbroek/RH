; -------- file: -------------------------- readinput.pro ------------ ;

; -------- begin -------------------------- readInput.pro ------------ ;

FUNCTION readInput, filename

;+
; NAME:
;	READINPUT
;
; PURPOSE:
;	This routine reads various input data from filee
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;       succes = READINPUT(fileName)
;
; INPUTS:
;	fileName:  Name of the file to be read.
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;	input:  (Communicated via inputcommon COMMON block).
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
;       inputcommon
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Tue Feb 13 15:18:57 2007 --
;-

@input.common

 inputData = {magneto_optical: 0L,  PRD_angle_dep: 0L,  XRD: 0L, $
              start_solution: 0L, Stokes_mode: 0L,  metallicity: 0.0D0, $
              backgr_pol: 0L, big_endian: 1L}

  WHILE (NOT existFile(filename, UNIT=inputUnit, /XDR)) DO BEGIN
    answer = dialog_message(/QUESTION, "Find new inputData file?")
    IF (answer EQ 'Yes') THEN BEGIN
      fileName = dialog_pickfile(FILTER='*.out', $
                                 TITLE='Input data file', $
                                 /MUST_EXIST, /READ, FILE=fileName)
      IF (fileName EQ '') THEN return, 0
    ENDIF ELSE $
     return, 0
  ENDWHILE

  on_ioerror, insufficient_data
  readu, inputUnit, inputData
  free_lun, inputUnit
  return, 1

insufficient_data:

  ;; --- Use Endianness of current machine in case of doubt -- ----- ;;

  inputData.big_endian = is_big_endian()

  print, "WARNING: Insufficient input data. Old version of writeinput_xdr.c ?"
  print, "Using:"
  help, /STRUCTURE, inputData
  return, 1

END
; -------- end ---------------------------- readInput.pro ------------ ;
