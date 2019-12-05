; -------- file: -------------------------- readflow.pro ------------- ;

; -------- begin -------------------------- readFlow.pro ------------- ;

FUNCTION readFlow, fileName

@geometry.common
@spectrum.common

  WHILE (NOT existFile(fileName, UNIT=flowUnit, /XDR)) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new flow file file?")
   IF (answer EQ 'Yes') THEN BEGIN
     fileName = dialog_pickfile(FILTER='flow_*', $
                                TITLE='Spectroscopic data file', $
                                /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '') THEN return, 0
   ENDIF ELSE $
    return, 0
  ENDWHILE

  IF (geometryType NE "TWO_D_PLANE") THEN return, 0

  Nx = geometry.Nx  &  Nz = geometry.Nz
  flow = {nspect: 0L,  Fx: dblarr(Nx, Nz),  Fz: dblarr(Nx, Nz)}

  readu, flowUnit, flow
  free_lun, flowUnit

  return, flow
END
; -------- end ---------------------------- readFlow.pro ------------- ;
