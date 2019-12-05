; -------- file: -------------------------- existfile.pro ------------ ;

; -------- begin -------------------------- existfile.pro ------------ ;

FUNCTION existFile, fileName, UNIT=unit, XDR=xdr, SWAP_ENDIAN=swap_endian

  ;; --- Check whether file exists, and if so open it for read only -- ;

  ON_IOERROR, error
  openr, unit, fileName, /GET_LUN, XDR=keyword_set(XDR), $
         SWAP_ENDIAN=keyword_set(SWAP_ENDIAN)
  IF (NOT keyword_set(UNIT)) THEN free_lun, unit
  return, 1

error:
  return, 0
END
; -------- end ---------------------------- existfile.pro ------------ ;
