FUNCTION tag_present, s, tag, POSITION=position

  ;; --- Determine if tag is present in structure s --  ---------;;

  IF (n_tags(s) LE 0) THEN BEGIN
    print, "Expression must be structure: S"
stop
    return, 0;
  ENDIF

  n = where(tag_names(s) EQ tag, count)
  IF (count GT 0) THEN BEGIN
    position = n[0]
    return, 1 
  ENDIF ELSE $
   return, 0
END
