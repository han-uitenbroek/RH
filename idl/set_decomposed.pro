PRO set_decomposed, COLOR_TABLE=color_table

  IF (!D.NAME EQ 'X') THEN BEGIN
    device, GET_VISUAL_DEPTH=depth
    IF (depth GT 8) THEN device, DECOMPOSED=0
  ENDIF ELSE $
   device, DECOMPOSED=1

  IF (keyword_set(COLOR_TABLE)) THEN loadct, color_table
END

