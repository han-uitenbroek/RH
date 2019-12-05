FUNCTION strsplit, string, pattern, _ref_extra=extra

  release = float(strmid(!VERSION.RELEASE, 0, 3))
  IF (release LT 5.3) THEN BEGIN
    IF (n_params() EQ 1) THEN pattern = ' '
    items = str_sep(string, pattern)
    nonempty = where(items NE '', count)
    IF (count GT 0) THEN items = items(nonempty)
    return, items
  ENDIF ELSE BEGIN
    forward_FUNCTION strtok
    IF (n_params() EQ 1) THEN return, strtok(string, _extra=extra) $
    ELSE return, strtok(string, pattern, _extra=extra)
  ENDELSE
END

