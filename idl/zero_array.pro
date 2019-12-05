FUNCTION zero_array, a, INDEFINITE=indefinite

  IF (keyword_set(INDEFINITE)) THEN BEGIN
    s  = size(a)
    Na = n_elements(s)
    return, reform(make_array(s[Na-1], VALUE=0.0, TYPE=s[Na-2]), $
                   s[1:s[0]], /OVERWRITE)
  ENDIF ELSE $
   return, a - a
END
