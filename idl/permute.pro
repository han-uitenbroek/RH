pro PERMUTE, a, CLOCKWISE=clockwise

  sz = size(a)
  ndim = sz[0]  &  n1 = sz[1]  &  n2 = sz[2]  &  n3 = sz[3]

  IF (keyword_set(CLOCKWISE)) THEN BEGIN
    a = reform(a, n1*n2, n3)
    a = reform(transpose(a), n3, n1, n2)
  ENDIF ELSE begin
    a = reform(a, n1, n2*n3)
    a = reform(transpose(a), n2, n3, n1)
  ENDELSE
END
