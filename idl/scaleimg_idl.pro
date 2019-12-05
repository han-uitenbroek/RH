FUNCTION scaleimg_idl, img, Ncol, Nrow, CUBIC=cubic

  imgSz = size(img)
  d1 = float(imgSz(1))  &  d2 = float(imgSz(2))

  IF (keyword_set(CUBIC)) THEN $
   return, poly_2d(img, [0, 0, (d1 - 1.0)/(Ncol-1), 0], $
                   [0, (d2 - 1.0)/(Nrow-1), 0, 0], 2, Ncol, Nrow, /CUBIC) $
  ELSE $
   return, poly_2d(img,[0, 0, (d1 - 1.0)/(Ncol-1), 0], $
                   [0, (d2 - 1.0)/(Nrow-1), 0, 0], 1, Ncol, Nrow)
END
