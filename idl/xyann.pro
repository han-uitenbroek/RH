; ----------------------------------------- xyann.pro ---------------- ;

; -------- begin -------------------------- xann.pro ----------------- ;

FUNCTION xann, fraction

  xann = !X.CRANGE(0) + fraction*(!X.CRANGE(1) - !X.CRANGE(0))
  IF (!X.TYPE) THEN xann = 10^xann

  return, xann
END
; -------- end ---------------------------- xann.pro ----------------- ;

; -------- begin -------------------------- yann.pro ----------------- ;

FUNCTION yann, fraction, LOGARITHM=logarithm

  yann = !Y.CRANGE(0) + fraction*(!Y.CRANGE(1) - !Y.CRANGE(0))
  IF (!Y.TYPE) THEN yann = 10^yann

  return, yann
END
; -------- end ---------------------------- yann.pro ----------------- ;
