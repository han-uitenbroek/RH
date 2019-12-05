PRO plot_pf, atom, T, PF=pf

;+
; NAME:
;	PLOT_PF
;
; PURPOSE:
;	This procedure plots the partition function of the different
;       ionization stages of atom as a funcion of temperature T.
;
; CATEGORY:
;	Data analysis
;
; CALLING SEQUENCE:
;	PLOT_pf, Atom, T
;
; INPUTS:
;	Atom   -- Structure containing atomic model.
;       T      -- Array of temperatures
;
; OUTPUTS:
;	None
;
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Thu Sep  3 17:31:28 1998 --
;-

  KBOLTZMANN = 1.380658E-23

  s_max = max(atom.stage, MIN=s_min)
  pf = fltarr(n_elements(T), s_max - s_min + 1)

  FOR s=s_min, s_max DO BEGIN
    index = where(atom.stage EQ s, count)
    IF (count GT 1) THEN BEGIN 
      E = atom.E[index] - atom.E[index[0]]
      FOR k=0, n_elements(T)-1 DO BEGIN 
        pf[k, s-s_min] = $
         total(exp(-E/(KBOLTZMANN*T[k])) * atom.g[index])
      ENDFOR
    ENDIF ELSE BEGIN 
      pf[*, s-s_min] = atom.g[index[0]]
    ENDELSE
  ENDFOR

  pf_max = max(pf, MIN=pf_min)
  plot, /NODATA, T, pf[*, 0], YRANGE=[0.95*pf_min, 1.05*pf_max], YSTYLE=1
  FOR s=s_min, s_max DO BEGIN
    color = 192-32*(s-s_min)
    oplot, T, pf[*, s-s_min], COLOR=color
    rhannotate, TEXT="stage " + string(FORMAT='(I2)', s), $
     xann(0.05), yann(0.2-0.025*(s-s_min)), MARKCOLOR=color, MARKTYPE=0
  ENDFOR
END
