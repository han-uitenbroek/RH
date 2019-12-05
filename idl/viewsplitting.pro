; ----------------------------------------- viewsplitting.pro -------- ;

function parselabel, g, label, S, L, J, QUIET=quiet

  words = strsplit(label, ' ', /EXTRACT)
  words = words(where(words ne ''))

  resolved = 0
  FOR n=n_elements(words)-1, 0, -1 DO BEGIN
    parity_pos = max([strpos(words[n], 'E'), strpos(words[n], 'O')])

    IF (parity_pos GT 1) THEN BEGIN
      multiplicity = 0
      orbit = ''
      reads, strmid(words[n], parity_pos-2, 1), multiplicity, FORMAT='(I)'
      orbit = strmid(words[n], parity_pos-1, 1)
      resolved = 1
    ENDIF
    IF (resolved) THEN break;
  ENDFOR

  IF (NOT resolved) THEN BEGIN
    IF (NOT keyword_set(QUIET)) THEN $
     print, 'PARSELABEL: Could not find parity of label ', label
    return, 0
  ENDIF

  default_orbits = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
  L = where(orbit EQ default_orbits, count)
  IF (count LE 0) THEN BEGIN
    IF (NOT keyword_set(QUIET)) THEN $
     print, 'PARSELABEL: Could not determine orbit from label ', label
    return, 0
  ENDIF ELSE L = L[0]

  S = (multiplicity - 1) / 2.0;
  J = (g - 1.0) / 2.0
stop
  IF (J GT L + S) THEN BEGIN
    IF (NOT keyword_set(QUIET)) THEN $
     print, 'PARSELABEL: Label ', label, ' is a composite level'
    return, 0
  ENDIF

  return, 1
END
      
FUNCTION zeemanstrength, Ju, Mu, Jl, Ml

  q  = long(Ml - Mu);
  dJ = long(Ju - Jl);

  CASE (dJ) OF
    0: BEGIN
      CASE (q) OF
         0: s = 2.0 * Mu^2
        -1: s = (Ju + Mu) * (Ju - Mu + 1.0)
         1: s = (Ju - Mu) * (Ju + Mu + 1.0)
      ENDCASE
    END

    1: BEGIN 
      CASE (q) OF
         0: s = 2.0 * (Ju^2 - Mu^2)
        -1: s = (Ju + Mu) * (Ju + Mu - 1.0)
         1: s = (Ju - Mu) * (Ju - Mu - 1.0)
      ENDCASE
    END 

    -1: BEGIN
      CASE (q) OF
         0: s = 2.0 * (Ju + 1.0)^2 - Mu^2
        -1: s = (Ju - Mu + 1.0) * (Ju - Mu + 2.0)
         1: s = (Ju + Mu + 1.0) * (Ju + Mu + 2.0)
      ENDCASE
    END
  ENDCASE

  return, s
END

FUNCTION gLande, S, L, J

  IF (J EQ 0.0) THEN $
   gL = 0.0 $
  ELSE $
   gL = 1.5 + (S*(S + 1.0) - L*(L + 1)) / (2.0*J*(J + 1.0))

  return, gL
END

FUNCTION zeemanshift, Su, Lu, Ju, Mu, Sl, Ll, Jl, Ml

  gu = gLande(Su, Lu, Ju)
  gl = gLande(Sl, Ll, Jl)

  return, gl*Ml - gu*Mu
END

FUNCTION get_splittings, Su, Lu, Ju, Sl, Ll, Jl

  components = replicate({q: 0L,  shift: 0.0,  strength: 0.0}, $
                         (2*Ju+1) * (2*Jl+1))

  Ncomponent = 0
  FOR Ml=-Jl, Jl, 1 DO BEGIN
    FOR Mu=-Ju, Ju, 1 DO BEGIN
      IF (abs(Mu - Ml) LE 1) THEN BEGIN
        components[Ncomponent].q = Ml - Mu
        components[Ncomponent].shift = $
         zeemanshift(Su, Lu, Ju, Mu, Sl, Ll, Jl, Ml)
        components[Ncomponent].strength = zeemanstrength(Ju, Mu, Jl, Ml)
        Ncomponent = Ncomponent + 1
      ENDIF
    ENDFOR
  ENDFOR

  norm = fltarr(3)
  FOR n=0, Ncomponent-1 DO $
    norm[components[n].q + 1] = $
     norm[components[n].q + 1] + components[n].strength

  FOR n=0, Ncomponent-1 DO BEGIN
    IF (norm[components[n].q + 1] GT 0.0) THEN $
     components[n].strength = $
     components[n].strength / norm[components[n].q + 1]
  ENDFOR

  return, components[0:Ncomponent-1]
END

FUNCTION effectiveLande, g_u, J_u, g_l, J_l

  return,  0.5*(g_u + g_l) + 0.25*(g_u - g_l) * $
   (J_u*(J_u + 1.0) - J_l*(J_l + 1.0))
END

; -------- begin -------------------------- XViewSplitting_Event.pro - ;

PRO XViewSplitting_Event, Event

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = '/tmp/viewSplitting-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR, /LANDSCAPE
      displaySplitting, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewSplitting-'
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of Zeeman line splitting components", $
             "", $
             "Version 1.0, Feb 4,2000 ", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewSplitting_Event.pro - ;

; -------- begin -------------------------- displaySplitting.pro ----- ;

PRO displaySplitting, state

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  colors = [16B, 100B, 225B]
;;  colors = [16B, 200B, 255B]
  thick = 2

  orbits = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
  multiplicity_l = strtrim(string(FORMAT='(I2)', 2*state.Sl + 1), 2)
  multiplicity_u = strtrim(string(FORMAT='(I2)', 2*state.Su + 1), 2)
  IF (state.Ju NE fix(state.Ju)) THEN BEGIN
    Jstring_u = strtrim(string(FORMAT='(I2, "/2")', 2*state.Ju), 2)
    Jstring_l = strtrim(string(FORMAT='(I2, "/2")', 2*state.Jl), 2)
  ENDIF ELSE BEGIN
    Jstring_u = strtrim(string(FORMAT='(I2)', state.Ju), 2)
    Jstring_l = strtrim(string(FORMAT='(I2)', state.Jl), 2)
  ENDELSE
  title = '!U'+multiplicity_u+'!N'+orbits[state.Lu]+'!D'+Jstring_u+'!N - ' + $
   '!U'+multiplicity_l+'!N'+orbits[state.Ll]+'!D'+Jstring_l+'!N'

  c = get_splittings(state.Su, state.Lu, state.Ju, $
                     state.Sl, state.Ll, state.Jl)
  c_pi = c[where(c.q EQ 0)]
  c_sp = c[where(c.q EQ 1)]
  c_sm = c[where(c.q EQ -1)]

  s_max = max(c_sp.strength)
  s_min = max(c_pi.strength)
  bottom = -s_min - (s_min + s_max) / 2.0

  sh_max = max(c.shift)

  plot, /NODATA, c.shift, c.strength, $
   XRANGE=[-sh_max, sh_max]*1.5, XSTYLE=1, $
   YRANGE=[bottom, s_max*1.2], YSTYLE=1, $
   XTITLE='Wavelength shift [Larmor units]', YTITLE='Normalized strength'
  oplot, !X.CRANGE, [0.0, 0.0], /LINE

  FOR n=0, n_elements(c_pi)-1 DO $
   oplot, [c_pi[n].shift, c_pi[n].shift], [0.0, c_pi[n].strength] + bottom, $
   COLOR=colors[0], THICK=thick

   FOR n=0, n_elements(c_sp)-1 DO $
   oplot, [c_sp[n].shift, c_sp[n].shift], [0.0, c_sp[n].strength], $
   COLOR=colors[1], THICK=thick
   
  FOR n=0, n_elements(c_sm)-1 DO $
   oplot, [c_sm[n].shift, c_sm[n].shift], [0.0, -c_sm[n].strength], $
   COLOR=colors[2], THICK=thick

  arrow, state.gL_eff, -0.3*s_min, state.gL_eff, 0.0, $
   /DATA, /SOLID, COLOR=colors[1]
  rhannotate, TEXT=string(FORMAT='("gL!Deff!N = ", F6.3)', state.gL_eff), $
   state.gL_eff+0.05, -0.3*s_min, CHARSIZE=0.8

  rhannotate, xann(0.5), yann(0.9), ALIGNMENT=0.5, TEXT=title, CHARSIZE=1.2

  rhannotate, xann(0.05), yann(0.15), TEXT='!7p!x - component', MARKTYPE=0, $
   MARKCOLOR=colors[0], THICK=thick
  rhannotate, xann(0.05), yann(0.10), TEXT='!7r!U+!N!x', MARKTYPE=0, $
   MARKCOLOR=colors[1], THICK=thick
  rhannotate, xann(0.05), yann(0.05), TEXT='!7r!U-!N!x', MARKTYPE=0, $
   MARKCOLOR=colors[2], THICK=thick

END
; -------- end ---------------------------- displaySplitting.pro ----- ;

; -------- begin -------------------------- splitWidgetSetup.pro ----- ;

FUNCTION splitWidgetSetup, g_j, label_j, g_i, label_i

  IF (NOT parselabel(g_j, label_j, Su, Lu, Ju) OR $
      NOT parselabel(g_i, label_i, Sl, Ll, Jl)) THEN return, 0

  gL_eff = effectiveLande(gLande(Su, Lu, Ju), Ju, gLande(Sl, Ll, Jl), Jl)

  state = {baseWidget: 0L, drawWidget: 0L, gL_eff: gL_eff, $
           Su: Su,  Lu: Lu,  Ju: Ju,  Sl: Sl,  Ll: Ll,  Jl: Jl}

  state.baseWidget = widget_base(TITLE='XViewSplitting', /COLUMN, $
                                 RESOURCE_NAME='XViewSplitting', $
                                 MBAR=menuBar)

  Base = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')
                             
  ;; --- Print and help menus --

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewSplitting', $
                             UVALUE='INFORMATION')

  infoFrame =  widget_base(Base, /FRAME, /ROW)
  frame1 = widget_base(infoFrame, /COLUMN)
  label  = widget_label(frame1, VALUE='Transition: ', /ALIGN_RIGHT)
  label  = widget_label(frame1, VALUE='Effective Lande facor: ', /ALIGN_RIGHT)

  frame2 = widget_base(infoFrame, /COLUMN)
  label  = widget_label(frame2, $
                        VALUE=label_j + ' --> ' + label_i, /ALIGN_LEFT)
  label  = widget_label(frame2, $
                        VALUE=string(FORMAT='(F5.3)', gL_eff), /ALIGN_LEFT)

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame =  widget_base(Base, /FRAME, /COLUMN)
  state.drawWidget = widget_draw(drawFrame, XSIZE=500, YSIZE=350)

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- splitWidgetSetup.pro ----- ;

; -------- begin -------------------------- XViewSplitting.pro ------- ;

PRO XViewSplitting, GROUP_LEADER=group_leader, g_j, label_j, g_i, label_i

;+
; NAME:
;	XVIEWSPLITTING
;
; PURPOSE:
;
; CATEGORY:
;	Data reduction
;
; CALLING SEQUENCE:
;
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Thu Jul  9 09:55:38 2009 --
;-

@atmos.common
@files.common

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  state = splitWidgetSetup(g_j, label_j, g_i, label_i)
  IF (n_tags(state) EQ 0) THEN return

  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  displaySplitting, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewSplitting', state.baseWidget, $
   EVENT_HANDLER='XViewSplitting_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewSplitting.pro ------- ;
