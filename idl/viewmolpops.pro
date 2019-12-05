; ----------------------------------------- viewmolpops.pro ---------- ;

; -------- begin -------------------------- XViewMolPops_Event.pro --- ;

PRO XViewMolPops_Event, Event

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': BEGIN
      widget_control, Event.top, /DESTROY
    END

    'PRINT': BEGIN
      filename = '/tmp/viewMolPops-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR
      displayMolPops, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewSource-'
    END

    'LOG_ON': displayMolPops, setLog(stash, Event.select)

    'RATIO_ON': displayMolPops, setRatio(stash, Event.select)

    'INFORMATION': result = dialog_message( /INFORMATION, $
            ["Display of molecular energy level data", $
             "", $
             "Version 1.0, Feb 24 1999", $
             "Han Uitenbroek (huitenbroek@cfa.harvard.edu)"] )
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewMolPops_Event.pro --- ;

; -------- begin -------------------------- displayMolecule.pro ------ ;

PRO displayMolPops, state

@geometry.common

  KM_TO_M = 1.0E+03

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  ytitle = (state.ratio) ? 'nv / nv!U*!N' : 'nv / n!Dmolecule!N'

  non_zero = where(*(state.molecule.n_ptr) GT 0.0)
  IF (state.ratio) THEN $
   beta = (*state.molecule.nv_ptr)[non_zero, *] / $
   (*state.molecule.nvstar_ptr)[non_zero, *] $
  ELSE BEGIN
    beta = (*state.molecule.nv_ptr)[non_zero, *]
    FOR v=0, state.molecule.Nv-1 DO $
     beta[*, v] = beta[*, v] / (*state.molecule.n_ptr)[non_zero]
  ENDELSE
  beta_max = max(beta, MIN=beta_min)

  IF (!D.NAME EQ 'PS') THEN thick = 4  ELSE  thick = 2

  plot, geometry.height[non_zero]/KM_TO_M, beta[*, 0], YLOG=state.log, $
   XTITLE='Height [km]', YTITLE=ytitle, THICK=thick, $
   YRANGE=[beta_min, beta_max]
  rhannotate, xann(0.02), yann(0.9), THICK=thick, $
   TEXT=string(FORMAT='("v = ", I1)', 0), MARKTYPE=0, MARKCOLOR=color

  FOR v=1, state.molecule.Nv-1 DO BEGIN
    color = 255B - v*16B
    format = (v LT 10) ? '(I1)' : '(I2)'
    oplot, geometry.height[non_zero] / KM_TO_M, beta[*, v], $
     COLOR=color, THICK=thick
     rhannotate, xann(0.02), yann(0.9 - v*0.05), THICK=thick, $
     TEXT=string(FORMAT=format, v), MARKTYPE=0, MARKCOLOR=color
  ENDFOR

END
; -------- end ---------------------------- displayMolecule.pro ------ ;

; -------- begin -------------------------- moleculeWidgetSetup.pro -- ;

FUNCTION molPopsWidgetSetup, molecule

  state = {baseWidget: 0L, drawWidget: 0L, molecule: molecule, $
           log: 0, ratio: 1}

  state.baseWidget = widget_base(TITLE='XViewMolPops', /COLUMN, $
                                 RESOURCE_NAME='XViewMolPops', MBAR=menuBar)
  base = widget_base(state.baseWidget, /COLUMN)

  fileMenu = widget_button(menuBar, VALUE='File', /MENU)
  button = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                         RESOURCE_NAME='quitbutton')

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame = widget_base(base, /FRAME, /COLUMN)
  state.drawWidget = widget_draw(drawFrame, XSIZE=500, YSIZE=350)

  makeupFrame = widget_base(base, /ROW, /FRAME)
  linlogFrame  = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                               UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                               UVALUE='LOG_ON')
  widget_control, (state.log) ? logButton : linearButton, /SET_BUTTON

  ratioFrame = widget_base(makeupFrame, /ROW, /EXCLUSIVE)
  offButton = widget_button(ratioFrame, VALUE='nv', $
                            UVALUE='RATIO_OFF')
  onButton = widget_button(ratioFrame, VALUE='nv/nv*', $
                           UVALUE='RATIO_ON')
  widget_control, (state.ratio) ? onButton : offButton, /SET_BUTTON

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- moleculeWidgetSetup.pro -- ;

; -------- begin -------------------------- XViewMolecule.pro -------- ;

PRO XViewMolPops, molecule, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWMOLTERM
;
; PURPOSE:
;	This procedure shows molecular energy levels of molecule molecule
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XVIEWMOLPOPS, molecule, GROUP_LEADER=group_leader
;
; INPUTS:
;	molecule:  Molecular structure
;
; KEYWORD PARAMETERS:
;	GROUP_LEADER:    ID of parent widget in hierarchy
;
; COMMON BLOCKS:
;	GEOMETRYCOMMON:  geometry, geometryType
;
; SIDE EFFECTS:
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Jun 13 17:38:33 2001 --
;-
  IF (n_tags(molecule) LE 0) THEN return
  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0
  state = molPopsWidgetSetup(molecule)

  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayMolPops, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewMolPops', state.baseWidget, $
   EVENT_HANDLER='XViewMolPops_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewMolPops.pro --------- ;
