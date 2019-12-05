; ----------------------------------------- viewmolterm.pro ---------- ;

; -------- begin -------------------------- XViewMolTerm_Event.pro --- ;

PRO XViewMolTerm_Event, Event

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': BEGIN
      widget_control, Event.top, /DESTROY
    END

    'PRINT': BEGIN
      filename = '/tmp/viewMolTerm-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR
      displayMolTerm, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewMolTerm-'
    END

    'INFORMATION': result = dialog_message( /INFORMATION, $
            ["Display of molecular energy level data", $
             "", $
             "Version 1.0, Feb 24 1999", $
             "Han Uitenbroek (huitenbroek@cfa.harvard.edu)"] )
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewMolTerm_Event.pro --- ;

; -------- begin -------------------------- displayMolecule.pro ------ ;

PRO displayMolTerm, state

  EV = 1.60217733E-19

  IF (ptr_valid(state.molecule.E_ptr)) THEN $
   E = *(state.molecule.E_ptr) / EV $
  ELSE $
   return

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  E_max = max(E, MIN=E_min)
  plot, [0, state.molecule.Nv], [E_min, E_max], $
   /NODATA, XMINOR=1, XTICKS=state.molecule.Nv + 1, $
   XTICKV=indgen(state.molecule.Nv), $
   XRANGE=[-0.5, state.molecule.Nv - 0.5], YRANGE=[-0.2, E_max + 0.3], $
   XSTYLE=9, YSTYLE=9, $
   XTITLE='Vibrational level', YTITLE='Energy [eV]'

   usersym, [-2, 2], [0, 0], THICK=1, COLOR=200B
   FOR v=0, state.molecule.Nv-1 DO $
    oplot, fltarr(state.molecule.NJ) + v, E[*, v], PSYM=8
    
END
; -------- end ---------------------------- displayMolecule.pro ------ ;

; -------- begin -------------------------- moleculeWidgetSetup.pro -- ;

FUNCTION molTermWidgetSetup, molecule

  state = {baseWidget: 0L, drawWidget: 0L, molecule: molecule}

  state.baseWidget = widget_base(TITLE='XViewMolTerm', /COLUMN, $
                                 RESOURCE_NAME='XViewMolTerm', MBAR=menuBar)
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

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- moleculeWidgetSetup.pro -- ;

; -------- begin -------------------------- XViewMolecule.pro -------- ;

PRO XViewMolTerm, molecule, GROUP_LEADER=group_leader

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
;	XVIEWMOLTERM, molecule, GROUP_LEADER=group_leader
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
;   --- Last modified: Thu Aug 22 13:33:38 2002 --
;-
  IF (n_tags(molecule) LE 0) THEN return
  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0
  state = molTermWidgetSetup(molecule)

  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayMolTerm, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewMolTerm', state.baseWidget, $
   EVENT_HANDLER='XViewMolTerm_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewMolTerm.pro --------- ;
