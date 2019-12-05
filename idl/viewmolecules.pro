; ----------------------------------------- viewmolecules.pro -------- ;

; -------- begin -------------------------- XViewMolecule_Event.pro -- ;

PRO XViewMolecule_Event, Event

@geometry.common
@files.common

  COMMON screen_common, screenSize, scaleFactor

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state

  CASE Action OF

    'QUIT': BEGIN
      widget_control, Event.top, /DESTROY
    END

    'TERMDIAGRAM': BEGIN
      XViewMolTerm, state.molecule, GROUP_LEADER=Event.top
    END

    'VIBRATION_POPS': BEGIN
      XViewMolPops, state.molecule, GROUP_LEADER=Event.top
    END

    'ORIENT': orient, state, EVENT_HANDLER='XViewMolecule_Event', $
     GROUP_LEADER=Event.top, TITLE='Orientation Molecular Density'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      displayMolecule, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      displayMolecule, state
    END

    'SURFACE_TRACK': displayMolecule, state, /TRACK

    'PRINT': BEGIN
      filename = '/tmp/viewMolecule-' + timeStamp() + '.ps'
      IF ((geometryType EQ "TWO_D_PLANE" OR $
           geometryType EQ "THREE_D_PLANE") AND state.displayType GT 0) THEN $
       font = -1 ELSE font = 0
      psopen, FILENAME=filename, /COLOR, FONT=font
      scaleFactor = [!D.X_SIZE, !D.Y_SIZE] / float(screenSize)
      displayMolecule, state
      psclose
      ScaleFactor = [1.0, 1.0]
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewMolecule-'
    END

    'LOG_ON': IF (Event.select) THEN $
     displayMolecule, setLog(stash, On) $
    ELSE $
     displayMolecule, setLog(stash, Off)
    
    'PANEL':     IF (Event.select) THEN $
     displayMolecule, setDisplayType(stash, 0)
    'WIREGRID':  IF (Event.select) THEN $
     displayMolecule, setDisplayType(stash, 1)
    'SHADESURF': IF (Event.select) THEN $
     displayMolecule, setDisplayType(stash, 2)

    'INFORMATION': result = dialog_message( /INFORMATION, $
            ["Display of molecular data", $
             "", $
             "Version 1.0, Oct 21, 1997", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"] )
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewMolecule_Event.pr0 -- ;

; -------- begin -------------------------- displayMolecule.pro ------ ;

PRO displayMolecule, state, TRACK=track

@geometry.common

  COMMON screen_common, screenSize, scaleFactor

  KM_TO_M = 1.0E+03

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF
  molpops = *(state.molecule.n_ptr)

  IF (state.log) THEN BEGIN
    nonzeros = where(molpops GT 0.0)
    minpops = MIN(molpops[nonzeros])
    molpops = (molpops > minpops)
  ENDIF

  concentration_label = 'n!D' + state.molecule.ID + '!N [m!U-3!N]'

  CASE (geometryType) OF
  "TWO_D_PLANE": BEGIN
    IF (keyword_set(TRACK)) THEN BEGIN
      tracktitle = 'densities for molecule ' + state.molecule.ID
      IF (state.log) THEN tracktitle = 'log of ' + tracktitle
      surf_track, molpops, geometry.x/KM_TO_M, geometry.z/KM_TO_M, $
       ZLOG=state.log, TITLE=tracktitle
    ENDIF ELSE begin
      CASE (state.displayType) OF
        0: BEGIN
          erase
          IF (state.log) THEN $
           molpops = scaleimg_idl(molpops, screenSize[0]-165, $
                                  screenSize[1]-75) > minpops $
          ELSE $
           molpops = scaleimg_idl(molpops, screenSize[0]-165, $
                                  screenSize[1]-75)
          panel, molpops, $
           geometry.x/KM_TO_M, geometry.z/KM_TO_M, $
           XPOS=60, YPOS=50, XTITLE='x [km]', YTITLE='z [km]', $
           SCALETEXT=concentration_label, ZLOG=state.log, /ORTHOSCOPIC, /ORDER
        END

        1: surface, molpops, geometry.x/KM_TO_M, geometry.z/KM_TO_M, $
         XTITLE='x [km]', YTITLE='z [km]', ZTITLE=concentration_label, $
         ZLOG=state.log, /T3D

        2: shade_surf, molpops, geometry.x/KM_TO_M, geometry.z/KM_TO_M, $
         XTITLE='x [km]', YTITLE='z [km]', ZTITLE=concentration_label, $
         ZLOG=state.log, /T3D
      ENDCASE
    ENDELSE
  END
  "ONE_D_PLANE": BEGIN 
    plot, geometry.height/KM_TO_M, molpops, XMARGIN=[12,3], $
     XTITLE='height [km]', YTITLE=concentration_label, YLOG=state.log
  END
  "SPHERICAL_SYMMETRIC": BEGIN
    plot, geometry.r/KM_TO_M, molpops, XMARGIN=[12,3], $
     XTITLE='radius [km]', YTITLE=concentration_label, YLOG=state.log
  END
  ELSE: message, '3-D viewing not yet supported', /INFORMATIONAL
  ENDCASE

END
; -------- end ---------------------------- displayMolecule.pro ------ ;

; -------- begin -------------------------- moleculeWidgetSetup.pro -- ;

FUNCTION moleculeWidgetSetup, molecule

@geometry.common

  COMMON screen_common, screenSize, scaleFactor

  scaleFactor = [1.0, 1.0]
  screenSize  = fix([550, 350])

  EV = 1.60217733E-19

  state = {baseWidget: 0L, drawWidget: 0L, molecule: molecule, $
           log: 1, displayType: 1, xSlider: 0L, zSlider: 0L}

  state.baseWidget = widget_base(TITLE='XViewMolecule', /COLUMN, $
                                 RESOURCE_NAME='XViewMolecule', MBAR=menuBar)
  moleculeBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  structMenu = widget_button(menuBar, VALUE='Structure', /MENU)
  termButton = widget_button(structMenu, VALUE='Termdiagram', $
                             UVALUE='TERMDIAGRAM')
  IF (NOT ptr_valid(molecule.E_ptr)) THEN $
   widget_control, termButton, SENSITIVE=0

  levelButton = widget_button(structMenu, VALUE='Vibration pops', $
                              UVALUE='VIBRATION_POPS')
  IF (NOT ptr_valid(molecule.nv_ptr)) THEN $
   widget_control, levelButton, SENSITIVE=0

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    orientMenu = widget_button(menuBar, VALUE='Orientation', /MENU)
    orientButton = widget_button(orientMenu, VALUE='set orientation', $
                                  UVALUE='ORIENT')
    trackButton = widget_button(orientMenu, VALUE='track surface', $
                                  UVALUE='SURFACE_TRACK')
  ENDIF

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewMolecule', $
                             UVALUE='INFORMATION')

  quantityFrame = widget_base(moleculeBase, /COLUMN, /FRAME)

  IDFrame  = widget_base(quantityFrame, /ROW)
  IDLabel  = widget_label(IDFrame, VALUE='Molecule ID:')
  IDWidget = widget_label(IDFrame, /FRAME, VALUE=molecule.ID)

  molFrame    = widget_base(quantityFrame, /ROW)
  EdissLabel  = widget_label(molFrame, VALUE='E_diss: ')
  EdissText   = string(molecule.Ediss/EV, FORMAT='(F5.2)')
  EdissWidget = widget_label(molFrame, /FRAME, VALUE=EdissText)
  EdissLabel  = widget_label(molFrame, VALUE=' [eV]')

  TminLabel  = widget_label(molFrame, VALUE='   T_min: ')
  TminText   = string(molecule.Tmin, FORMAT='(F6.0)')
  TminWidget = widget_label(molFrame, /FRAME, VALUE=TminText)
  TminLabel  = widget_label(molFrame, VALUE=' [K]')

  TmaxLabel  = widget_label(molFrame, VALUE='   T_max: ')
  TmaxText   = string(molecule.Tmax, FORMAT='(F6.0)')
  TmaxWidget = widget_label(molFrame, /FRAME, VALUE=TmaxText)
  TmaxLabel  = widget_label(molFrame, VALUE=' [K]')

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame =  widget_base(moleculeBase, /FRAME, /COLUMN)
  state.drawWidget = widget_draw(drawFrame, $
                                 XSIZE=screenSize[0], YSIZE=screenSize[1])

  makeupFrame  = widget_base(moleculeBase, /ROW)
  linlogFrame  = widget_base(makeupFrame, /ROW, /EXCLUSIVE, /FRAME)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                               UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                               UVALUE='LOG_ON')
  widget_control, (state.log) ? logButton : linearButton, /SET_BUTTON
  
  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    shadeFrame  = widget_base(makeupFrame, /ROW, /EXCLUSIVE, /FRAME)
    panelButton = widget_button(shadeFrame, VALUE='Panel', $
                                UVALUE='PANEL')
    wireButton  = widget_button(shadeFrame, VALUE='Wire Grid', $
                                 UVALUE='WIREGRID')
    shadeButton = widget_button(shadeFrame, VALUE='Shade Surface', $
                                  UVALUE='SHADESURF')
    widget_control, wireButton, /SET_BUTTON
  ENDIF

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- moleculeWidgetSetup.pro -- ;

; -------- begin -------------------------- XViewMolecule.pro -------- ;

PRO XViewMolecule, moleculeNo, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWMOLECULE
;
; PURPOSE:
;	This procedure shows molecular populations of molecule No MoleculeNo
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XVIEWMOLECULE, Moleculeno, GROUP_LEADER=group_leader
;
; INPUTS:
;	Moleculeno:  Index number of the molecule in list in ATMOSCOMMON
;
; KEYWORD PARAMETERS:
;	GROUP_LEADER:    ID of parent widget in hierarchy
;
; COMMON BLOCKS:
;	ATMOSCOMMON:    atmos, metals, molecules, nHmin
;
; SIDE EFFECTS:
;	The molecural populations are put in a heap variable. The pointer to
;       this variable is freed when the routine is exited.
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Mon Jan 29 09:07:26 2001 --
;-
@atmos.common

  IF (n_tags(molecules) LE 0) THEN return
  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0
  state = moleculeWidgetSetup(molecules[moleculeNo])

  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  surfr, AX=30, AZ=30
  displayMolecule, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewMolecule', state.baseWidget, $
   EVENT_HANDLER='XViewMolecule_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewMolecule.pro -------- ;
