 ; ----------------------------------------- viewpops.pro ------------- ;

FUNCTION levelMenu_Event_Func, event

  widget_control, event.id, GET_UVALUE=uvalue

  return, {ID:event.handler, TOP:event.top, HANDLER:0L, VALUE:uvalue}
END
; -------- begin -------------------------- setOplot.pro ------------- ;

FUNCTION setOplot, stash, toggle
  widget_control, stash, GET_UVALUE=state

  state.oplot = toggle
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setOplot.pro ------------- ;

; -------- begin -------------------------- setLevel.pro ------------- ;

FUNCTION setLevel, stash, level
  widget_control, stash, GET_UVALUE=state

  IF (NOT state.oplot) THEN BEGIN
    widget_control, state.labelNo,   SET_VALUE=string(FORMAT='(I2)', level)
    widget_control, state.labelText, $
     SET_VALUE=(*state.atom_ptr).labels(level)
  ENDIF

  state.level = level
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setLevel.pro ------------- ;

; -------- begin -------------------------- XViewPops_Event.pro ------ ;

PRO XViewPops_Event, Event

@files.common
@geometry.common

  Off = 0  &  On = 1

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(Event.top, /CHILD)
  widget_control, stash, GET_UVALUE=state
     
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY
   
    'PRINT': BEGIN
      filename = '/tmp/viewPops-' + timeStamp() + '.ps'
      IF ((geometryType EQ "TWO_D_PLANE" OR $
           geometryType EQ "THREE_D_PLANE")) THEN $
       font = -1 ELSE font = 0
      psopen, FILENAME=filename, /COLOR, FONT=font
      drawPops, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewPops-'
    END

    'TERMDIAGRAM': BEGIN
      XViewTermDiag, state.atom_ptr, GROUP_LEADER=Event.top
    END

    'ORIENT': orient, state, EVENT_HANDLER='XViewPops_Event', $
     GROUP_LEADER=Event.top, TITLE='Orientation pops'
  
    'ORIENT_UPDATE': BEGIN
      OrientationUpdate, state
      drawPops, state
    END

    'ORIENT_RESET': BEGIN
      OrientationUpdate, state, /SET
      drawPops, state
    END

    'NEWPOPSFILE': BEGIN
      result = readPops(*atom_ptr, $
                        dialog_pickfile(FILTER='*.out', $
                                        TITLE='Populations File', $
                                        /MUST_EXIST, /READ))
      drawPops, setLevel(stash, 0)
    END
 
    'LEVELMENU': drawPops, setLevel(stash, Event.value)
   
    'LOG_ON': IF (Event.select) THEN $
       drawPops, setLog(stash, On)  ELSE  drawPops, setLog(stash, Off)

    'SHADE_ON': IF (Event.select) THEN $
     drawPops, setShade(stash, On)  ELSE  drawPops, setShade(stash, Off)

    'OPLOT_ON': IF (Event.select) THEN $
     r = setOplot(stash, On)  ELSE BEGIN
      r = setLevel(stash, state.level)
      drawPops, setOplot(stash, Off)
    ENDELSE
  
    'TYPE_0': IF (Event.select) THEN drawPops, setType(stash, 0)
    'TYPE_1': IF (Event.select) THEN drawPops, setType(stash, 1)
    'TYPE_2': IF (Event.select) THEN drawPops, setType(stash, 2)
    'TYPE_3': IF (Event.select) THEN drawPops, setType(stash, 3)
   
    'INFORMATION': result = dialog_message(/INFORMATION, $
                       ["Display of population numbers", $
                        "", $
                        "", $
                        "Version 1.0, Jun 2, 1995", $
                        "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
    ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewpops_Event.pro ------ ;

; -------- begin -------------------------- popWidgetSetup.pro ------- ;

FUNCTION popWidgetSetup, theatom

@geometry.common

  state = {atom_ptr: theatom, baseWidget: 0L, drawWidget: 0L, $
           labelNo: 0L, labelText: 0L, xSlider: 0L, zSlider: 0L, $
           level: 0, log: 1, shade: 0, type: 0, oplot: 0}

  state.baseWidget = widget_base(TITLE='XViewPops', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='XViewPops')

  popsBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu  = widget_button(menuBar, VALUE='File', /MENU)
  openAtmosButton = widget_button(fileMenu, VALUE='new population', $
                                   UVALUE='NEWPOPSFILE')

  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewPops', $
                              UVALUE='INFORMATION')

  levelMenu = widget_button(menuBar, Value='Levels', UVALUE='LEVELMENU', $
                             /MENU, EVENT_FUNC='levelMenu_Event_Func')
  FOR i=0,n_elements((*theatom).labels)-1 DO BEGIN
    menuItem = widget_button(levelMenu, VALUE=(*theatom).labels(i), UVALUE=i)
  ENDFOR

  termMenu   = widget_button(menuBar, VALUE='Termdiagram', /MENU)
  termButton = widget_button(termMenu, VALUE='termdiagram', $
                              UVALUE='TERMDIAGRAM')

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    orientMenu = widget_button(menuBar, VALUE='Orientation', /MENU)
    orientButton = widget_button(orientMenu, VALUE='set orientation', $
                                  UVALUE='ORIENT')
  ENDIF

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')


  ;; --- Draw frame base widget --                      ------------- ;;

  drawFrame = widget_base(popsBase, /FRAME, /COLUMN)

  labelFrame      = widget_base(drawFrame, /ROW)
  nolabel         = widget_label(labelFrame, VALUE='Level No:')
  state.labelNo   = widget_label(labelFrame, /FRAME)
  labelLabel      = widget_label(labelFrame, VALUE='      Label:')
  state.labelText = widget_label(labelFrame, /FRAME, $
                                  VALUE='', /DYNAMIC_RESIZE)

  state.drawWidget = widget_draw(drawFrame, XSIZE=500, YSIZE=400)

  makeupFrame  = widget_base(popsBase, /ROW, /FRAME)
  linlogFrame  = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE)
  linearButton = widget_button(linlogFrame, VALUE='Linear', $
                                UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Logarithmic', $
                                UVALUE='LOG_ON')
  widget_control, logButton, /SET_BUTTON
  
  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    shadeFrame = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE)
    wireButton = widget_button(shadeFrame, VALUE='Wire Grid', $
                                UVALUE='SHADE_OFF')
    shadeButton  = widget_button(shadeFrame, VALUE='Shade Surface', $
                                  UVALUE='SHADE_ON')
    widget_control, wireButton, /SET_BUTTON
  ENDIF ELSE BEGIN
    oplotFrame = widget_base(makeupFrame, /COLUMN, /EXCLUSIVE)
    plotButton = widget_button(oplotFrame, VALUE='Plot', $
                                  UVALUE='OPLOT_OFF')
    oplotButton = widget_button(oplotFrame, VALUE='Oplot', $
                                UVALUE='OPLOT_ON')
    widget_control, plotButton, /SET_BUTTON
  ENDELSE

  typeFrame   = widget_base(makeupFrame, COLUMN=2, /EXCLUSIVE)
  typeButton0 = widget_button(typeFrame, VALUE='n',    UVALUE='TYPE_0')
  typeButton1 = widget_button(typeFrame, VALUE='n*',   UVALUE='TYPE_1')
  typeButton2 = widget_button(typeFrame, VALUE='n/n*', UVALUE='TYPE_2')
  typeButton3 = widget_button(typeFrame, VALUE='n/n0', UVALUE='TYPE_3')
  widget_control, typeButton0, /SET_BUTTON

  ;; --- Control frame base widget --                -------------- ;;

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- popWidgetSetup.pro ------- ;

; -------- begin -------------------------- drawPops.pro ------------- ;

PRO drawPops, state

@geometry.common
@atmos.common

  COMMON oplot_save_thePops, thePops

  Off = 0  &  On = 1

  theatom = state.atom_ptr
  n       = (*theatom).n_ptr
  nstar   = (*theatom).nstar_ptr

  IF ((state.level LT 0)  OR  (state.level GE (*theatom).Nlevel)) THEN return

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
    thick = 1
  ENDIF ELSE $
   thick = 3

  IF (state.type LE 1) THEN $
   poptitle = 'Population number m!U-3!N' $
  ELSE $
   poptitle = 'Relative populations'

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    CASE (state.type) OF
      0: nv = (*n)[*, *, state.level]
      1: nv = (*nstar)[*, *, state.level]
      2: nv = (*n)[*, *, state.level] / (*nstar)[*, *, state.level]
      3: nv = (*n)[*, *, state.level] / (*n)[*, *, 0]
    ENDCASE
    IF (state.log) THEN  zlog = 1  ELSE  zlog = 0

    IF (state.shade) THEN BEGIN
      shade_surf, nv, geometry.x/1.0E3, geometry.z/1.0E3, /T3D, $
       CHARSIZE=1.4, XTITLE='x', YTITLE='z', ZTITLE=poptitle, ZLOG=zlog
    ENDIF ELSE BEGIN
      surface, nv, geometry.x/1.0E3, geometry.z/1.0E3, /T3D, CHARSIZE=1.4, $
       XTITLE='x', YTITLE='z', ZTITLE=poptitle, BOTTOM=200B, ZLOG=zlog
    ENDELSE
  ENDIF ELSE BEGIN
    IF (state.log) THEN  ylog = 1  ELSE  ylog = 0
    IF (state.oplot NE On) THEN $
     thePops = [state.level] $
    ELSE BEGIN 
      idx = where(thePops EQ state.level, count)
      IF (count EQ 0) THEN thePops = [thePops, state.level]
    ENDELSE

    FOR i=0, n_elements(thePops)-1 DO BEGIN
      CASE (state.type) OF
        0: nv = (*n)[*, thePops[i]]
        1: nv = (*nstar)[*, thePops[i]]
        2: nv = (*n)[*, thePops[i]] / (*nstar)[*, thePops[i]]
        3: nv = (*n)[*, thePops[i]] / (*n)[*, 0]
      ENDCASE
      IF (i EQ 0) THEN $      
       plot, geometry.cmass, nv, XTITLE='Column Mass [kg m!U-2!N]', $
       YTITLE=poptitle, YLOG=ylog, /XLOG, THICK=thick $
      ELSE BEGIN
        oplot, geometry.cmass, nv, COLOR=200B - i*32B, THICK=thick
        rhannotate, xann(0.7), yann(i*0.05), $
         TEXT=(*theatom).labels[thePops[i]], $
         MARKTYPE=0, MARKCOLOR=(200B - i*32B), THICK=thick
      ENDELSE
    ENDFOR
  ENDELSE
END
; -------- end ---------------------------- DrawPops.pro ------------- ;

; -------- begin -------------------------- XViewPops.pro ------------ ;

PRO XViewPops, atom_ptr, LEVEL=level, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWPOPS
;
; PURPOSE:
;	This routine displays the poulation densities of the atomic model
;       pointed to by Atom_ptr
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XVIEWPOPS, Atom_ptr, LEVEL=level, GROUP_LEADER=group_leader
;
; INPUTS:
;	Atom_ptr:  Pointer to the heap variable containing atomic data
;                  structure
;
; KEYWORD PARAMETERS:
;	LEVEL:  Index of initial atomic level of which populations are to
;               be displayed.
;
;	GROUP_LEADER:  ID of parent widget in hierarchy.
;
; COMMON BLOCKS:
;	atmos.common
;       files.common
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Mon Mar 12 16:26:47 2001 --
;-


@atmos.common
@files.common



  IF (NOT keyword_set(GROUP_LEADER)) THEN BEGIN
    group_leader=0
    surfr, AX=30, AZ=30
  ENDIF

  IF (n_elements(LEVEL) EQ 0) THEN level = 0

  IF (NOT ptr_valid(atom_ptr)) THEN return
  IF (n_elements(atmos) EQ 0) THEN return

  state = popWidgetSetup(atom_ptr)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  drawPops, setLevel(widget_info(state.baseWidget, /CHILD), level)

  ;; --- Register with the XManager --                 --------------- ;

  xmanager, 'XViewPops', state.baseWidget, $
   EVENT_HANDLER='XViewPops_Event', GROUP_LEADER=group_leader

END
; -------- end ---------------------------- XViewPops.pro ------------ ;
