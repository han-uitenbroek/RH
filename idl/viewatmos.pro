; ----------------------------------------- viewatmos.pro ------------ ;

; -------- begin -------------------------- setDisplayType.pro ------- ;

FUNCTION setDisplayType, stash, type

  widget_control, stash, GET_UVALUE=state

  state.displayType = type
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setDisplayType.pro ------- ;

; -------- begin -------------------------- XViewAtmos_Event.pro ----- ;

PRO XViewAtmos_Event, Event

@files.common

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'GRID':     XViewGrid, GROUP_LEADER=Event.top
    'GEOMETRY': XViewAngles, GROUP_LEADER=Event.top

    'NEWATMOSFILE': BEGIN
      IF (result = readAtmos(dialog_pickfile(FILTER='*.out', $
                                             TITLE='Atmos File', $
                                             /MUST_EXIST, /READ, $
                                             FILE=atmosFile))) THEN $
       displayAtmos, state $
      ELSE $
       widget_control, Event.top, /DESTROY
    END
    'PANEL':     IF (Event.select) THEN displayAtmos, setDisplayType(stash, 0)
    'WIREGRID':  IF (Event.select) THEN displayAtmos, setDisplayType(stash, 1)
    'SHADESURF': IF (Event.select) THEN displayAtmos, setDisplayType(stash, 2)

    'INFORMATION': result = dialog_message(/INFORMATION, $
          ["Display of atmospheric data", $
           "", $
           "Temperature, electron density, total hydrogen density,", $
           "and microturbulent velocities can be view as image panels,", $
           "and wire grid or shaded surface plots.", $
           "", $
           "Version 1.0, Jun 2, 1995", $
           "Han Uitenbroek (huitenbroek@cfa.harvard.edu)"])
    ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewAtm_Event.pro ------- ;

; -------- begin -------------------------- displayAtmos.pro --------- ;

PRO displayAtmos, state

@geometry.common
  COMMON screen_common, ScreenSize, ScaleFactor
@atmos.common

  widget_control, state.IDtext, SET_VALUE=atmos.ID

  CASE geometryType OF 
    "ONE_D_PLANE": BEGIN
      widget_control, state.NzText, $
       SET_VALUE=string(FORMAT='(I3)', geometry.Ndep)
      xtitle = 'height [km]'
    END
    "TWO_D_PLANE": BEGIN
      widget_control, state.NxText, $
       SET_VALUE=string(FORMAT='(I3)', geometry.Nx)
      widget_control, state.NzText, $
       SET_VALUE=string(FORMAT='(I3)', geometry.Nz)
      xtitle = 'x [km]'  &  ytitle = 'z [km]'
    END 
    "SPHERICAL_SYMMETRIC": BEGIN 
      widget_control, state.NzText, $
       SET_VALUE=string(FORMAT='(I3)', geometry.Nradius)
      xtitle = 'r [km]'
    END 
  ENDCASE

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    erase
    xkm = geometry.x / 1.0E+3
    zkm = geometry.z / 1.0E+3
    CASE (state.displayType) OF
      0: BEGIN
        ScreenSize = state.ScreenSize  &  ScaleFactor = state.ScaleFactor

        panel, scaleimg_idl(atmos.T, 200, 175), XPOS=75, YPOS=50, $
         xkm, zkm, XTITLE=xtitle, YTITLE=ytitle, $
         SCALETEXT='T [K]', /ORDER, /ORTHOSCOPIC
        
        panel, scaleimg_idl(alog10(atmos.n_elec), 200, 175), $
         XPOS=400, YPOS=50, xkm, zkm, XTITLE=xtitle, $
         SCALETEXT='log Ne [m!U-3!N]', /ORDER, /ORTHOSCOPIC
        
        panel, scaleimg_idl(alog10(total(atmos.nH, 3)), 200, 175), $
         XPOS=400, YPOS=275, xkm, zkm, XTITLE=xtitle, $
         SCALETEXT='log NHydr [m!U-3!N]', /ORDER, /ORTHOSCOPIC
        
        panel, scaleimg_idl(atmos.vturb/1.0E3, 200, 175), XPOS=75, YPOS=275, $
         xkm, zkm, XTITLE=xtitle, YTITLE=ytitle, $
         SCALETEXT='Vturb [km/s]', /ORDER, /ORTHOSCOPIC
      END

      1: BEGIN
        !P.NOERASE = 1

        t3d, SCALE=[0.57, 0.57, 0.57, 0.57]
        surface, atmos.T, xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, $
         /T3D, BOTTOM=200B, ZTITLE='T [K]', CHARSIZE=3.0

        t3d, TRANSLATE=[0.47, 0.0, 0.0]
        surface, alog10(atmos.n_elec), xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, CHARSIZE=3.0, $
         /T3D, BOTTOM=200B, ZTITLE='log Ne [m!U-3!N]'

        t3d, TRANSLATE=[0.0, 0.47, 0.0]
        surface, alog10(total(atmos.nH, 3)), xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, CHARSIZE=3.0, $
         /T3D, BOTTOM=200B, ZTITLE='log NHydr [m!U-3!N]'

        t3d, TRANSLATE=[-0.47, 0.0, 0.0]
        surface, atmos.vturb/1.0E3, xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, $
         /T3D, BOTTOM=200B, ZTITLE='Vturb [km/s]', CHARSIZE=3.0

        surfr, AX=30, AZ=30  &  !P.NOERASE = 0
      END
      2: BEGIN
        !P.NOERASE = 1

        t3d, SCALE=[0.57, 0.57, 0.57, 0.57]
        shade_surf, atmos.T, xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, $
         /T3D, ZTITLE='T [K]', CHARSIZE=3.0

        t3d, TRANSLATE=[0.47, 0.0, 0.0]
        shade_surf, alog10(atmos.n_elec), xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, CHARSIZE=3.0, $
         /T3D, ZTITLE='log Ne [m!U-3!N]'

        t3d, TRANSLATE=[0.0, 0.47, 0.0]
        shade_surf, alog10(total(atmos.nH, 3)), xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, CHARSIZE=3.0, $
         /T3D, ZTITLE='log NHydr [m!U-3!N]'

        t3d, TRANSLATE=[-0.47, 0.0, 0.0]
        shade_surf, atmos.vturb/1.0E3, xkm, zkm, $
         XTITLE=XTITLE, YTITLE=ytitle, $
         /T3D, ZTITLE='Vturb [km/s]', CHARSIZE=3.0

        surfr, AX=30, AZ=30 &  !P.NOERASE = 0
      END
    ENDCASE
  ENDIF ELSE BEGIN
    symbolColor = 208B
    !P.MULTI = [0, 2, 3]
    IF (geometryType EQ "ONE_D_PLANE") THEN BEGIN 
      x = geometry.height / 1.0E3
      xtitle = 'Height [km]'
      vel_km = geometry.vz / 1.0E3
    ENDIF ELSE BEGIN 
      x = geometry.r / 1.0E3
      xtitle = 'Radius [km]'
      vel_km = geometry.vr / 1.0E3
    ENDELSE
    plot, x, geometry.cmass, XTITLE=xtitle, $
     YTITLE='Column Mass [kg m!U-2!N]', /YLOG, CHARSIZE=1.4
    oplot, x, geometry.cmass, PSYM=5, COLOR=symbolColor
    plot, x, geometry.tau500, XTITLE=xtitle, $
     YTITLE='tau_500', /YLOG, CHARSIZE=1.4
    plot, x, atmos.T, XTITLE=xtitle, YTITLE='T [K]', CHARSIZE=1.4
    plot, x, atmos.n_elec, XTITLE=xtitle, YTITLE='Ne [m!U-3!N]', $
     /YLOG, CHARSIZE=1.4
    plot, x, total(atmos.nH, 2), XTITLE=xtitle, YTITLE='NHydr [m!U-3!N]', $
     /YLOG, CHARSIZE=1.4

    v_max = max([atmos.vturb/1.0E3, vel_km], MIN=v_min)
    plot, x, atmos.vturb/1.0E3, XTITLE=xtitle, YTITLE='Vturb [km/s]', $
     CHARSIZE=1.4, YSTYLE=11, YRANGE=1.05*[v_min, v_max]
    oplot, x, vel_km, COLOR=208B
    axis, /YAXIS, YTITLE='Velocity [km/s]', CHARSIZE=1.4
    
    !P.MULTI = 0
  ENDELSE  
END
; -------- end ---------------------------- displayAtmos.pro --------- ;

; -------- begin -------------------------- atmosWidgetSetup.pro ----- ;

FUNCTION atmosWidgetSetup, atmosFile

@geometry.common
@atmos.common

  IF (geometryType EQ "TWO_D_PLANE") THEN $
   ScreenSize = fix([700, 500]) $
  ELSE $
   ScreenSize = fix([700, 700])

  state = {baseWidget: 0L, drawWidget: 0L, $
           ScreenSize: ScreenSize, ScaleFactor: float([1.0, 1.0]), $
           IDtext: 0L, NxText: 0L,  NzText: 0L, displayType: 0}

  state.baseWidget = widget_base(TITLE='XViewAtmos', /COLUMN, MBAR=menuBar, $
                                 RESOURCE_NAME='XViewAtmos')
  atmosBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu  = widget_button(menuBar, VALUE='File', /MENU)
  openAtmosButton = widget_button(fileMenu, VALUE='open Atmos file', $
                                  UVALUE='NEWATMOSFILE')
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewAtmos', $
                             UVALUE='INFORMATION')

  gridMenu = widget_button(menuBar, VALUE='Grids', /MENU)
  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN 
    gridButton = widget_button(gridMenu, VALUE='spatial grid', $
                               UVALUE='GRID')
  ENDIF
  geomButton = widget_button(gridMenu, VALUE='angular grid', $
                             UVALUE='GEOMETRY')

  ;; --- Draw frame base widget --                      ------------- ;;

  drawFrame = widget_base(atmosBase, /FRAME, /COLUMN)
  labelFrame  = widget_base(drawFrame, /ROW)

  IDlabel = widget_label(labelFrame, VALUE='atmos_ID:')
  state.IDtext = widget_label(labelFrame, /FRAME, /DYNAMIC_RESIZE)

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      NdepLabel    = widget_label(labelFrame, $
                                  VALUE='     Spatial dimension:  Ndep')
      state.Nztext = widget_label(labelFrame, /FRAME, /DYNAMIC_RESIZE)
    END
    "TWO_D_PLANE": BEGIN
      NxLabel      = widget_label(labelFrame, $
                                  VALUE='     Spatial dimensions:  Nx')
      state.Nxtext = widget_label(labelFrame, /FRAME, /DYNAMIC_RESIZE)
      NzLabel      = widget_label(labelFrame, VALUE='  Nz')
      state.Nztext = widget_label(labelFrame, /FRAME, /DYNAMIC_RESIZE)
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      NradiusLabel = widget_label(labelFrame, $
                                  VALUE='     Spatial dimension:  Nradius')
      state.Nztext = widget_label(labelFrame, /FRAME, /DYNAMIC_RESIZE)
      RadiusLabel  = widget_label(labelFrame, VALUE='     Radius')
      RadiusLabel  = widget_label(labelFrame, /FRAME, $
                                  VALUE=string(geometry.Radius/1.0E3, $
                                               FORMAT='(E9.3)'))
      kmLabel  = widget_label(labelFrame, VALUE='[km]')
    END
  ENDCASE
  state.drawWidget = widget_draw(drawFrame, XSIZE=state.ScreenSize(0), $
                                 YSIZE=state.ScreenSize(1))

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    panelFrame  = widget_base(atmosBase, /FRAME, /EXCLUSIVE, /ROW)
    panelButton = widget_button(panelFrame, VALUE='Panel', UVALUE='PANEL')
    wireButton  = widget_button(panelFrame, VALUE='WireGrid', $
                                UVALUE='WIREGRID')
    shadeButton = widget_button(panelFrame, VALUE='ShadeSurface', $
                                UVALUE='SHADESURF')
    widget_control, panelButton, /SET_BUTTON
  ENDIF

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- atmosWidgetSetup.pro ----- ;

; -------- begin -------------------------- XViewAtmos.pro ----------- ;

PRO XViewAtmos, atmosFile, GROUP_LEADER=group_leader

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0
  IF (NOT readAtmos(atmosFile)) THEN return

  state = atmosWidgetSetup(atmosFile)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  surfr, AX=30, AZ=30

  displayAtmos, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAtmos', state.baseWidget, $
   EVENT_HANDLER='XViewAtmos_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewAtmos.pro ----------- ;
