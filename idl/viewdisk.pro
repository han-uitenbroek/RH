; ----------------------------------------- viewdisk.pro ------------- ;

; -------- begin -------------------------- XViewDisk_Event.pro ------ ;

PRO XViewDisk_Event, Event

@files.common



  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  widget_control, Event.id, GET_UVALUE=Action
  stash = widget_info(  Event.top, /CHILD )
  widget_control, stash, GET_UVALUE=state

  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'SLIDER':  drawDisk, setLambdaNo(stash, Event.value)

    'NEWSPECTRUMFILE': BEGIN
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of disk intensity", $
             "", $
             "", $
             "Version 1.0, Jun 2, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewJ_Event.pro --------- ;

; -------- begin -------------------------- WidgetSetup_J.pro -------- ;

FUNCTION diskWidgetSetup, lambda, lambdaDisplay

  state = {baseWidget: 0L, drawWidget: 0L, $
           ScreenSize: fix([450, 350]), ScaleFactor: float([1.0, 1.0]), $
           radius: 200, $
           lambdaText: 0L, lambdaSlider: 0L, lambdaNo: lambdaDisplay}

  state.baseWidget = widget_base( TITLE='XViewDisk', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='XViewDisk' )

  diskBase = widget_base( state.baseWidget, /COLUMN )

  fileMenu  = widget_button( menuBar, VALUE='File', /MENU)
  openSpectrumButton = widget_button( fileMenu, $
                                      VALUE='open spectrum file', $
                                   UVALUE='NEWSPECTRUMFILE' )
  quitButton = widget_button( fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button( menuBar, VALUE='Help', /MENU, /HELP )
  infoButton = widget_button( helpMenu, VALUE='XViewDisk', $
                              UVALUE='INFORMATION' )

  ;; --- Draw frame base widget --                      ------------- ;;

  drawFrame = widget_base( diskBase, /FRAME, /COLUMN )

  waveFrame = widget_base( drawFrame, /ROW )
  waveLabel = widget_label( waveFrame, VALUE='Wavelength: ')
  state.lambdaText = widget_label( waveFrame, /FRAME, $
                                   VALUE=string(FORMAT='(F10.3,X)','') )
  nmLabel = widget_label( waveFrame, VALUE='[nm]')

  state.drawWidget = widget_draw( drawFrame, XSIZE=state.ScreenSize(0), $
                                  YSIZE=state.ScreenSize(1) )

  waveFrame = widget_base( diskBase, /FRAME )
  state.lambdaSlider = widget_slider( waveFrame, TITLE='Wavelength No:', $
                                      UVALUE='SLIDER', XSIZE=255, $
                                      MAXIMUM=n_elements(lambda) - 1, $
                                      VALUE=lambdaDisplay )

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- WidgetSetup_Disk.pro --- ;

; -------- begin -------------------------- drawDisk.pro ----------- ;

PRO drawDisk, state

  COMMON screen_common, ScreenSize, ScaleFactor
@geometry.common
@spectrum.common

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo
  ScreenSize = state.ScreenSize  &  ScaleFactor = state.ScaleFactor
  erase

  lambdaDisplay = state.lambdaNo
  Nr = long(state.Radius)
  d  = shift(dist(Nr, Nr), Nr/2, Nr/2)
  mu = sqrt((1.0 - (d/float(Nr/2))^2) > 0.0)
  zeros = where(mu EQ 0.0)
  mu = reform(mu, Nr*Nr, /OVERWRITE)

  tabinv, geometry.xmu, mu, mueff

  imlambda = reform(linear(cmprss(spectrum.I(state.lambdaNo, *)), mueff), $
                    Nr, Nr)
  imlambda[zeros] = 0.0
  panel, scaleimg_idl(imlambda, ScreenSize[1] - 50, ScreenSize[1] - 50), $
   SCALETEXT='Intensity [J m!U-2!N Hz!U-1!N s!U-1!N sr!U-1!N]'

END
; -------- end ---------------------------- drawDisk.pro ------------- ;

; -------- begin -------------------------- XViewDisk.pro ------------ ;

PRO XViewDisk, GROUP_LEADER=group_leader

@spectrum.common
@geometry.common
@files.common



  IF (n_elements(spectrum) EQ 0) THEN $
   IF (NOT readSpectrum(spectrumFile)) THEN return

  IF (NOT keyword_set( GROUP_LEADER)) THEN BEGIN
    group_leader  = 0
    lambdaDisplay = 0
  ENDIF

  state = diskWidgetSetup( spectrum.lambda, lambdaDisplay )
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  drawDisk, setLambdaNo( widget_info(state.baseWidget, /CHILD), lambdaDisplay )

  ;; --- Register with the XManager --               -------------- ;;

  xmanager, 'XViewDisk', state.baseWidget, $
   EVENT_HANDLER='XViewDisk_Event', GROUP_LEADER=group_leader

END
; -------- end ---------------------------- XViewDisk.pro ------------ ;
