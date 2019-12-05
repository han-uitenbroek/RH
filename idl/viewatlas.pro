; ----------------------------------------- viewatlas.pro ------------ ;

; -------- begin -------------------------- setAtlas.pro ------------- ;

FUNCTION setAtlas, state, KPNO=kpno, HAWAII=hawaii, KPK=kpk, $
                   SUMER=sumer, ATMOS=atmos, SOLFLUX=SolFlux, $
                   KPSPOT=kpspot, BRAULTNECKEL=braultneckel


  ;; --- Just to be on the safe side with memory requirements -- ----- ;

  state.sample = 0.01

  widget_control, state.minField, $
   SET_VALUE=string(state.lambdaMin, FORMAT='(F8.3)')
  widget_control, state.maxField, $
   SET_VALUE=string(state.lambdaMax, FORMAT='(F8.3)')
  widget_control, state.sampleField, $
   SET_VALUE=string(state.sample, FORMAT='(F8.5)')

  IF (keyword_set(KPNO))         THEN state.atlas = "KPNO"
  IF (keyword_set(HAWAII))       THEN state.atlas = "Hawaii_UV"
  IF (keyword_set(KPK))          THEN state.atlas = "Harvard_UV"
  IF (keyword_set(SUMER))        THEN state.atlas = "SUMER"
  IF (keyword_set(ATMOS))        THEN state.atlas = "ATMOS"
  IF (keyword_set(SOLFLUX))      THEN state.atlas = "SOLFLUX"
  IF (keyword_set(KPSPOT))       THEN state.atlas = "KPSPOT"
  IF (keyword_set(BRAULTNECKEL)) THEN state.atlas = "BRAULTNECKEL"

  widget_control, state.atlasText, SET_VALUE=state.atlas
  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setAtlas.pro ------------- ;

; -------- begin -------------------------- setAtlasParams.pro ------- ;

FUNCTION setAtlasParams, state
 
  widget_control, state.minField, GET_VALUE=lambdaMin
  widget_control, state.maxField, GET_VALUE=lambdaMax
  state.lambdaMin = lambdaMin(0)  &  state.lambdaMax = lambdaMax(0)
 
  widget_control, state.sampleField, GET_VALUE=sample
  state.sample = sample(0)  
 
  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setAtlasParams.pro ------- ;

; -------- begin -------------------------- XViewAtlas_Event.pro ----- ;

PRO XViewAtlas_Event, Event

  On = 1  &  Off = 0

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info( Event.handler, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'SET_ATLAS_PARAMS':  displayAtlas, setAtlasParams(state)

    'KPNO_ATLAS':   displayAtlas, setAtlas(state, /KPNO)
    'HAWAII_ATLAS': displayAtlas, setAtlas(state, /HAWAII)
    'HARVARD_UV':   displayAtlas, setAtlas(state, /KPK)
    'SUMER':        displayAtlas, setAtlas(state, /SUMER)
    'ATMOS':        displayAtlas, setAtlas(state, /ATMOS)
    'SOLFLUX':      displayAtlas, setAtlas(state, /SOLFLUX)
    'BRAULTNECKEL': displayAtlas, setAtlas(state, /BRAULTNECKEL)
    'KPSPOT':       displayAtlas, setAtlas(state, /KPSPOT)

    'PRINT': BEGIN
      filename = '/tmp/viewAtlas-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR
      displayAtlas, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewAtlas-'
    END

    'XLOADCT':     XLoadct, GROUP=Event.top

    'LOG_ON': IF (Event.select) THEN $
     displayAtlas, setLog(stash, On)  ELSE  displayAtlas, setLog(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Interactive display of solar spectrum atlases", $
             "", $
             "Enter minimum and maximum wavelength at the top, just", $
             "below the menubar. Third entry is the wavelength sample", $
             "interval. The true sampling interval will be an integer", $
             "times the atlas' maximum resolution", $
             "", $
             "Version 1.0, Dec 22, 1995", $
             "Han Uitenbroek (HUitenbroek@nso.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewAtlas_Event.pro ----- ;

; -------- begin -------------------------- displayAtlas.pro --------- ;

PRO displayAtlas, state

  limbColor = 150B

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  IF ( (state.lambdaMin EQ 0.0) AND (state.lambdaMax EQ 0.0) ) THEN BEGIN
    erase
    return
  ENDIF

  CASE (state.atlas) OF
    "KPNO": BEGIN
      KPNOatlas, state.lambdaMin, state.lambdaMax, SAMPLE=state.sample, $
       Iatlas_c, latlas_c
      KPNOatlas, state.lambdaMin, state.lambdaMax, SAMPLE=state.sample, $
       Iatlas_l, latlas_l, /LIMB
      ytitle = 'Relative Intensity'
      xmargin = [12, 2]
    END
    "Hawaii_UV": BEGIN
      Hawaiiatlas, state.lambdaMin, state.lambdaMax, SAMPLE=state.sample, $
       Iatlas_c, latlas_c
      ytitle = 'Intensity [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'
      xmargin = [13, 2]
    END
    "Harvard_UV": BEGIN
      KPKatlas, state.lambdaMin, state.lambdaMax, SAMPLE=state.sample, $
       Iatlas_c, latlas_c
      ytitle = 'Intensity [J m!U-2!N s!U-1!N Hz!U-1!N sr!U-1!N]'
      xmargin = [13, 2]
    END
    "SUMER": BEGIN
      sumeratlas, state.lambdaMin, state.lambdaMax, SAMPLE=state.sample, $
       Iatlas_c, latlas_c
      sumeratlas, state.lambdaMin, state.lambdaMax, SAMPLE=state.sample, $
       Iatlas_l, latlas_l, /NETWORK
      ytitle = 'Counts/sec/pixel'
      xmargin = [13, 2]
    END 
    "ATMOS": BEGIN
      ATMOSatlas, state.lambdaMin, state.lambdaMax, Iatlas_c, latlas_c
      ytitle = 'Relative Intensity'
      xmargin = [12, 2]
    END
    "SOLFLUX": BEGIN
      SolFluxatlas, state.lambdaMin, state.lambdaMax, Iatlas_c, latlas_c
      ytitle = 'Flux [J m!U-2!N s!U-1!N Hz!U-1!N]'
      xmargin = [13, 2]
    END
    "BRAULTNECKEL": BEGIN
      braultneckel, state.lambdaMin, state.lambdaMax, Iatlas_c, latlas_c
      ytitle = 'Intensity [J m!U-2!N s!U-1!N Hz!U-1!N]'
      xmargin = [13, 2]
    END
    "KPSPOT": BEGIN
      KPspotatlas, state.lambdaMin, state.lambdaMax, Iatlas_c, latlas_c
      ytitle = 'Relative Intensity'
      xmargin = [13, 2]
    END
  ENDCASE
  IF (n_elements(Iatlas_c) LE 1) THEN BEGIN
    r = dialog_message("Unable to sample atlas with specified parameters", $
                       /ERROR)
    erase
    return
  ENDIF
  IF (n_elements(Iatlas_l) GT 1) THEN limb = 1  ELSE  limb = 0
  IF (limb) THEN $
   Imax = max([Iatlas_c, Iatlas_l], MIN=Imin) $
  ELSE $
   Imax = max(Iatlas_c, MIN=Imin)

  Natlas = n_elements(latlas_c)
  trueSample = (latlas_c(Natlas-1) - latlas_c(0)) / (Natlas - 1)
  widget_control, state.sampleField, $
   SET_VALUE=string(trueSample, FORMAT='(F8.5)')

  plot, latlas_c, Iatlas_c, XTITLE='Wavelength [nm]', $
   YTITLE=ytitle, /NODATA, YRANGE=[Imin, Imax], XMARGIN=xmargin, $
   YLOG=state.log
  IF (limb) THEN oplot, latlas_l, Iatlas_l, COLOR=limbColor
  oplot, latlas_c, Iatlas_c


  IF (state.atlas EQ 'KPNO') THEN BEGIN
    rhannotate, xann(0.05), yann(0.1), TEXT='Center', MARKTYPE=0, CHARSIZE=1.2
    IF (limb) THEN $
     rhannotate, xann(0.05), yann(0.05), $
     TEXT='Limb [!7l!x=0.2]', MARKCOLOR=limbColor, MARKTYPE=0, CHARSIZE=1.2
  ENDIF

  IF (state.atlas EQ 'SUMER') THEN BEGIN
    rhannotate, xann(0.8), yann(0.9), TEXT='Quiet', MARKTYPE=0, CHARSIZE=1.2
    IF (limb) THEN $
     rhannotate, xann(0.8), yann(0.85), $
     TEXT='Network', MARKCOLOR=limbColor, MARKTYPE=0, CHARSIZE=1.2
  ENDIF
END
; -------- end ---------------------------- displayAvg.pro ----------- ;

; -------- begin -------------------------- atlasWidgetSetup.pro ----- ;

FUNCTION atlasWidgetSetup, lambdaMin, lambdaMax, sample, ATLAS=atlas

  state = {log: 0, baseWidget: 0L, drawWidget: 0L, $
           atlas: atlas,  atlasText: 0L, $
           lambdaMin: lambdaMin, lambdaMax: lambdaMax, sample: sample, $
           minField: 0L, maxField: 0L, sampleField: 0L, absolute: 0}

  state.baseWidget = widget_base( TITLE='XViewAtlas', /COLUMN, $
                                  RESOURCE_NAME='XViewAtlas', $
                                  MBAR=menuBar )

  atlasBase = widget_base( state.baseWidget, /COLUMN )

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  atlasMenu = widget_button( menuBar, VALUE='Atlas', /MENU)
  KPNObutton = widget_button(atlasMenu, VALUE='KPNO atlas (optical)', $
                             UVALUE='KPNO_ATLAS' )
  Hawaiibutton = widget_button(atlasMenu, VALUE='Hawaii (UV)', $
                               UVALUE='HAWAII_ATLAS' )
  KPKbutton = widget_button(atlasMenu, VALUE='Harvard (UV)', $
                               UVALUE='HARVARD_UV' )
  SUMERbutton = widget_button(atlasMenu, VALUE='SUMER/SOHO', $
                               UVALUE='SUMER' )
  ATLASbutton = widget_button(atlasMenu, VALUE='ATMOS (IR)', $
                               UVALUE='ATMOS' )
  ATLASbutton = widget_button(atlasMenu, VALUE='Solar Flux', $
                               UVALUE='SOLFLUX' )
  LNbutton = widget_button(atlasMenu, VALUE='Brault & Neckel', $
                             UVALUE='BRAULTNECKEL')
  KPbutton = widget_button(atlasMenu, VALUE='Kitt Peak spot atlas', $
                             UVALUE='KPSPOT')


  toolMenu = widget_button( menuBar, VALUE='Tools', /MENU )
  xloadctButton = widget_button( toolMenu, VALUE='XLoadct', $
                                 UVALUE='XLOADCT' )

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button( menuBar, VALUE='Help', /MENU, /HELP )
  infoButton = widget_button( helpMenu, VALUE='XViewAtlas', $
                              UVALUE='INFORMATION' )

  fieldFrame = widget_base( atlasBase, /FRAME, /ROW )
  blueLabel  = widget_label( fieldFrame, VALUE='Wavelength  min:' )
  state.minField = widget_text(fieldFrame, $
                               UVALUE='SET_ATLAS_PARAMS', XSIZE=8, YSIZE=1, $
                               VALUE=string(FORMAT='(F8.3)', state.lambdaMin),$
                               /EDITABLE, RESOURCE_NAME='text' )
  redLabel = widget_label( fieldFrame, VALUE='  max:' )
  state.maxField = widget_text(fieldFrame, $
                               UVALUE='SET_ATLAS_PARAMS', XSIZE=8, YSIZE=1, $
                               VALUE=string(FORMAT='(F8.3)', state.lambdaMax),$
                               /EDITABLE, RESOURCE_NAME='text' )
  sampleLabel = widget_label( fieldFrame, VALUE='    sample:' )
  state.sampleField = widget_text(fieldFrame, $
                               UVALUE='SET_ATLAS_PARAMS', XSIZE=8, YSIZE=1, $
                               VALUE=string(FORMAT='(F8.5)', state.sample),$
                               /EDITABLE, RESOURCE_NAME='text' )
  nmLabel = widget_label( fieldFrame, VALUE='[nm]' )

  atlasLabel = widget_label( fieldFrame, VALUE='   Atlas:' )
  state.atlasText = widget_label( fieldFrame, /FRAME, /DYNAMIC_RESIZE, $
                                  VALUE=state.atlas )

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame  = widget_base( atlasBase, /FRAME, /COLUMN )
  state.drawWidget = widget_draw( drawFrame, XSIZE=700, YSIZE=350 )

  fieldFrame = widget_base(atlasBase, /FRAME, /ROW, /EXCLUSIVE)
  linearButton = widget_button( fieldFrame, VALUE='Lin', $
                                UVALUE='LOG_OFF' )
  logButton    = widget_button( fieldFrame, VALUE='Log', $
                                UVALUE='LOG_ON' )
  widget_control, linearButton, /SET_BUTTON

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- atlasWidgetSetup.pro ----- ;

; -------- begin -------------------------- XViewAtlas.pro ----------- ;

PRO XViewAtlas, lambdaMin, lambdaMax, SAMPLE=sample, $
                GROUP_LEADER=group_leader, $
                HAWAII=hawaii, KPK=kpk, SUMER=sumer, ATMOS=atmos, $
                SOLFLUX=SolFlux
;+
; NAME:
;	XVIEWATLAS
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
;   --- Last modified: Tue Jun  7 16:29:24 2011 --
;-

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0
  IF (NOT keyword_set(SAMPLE)) THEN sample = 0.01
  IF (keyword_set(HAWAII)) THEN $
   atlas = "Hawaii_UV" $
  ELSE IF (keyword_set(KPK)) THEN $
   atlas = "Harvard_UV" $
  ELSE IF (keyword_set(SUMER)) THEN $
   atlas = "SUMER" $
  ELSE IF (keyword_set(ATMOS)) THEN $
   atlas = "ATMOS" $
  ELSE IF (keyword_set(SOLFLUX)) THEN $
   atlas = "SOLFLUX" $
  ELSE $
   atlas = "KPNO"

  IF (n_params(0) LT 2) THEN BEGIN
    lambdaMin = 0.0
    lambdaMax = 0.0
  ENDIF

  state = atlasWidgetSetup(lambdaMin, lambdaMax, sample, ATLAS=atlas)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayAtlas, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAtlas', state.baseWidget, $
   EVENT_HANDLER='XViewAtlas_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewAtlas.pro ----------- ;
