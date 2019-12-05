; ----------------------------------------- viewaalpha.pro ----------- ;

; -------- begin -------------------------- XViewAlpha_Event.pro ----- ;

PRO XViewAlpha_Event, Event

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info( Event.handler, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': BEGIN
      ptr_free, state.cont_ptr
      widget_control, Event.top, /DESTROY
    END
 
    'LOG_ON': IF (Event.select) THEN $
     displayAlpha, setLog(stash, On)  ELSE  displayAlpha, setLog(stash, Off)

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of ionization cross-section", $
             "", $
             "Version 1.0, Feb 9, 1996", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewAlpha_Event.pro ----- ;

; -------- begin -------------------------- displayAlpha.pro --------- ;

PRO displayAlpha, state

  MEGABARN_TO_M2 = 1.0E-22
  HYDROGENIC = 3
  EXPLICIT   = 4

  thecontinuum = state.cont_ptr

  widget_control, state.NlambText, SET_VALUE=string(FORMAT='(I3)', $
                                                    (*thecontinuum).Nlambda)
  CASE ((*thecontinuum).shape) OF
    HYDROGENIC: shapeStr = "BF_HYDROGENIC"
    EXPLICIT:   shapeStr = "BF_EXPLICIT"
  ENDCASE
  widget_control, state.shapeText, SET_VALUE=shapeStr
  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo

  alpha0 = (*thecontinuum).strength / MEGABARN_TO_M2
  lambda = *((*thecontinuum).lambda_ptr)
  alpha  = *((*thecontinuum).alpha_ptr)

  IF ((*thecontinuum).shape EQ EXPLICIT) THEN BEGIN
    alpha = alpha / MEGABARN_TO_M2

    plot, lambda, alpha, XTITLE='Wavelength [nm]', $
     YTITLE='Cross-section [MegaBarn]', YLOG=state.log
    oplot, lambda, alpha0 * (lambda/(*thecontinuum).lambda0)^3, COLOR=175B
  ENDIF ELSE BEGIN
    lambda = lambda[0] + ((*thecontinuum).lambda0 - lambda[0]) * $
     findgen((*thecontinuum).Nlambda) / ((*thecontinuum).Nlambda - 1)
    alpha = alpha0 * (lambda/(*thecontinuum).lambda0)^3

    plot, lambda, alpha, XTITLE='Wavelength [nm]', $
     YTITLE='Cross-section [MegaBarn]', YLOG=state.log
  ENDELSE
  oplot, lambda, alpha, PSYM=5, COLOR=200B

  usersym, [0.0, 0.0], [-3.0, 3.0], THICK=2, COLOR=100B
  oplot, [(*thecontinuum).lambda0], [alpha0], PSYM=8

  xyouts, 0.15, 0.9, $
   'transition ' + string(FORMAT='(I2, "-", I2)', $
                          (*thecontinuum).j, (*thecontinuum).i), $
   /NORM, CHARSIZE=1.2, COLOR=100B
END
; -------- end ---------------------------- displayAlpha.pro --------- ;

; -------- begin -------------------------- alphaWidgetSetup.pro ----- ;

FUNCTION alphaWidgetSetup, cont

  state = {baseWidget: 0L, drawWidget: 0L, $
           NlambText: 0L,  shapeText: 0L,  log: 1, $
           cont_ptr: ptr_new(cont)}

  state.baseWidget = widget_base(TITLE='XViewAlpha', /COLUMN, $
                                 RESOURCE_NAME='XViewAlpha', $
                                 MBAR=menuBar)

  alphaBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewalpha', $
                             UVALUE='INFORMATION')

  labelFrame = widget_base(alphaBase, /ROW, /FRAME)
  NlambLabel = widget_label(labelFrame, VALUE='  Nlamb:')
  state.NlambText = widget_label(labelFrame, /FRAME)
  shapeLabel = widget_label(labelFrame, VALUE='    shape:')
  state.shapeText = widget_label(labelFrame, /FRAME, /DYNAMIC_RESIZE)

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame =  widget_base(alphaBase, /FRAME, /COLUMN)
  state.drawWidget = widget_draw( drawFrame, XSIZE=500, YSIZE=375)

  linlogFrame  = widget_base(drawFrame, /ROW, /EXCLUSIVE, /FRAME)
  linearButton = widget_button(linlogFrame, VALUE='Lin', $
                               UVALUE='LOG_OFF')
  logButton    = widget_button(linlogFrame, VALUE='Log', $
                                UVALUE='LOG_ON')
  widget_control, logButton, /SET_BUTTON

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- alphaWidgetSetup.pro ----- ;

; -------- begin -------------------------- XViewAlpha.pro ----------- ;

PRO XViewAlpha, transition, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWALPHA
;
; PURPOSE:
;       Displays bound-free radiative cross-section of specified transition
;       as function of wavelength.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XVIEWALPHA, Transition, GROUP_LEADER=group_leader
;
; INPUTS:
;	Transition: Structure containing bound-free transition.
;
; KEYWORD PARAMETERS:
;	GROUP_LEADER:  Handle of parent widget in hierarchy
;
; SIDE EFFECTS:
;	The input structure is placed in a heap variable. At exit the pointer
;       to this variable is freed.
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Jan 10 15:01:37 2018 --
;-

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  state = alphaWidgetSetup( transition )
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  displayAlpha, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAlpha', state.baseWidget, $
   EVENT_HANDLER='XViewAlpha_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewAlpha.pro ----------- ;
