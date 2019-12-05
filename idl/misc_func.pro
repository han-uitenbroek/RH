
; -------- begin -------------------------- OrientationUpdate.pro ---- ;

PRO OrientationUpdate, state, SET=set

  IF (keyword_set(SET)) THEN BEGIN
    ax = 30
    az = 30
    WIDGET_CONTROL, state.xSlider, SET_VALUE=ax
    WIDGET_CONTROL, state.zSlider, SET_VALUE=az
  ENDIF ELSE BEGIN
    WIDGET_CONTROL, state.xSlider, GET_VALUE=ax
    WIDGET_CONTROL, state.zSlider, GET_VALUE=az
  ENDELSE

  surfr, AX=ax, AZ=az
END
; -------- end ---------------------------- OrientationUpdate.pro ---- ;


; -------- begin -------------------------- setLog.pro --------------- ;

FUNCTION setLog, stash, toggle
  widget_control, stash, GET_UVALUE=state

  state.log = toggle
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setLog.pro --------------- ;

; -------- begin -------------------------- setShade.pro ------------- ;

FUNCTION setShade, stash, toggle
  widget_control, stash, GET_UVALUE=state

  state.shade = toggle
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setShade.pro ------------- ;

; -------- begin -------------------------- setRelative.pro ---------- ;

FUNCTION setRelative, stash, toggle
  widget_control, stash, GET_UVALUE=state

  state.relative = toggle
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setRelative.pro ---------- ;

; -------- begin -------------------------- setWire.pro -------------- ;

FUNCTION setWire, stash, toggle
  widget_control, stash, GET_UVALUE=state

  state.wire = toggle
  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setWire.pro -------------- ;

; -------- begin -------------------------- setType.pro -------------- ;

FUNCTION setType, stash, type
  widget_control, stash, GET_UVALUE=state

  state.type = type
  widget_control, stash, SET_UVALUE=state
  return, state   
END
; -------- end ---------------------------- setType.pro -------------- ;

