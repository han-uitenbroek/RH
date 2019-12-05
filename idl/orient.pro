PRO Orient_Event, Event

  stash = widget_info( Event.top, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'ORIENT_DISMISS': BEGIN
      widget_control, state.baseWidget, /SENSITIVE
      widget_control, Event.top, /DESTROY
    END

    'ORIENT_UPDATE': BEGIN
      WIDGET_CONTROL, state.xSlider, GET_VALUE=ax
      WIDGET_CONTROL, state.zSlider, GET_VALUE=az
      surfr, AX=ax, AZ=az
    END

    'ORIENT_RESET': BEGIN
      ax = 30  &  az = 30
      WIDGET_CONTROL, state.xSlider, SET_VALUE=ax
      WIDGET_CONTROL, state.xSlider, SET_VALUE=az
      surfr, AX=ax, AZ=az
    END

  ENDCASE
END

FUNCTION orientSetup, state, event_handler, group_leader, title

  ax = 30  &  az = 30

  baseWidget = widget_base( TITLE=title, /ROW, $
                            RESOURCE_NAME='Orientation' )

  orientBase = widget_base( baseWidget, /COLUMN )
  quitButton = widget_button( orientBase, VALUE='Dismiss', $
                              UVALUE='ORIENT_DISMISS' )

  sliderFrame = widget_base( orientBase, /COLUMN, /FRAME )
  state.zSlider = WIDGET_SLIDER( sliderFrame, MIN=-180, MAX=180, $
                                 VALUE=az, TITLE="Z angle", XSIZE=200, $
                                 UVALUE='ORIENT_UPDATE', $
                                 EVENT_PRO=Event_Handler, $
                                 RESOURCE_NAME='slider' )
  state.xSlider = WIDGET_SLIDER( sliderFrame, MIN=-180, MAX=180, $
                                 VALUE=ax, TITLE="X angle", XSIZE=200, $
                                 UVALUE='ORIENT_UPDATE', $
                                 EVENT_PRO=Event_Handler, $
                                 RESOURCE_NAME='slider' )

  resetButton = widget_button( sliderFrame, VALUE='Reset', $
                               UVALUE='ORIENT_RESET', $
                               EVENT_PRO=Event_Handler )

  widget_control, baseWidget, /REALIZE, GROUP_LEADER=group_leader
  return, baseWidget
END

PRO Orient, state, GROUP_LEADER=group_leader, EVENT_HANDLER=event_handler, $
  TITLE=title

  IF (NOT keyword_set(GROUP_LEADER)) THEN $
   group_leader = 0 $
  ELSE $
   widget_control, group_leader, SENSITIVE=0
   
  IF (NOT keyword_set(EVENT_HANDLER)) THEN BEGIN
    state = { xSlider: 0L, zSlider: 0L }
    event_handler = 'Orient_Event'
  ENDIF
  IF (NOT keyword_set(TITLE))  THEN title = 'Orientation'

  baseWidget = orientSetup(state, event_handler, group_leader, title)
  widget_control, widget_info(baseWidget, /CHILD), SET_UVALUE=state
  xmanager, 'Orient', baseWidget, GROUP_LEADER=group_leader
END
