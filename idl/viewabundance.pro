; ----------------------------------------- viewabundance.pro -------- ;

; -------- begin -------------------------- XViewAbund_Event.pro ----- ;

PRO XViewAbund_Event, Event

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info( Event.handler, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of abundances", $
             "", $
             "Version 1.0, Aug 3, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewAbund_Event.pro ----- ;

; -------- begin -------------------------- displayAbund.pro --------- ;

PRO displayAbund, stash

@atmos.common

  widget_control, stash, GET_UVALUE=state
  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo

  hydrogen = where(atmos.elements.id EQ "H ")
  IF (atmos.elements(hydrogen).abund EQ 12.0) THEN $
   abundValues = atmos.elements.abund $
  ELSE $
   abundValues = alog10(atmos.elements.abund) + 12.0

  bar_plot, abundValues[0:29], BACKGROUND=!P.BACKGROUND, $
   BARNAMES=strtrim(atmos.elements[0:29].ID), $
   COLORS=byte(25 + 200.0 * findgen(30)/29.0), $
   /OUTLINE, YTITLE='!U10!Nlog(A) + 12'
END
; -------- end ---------------------------- displayAbund.pro --------- ;

; -------- begin -------------------------- abundWidgetSetup.pro ----- ;

FUNCTION abundWidgetSetup

  state = {baseWidget: 0L, drawWidget: 0L}

  state.baseWidget = widget_base(TITLE='XViewAbund', /COLUMN, $
                                 RESOURCE_NAME='XViewAbund', $
                                 MBAR=menuBar)

  gridBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')
                             
  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewAbund', $
                             UVALUE='INFORMATION')

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame =  widget_base(gridBase, /FRAME, /COLUMN)
  state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=400)

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- abundWidgetSetup.pro ----- ;

; -------- begin -------------------------- XViewAbund.pro ----------- ;

PRO XViewAbund, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWABUND
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
;   --- Last modified: Tue Oct 30 04:34:02 2007 --
;-

@atmos.common
@files.common



  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  IF (n_elements(atmos) EQ 0) THEN $
   IF ( (NOT readGeometry(geometryFile)) OR $
        (NOT readAtmos(atmosFile)) ) THEN return

  state = abundWidgetSetup( )
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader

  displayAbund, widget_info( state.baseWidget, /CHILD )

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAbund', state.baseWidget, $
   EVENT_HANDLER='XViewAbund_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewAbund.pro ----------- ;
