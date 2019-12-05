; ----------------------------------------- viewatom.pro ------------- ;

; -------- begin -------------------------- displayTrans.pro --------- ;

PRO displayTrans, state

  theatom = state.atom_ptr

  fmtbb = '(I2, " - ", I2, 4X, F9.3, 5X, E8.2, 2X, A5)'
  fmtbf = '(I2, " - ", I2, 2X, F9.3, 5X, E8.2, X)'
  fmtfx = '(I2, " - ", I2, 2X, F9.3, 4X, E8.2, 2X, A18, 2X, F7.1, 2X, A11)'

  statusStr = ((*theatom).active) ?  'ACTIVE' : 'PASSIVE'
  widget_control, state.statusWidget, $
   SET_VALUE=string(FORMAT='(A)', statusStr)

  case (state.type) of
    'BOUNDBOUND': BEGIN
      widget_control, state.lineList, SET_VALUE=strarr(10)
      IF ((*theatom).Nline GT 0) THEN lines = strarr((*theatom).Nline)
      FOR kr=0, (*theatom).Nline-1 DO BEGIN
        CASE ((*theatom).transition(kr).shape) OF
          0:  shapeStr = 'GAUSS'
          1:  shapeStr = 'VOIGT'
          2:  shapeStr = ' PRD '
        ENDCASE
        wavelength = (*theatom).transition(kr).lambda0
        lines(kr) = string(FORMAT=fmtbb, $
                           (*theatom).transition(kr).j, $
                           (*theatom).transition(kr).i, $
                           wavelength, (*theatom).transition(kr).strength, $
                           shapeStr)
      ENDFOR
      IF ((*theatom).Nline GT 0) THEN $
       widget_control, state.lineList, SET_VALUE=lines
    END
    'BOUNDFREE': BEGIN
      widget_control, state.contList, SET_VALUE=strarr(10)
      IF ((*theatom).Ncont GT 0) THEN continua = strarr((*theatom).Ncont)
      FOR krc=0, (*theatom).Ncont-1 DO BEGIN
        kr = (*theatom).Nline + krc
        wavelength = (*theatom).transition(kr).lambda0
        continua(krc) = string(FORMAT=fmtbf, $
                               (*theatom).transition(kr).j, $
                               (*theatom).transition(kr).i, $
                               wavelength, (*theatom).transition(kr).strength)
      ENDFOR
      IF ((*theatom).Ncont GT 0) THEN $
       widget_control, state.contList, SET_VALUE=continua
    END
    'FIXED': BEGIN
      widget_control, state.fixedList, SET_VALUE=strarr(10)
      IF ((*theatom).Nfixed GT 0) THEN fixtrans = strarr((*theatom).Nfixed)
      FOR kf=0, (*theatom).Nfixed-1 DO BEGIN
        CASE ((*theatom).fixed[kf].option) OF
          0:  optionStr = 'TRAD_ATMOSPHERIC  '
          1:  optionStr = 'TRAD_PHOTOSPHERIC '
          2:  optionStr = 'TRAD_CHROMOSPHERIC'
        ENDCASE
        CASE ((*theatom).fixed[kf].type) OF
          0:  typeStr = 'BOUND_BOUND'
          1:  typeStr = 'BOUND_FREE '
        ENDCASE
        fixtrans[kf] = string(FORMAT=fmtfx, $
                              (*theatom).fixed[kf].j, $
                              (*theatom).fixed[kf].i, $
                              (*theatom).fixed[kf].lambda0, $
                              (*theatom).fixed[kf].strength, $
                              optionStr, (*theatom).fixed[kf].Trad, typeStr)
      ENDFOR
      IF ((*theatom).Nfixed GT 0) THEN $
       widget_control, state.fixedList, SET_VALUE=fixtrans
    END
  END

END
; -------- begin -------------------------- XViewAtom_Event.pro ------ ;

PRO XViewTrans_Event, Event

@files.common
@atmos.common

  Line = 0  &  Cont = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info( Event.handler, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF

    'QUIT': BEGIN
      IF (state.free) THEN free_atom, state.atom_ptr
      widget_control, Event.top, /DESTROY
    END

    'TERMDIAGRAM': BEGIN
      XViewTermDiag, state.atom_ptr, GROUP_LEADER=Event.top
    END

    'LINELIST': BEGIN
      theatom = state.atom_ptr
      kr = Event.index
      IF ((*theatom).active) THEN $
       XViewIe, BLUE=(*theatom).transition(kr).Nblue, $
                RED=(*theatom).transition[kr].Nblue + $
                (*theatom).transition[kr].Nlambda - 1, GROUP_LEADER=Event.top
    END
    'CONTLIST': BEGIN
      theatom = state.atom_ptr
      kr = Event.index + (*theatom).Nline
      IF ((*theatom).active) THEN $
       XViewIe, BLUE=(*theatom).transition(kr).Nblue, $
       RED=(*theatom).transition[kr].Nblue + $
       (*theatom).transition[kr].Nlambda - 1, GROUP_LEADER=Event.top $
      ELSE $
       XViewAlpha, (*theatom).transition(kr), GROUP_LEADER=Event.top
    END

    'INFORMATION': result = dialog_message( /INFORMATION, $
            ["Display of atomic transition data", $
             "", $
             "-- Click on transition to display corresponding", $
             "intensities in graph", $
             "", $
             "Version 1.0, Jun 2, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"] )
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewTrans_Event.pro ----- ;

; -------- begin -------------------------- transWidgetSetup.pro ----- ;

FUNCTION transWidgetSetup, atom_ptr, TYPE=type, FREE=free

  IF (keyword_set(FREE)) THEN free = 1  ELSE  free = 0
  IF (keyword_set(TYPE)) THEN transtype = type  ELSE  transtype = 'BOUNDBOUND'

  state = {atom_ptr: atom_ptr, free: free, $
           baseWidget: 0L,  lineList: 0L,  contList: 0L,  fixedList: 0L, $
           type: transtype,  statusWidget: 0L}

  state.baseWidget = widget_base(TITLE='XViewTrans', /COLUMN, $
                                 RESOURCE_NAME='XViewTrans', MBAR=menuBar)
  transBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button( menuBar, VALUE='File', /MENU)
  quitButton = widget_button( fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button( menuBar, VALUE='Help', /MENU, /HELP )
  infoButton = widget_button( helpMenu, VALUE='XViewTrans', $
                              UVALUE='INFORMATION' )

  termMenu   = widget_button( menuBar, VALUE='Termdiagram', /MENU )
  termButton = widget_button( termMenu, VALUE='termdiagram', $
                              UVALUE='TERMDIAGRAM' )

  statusFrame = widget_base(transBase, /ROW, /ALIGN_RIGHT)
  statusLabel = widget_label(statusFrame, VALUE='Status: ')
  state.statusWidget = widget_label(statusFrame, /FRAME, /DYNAMIC_RESIZE)

  ; --- Establish the list widgets for lines, continua and fixed transitions
  ;     Remember that the font has to be fixed size! --        ------- ;

  transFrame = widget_base( transBase, /ROW, /FRAME )

  CASE (state.type) OF
    'BOUNDBOUND': BEGIN
      lineFormat = ' j -  i  lambda [nm]   Aji [s^-1]  shape   '
      lineFrame  = widget_base(transFrame, /COLUMN)
      lineLabel  = widget_label(lineFrame, VALUE='Bound-Bound Transitions')
      lineLabel  = widget_label(lineFrame, VALUE=lineFormat, $
				RESOURCE_NAME='listheader')
      state.lineList = widget_list(lineFrame, UVALUE='LINELIST', $
                                   RESOURCE_NAME='list', $
                                   YSIZE=10, VALUE=lineFormat)
    END
    'BOUNDFREE': BEGIN
      contFormat = ' j -  i  edge [nm]  alpha [m^2]   '
      contFrame  = widget_base(transFrame, /COLUMN)
      contLabel  = widget_label(contFrame, VALUE='Bound-Free Transitions')
      contLabel  = widget_label(contFrame, VALUE=contFormat, $
				RESOURCE_NAME='listheader')
      state.contList = widget_list(contFrame, UVALUE='CONTLIST', $
                                   RESOURCE_NAME='list', $
                                   YSIZE=10, VALUE=contFormat)
    END
    'FIXED': BEGIN
      fixedFormat = ' j -  i  lambda [nm]  strength         option' + $
                    '         Trad      type     '
      fixedFrame  = widget_base(transFrame, /COLUMN)
      fixedLabel  = widget_label(fixedFrame, VALUE='Fixed Transitions')
      fixedLabel  = widget_label(fixedFrame, VALUE=fixedFormat, $
				 RESOURCE_NAME='listheader')
      state.fixedList = widget_list(fixedFrame, UVALUE='FIXEDLIST', $
                                    RESOURCE_NAME='list', $
                                    YSIZE=10, VALUE=fixedFormat)
    END
  END

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- transWidgetSetup.pro ----- ;

; -------- begin -------------------------- XViewTrans.pro ----------- ;

PRO XViewTrans, atom_ptr, type, GROUP_LEADER=group_leader, FREE=free

;+
; NAME:
;	XVIEWTRANS
;
; PURPOSE:
;	This procedure displays the specified transitions (either of the three
;       possible cases: bound-bound, bound-free, fixed) of the atomic model
;       pointed to by atom_ptr.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XVIEWTRANS, Atom_ptr, type, GROUP_LEADER=group_leader, FREE=free
;
; INPUTS:
;	Atom_ptr:  Pointer to the heap variable containing atomic data
;                  structure
;
;       type:      Type of transition. Has to be string of either:
;                  "BOUNDBOUND", "BOUNDFREE", "FIXED".
;
; KEYWORD PARAMETERS:
;	GROUP_LEADER:   ID of parent widget in hierarchy.
;       FREE:           Free the heap variables on exit
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Tue Apr 21 16:58:13 2009 --
;-


  IF (NOT ptr_valid(atom_ptr)) THEN return
  IF (n_params(0) LT 2) THEN return
  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  CASE (type) OF
    'BOUNDBOUND': $
     state = transWidgetSetup(atom_ptr, TYPE='BOUNDBOUND', $
                              FREE=keyword_set(FREE))
    'BOUNDFREE': $
     state = transWidgetSetup(atom_ptr, TYPE='BOUNDFREE', $
                              FREE=keyword_set(FREE))
    'FIXED': $
     state = transWidgetSetup(atom_ptr, TYPE='FIXED', $
                              FREE=keyword_set(FREE))
    ELSE: BEGIN
      result = dialog_message("Not an allowed type of transition: " + type, $
                              /ERROR)
      return
    END
  END

  displayTrans, state
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader


  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewTrans', state.baseWidget, $
   EVENT_HANDLER='XViewTrans_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewTrans.pro ----------- ;
