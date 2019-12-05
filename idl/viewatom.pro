; ----------------------------------------- viewatom.pro ------------- ;

; -------- begin -------------------------- setAtom.pro -------------- ;

FUNCTION setAtom, stash, atom

  widget_control, stash, GET_UVALUE=state

  IF (ptr_valid(state.atom_ptr)) THEN free_atom, state.atom_ptr
  state.atom_ptr = ptr_new(atom)

  widget_control, stash, SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- setAtom.pro -------------- ;

; -------- begin -------------------------- XViewAtom_Event.pro ------ ;

PRO XViewAtom_Event, Event

@files.common

  Line = 0  &  Cont = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF

    'QUIT': BEGIN
      free_atom, state.atom_ptr
      widget_control, Event.top, /DESTROY
    END
      
    'TERMDIAGRAM': BEGIN
      XViewTermDiag, state.atom_ptr, GROUP_LEADER=Event.top
    END

    'LEVELLIST': BEGIN
      IF (ptr_valid((*state.atom_ptr).n_ptr)) THEN $
       XViewPops, state.atom_ptr, LEVEL=Event.index, GROUP_LEADER=Event.top $
      ELSE $
       messageresult = dialog_message('No populations available for this atom')
    END

    'BOUNDBOUND': BEGIN
      XViewTrans, state.atom_ptr, 'BOUNDBOUND', GROUP_LEADER=Event.top
    END
    'BOUNDFREE': BEGIN
      XViewTrans, state.atom_ptr, 'BOUNDFREE', GROUP_LEADER=Event.top
    END
    'FIXED': BEGIN
      XViewTrans, state.atom_ptr, 'FIXED', GROUP_LEADER=Event.top
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of atomic data", $
             "", $
             "-- Click on level to display corresponding population numbers", $
             "", $
             "Version 1.0, Jun 2, 1995", $
             "Han Uitenbroek (huitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewAtm_Event.pro ------- ;

; -------- begin -------------------------- displayAtom.pro ---------- ;

PRO displayAtom, state

  one_Ev = 1.60219E-19
  theatom = state.atom_ptr

  statusStr = ((*theatom).active) ? 'ACTIVE' : 'PASSIVE'
  widget_control, state.statusWidget, $
   SET_VALUE=string(FORMAT='(A)', statusStr)

  widget_control, state.levelWidget, $
   SET_VALUE=strtrim(string(FORMAT='(I3)', (*theatom).Nlevel), 2)
  widget_control, state.lineWidget, $
   SET_VALUE=strtrim(string(FORMAT='(I3)', (*theatom).Nline), 2)
  widget_control, state.contWidget, $
   SET_VALUE=strtrim(string(FORMAT='(I3)', (*theatom).Ncont), 2)
  widget_control, state.fixedWidget, $
   SET_VALUE=strtrim(string(FORMAT='(I3)', (*theatom).Nfixed), 2)

  widget_control, state.abundWidget, $
   SET_VALUE=strtrim(string(FORMAT='(E8.2)', (*theatom).abundance), 2)
  widget_control, state.weightWidget, $
   SET_VALUE=strtrim(string(FORMAT='(F6.2)', (*theatom).weight), 2)

  fmt1 = '(I2, ":  ", A20, 2X, F5.1, 5X, F6.3, X)'

  levels = strarr((*theatom).Nlevel)
  FOR i=0, (*theatom).Nlevel-1 DO $
   levels(i) = string(FORMAT=fmt1, i, (*theatom).labels(i), $
                      (*theatom).g(i), $
                      (*theatom).E(i) / one_eV)
  widget_control, state.levelList, SET_VALUE=levels
END
; -------- end ---------------------------- displayAtom.pro ---------- ;

; -------- begin -------------------------- atomWidgetSetup.pro ------ ;

FUNCTION atomWidgetSetup, atom, $
                          NO_INTENSITIES=no_intensities, $
                          GROUP_LEADER=group_leader

  state = {atom_ptr: ptr_new(atom), $
           baseWidget: 0L,  abundWidget: 0L,  weightWidget: 0L, $
           levelWidget: 0L,  lineWidget: 0L,  contWidget: 0L, $
           fixedWidget: 0L,  levelList: 0L,  $
           no_intensities: keyword_set(NO_INTENSITIES),  statusWidget: 0L}

  state.baseWidget = widget_base(TITLE='XViewAtom', /COLUMN, $
                                  RESOURCE_NAME='XViewAtom', MBAR=menuBar)
  atomBase = widget_base(state.baseWidget, /COLUMN)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  transitionMenu =  widget_button(menuBar, VALUE='Transitions', /MENU)
  boundboundButton = widget_button(transitionMenu, $
                                   VALUE='bound-bound', UVALUE='BOUNDBOUND')
  IF (atom.Nline EQ 0) THEN widget_control, boundboundButton, SENSITIVE=0

  boundfreeButton = widget_button(transitionMenu, $
                                  VALUE='bound-free', UVALUE='BOUNDFREE')
  IF (atom.Ncont EQ 0) THEN widget_control, boundfreeButton, SENSITIVE=0

  fixedButton = widget_button(transitionMenu, $
                              VALUE='fixed', UVALUE='FIXED')
  IF (atom.Nfixed EQ 0) THEN widget_control, fixedButton, SENSITIVE=0


  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='XViewAtom', $
                              UVALUE='INFORMATION')

  termMenu   = widget_button(menuBar, VALUE='Termdiagram', /MENU)
  termButton = widget_button(termMenu, VALUE='termdiagram', $
                              UVALUE='TERMDIAGRAM')

  quantityFrame = widget_base(atomBase, /ROW, /FRAME)
  keywordFrame = widget_base(quantityFrame, /COLUMN)
  valueFrame   = widget_base(quantityFrame, /COLUMN)

  keyword = widget_label(keywordFrame, VALUE='Status:', /ALIGN_RIGHT)
  state.statusWidget = widget_label(valueFrame, $
                                    /DYNAMIC_RESIZE, /ALIGN_LEFT)

  levelLabel = widget_label(keywordFrame, VALUE='Nlevel:', /ALIGN_RIGHT)
  state.levelWidget = widget_label(valueFrame, $
                                 /DYNAMIC_RESIZE, /ALIGN_LEFT)

  lineLabel = widget_label(keywordFrame, VALUE='Nline:', /ALIGN_RIGHT)
  state.lineWidget = widget_label(valueFrame, $
                                  /DYNAMIC_RESIZE, /ALIGN_LEFT)

  lineLabel = widget_label(keywordFrame, VALUE='Ncont:', /ALIGN_RIGHT)
  state.contWidget = widget_label(valueFrame, $
                                  /DYNAMIC_RESIZE, /ALIGN_LEFT)

  lineLabel = widget_label(keywordFrame, VALUE='fixed:', /ALIGN_RIGHT)
  state.fixedWidget = widget_label(valueFrame, $
                                   /DYNAMIC_RESIZE, /ALIGN_LEFT)

  abundLabel = widget_label(keywordFrame, VALUE='Abundance:', /ALIGN_RIGHT)
  state.abundWidget = widget_label(valueFrame, $
                                   /DYNAMIC_RESIZE, /ALIGN_LEFT)

  weightLabel = widget_label(keywordFrame, $
                             VALUE='Atomic Weight:', /ALIGN_RIGHT)
  state.weightWidget = widget_label(valueFrame, $
                                    /DYNAMIC_RESIZE, /ALIGN_LEFT)

  ; --- Establish the list widgets for levels, lines and continua.
  ;     Remember that the font has to be fixed size! The actual font is
  ;     set in the file Idl --                          --------------- ;

  levelFormat = $
   string(FORMAT='(A, 9X, "Label", 14X, "g", 6X, "E [eV]   ")', "No")

  atomFrame  = widget_base(atomBase, /COLUMN, /FRAME)
  levelFrame = widget_base(atomFrame, /COLUMN)
  levelLabel = widget_label(levelFrame, VALUE='Levels', /ALIGN_CENTER)
  levelLabel = widget_label(levelFrame, VALUE=levelFormat, $
                            RESOURCE_NAME='listheader', /ALIGN_CENTER)
  state.levelList = widget_list(levelFrame, UVALUE='LEVELLIST', $
                                RESOURCE_NAME='list', /ALIGN_CENTER, $
                                YSIZE=10, VALUE=levelFormat)

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- atomWidgetSetup.pro ------ ;

; -------- begin -------------------------- XViewAtom.pro ------------ ;

PRO XViewAtom, atom, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWATOM
;
; PURPOSE:
;	This routines provides a widget driven display of an atomic model.
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;       XVIEWATOM, Atom, GROUP_LEADER=group_leader
;
; INPUTS:
;	Atom:   Atomic data structure to be viewed.
;
; KEYWORD PARAMETERS:
;	GROUP_LEADER:    ID of parent widget in hierarchy.
;	
;
; SIDE EFFECTS:
;	The input atomic structure is put in a heap variable. The pointer to
;       this variable is freed when the routine is exited.
;
; EXAMPLE:
;	When you want to look at an input model atom for calcium use:
;
;         XVIEWATOM, RAWATOM('CaII.atom', /VACUUM_TO_AIR)
;
;	When you want to look at an output model atom for calclium use:
;
;         XVIEWATOM, READATOM('atom.out')
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Tue Apr 21 16:47:21 2009 --
;-

  IF (n_tags(atom) LE 0) THEN return
  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader=0

  IF (keyword_set(NO_INTENSITIES)) THEN $
   state = atomWidgetSetup(atom, /NO_INTENSITIES, GROUP_LEADER=group_leader) $
  ELSE $
   state = atomWidgetSetup(atom, GROUP_LEADER=group_leader)

  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayAtom, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewAtom', state.baseWidget, $
   EVENT_HANDLER='XViewAtom_Event', GROUP_LEADER=group_leader
END
; -------- end ---------------------------- XViewAtom.pro ------------ ;
