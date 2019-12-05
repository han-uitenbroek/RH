; -------- file: -------------------------- viewtermdiag.pro --------- ;

; -------- begin -------------------------- setBox.pro --------------- ;

PRO setBox, stash, x, y
  widget_control, stash, GET_UVALUE=state

  state.box_x = x
  state.box_y = y
  widget_control, stash, SET_UVALUE=state
END
; -------- end ---------------------------- setBox.pro --------------- ;

; -------- begin -------------------------- drawBox.pro -------------- ;

PRO drawBox, stash, x, y
  widget_control, stash, GET_UVALUE=state

  widget_control, state.drawWidget, GET_VALUE=WindowNo
  wset, WindowNo

  bx = state.box_x  &  by = state.box_y
  plots, /DEVICE, COLOR=!P.BACKGROUND, $
   [bx(0), bx(1), bx(1), bx(0), bx(0)], [by(0), by(0), by(1), by(1), by(0)]

  bx(1) = x  &  by(1) = y
  plots, /DEVICE, COLOR=!P.COLOR, $
   [bx(0), bx(1), bx(1), bx(0), bx(0)], [by(0), by(0), by(1), by(1), by(0)]
  state.box_x = bx  &  state.box_y = by

  widget_control, stash, SET_UVALUE=state
END
; -------- end ---------------------------- drawBox.pro -------------- ;

; -------- begin -------------------------- setZoomState.pro --------- ;

PRO setZoomState, stash, value
  widget_control, stash, GET_UVALUE=state

  state.zoomActive = value
  widget_control, stash, SET_UVALUE=state
END
; -------- end ---------------------------- setZoomState.pro --------- ;

; -------- begin -------------------------- setRange.pro ------------- ;

FUNCTION setRange, stash, xrange, yrange
  widget_control, stash, GET_UVALUE=state

  state.xrange = cmprss(xrange)
  state.yrange = cmprss(yrange)
  widget_control, stash, SET_UVALUE=state

  return, state
END
; -------- end ---------------------------- setRange.pro ------------- ;

; -------- begin -------------------------- getTerms.pro ------------- ;

FUNCTION getTerms, theatom

  ;; --- This routine sorts out the structure of the atomic model and
  ;;     returns it as a structure --                     -------------;

  Nlevel = (*theatom).Nlevel
  IF ((*theatom).stage(Nlevel-1) GT (*theatom).stage(Nlevel-2)) THEN $
   Nlevel = Nlevel - 1

  config  = strarr(Nlevel)  &  term  = config   &  xlabel = config
  iconfig = intarr(Nlevel)  &  iterm = iconfig  &  iorbit = iconfig
  multiplet = strarr(Nlevel)  &  imultiplet = intarr(Nlevel)

  Nline = (*theatom).Nline
  IF (Nline GT 0) THEN BEGIN
    linegr  = intarr(2, Nline)
    ilinegr = intarr(Nline)
  ENDIF ELSE BEGIN
    linegr  = 0
    ilinegr = 0
  ENDELSE

  orbits = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L']
  Norbit = n_elements(orbits)
  iorbit = intarr(Nlevel)

  Nconfig = 0  &  Nterm = 0  &  Nmultiplet = 0

  FOR i=0, Nlevel-1 DO BEGIN
    foundConfig = 0  &  foundTerm = 0  &  foundMultiplet = 0

;;    words = str_sep((*theatom).labels(i), ' ')
    words = strsplit((*theatom).labels(i), ' ', /EXTRACT)
    words = words(where(words ne ''))  &  Nwords = n_elements(words)

    ;; --- Find configuration and term for each level --  ------------ ;
    ;; Rules:
    ;;
    ;; - label = 'O I 2P4 3PE 2' --> config = 'O I 2P4' and term = '3PE'
    ;;
    ;; - label = 'O I 2P4 3PE'   --> config = 'O I 2P4' and term = '3PE'

    nw = Nwords - 1
    WHILE ((strpos(words(nw), 'E') LT 0) AND  $
           (strpos(words(nw), 'O') LT 0)) DO nw = nw - 1
    config(i) = words(0)
    FOR nww=1,nw-1 DO config(i) = config(i) + ' ' + words(nww)
    term(i)   = words(nw)
    multiplet(i) = config(i) + ' ' + term(i)

    ;; --- Number different configurations --             ------------ ;

    FOR ii=0, i-1 DO IF (config(i) eq config(ii)) THEN BEGIN
      iconfig(i)  = iconfig(ii)
      foundConfig = 1
    ENDIF
    IF (NOT foundConfig) THEN BEGIN
      iconfig(i) = Nconfig
      Nconfig    = Nconfig + 1
    ENDIF

    ;; --- Number different terms. These also label the x-axis --  --- ;

    FOR ii=0, i-1 DO  IF (term(i) eq term(ii)) THEN BEGIN
      iterm(i)  = iterm(ii)
      foundTerm =1
    ENDIF
    IF (NOT foundTerm) THEN BEGIN
      iterm(i) = Nterm
      xlabel(Nterm) = term(i) 
      Nterm = Nterm + 1
    ENDIF

    ;; --- Number multiplets with same config and term --   ---------- ;

    FOR ii=0, i-1 DO  IF (multiplet(i) eq multiplet(ii)) THEN BEGIN
      imultiplet(i)  = imultiplet(ii)
      foundMultiplet = 1
    ENDIF
    IF (NOT foundMultiplet) THEN BEGIN
      imultiplet(i) = Nmultiplet
      Nmultiplet = Nmultiplet + 1
    ENDIF

    theOrbit = strmid(words(nw), 1, 1)
    FOR ii=0, Norbit-1 DO IF (theOrbit EQ orbits(ii)) THEN iorbit(i) = ii
  ENDFOR

  ;; --- Number lines between each multiplet --           ------------ ;

  Nlinegr = 0
  FOR kr=0, (*theatom).Nline-1 DO BEGIN
    foundLinegr = 0
    linegr(0, kr) = imultiplet((*theatom).transition(kr).i)
    linegr(1, kr) = imultiplet((*theatom).transition(kr).j)
    
    FOR ikr=0, kr-1 DO $
     IF (linegr(0, kr) EQ linegr(0, ikr) AND $
         linegr(1, kr) EQ linegr(1, ikr)) THEN BEGIN
      ilinegr(kr) = ilinegr(ikr)
      foundLinegr = 1
    ENDIF
    IF (NOT foundLinegr) THEN BEGIN
      ilinegr(kr) = Nlinegr
      Nlinegr = Nlinegr + 1
    ENDIF
  ENDFOR

  return, { config: config, Nconfig: Nconfig, iconfig: iconfig, $
            term: term, Nterm: Nterm, iterm: iterm, $
            multiplet: multiplet, Nmultiplet: Nmultiplet, $
            imultiplet: imultiplet, $
            linegr: linegr, ilinegr: ilinegr, Nlinegr: Nlinegr, $
            iorbit: iorbit, xlabel: xlabel(0:Nterm-1) }
END
; -------- end ---------------------------- getTerms.pro ------------- ;

; -------- begin -------------------------- orderMultiplet.pro ------- ;

PRO orderMultiplet, td

  ;; --- Order the different terms and x-labels according to their
  ;;     multiplet system (singlet, doublet, ..., up to septet) 
  ;;     and orbital --                                          ----- ;

  Nacc   = 0
  Nlabel = n_elements(td.xlabel)  &  index = intarr(Nlabel)

  ;; --- First sort according to multiplet system --    -------------- ;

  FOR n=1, 7 DO BEGIN
    mstr = string(FORMAT='(I1)', n)
    multiplet = where(strmid(td.xlabel, 0, 1) EQ mstr, Nmultiplet)
    IF (Nmultiplet GT 0) THEN BEGIN
      IF (Nmultiplet GT 1) THEN BEGIN
        iorbit = intarr(Nmultiplet)
        FOR i=0,Nmultiplet-1 DO BEGIN
          io = where(td.xlabel(multiplet(i)) EQ td.term)
          iorbit(i) = td.iorbit(io(0))
        ENDFOR

        ;; --- Then sort according to orbital --        -------------- ;

        index(Nacc: Nacc+Nmultiplet-1) = multiplet(sort(iorbit))
      ENDIF ELSE index(Nacc) = multiplet
      Nacc = Nacc + Nmultiplet
    ENDIF
  ENDFOR

  iterm = td.iterm
  td.xlabel = td.xlabel(index)
  FOR i=0, Nlabel-1 DO $
   iterm(where(td.iterm EQ index(i))) = i
  td.iterm = iterm
END
; -------- end ---------------------------- orderMultiplet.pro ------- ;

; -------- begin -------------------------- getOffset.pro ------------ ;

PRO getOffset, td, eV, xOffset, yOffset

  ;; --- Find horizontal and vertical offsets to prevent levels
  ;;     and transitions from overlapping --             ------------- ;

  dx = 0.025  &  frac_y = 0.01

  NeV = n_elements(eV)
  Nlevel  = n_elements(td.config)
  xOffset = fltarr(Nlevel)  &  yOffset =  xOffset

  ;; First go through levels with the same configuration

  FOR i=0, td.Nmultiplet-1 DO BEGIN
    index = where(i EQ td.imultiplet, Nindx)

    IF (Nindx GT 1) THEN BEGIN
      xOffset(index) = xOffset(index) - dx * findgen(Nindx)

      dE = eV(index) - total(eV(index))/Nindx
      deltaE = max(dE) - min(dE)
      IF (deltaE EQ 0.0) THEN BEGIN
        yOffset[index] = 0.0
      ENDIF ELSE BEGIN
        dE = frac_y / min([deltaE/eV(NeV-1), frac_y]) * dE
        yOffset(index) = yOffset(index) + dE
      ENDELSE
    ENDIF
  ENDFOR

  ;; Then go through the levels with the same terms

  FOR i=0, td.Nterm DO BEGIN
    index = where(i EQ td.iterm, Nindx)
    IF (Nindx GT 1) THEN BEGIN

      ;; Shift only if not in the same configuration

      FOR ii=1, Nindx-1 DO $
       IF (td.iconfig(index(ii)) NE td.iconfig(index(ii-1))) THEN $
       xOffset(index(ii)) = xOffset(index(ii)) - dx * ii
    ENDIF
  ENDFOR
END
; -------- end ---------------------------- getOffset.pro ------------ ;

; -------- begin -------------------------- drawTermDiag.pro --------- ;

PRO drawTermDiag, state

  theatom = state.atom_ptr

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  levelwidth = 0.1  &  dline = 0.3
  headSize = !D.X_SIZE/128.

  GAUSS = 0  &  VOIGT = 1  &  PRD = 2
  levelColor = 100B  &  lineColor = 135B  &  contColor = 200B
  PRDlineColor = 16B

  one_eV = 1.60219E-19
  eV = (*theatom).E / one_eV

  td = getTerms(state.atom_ptr)
  orderMultiplet, td
  getOffset, td, eV, xOff, yOff

  xtickv = indgen(td.Nterm)  &  xticks = td.Nterm-1

  ;; --- In case the zoom function has been activated --  ------------ ;

  IF ( (state.xrange(0) NE 0.0)  AND  (state.xrange(1) NE 0.0) ) THEN BEGIN
   xtickv = where( (xtickv GE state.xrange(0)) AND $
                   (xtickv LE state.xrange(1)) )
   xticks = n_elements(xtickv) - 1
  ENDIF

  ;; --- Establish the frame --                          ------------- ;

  plot, /NODATA, [-0.2, td.Nterm-0.5], [-0.4, eV((*theatom).Nlevel-1)+0.4], $
   XTICKS=xticks, XTICKLEN=-0.02, XTICKV=xtickv, $
   XTICKNAME=td.xlabel(xtickv), XRANGE=state.xrange, XSTYLE=9, $
   YTICKLEN=-0.02, YTITLE='Energy [eV]', YRANGE=state.yrange, YSTYLE=9

  levelsymmx = [-1, 1] * levelwidth
  levelsymmy = [0.0, 0.0]

  ;; --- Draw the levels first using both x and y-offset --  --------- ;

  oplot, [-levelWidth, td.Nterm-(1.0-levelWidth)], $
   [1, 1]*eV((*theatom).Nlevel-1), THICK=2
  FOR i=0, (*theatom).Nlevel-2 DO $
   oplot, td.iterm(i) + levelsymmx + xOff(i), $
   eV(i) + levelsymmy + yOff(i), THICK=2, COLOR=levelColor

  ;; --- Add labels for each multiplet --                ------------- ;

  FOR i=0, td.Nmultiplet-1 DO BEGIN
    index = where(i EQ td.imultiplet, Nindx)
    xyouts, td.iterm(index(0))+ 1.2*levelWidth, total(eV(index))/Nindx, $
     td.config(index(0)), CHARSIZE=max([0.8, 2.5/td.Nterm^0.7])
  ENDFOR
  xyouts, 0.0, (1.01)*eV((*theatom).Nlevel-1), $
   (*theatom).labels((*theatom).Nlevel-1)

  ;; --- Then the line transitions --                    ------------- ;

  FOR kr=0, (*theatom).Nline-1 DO BEGIN
    i = (*theatom).transition(kr).i
    j = (*theatom).transition(kr).j
    x = [td.iterm(j) + xOff(j), td.iterm(i) + xOff(i)]
    y = [eV(j) + yOff(j), eV(i) + yOff(i)]
    IF ((*theatom).transition(kr).shape EQ PRD) THEN $
     color = PRDlineColor $
    ELSE $
     color = lineColor

    ;; --- Check for inter-system lines --               ------------- ;

    IF (strmid(td.term(j), 0, 1) NE strmid(td.term(i), 0, 1)) THEN BEGIN
      oplot, x, y, COLOR=color, LINE=2
    ENDIF ELSE BEGIN
      arrow, x(0), y(0), x(1), y(1), /SOLID, /DATA, $
       HSIZE=headSize, COLOR=color
    ENDELSE
  ENDFOR

  ;; --- The continua --                                 ------------- ;

  FOR kr=(*theatom).Nline, (*theatom).Nline+(*theatom).Ncont-1 DO BEGIN
    i = (*theatom).transition(kr).i
    j = (*theatom).transition(kr).j
    arrow, td.iterm(i) + xOff(i), eV(i) + yOff(i), $
     td.iterm(i) + xOff(i), eV(j), /SOLID, /DATA, $
     HSIZE=headSize, COLOR=contColor
  ENDFOR

  ;; --- And finally draw the groups of lines between multiplets ----- ;

  IF (td.Nlinegr GT 0) THEN lineOffs = fltarr(td.Nlinegr)
  FOR kr=0, (*theatom).Nline-1 DO BEGIN
    i = (*theatom).transition(kr).i  &  j = (*theatom).transition(kr).j
    x = [td.iterm(i) + xOff(i), td.iterm(j) + xOff(j)]
    y = [eV(i) + yOff(i), eV(j) + yOff(j)]

    lambdalabel = $
     strtrim(string(FORMAT='(F9.2)', (*theatom).transition(kr).lambda0), 2)

    dc = convert_coord(x, y, /TO_DEVICE)
    orient = !RADEG*atan((dc(1,1) - dc(1,0)) / (dc(0,1) - dc(0,0)))
    IF ((x(1) - x(0)) GT 0.0) THEN $
      ff = 0.2 + lineOffs(td.ilinegr(kr)) $
    ELSE $
      ff = 0.8 - lineOffs(td.ilinegr(kr))
    xpos = x(0) + ff*(x(1) - x(0))
    ypos = y(0) + ff*(y(1) - y(0))

    IF (!D.NAME EQ 'PS') THEN BEGIN
      charsize = 0.7 
      color = ((*theatom).transition(kr).shape EQ PRD) ? $
       PRDlineColor : lineColor
    ENDIF ELSE BEGIN
      charsize = max([1.0, 2.3/td.Nterm^0.7])
      color = !P.COLOR
    ENDELSE

    xyouts, xpos, ypos, lambdalabel, ORIENTATION=orient, $
     CHARSIZE=charsize, ALIGNMENT=0.5, COLOR=color

    lineOffs(td.ilinegr(kr)) = lineOffs(td.ilinegr(kr)) + dline
  ENDFOR
END
; -------- end ---------------------------- drawTermDiag.pro --------- ;

; -------- begin -------------------------- XViewTerm_Event.pro  ----- ;

PRO XViewTerm_Event, Event

@files.common

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info( Event.handler, /CHILD )
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = '/tmp/viewTerm-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR
      drawTermDiag, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewTerm-'
    END

    'ZOOM': setZoomState, stash, 1

    'RESET': BEGIN
      setZoomState, stash, 0
      drawTermDiag, setRange(stash, [0, 0], [0, 0])
    END

    'BUTTON_EVENT': BEGIN
      IF (state.zoomActive GT 0) THEN BEGIN
        CASE (Event.type) OF
          0: BEGIN
            setBox, stash, [Event.x, 0], [Event.y, 0]
            setZoomState, stash, 2
          END
          1: BEGIN
            setBox, stash, [state.box_x(0), Event.x], $
                           [state.box_y(0), Event.y]
            xy = convert_coord( /DEVICE, /TO_DATA, $
                                state.box_x, state.box_y )
            drawTermDiag, setRange( stash, xy(0,*), xy(1, *) )
            setZoomState, stash, 0
          END
          2:  IF (state.zoomActive EQ 2) THEN $
           drawBox, stash, Event.x, Event.y
        ENDCASE
      ENDIF
    END

    'XLOADCT':  XLoadct, GROUP=Event.top
    'XPALETTE': XPalette, GROUP=Event.top

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of term diagram", $
             "", $
             "Version 1.0, Jun 19, 1995", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE

END
; -------- end ---------------------------- XViewTerm_Event.pro ------ ;

; -------- begin -------------------------- termWidgetSetup.pro ------ ;

FUNCTION termWidgetSetup, atom_ptr

  xsize = 640  &  ysize = 512

  state = {atom_ptr: atom_ptr, baseWidget: 0L, drawWidget: 0L, $
           xrange: [0.0, 0.0], yrange: [0.0, 0.0], zoomActive: 0, $
           box_x: [0L, 0L], box_y: [0L, 0L]}

  state.baseWidget = widget_base( TITLE='XViewTermDiag', /COLUMN, $
                                  RESOURCE_NAME='XViewTermDiag', $
                                  MBAR=menuBar )

  termdiagBase = widget_base( state.baseWidget, /COLUMN )

  fileMenu    = widget_button(menuBar, VALUE='File', /MENU)
  quitButton  = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                              RESOURCE_NAME='quitbutton')

  zoomMenu    = widget_button( menuBar, VALUE='Zoom', /MENU)
  zoomButton  = widget_button( zoomMenu, VALUE='zoom in', UVALUE='ZOOM' )
  resetButton = widget_button( zoomMenu, VALUE='reset', UVALUE='RESET' )

  toolMenu = widget_button( menuBar, VALUE='Tools', /MENU )
  xloadctButton = widget_button( toolMenu, VALUE='XLoadct', $
                                 UVALUE='XLOADCT' )
  xpaletteButton = widget_button( toolMenu, VALUE='XPalette', $
                                 UVALUE='XPALETTE' )

  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu   = widget_button( menuBar, VALUE='Help', /MENU, /HELP )
  infoButton = widget_button( helpMenu, VALUE='XViewTermDiagram', $
                              UVALUE='INFORMATION' )

  ;; --- Draw frame base widget --                      ------------- ;;

  state.drawWidget = widget_draw( termdiagBase, XSIZE=xsize, YSIZE=ysize, $
                                  /BUTTON_EVENTS, /MOTION_EVENTS, $
                                  UVALUE='BUTTON_EVENT' )

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- termWidgetSetup.pro ------ ;

; -------- begin --------------------------             -------------- ;

PRO XViewTermDiag, atom_ptr, GROUP_LEADER=group_leader

;+
; NAME:
;	XVIEWTERMDIAG
;
; PURPOSE:
;	This routine displays a Grottrian termdiagram associated to the
;       atomic model pointed to by Atom_ptr
;
; CATEGORY:
;	Widgets.
;
; CALLING SEQUENCE:
;	XVIEWTERMDIAG, Atom_ptr, GROUP_LEADER=group_leader
;
; INPUTS:
;	Atom_ptr:  pointer to heap variable containing atomic structure.
;
; KEYWORD PARAMETERS:
;	GROUP_LEADER:  ID of parent widget in hierarchy.
;
;
; EXAMPLE:
;	XVIEWTERMDIAG, ptr_new(rawatom('CaII.atom', /VACUUM_TO_AIR))
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Jun  9 16:52:34 2010 --
;-

  IF (NOT ptr_valid(atom_ptr)) THEN return
  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader = 0

  state = termWidgetSetup(atom_ptr)

  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  drawTermDiag, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewTermDiag', state.baseWidget, $
   EVENT_HANDLER='XViewTerm_Event', GROUP_LEADER=group_leader

END
; -------- end ----------------------------             -------------- ;
