; -------- begin -------------------------- readContrib.pro ---------- ;

FUNCTION readContrib, fileName, LAMBDA=lambda

@geometry.common

  CASE geometryType OF
    "ONE_D_PLANE": $
     contr_str = {label: "", type: 0L, opac: dblarr(geometry.Ndep)}

    "TWO_D_PLANE": $
     contr_str = {label: "", type: 0L, $
                  opac: dblarr(geometry.Nx, geometry.Nz)}

    "THREE_D_PLANE": $
     contr_str = {label: "", type: 0L, $
                  opac: dblarr(geometry.Nx, geometry.Ny, geometry.Nz)}

    "SPHERICAL_SYMMETRIC": $
     contr_str = {label: "", type: 0L, opac: dblarr(geometry.Nradius)}
  ENDCASE

  openr, lun, fileName, /GET_LUN, /XDR
  lambda = 0.0D+0
  readu, lun, lambda
  Ncontr = 0L
  WHILE (NOT EOF(lun)) DO BEGIN
    readu, lun, contr_str
    Ncontr = Ncontr + 1
  ENDWHILE
  point_lun, lun, 0
  contrib = replicate(contr_str, Ncontr)
  readu, lun, lambda, contrib
  free_lun, lun

  return, contrib
END
; -------- end ---------------------------- readContrib.pro ---------- ;

; -------- begin -------------------------- XViewBopac_Event.pro ----- ;

PRO XViewBopac_Event, Event

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT': widget_control, Event.top, /DESTROY

    'PRINT': BEGIN
      filename = '/tmp/viewBopac-' + timeStamp() + '.ps'
      psopen, FILENAME=filename, /COLOR
      displayBopac, state
      psclose
      r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
    END

    'PNG': BEGIN
      widget_control, state.drawWidget, GET_VALUE=WindowNo
      rhwritepng, WindowNo, '/tmp/viewBopac-'
    END

    'SETTRESHOLD': BEGIN
      widget_control, state.tresholdField, GET_VALUE=treshold
      state.treshold = (0.01 > treshold < 0.99)
      widget_control, state.tresholdField, SET_VALUE=state.treshold
      widget_control, stash, SET_UVALUE=state
      displayBopac, state
    END

    'SLIDER': BEGIN
      displayBopac, state
    END

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["Display of background opacity contributions", $
             "", $
             "Version 1.0, Nov 6, 1998", $
             "Han Uitenbroek (HUitenbroek@cfa.harvard.edu)"])
  ELSE:
  ENDCASE
END
; -------- end ---------------------------- XViewBopac_Event.pro ----- ;

; -------- begin -------------------------- displayBopac.pro --------- ;

PRO displayBopac, state

@geometry.common

  MEGAMETER_TO_M = 1.0E6

  IF (!D.NAME EQ 'X') THEN BEGIN
    widget_control, state.drawWidget, GET_VALUE=WindowNo
    wset, WindowNo
  ENDIF

  c = readContrib(state.file, LAMBDA=lambda)

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      z = geometry.cmass
      xlog = 1
      xtitle = 'column mass [kg m!U-2!N]'
      opac = c.opac
    END

    "TWO_D_PLANE": BEGIN
      z = geometry.z / MEGAMETER_TO_M
      xlog = 0
      xtitle = 'Height [Mm]'
      widget_control, state.slider, GET_VALUE=x_index
      widget_control, state.xlabel, $
       SET_VALUE=string(geometry.x[x_index]/MEGAMETER_TO_M, FORMAT='(F7.3)')
      opac = cmprss(c.opac[x_index, *])
    END
  ENDCASE
  opactotal = total(opac, 2)

  IF (!D.NAME EQ 'PS') THEN BEGIN
    title = string(FORMAT='("display treshold: ", I3, "%")', $
                   fix(state.treshold*100))
    IF (geometryType EQ "TWO_D_PLANE") THEN $
     title = title + string(FORMAT='(", x = ", F7.2, " [Mm]")', $
                            geometry.x[x_index]/MEGAMETER_TO_M)
    thick = 4
  ENDIF ELSE BEGIN
    thick = 2
    title = ''
  ENDELSE

  Nz = n_elements(z)
  plot, z, dblarr(Nz) + 1.0, XLOG=xlog, $
   YRANGE=[-0.02, 1.1], YSTYLE=1, /NODATA, TITLE=title, $
   XTITLE=xtitle, YTITLE='Relative background opacity'
  IF (xlog) THEN $
   oplot, 10^!X.CRANGE, [1.0, 1.0], /LINE $
  ELSE $
   oplot, !X.CRANGE, [1.0, 1.0], /LINE

  color = 32B
  IF (!D.NAME EQ 'PS') THEN charsize = 0.8 ELSE charsize = 1.0
  FOR n=0, n_elements(c)-1 DO BEGIN
    linetype = (c[n].type) ? 2 : 0
    relopac = opac[*, n]/opactotal

    maxopac = max(relopac, max_position)
    IF (maxopac GE state.treshold) THEN BEGIN
      oplot, z, relopac, LINE=linetype, COLOR=color, THICK=thick
      align = $
       convert_coord([z[max_position], 0.0, 0.0], /TO_NORMAL)
      xyouts, z[max_position], maxopac+0.01, $
       c[n].label, CHARSIZE=charsize, COLOR=color, $
       ALIGNMENT=align[0] / (!x.window[1] - !x.window[0])
      color = color + 32B
    ENDIF 
  ENDFOR
  IF (!D.NAME EQ 'PS') THEN $
   rhannotate, xann(0.75), 1.06, $
   TEXT=string(FORMAT='("!7l = !x", F9.3, " [nm]")', lambda), $
   CHARSIZE=0.8, CHARCOLOR=16B
END

; -------- end ---------------------------- displayBopac.pro --------- ;

; -------- begin -------------------------- bopacWidgetSetup.pro ----- ;

FUNCTION bopacWidgetSetup, fileName

@geometry.common

  state = {baseWidget: 0L, drawWidget: 0L, log: 0, file: fileName, $
           treshold: 0.1, tresholdField: 0L, slider: 0L, xlabel: 0L}

  openr, lun, fileName, /XDR, /GET_LUN
  lambda = 0.0D+0
  readu, lun, lambda
  free_lun, lun

  state.baseWidget = widget_base(TITLE='XViewBopac', /COLUMN, $
				 RESOURCE_NAME='XViewBopac', $
				 MBAR=menuBar)

  base = widget_base(state.baseWidget, /COLUMN)

  fileMenu = widget_button(menuBar, VALUE='File', /MENU)
  button = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                         RESOURCE_NAME='quitbutton')


  printMenu= widget_button(menuBar, VALUE='Print', /MENU)
  printButton = widget_button(printMenu, VALUE='PostScript', UVALUE='PRINT')
  pngButton   = widget_button(printMenu, VALUE='png', UVALUE='PNG')

  helpMenu = widget_button(menuBar, VALUE='Help', /HELP)
  button = widget_button(helpMenu, UVALUE='INFORMATION', $
                         VALUE='about XViewBopac')

  frame = widget_base(base, /ROW, /FRAME)
  label = widget_label(frame, VALUE='wavelength: ')
  label = widget_label(frame, /FRAME, VALUE=string(lambda, FORMAT='(F9.3)'))
  label = widget_label(frame, VALUE='[nm]')

  ;; --- Draw frame --                                  ------------- ;;

  drawFrame = widget_base(base, /FRAME, /COLUMN)
  label = widget_label(drawFrame, VALUE="Background opacity contributions", $
                       /ALIGN_CENTER)
  state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=450)

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    state.Slider = widget_slider(base, TITLE='x index:', /ALIGN_CENTER, $
                                 UVALUE='SLIDER', MAXIMUM=geometry.Nx - 1, $
                                 VALUE=geometry.Nx/2, XSIZE=300)
    label = widget_label(frame, VALUE='   x: ')
    state.xlabel = widget_label(frame, /FRAME, $
                                VALUE=string(geometry.x[geometry.Nx/2]/1.0E6, $
                                             FORMAT='(F7.3)'))
    label = widget_label(frame, VALUE='[Mm]')
  ENDIF

  state.tresholdField = cw_field(frame, TITLE='      Treshold: ', /FLOATING, $
                                 UVALUE='SETTRESHOLD', VALUE=state.treshold, $
                                 XSIZE=12, /RETURN_EVENT)

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
end
; -------- end ---------------------------- avgWidgetSetup.pro ------- ;


PRO XViewBopac, fileName, GROUP_LEADER=group_leader

  IF (NOT keyword_set(GROUP_LEADER)) THEN group_leader = 0

  state = bopacWidgetSetup(fileName)
  widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=group_leader
  displayBopac, state

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, 'XViewBopac', state.baseWidget, $
   EVENT_HANDLER='XViewBopac_Event', GROUP_LEADER=group_leader
END
