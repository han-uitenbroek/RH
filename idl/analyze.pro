; -------- file: -------------------------- analyze.pro -------------- ;

; -------- begin -------------------------- metalMenu_Event_func.pro - ;

FUNCTION metalMenu_Event_Func, event

  widget_control, event.id, GET_UVALUE=uvalue

  return, {ID:event.handler, TOP:event.top, HANDLER:0L, VALUE:uvalue}
END
; -------- end ---------------------------- metalMenu_Event_func.pro - ;

; -------- begin -------------------------- moleculeMenu_Event_func.pro;

FUNCTION moleculeMenu_Event_Func, event

  widget_control, event.id, GET_UVALUE=uvalue

  return, {ID:event.handler, TOP:event.top, HANDLER:0L, VALUE:uvalue}
END
; -------- end ---------------------------- moleculeMenu_Event_func.pro;

; -------- begin -------------------------- moleculeMenu_Event_func.pro;

FUNCTION contributions_Event_Func, event

  widget_control, event.id, GET_UVALUE=uvalue

  return, {ID:event.handler, TOP:event.top, HANDLER:0L, VALUE:uvalue}
END
; -------- end ---------------------------- moleculeMenu_Event_func.pro;

; -------- begin -------------------------- atoms_Event_func.pro -- --- ;

FUNCTION atoms_Event_Func, event

  widget_control, event.id, GET_UVALUE=uvalue

  return, {ID:event.handler, TOP:event.top, HANDLER:0L, VALUE:uvalue}
END
; -------- end ---------------------------- atoms_Event_func.pro -- --- ;

; -------- begin -------------------------- putLogo.pro -------------- ;

PRO putLogo, state, geometryType

  rh_IDL_path = getenv('RH_IDL_PATH')
  IF (strlen(rh_IDL_path) GT 0) THEN BEGIN
    widget_control, state.logoWidget, GET_VALUE=WindowNo
    wset, WindowNo

    CASE geometryType OF
      "ONE_D_PLANE":         logo_name = '1-D'
      "TWO_D_PLANE":         logo_name = '2-D'
      "THREE_D_PLANE":       logo_name = '3-D'
      "SPHERICAL_SYMMETRIC": logo_name = 'Sphere'
    ENDCASE

    xyouts, /DEVICE, 25, 50, 'RH', FONT=1, COLOR=100, CHARSIZE=6.0
    xyouts, /DEVICE, 115, 40, 'Analyze', FONT=1, COLOR=100, CHARSIZE=6.0
    xyouts, /DEVICE, 315, 20, logo_name, FONT=1, COLOR=120, CHARSIZE=6.0
  ENDIF ELSE $
   message, "Could not find logo. Check environment RH_IDL_PATH", $
   /INFORMATIONAL

END
; -------- end ---------------------------- putLogo.pro -------------- ;

; -------- begin -------------------------- Analyze_Event.pro -------- ;

PRO Analyze_Event, Event

@files.common
@atmos.common
@geometry.common

  Off = 0  &  On = 1

  ;; --- Main event handler --                          -------------- ;

  stash = widget_info(Event.handler, /CHILD)
  widget_control, stash, GET_UVALUE=state

  widget_control, Event.id, GET_UVALUE=Action
  CASE Action OF
    'QUIT':  BEGIN
      FOR kr=H.Nline, H.Nline+H.Ncont-1 DO BEGIN
        ptr_free, H.transition[kr].lambda_ptr
        ptr_free, H.transition[kr].alpha_ptr
      ENDFOR

      IF (ptr_valid(metals[0])) THEN $
       FOR n=0, n_elements(metals)-1 DO free_atom, metals[n]

      free_lun, Junit
      widget_control, Event.top, /DESTROY
    END

    'ATMOSPHERE':   XViewAtmos, atmosFile, GROUP_LEADER=Event.top
    'GRID':         XViewGrid, GROUP_LEADER=Event.top
    'GEOMETRY':     XViewAngles, GROUP_LEADER=Event.top
    'ABUNDANCE':    XViewAbund, GROUP_LEADER=Event.top
    'HYDROGEN':     XViewAtom, H, GROUP_LEADER=Event.top
    'METALMENU':    XViewAtom, *metals[Event.value], GROUP_LEADER=Event.top
    'CONTRIBUTIONS': XViewBopac, Event.value, GROUP_LEADER=Event.top
    'MOLECULEMENU': XViewMolecule, Event.value, GROUP_LEADER=Event.top
    'BACKGR_EDGES': XDisplayFile, TEXT=getEdges(/TEXT), GROUP=Event.top, $
                     FONT='7x13bold', WIDTH=35, $
                     TITLE='Edges due to background atoms'
    'BACKGR_LINES': XDisplayFile, TEXT=getBackgrLines(/TEXT), $
                     GROUP=Event.top, FONT='7x13bold', WIDTH=60, $
                     TITLE='Lines from background atoms'
    'ATOM':         XViewAtom, readAtom(Event.value), GROUP_LEADER=Event.top

    'EMERGENT':     XViewIe, GROUP_LEADER=Event.top
    'MEAN_INTENS':  XViewJ, GROUP_LEADER=Event.top
    'SOURCE_FNCT':  XViewSource, GROUP_LEADER=Event.top
    'FLUX':         XViewFlux, GROUP_LEADER=Event.top
    'SPATIALAVG':   XViewAvg, GROUP_LEADER=Event.top
    'DISKIMAGE':    XViewDisk, GROUP_LEADER=Event.top

    'VIEWRAY':      XViewRay, atmos.T, 0, geometry.Nx/2, $
     ZTITLE='Temperature [K]', GROUP_LEADER=Event.top

    'KPNO_ATLAS':   XViewAtlas, GROUP_LEADER=Event.top
    'HAWAII_ATLAS': XViewAtlas, GROUP_LEADER=Event.top, /HAWAII
    'KPK_ATLAS':    XViewAtlas, GROUP_LEADER=Event.top, /KPK
    'SUMER':        XViewAtlas, GROUP_LEADER=Event.top, /SUMER
    'ATMOS':        XViewAtlas, GROUP_LEADER=Event.top, /ATMOS
    'SOLFLUX':      XViewAtlas, GROUP_LEADER=Event.top, /SOLFLUX
    'LABS&NECKEL':  result = dialog_message('Not yet implemented', /INFO)

    'XLOADCT':     XLoadct, GROUP=Event.top
    'XPALETTE':    XPalette, GROUP=Event.top

    'SHOWTIMES':   XDisplayFile, 'time.out', GROUP=Event.top, WIDTH=60, $
     FONT='-schumacher-clean-bold-r-normal--13-130-72-*-*-80-iso8859-1'

    'INFORMATION': result = dialog_message(/INFORMATION, $
            ["This is an IDL-widget program to analyze the results from", $
             "calculations with the family of radiative transfer codes", $
             "based on the Rybicki&Hummer ALI method (see A&A 245, 171", $
             "and A&A 262, 209)", $
             " ", $
             "Version 1.0, Dec 7, 1995", $
             "Han Uitenbroek (huitenbroek@cfa.harvard.edu)"])
    'IDLHELP':  spawn, /NOSHELL, 'idlhelp'  
    ELSE:
  ENDCASE

END
; -------- end ---------------------------- Analyze_Event.pro -------- ;
; -------- begin -------------------------- analyzeWidgetSetup.pro --- ;

FUNCTION analyzeWidgetSetup, geometryType, metals, molecules

@files.common

  state = {baseWidget: 0L, logoWidget: 0L, directory: ''}

  CASE geometryType OF
    "ONE_D_PLANE":         title = 'Analyze_1D'
    "TWO_D_PLANE":         title = 'Analyze_2D'
    "THREE_D_PLANE":       title = 'Analyze_3D'
    "SPHERICAL_SYMMETRIC": title = 'Analyze_Sphere'
  ENDCASE

  state.baseWidget = widget_base(TITLE=title, MBAR=menuBar)

  analyzeBase = widget_base(state.baseWidget, /COLUMN)
  state.logoWidget = widget_draw(analyzeBase, XSIZE=450, YSIZE=100, /FRAME)

  fileMenu   = widget_button(menuBar, VALUE='File', /MENU)
  quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

  helpMenu   = widget_button(menuBar, VALUE='Help', /MENU, /HELP)
  infoButton = widget_button(helpMenu, VALUE='about Analyze', $
                             UVALUE='INFORMATION')
  IDLbutton  = widget_button(helpMenu, VALUE='IDL Help', $
                             UVALUE='IDLHELP')

  atmosMenu   = widget_button(menuBar, VALUE='Atmosphere', /MENU)
  atmosButton = widget_button(atmosMenu, VALUE='model', $
                              UVALUE='ATMOSPHERE')

  IF (geometryType EQ "TWO_D_PLANE") THEN BEGIN
    gridButton   = widget_button(atmosMenu, VALUE='spatial grid', $
                                  UVALUE='GRID')
  ENDIF
  anglesButton = widget_button(atmosMenu, VALUE='angular grid', $
                                UVALUE='GEOMETRY')

  backgrMenu = widget_button(atmosMenu, Value='background', /MENU)
  abundButton  = widget_button(backgrMenu, VALUE='abundances', $
                                UVALUE='ABUNDANCE')
  hydrogenButton = widget_button(backgrMenu, VALUE='hydrogen', $
                                 UVALUE='HYDROGEN')
  IF (ptr_valid(metals[0])) THEN BEGIN
    metalMenu = widget_button(backgrMenu, VALUE='metals', $
                               UVALUE='METALMENU', $
                               /MENU, EVENT_FUNC='metalMenu_Event_Func')
    FOR n=0, n_elements(metals)-1 DO BEGIN
     ID = strupcase(strmid((*metals[n]).labels[0], 0, 2))
     metalButton = widget_button(metalMenu, VALUE=ID, UVALUE=n)
    ENDFOR
  ENDIF
  IF (n_tags(molecules) GT 0) THEN BEGIN
    moleculeMenu = widget_button(backgrMenu, VALUE='molecules', $
                               UVALUE='MOLECULEMENU', $
                               /MENU, EVENT_FUNC='moleculeMenu_Event_Func')
    FOR n=0,n_elements(molecules)-1 DO $
     moleculeButton = widget_button(moleculeMenu, $
                                    VALUE=molecules[n].id, UVALUE=n)
  ENDIF

  edgesButton = widget_button(backgrMenu, VALUE='edges', $
                              UVALUE='BACKGR_EDGES')
  linesButton = widget_button(backgrMenu, VALUE='lines', $
                              UVALUE='BACKGR_LINES')


  contribMenu = widget_button(backgrMenu, VALUE='contributions', $
			      UVALUE='CONTRIBUTIONS', $
			      /MENU, EVENT_FUNC='contributions_Event_Func')
  files = file_search('bopac*.???', COUNT=count)
  lambda = float(strmid(files, 5))
  files = files(sort(lambda))
  IF (count EQ 0) THEN $
    widget_control, contribMenu, SENSITIVE=0 $
  ELSE BEGIN
    FOR n=0, count-1 DO $
      button = widget_button(contribMenu, VALUE=strmid(files[n], 5), $
			     UVALUE=files[n])
  ENDELSE
  IF (geometryType EQ "TWO_D_PLANE") THEN $
   viewrayButton = widget_button(atmosMenu, VALUE='view ray', $
                                 UVALUE='VIEWRAY')

  atomMenu   = widget_button(menuBar, VALUE='Atoms', /MENU, $
                              EVENT_FUNC='atoms_Event_Func', UVALUE='ATOM')

  atom_files = file_search("atom.*.out", COUNT=count)
  IF (count EQ 0) THEN $
   widget_control, atomMenu, SENSITIVE=0 $
  ELSE BEGIN
    
    FOR n=0, count-1 DO BEGIN
      atomID = (strsplit(atom_files[n], ".", /EXTRACT))[1]
      atomButton = widget_button(atomMenu, VALUE=atomID, $
                                 UVALUE=atom_files[n])
    ENDFOR
  ENDELSE

  intensityMenu = widget_button(menuBar, VALUE='Intensity', /MENU)
  meanButton    = widget_button(intensityMenu, VALUE='angular mean', $
                                UVALUE='MEAN_INTENS')
  sourceButton  = widget_button(intensityMenu, VALUE='source function', $
                                UVALUE='SOURCE_FNCT')
  emergeButton  = widget_button(intensityMenu, VALUE='emergent', $
                                UVALUE='EMERGENT')

  spectrButton = widget_button(intensityMenu, VALUE='flux', $
                               UVALUE='FLUX')

  If (geometryType EQ "TWO_D_PLANE" OR geometryType EQ "THREE_D_PLANE") THEN $
   spectrButton = widget_button(intensityMenu, VALUE='spatial average', $
                                UVALUE='SPATIALAVG')
  IF (geometryType EQ "SPHERICAL_SYMMETRIC") THEN $
   diskButton = widget_button(intensityMenu, VALUE='disk image', $
                              UVALUE='DISKIMAGE')

  runMenu = widget_button(menuBar, VALUE='Statistics', /MENU)
  timesButton = widget_button(runMenu, VALUE='Show run time', $
                              UVALUE='SHOWTIMES')

  atlasMenu  = widget_button(menuBar, VALUE='Atlas', /MENU)
  IF (strlen(getenv('RH_ATLAS_PATH')) GT 0) THEN BEGIN
    button = widget_button(atlasMenu, VALUE='KPNO atlas (optical)', $
                           UVALUE='KPNO_ATLAS')
    button = widget_button(atlasMenu, VALUE='Hawaii (UV)', $
                           UVALUE='HAWAII_ATLAS')
    button = widget_button(atlasMenu, VALUE='Harvard (UV)', $
                           UVALUE='KPK_ATLAS')
    button = widget_button(atlasMenu, VALUE='SUMER/SOHO', $
                           UVALUE='SUMER')
    button = widget_button(atlasMenu, VALUE='ATMOS (IR)', $
                           UVALUE='ATMOS')
    button = widget_button(atlasMenu, VALUE='Solar Flux', $
                           UVALUE='SOLFLUX')
    button = widget_button(atlasMenu, VALUE='Labs & Neckel', $
                           UVALUE='LABS&NECKEL')
  ENDIF ELSE $
   widget_control, atlasMenu, SENSITIVE=0

  toolMenu = widget_button(menuBar, VALUE='Tools', /MENU)
  xloadctButton = widget_button(toolMenu, VALUE='XLoadct', $
                                UVALUE='XLOADCT')
  xpaletteButton = widget_button(toolMenu, VALUE='XPalette', $
                                 UVALUE='XPALETTE')

  widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
  return, state
END
; -------- end ---------------------------- analyzeWidegtSetup.pro --- ;

; -------- begin -------------------------- analyze.pro -------------- ;

PRO analyze, GEOMETRYFILENAME=geometryFileName, $
             ATMOSFILENAME=atmosFileName, $
             METALFILENAME=metalFileName, $
             MOLECULEFILENAME=moleculeFileName, HELP=help

;+
; NAME:
;	ANALYZE
;
; PURPOSE:
;       Main section of IDL-widget based analysis package.
;
; CATEGORY:
;	Data analysis, Widgets.
;
; CALLING SEQUENCE:
;       ANALYZE
;
; INPUTS:
;
; OPTIONAL INPUTS:
;	
; KEYWORD PARAMETERS:
;       GEOMETRYFILENAME:   File name for geometry data output file if
;                           different from "geometry.out"
;
;       ATMOSFILENAME:      File name for geometry data output file if
;                           different from "atmos.out"
;
;       METALFILENAME:      File name for geometry data output file if
;                           different from "metals.out"
;
;       ATOMFILENAME:       File name for geometry data output file if
;                           different from "atom.out"
;
;       MOLECULEFILENAME:   File name for geometry data output file if
;                           different from "molecules.out"
;
;       HELP:               Show calling sequence.
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;     @geometry.common
;     @atmos.common
;     @spectrum.common
;     @opacity.common
;     @files.common
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Tue Apr 21 16:59:46 2009 --
;-

  ;; --- Define the common blocks --                    -------------- ;

@geometry.common
@atmos.common
@spectrum.common
@opacity.common
@files.common

  ProgramName = 'Analyze'

  IF (xregistered(ProgramName)) THEN return

  IF (!D.NAME NE 'X') THEN BEGIN
    print, 'analyze.pro only works with the X-windows device'
    return
  ENDIF

  IF (fix(strmid(!VERSION.RELEASE, 0, 1)) LT 5) THEN BEGIN
    print, 'analyze.pro requires IDL release 5 or higher'
    return
  ENDIF

  set_decomposed, COLOR_TABLE=15

  !P.COLOR=0B
  !P.BACKGROUND=255B
  !P.FONT=-1

  IF (keyword_set(HELP)) THEN BEGIN
    print, "analyze, GEOMETRYFILENAME=geometryFileName, "
    print, "         ATMOSFILENAME=atmosFileName, "
    print, "         METALFILENAME=metalFileName, "
    print, "         MOLECULEFILENAME=moleculeFileName, HELP=help"
  ENDIF

  success = readInput('input.out')

  ;; --- Load the default files for geometry, atmosphere, and atom --- ;

  IF (NOT keyword_set(GEOMETRYFILENAME)) THEN $
   geometryFile = 'geometry.out' $
  ELSE $
   geometryFile = geometryFileName
  IF (NOT readGeometry(geometryFile)) THEN return

  IF (NOT keyword_set(METALFILENAME)) THEN $
   metalFile = 'metals.out' $
  ELSE $
   metalFile = metalFileName
  IF (NOT keyword_set(MOLECULEFILENAME)) THEN $
   moleculeFile = 'molecules.out' $
  ELSE $
   moleculeFile = moleculeFileName
  IF (NOT keyword_set(ATMOSFILENAME)) THEN $
   atmosFile = 'atmos.out' $
  ELSE $
   atmosFile = atmosFileName
  IF (NOT readAtmos(atmosFile)) THEN return

  spectrumFile  = 'spectrum.out'
  r = readSpectrum(spectrumFile)

  ;; --- Other default file names --                    -------------- ;

  JFile = 'J.dat'
  IF (NOT openJ(JFile)) THEN return
  opacFile = 'opacity.out'
  backgroundFile = 'background.dat'
  lambdaDisplay = 0

  ;; Set up default value of T3D projection matrix

  surfr, AX=30, AZ=30

  state = analyzeWidgetSetup(geometryType, metals, molecules)
  widget_control, state.baseWidget, /REALIZE
  putLogo, state, geometryType

  ;; --- Register with the XManager --                 -------------- ;;

  xmanager, ProgramName, state.baseWidget, EVENT_HANDLER='Analyze_Event', $
   /NO_BLOCK
END
; -------- end ---------------------------- analyze.pro -------------- ;
