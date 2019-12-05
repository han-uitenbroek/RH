FUNCTION rawAtom, fileName, VACUUM_TO_AIR=vacuum_to_air

;+
; NAME:
;	RAWATOM
;
; PURPOSE:
;	This function reads an atomic data input file from the RH distribution.
;
; CATEGORY:
;	I/O
;
; CALLING SEQUENCE:
;       Result = RAWATOM(Filename)
;
; INPUTS:
;	Filename: File name containing atomic input data.
;
; KEYWORD PARAMETERS:
;       VACUUM_TO_AIR:  If set transition wavelengths are translated to air
;                       
; OUTPUTS:
;	The function returns a structure containing the requested atomic model.
;
; SIDE EFFECTS:
;	The returned atomic structure contains pointers to heap variables.
;
; EXAMPLE:
;		atom = RAWATOM( '~/src/rh/Atoms/CaII.atom' )
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Jan 10 15:03:51 2018 --
;-

  LABEL_WIDTH    = 20
  MAX_LINE_WIDTH = 132

  WHILE (NOT existFile( fileName, UNIT=atomUnit )) DO BEGIN
   answer = dialog_message(/QUESTION, "Find atom file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.dat', TITLE='Atomic data file', $
                               /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '' ) THEN return, 0
   ENDIF ELSE $
    return, 0
  ENDWHILE

  cLight  = 2.99792458E+08
  hPlanck = 6.626176E-34 
  mElect  = 9.109534E-31
  qElect  = 1.60219E-19
  epsil0  = 8.85419E-12

  CM_TO_M = 1.0E-02
  NM_TO_M = 1.0E-09

  atomID = ''
  inputLine = getline(atomUnit, /EXIT_ON_EOF)
  reads, inputLine, FORMAT='(A2)', atomID

;;  inputLine = getline(atomUnit, /EXIT_ON_EOF)
;;  reads, inputLine, abundance, weight
;;  abundance = 10^(abundance - 12.0)

  Nlevel = 0L  &  Nline = 0L  &  Ncont = 0L  &  Nfixed = 0L
  inputLine = getline(atomUnit, /EXIT_ON_EOF)
  reads, inputLine, Nlevel, Nline, Ncont, Nfixed

  ;; --- Structure definitions --                       -------------- ;

  Nrad = Nline + Ncont
  IF (Nrad GT 0) THEN BEGIN
    trans = replicate({type: 0L,  i: 0L,  j: 0L, Nlambda: 0L, $
                       blue: 0L,  red: 0L, $
                       lambda0: 0.0,  lambdamin: 0.0, $
                       shape: 0L,  strength: 0.0, $
                       lambda_ptr: ptr_new(), alpha_ptr: ptr_new()}, Nrad)
  ENDIF ELSE $
   trans = 0

  IF (Nfixed GT 0) THEN BEGIN
    fixed = replicate({type: 0L,  option: 0L,  i: 0L,  j: 0L, $
                       lambda0: 0.0, strength: 0.0,  Trad: 0.0}, Nfixed)
  ENDIF ELSE $
   fixed = 0

  atom = {Nlevel: Nlevel, Nline: Nline, Ncont: Ncont, Nfixed: Nfixed, $
          abundance: 0.0,  weight: 0.0, active: 0, $
          labels: replicate(string(FORMAT='(A20)', ''), Nlevel),  $
          g: fltarr(Nlevel), E: fltarr(Nlevel), stage: lonarr(Nlevel), $
          transition: trans, fixed: fixed, $
          n_ptr: ptr_new(), nstar_ptr: ptr_new()}

  FOR i=0, Nlevel-1 DO BEGIN
    inputLine = getline(atomUnit, /EXIT_ON_EOF)
    reads, inputLine, E, g
    atom.E(i) = E * (hPlanck * cLight) / CM_TO_M
    atom.g(i) = g

    atom.labels(i) = strmid(inputLine, strpos(inputLine, "'")+1, LABEL_WIDTH)
    atom.stage(i)  = strmid(inputLine, rstrpos(inputLine, "'")+1, $
                            MAX_LINE_WIDTH)
  ENDFOR

  C = 2 * !PI * (qElect/epsil0) * (qElect/mElect) / cLight;

  line = atom.transition(0)
  FOR kr=0, Nline-1 DO BEGIN
    line.type = 1

    inputLine = getline(atomUnit, /EXIT_ON_EOF)
;;    items = str_sep(inputLine, ' ')
    items = strsplit(inputLine, ' ', /EXTRACT)
    items = items(where(items NE ''))

    line.i = long(items(1))  &  line.j = long(items(0))
    line.Nlambda = long(items(4))
    line.lambda0 = (hPlanck * cLight) / $
     (atom.E(line.j) - atom.E(line.i))
    line.strength = C / line.lambda0^2 * $
     (atom.g(line.i) / atom.g(line.j)) * float(items(2))
    line.lambda0 = line.lambda0 / NM_TO_M
    IF (keyword_set(VACUUM_TO_AIR)) THEN $
      line.lambda0 = vacuumtoair(line.lambda0)

    CASE (items(3)) OF
      'GAUSS': line.shape = 0
      'VOIGT': line.shape = 1
      'PRD'  : line.shape = 2
    ENDCASE

    atom.transition(kr) = line
  ENDFOR

  IF (Ncont GT 0) THEN BEGIN
    cont = atom.transition(Nline)
    FOR kr=Nline, Nline+Ncont-1 DO BEGIN
      cont.type = 2

      inputLine = getline(atomUnit, /EXIT_ON_EOF)
      items = strsplit(inputLine, ' ', /EXTRACT)
      items = items(where(items NE ''))

      cont.i = long(items(1))  &  cont.j = long(items(0))
      cont.Nlambda = long(items(3))
      cont.strength = float(items(2))
      cont.lambda0 = (hPlanck * cLight) / $
       (atom.E(cont.j) - atom.E(cont.i)) / NM_TO_M

      CASE (items(4)) OF
        'HYDROGENIC': cont.shape = 3
        'EXPLICIT'  : cont.shape = 4
      ENDCASE
      cont.lambdaMin = float(items(5))

      lambda = fltarr(cont.Nlambda)  &  alpha = lambda
      dummy = fltarr(2)
      IF (cont.shape EQ 4) THEN BEGIN
        FOR la=0, cont.Nlambda-1 DO BEGIN
          inputLine = getline(atomUnit, /EXIT_ON_EOF)
          reads, inputLine, dummy
          lambda[la] = dummy[0]  &  alpha[la] = dummy[1]
        ENDFOR
      ENDIF ELSE BEGIN
        lambda = [cont.lambdaMin]
        alpha  = [cont.strength * (cont.lambda0/lambda[0])^3]
      ENDELSE
      IF (keyword_set(VACUUM_TO_AIR)) THEN BEGIN
        cont.lambda0 = vacuumtoair(cont.lambda0)
        lambda = vacuumtoair(lambda)
      ENDIF 

      cont.lambda_ptr = ptr_new(lambda)
      cont.alpha_ptr  = ptr_new(alpha)

      atom.transition(kr) = cont
    ENDFOR
  ENDIF

  IF (Nfixed GT 0) THEN BEGIN
    fixed = atom.fixed(0)
    FOR kf=0, Nfixed-1 DO BEGIN
      inputLine = getline(atomUnit, /EXIT_ON_EOF)
      items = strsplit(inputLine, ' ', /EXTRACT)
      items = items(where(items NE ''))

      fixed.i = long(items(1))  &  fixed.j = long(items(0))
      fixed.lambda0 = (hPlanck * cLight) / $
       (atom.E(fixed.j) - atom.E(fixed.i)) / NM_TO_M
      IF (keyword_set(VACUUM_TO_AIR)) THEN $
       fixed.lambda0 = vacuumtoair(fixed.lambda0)

      fixed.strength = float(items(2))
      fixed.Trad     = float(items(3))

      CASE (items(4)) OF
        'TRAD_ATMOSPHERIC':   fixed.option = 0
        'TRAD_PHOTOSPHERIC':  fixed.option = 1
        'TRAD_CHROMOSPHERIC': fixed.option = 2
      ENDCASE

      IF (atom.stage(fixed.j) GT atom.stage(fixed.i)) THEN $
       fixed.type = 2 $
      ELSE $
       fixed.type = 1

      atom.fixed(kf) = fixed
    ENDFOR
  ENDIF

  free_lun, atomUnit
  return, atom
END

