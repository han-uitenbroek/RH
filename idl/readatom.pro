; -------- file: -------------------------- readatom.pro ------------- ;

; -------- begin -------------------------- readPops.pro ------------- ;

FUNCTION readPops, atom, fileName

@geometry.common
@files.common

  IF (n_tags(atom) LE 0) THEN return, 0

  openr, popsUnit, fileName, /XDR, /GET_LUN

  ;; --- Store name in common block Files_Common --     -------------- ;

  CASE geometryType OF
    "ONE_D_PLANE":         nstar = dblarr(geometry.Ndep, atom.Nlevel)
    "SPHERICAL_SYMMETRIC": nstar = dblarr(geometry.Nradius, atom.Nlevel)
    "TWO_D_PLANE":         nstar = dblarr(geometry.Nx, geometry.Nz, $
                                          atom.Nlevel)
    "THREE_D_PLANE":       nstar = dblarr(geometry.Nx, geometry.Ny, $
                                          geometry.Nz, atom.Nlevel)
  ENDCASE
  n = nstar

  atmosID = ''
  Nlevel = 0L  &  Nspace = 0L
  readu, popsUnit, atmosID, Nlevel, Nspace, n, nstar
  IF (ptr_valid(atom.nstar_ptr)) THEN  ptr_free, atom.nstar_ptr
  atom.nstar_ptr = ptr_new(nstar)

  IF (ptr_valid(atom.n_ptr)) THEN  ptr_free, atom.n_ptr
  atom.n_ptr = ptr_new(n)

  free_lun, popsUnit
  return, 1
END
; -------- end ---------------------------- readPops.pro ------------- ;

; -------- begin -------------------------- readRates.pro ------------ ;

FUNCTION readRates, atom, fileName

@geometry.common

  IF (n_tags(atom) LE 0) THEN return, 0

  openr, ratesUnit, fileName, /XDR, /GET_LUN

  CASE geometryType OF
    "ONE_D_PLANE":         Rij = dblarr(geometry.Ndep)
    "SPHERICAL_SYMMETRIC": Rij = dblarr(geometry.Nradius)
    "TWO_D_PLANE":         Rij = dblarr(geometry.Nx, geometry.Nz)
    "THREE_D_PLANE":       Rij = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)
  ENDCASE
  Rji = Rij

  FOR kr=0, (atom.Nline + atom.Ncont)-1 DO BEGIN
    readu, ratesUnit, Rij, Rji
    atom.transition[kr].Rij_ptr = ptr_new(Rij)
    atom.transition[kr].Rji_ptr = ptr_new(Rji)
  ENDFOR

  free_lun, ratesUnit
  return, 1
END
; -------- end ---------------------------- readRates.pro ------------ ;

; -------- begin -------------------------- readCollisions.pro ------- ;

FUNCTION readCollisions, atom, fileName

@geometry.common

  ;; --- The following convention is used for the collision matrix:
  ;;
  ;;     C[k, i, j] is the collision rate of transition i -> j at depth k
  ;; 

  IF (n_tags(atom) LE 0) THEN return, 0

  openr, ratesUnit, fileName, /XDR, /GET_LUN

  CASE geometryType OF
    "ONE_D_PLANE": $
     Cij = dblarr(geometry.Ndep, atom.Nlevel, atom.Nlevel)
    "SPHERICAL_SYMMETRIC": $
     Cij = dblarr(geometry.Nradius, atom.Nlevel, atom.Nlevel)
    "TWO_D_PLANE": $
     Cij = dblarr(geometry.Nx, geometry.Nz, atom.Nlevel, atom.Nlevel)
    "THREE_D_PLANE":$
     Cij = dblarr(geometry.Nx, geometry.Ny, geometry.Nz, $
                  atom.Nlevel, atom.Nlevel)
  ENDCASE

  readu, ratesUnit, Cij
  atom.Cij_ptr = ptr_new(Cij)

  free_lun, ratesUnit
  return, 1
END
; -------- end ---------------------------- readCollisions.pro ------- ;

; -------- begin -------------------------- readDamping.pro ---------- ;

FUNCTION readDamping, atom, fileName

@geometry.common

  IF (n_tags(atom) LE 0) THEN return, 0

  openr, ratesUnit, fileName, /XDR, /GET_LUN

  CASE geometryType OF
    "ONE_D_PLANE":         Adamp = dblarr(geometry.Ndep)
    "SPHERICAL_SYMMETRIC": Adamp = dblarr(geometry.Nradius)
    "TWO_D_PLANE":         Adamp = dblarr(geometry.Nx, geometry.Nz)
    "THREE_D_PLANE":$
     Adamp = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)
  ENDCASE
  vbroad = Adamp

  readu, ratesUnit, vbroad
  atom.vbroad_ptr = ptr_new(vbroad)
  FOR kr=0, atom.Nline-1 DO BEGIN
    readu, ratesUnit, Adamp
    atom.transition[kr].Adamp_ptr = ptr_new(Adamp)
  ENDFOR

  free_lun, ratesUnit
  return, 1
END
; -------- end ---------------------------- readDamping.pro ---------- ;

; -------- begin -------------------------- read_xdr_atom.pro -------- ;

FUNCTION read_xdr_atom, unit, ESSENTIALS=essentials

;+
; NAME:
;	READ_XDR_ATOM
;
; PURPOSE:
;       This routine reads atomic output data in RH format from the specified
;       logical unit.
;
; CATEGORY:
;       I/O
;
; CALLING SEQUENCE:
;       Result = READ_XDR_ATOM( Unit )
;
; INPUTS:
;	Unit:   Logical unit for output file.
;
; OUTPUTS:
;	The function returns a structure containing the requested atomic model.
;
; SIDE EFFECTS:
;       The returned atomic structure contains pointers to heap variables.
;
; RESTRICTIONS:
;	The logical unit has to be opened with the /XDR keyword.
;
; EXAMPLE:
;       Atom = READ_XDR_ATOM( Unit )
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Apr 22 10:43:55 2009 --
;-

  point_lun, -unit, startPosition

  HYDROGENIC = 3
  EXPLICIT   = 4

  active = 0  &  Nlevel = 0L  &  Nline = 0L  &  Ncont = 0L  &  Nfixed = 0L
  readu, unit, active, Nlevel, Nline, Ncont, Nfixed

  ;; --- Structure definitions --                       -------------- ;

  Nrad = Nline + Ncont
  IF (Nrad GT 0) THEN BEGIN
    trans = replicate({type: 0L,  i: 0L,  j: 0L, Nlambda: 0L, Nblue: 0L, $
                       lambda0: 0.0D+0,  shape: 0L,  strength: 0.0D+0, $
                       lambda_ptr: ptr_new(), alpha_ptr: ptr_new(), $
                       Rij_ptr: ptr_new(), Rji_ptr: ptr_new(), $
                       Adamp_ptr: ptr_new()}, Nrad)
  ENDIF ELSE $
   trans = 0

  IF (Nfixed GT 0) THEN BEGIN
    fixed = replicate({type: 0L,  option: 0L,  i: 0L,  j: 0L, $
                       lambda0: 0.0D+0, strength: 0.0D+0,  Trad: 0.0D+0}, $
                      Nfixed)
  ENDIF ELSE $
   fixed = 0

  atom = {active: active, $
          Nlevel: Nlevel, Nline: Nline, Ncont: Ncont, Nfixed: Nfixed, $
          abundance: 0.0D+0,  weight: 0.0D+0, $
          labels: replicate(string(FORMAT='(A20)', ''), Nlevel),  $
          g: dblarr(Nlevel), E: dblarr(Nlevel), stage: lonarr(Nlevel), $
          transition: trans, fixed: fixed, $
          n_ptr: ptr_new(), nstar_ptr: ptr_new(), Cij_ptr: ptr_new(), $
          vbroad_ptr: ptr_new()}

  ;; --- Do the actual reading here --                  -------------- ;

  point_lun, unit, startPosition
  readstructure, atom, unit, N_MAX_TAGS=11
  trans = atom.transition[0]
  FOR kr=0, Nrad-1 DO BEGIN
    readstructure, trans, unit, N_MAX_TAGS=8
    atom.transition[kr] = trans
  ENDFOR

  ;; --- Read the bound-free wavelengths and cross-sections -- ------- ;

  FOR kr=Nline, Nrad-1 DO BEGIN
    IF (atom.transition[kr].shape EQ HYDROGENIC) THEN BEGIN
      lambda = [0.0D+0]  &  alpha = lambda
      readu, unit, lambda
    ENDIF ELSE BEGIN
      lambda = dblarr(atom.transition[kr].Nlambda)  &  alpha = lambda
      readu, unit, lambda, alpha
    ENDELSE
    atom.transition[kr].lambda_ptr = ptr_new(lambda)
    atom.transition[kr].alpha_ptr  = ptr_new(alpha)
  ENDFOR

  ;; --- Read the fixed transitions --                  -------------- ;

  IF (Nfixed GT 0) THEN BEGIN 
    fixed_tmp = atom.fixed[0]
    FOR kf=0, Nfixed-1 DO BEGIN
      readstructure, fixed_tmp, unit
      atom.fixed[kf] = fixed_tmp
    ENDFOR
  ENDIF

  IF (NOT keyword_set(ESSENTIALS)) THEN BEGIN

    atomid = strtrim(strupcase(strmid(atom.labels[0], 0, 2)), 2)

    ;; --- Read the population numbers if available --  -------------- ;

    popsFile = string(FORMAT='("pops.", A, ".out")', atomid)
    r = file_search(popsFile, COUNT=count)
    IF (count GT 0) THEN r = readPops(atom, popsFile)

    ;; --- Read the radiative rates if available --     -------------- ;

    ratesfile = string(FORMAT='("radrate.", A, ".out")', atomid)
    r = file_search(ratesfile, COUNT=count)
    IF (count GT 0) THEN r = readRates(atom, ratesfile)
    
    ;; --- Read the radiative rates if available --     -------------- ;

    ratesfile = string(FORMAT='("radrate.", A, ".out")', atomid)
    r = file_search(ratesfile, COUNT=count)
    IF (count GT 0) THEN r = readRates(atom, ratesfile)

    ;; --- Read the collisional rates if available --   -------------- ;

    ratesfile = string(FORMAT='("collrate.", A, ".out")', atomid)
    r = file_search(ratesfile, COUNT=count)
    IF (count GT 0) THEN r = readCollisions(atom, ratesfile)

    ;; --- Read the broadening velocity and damping paramaters
    ;;     if available --                              -------------- ;

    ratesfile = string(FORMAT='("damping.", A, ".out")', atomid)
    r = file_search(ratesfile, COUNT=count)
    IF (count GT 0) THEN r = readDamping(atom, ratesfile)
  ENDIF

  return, atom
END
; -------- end ---------------------------- read_xdr_atom.c ---------- ;

; -------- begin -------------------------- readAtom.pro ------------- ;

FUNCTION readAtom, fileName, ESSENTIALS=essentials

;+
; NAME:
;	READATOM
;
; PURPOSE:
;	This function reads an atomic data output file in the RH format
;
; CATEGORY:
;	I/O
;
; CALLING SEQUENCE:
;       Result = READATOM( Filename )
;
; INPUTS:
;	Filename: File name containing atomic output data.
;
; OUTPUTS:
;	The function returns a structure containing the requested atomic model.
;
; COMMON BLOCKS:
;	files.common
;
; SIDE EFFECTS:
;	The returned atomic structure contains pointers to heap variables.
;
; PROCEDURE:
;	This procedure uses READ_XDR_ATOM to do the actual reading
;
; EXAMPLE:
;       Atom = READTOM( 'atom.out' )
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Wed Oct 29 16:15:13 1997 --
;-


@files.common

  WHILE (NOT existFile( fileName, UNIT=atomUnit, /XDR )) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new atom file?")
   IF (answer EQ 'Yes') THEN BEGIN
    fileName = dialog_pickfile(FILTER='*.out', TITLE='Atomic data file', $
                          /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '' ) THEN return, 0
   ENDIF ELSE $
    return, 0

  ENDWHILE

  ;; --- Store name in common block Files_Common --     -------------- ;

  atom = read_xdr_atom(atomUnit, ESSENTIALS=keyword_set(ESSENTIALS))
  free_lun, atomUnit

  return, atom
END
; -------- end ---------------------------- readAtom.pro ------------- ;

; -------- begin -------------------------- readMetals.pro ----------- ;

PRO readMetals, fileName

@files.common
@atmos.common

  WHILE (NOT existFile( fileName, UNIT=metalUnit, /XDR )) DO BEGIN
   answer = dialog_message(/QUESTION, "Find new metals file?")
   IF (answer EQ 'Yes') THEN BEGIN
     fileName = dialog_pickfile(FILTER='*.out', TITLE='Atomic data file', $
                                /MUST_EXIST, /READ, FILE=fileName)
    IF (fileName EQ '' ) THEN return
   ENDIF ELSE $
    return

  ENDWHILE

  ;; --- Store name in common block Files_Common --     -------------- ;

  metalsFile = fileName
  H = read_xdr_atom(metalUnit)

  Nmetal = 0L;
  readu, metalUnit, Nmetal;
  IF (Nmetal GT 0) THEN BEGIN
    metals = ptrarr(Nmetal)
    FOR n=0, Nmetal-1 DO BEGIN
      metals[n] = ptr_new(read_xdr_atom(metalUnit))
    ENDFOR
  ENDIF ELSE BEGIN
    metals = ptr_new()
  ENDELSE

  free_lun, metalUnit
  return
END
; -------- end ---------------------------- readMetals.pro ----------- ;

; -------- begin -------------------------- getEdges.pro ------------- ;

FUNCTION getEdges, TEXT=text, LABELS=labels

@atmos.common

  IF (NOT H.active) THEN BEGIN
    index  = H.Nline + indgen(H.Ncont)
    edges  = H.transition[index].lambda0
    labels = H.labels[H.transition[index].i]
  ENDIF

  FOR n=0, n_elements(metals)-1 DO BEGIN
    IF ((*metals[n]).Ncont GT 0) THEN BEGIN
      index = (*metals[n]).Nline + indgen((*metals[n]).Ncont)
      IF (n_elements(edges) GT 0) THEN BEGIN
        edges  = [edges, $
                  (*metals[n]).transition[index].lambda0]
        labels = [labels, $
                  (*metals[n]).labels[(*metals[n]).transition[index].i]]
      ENDIF ELSE BEGIN
        edges  = (*metals[n]).transition[index].lambda0
        labels = (*metals[n]).labels[(*metals[n]).transition[index].i]
      ENDELSE
    ENDIF
  ENDFOR

  s = sort(edges)
  edges = edges[s]
  labels = labels[s]

  IF (keyword_set(TEXT)) THEN BEGIN
    printout = strarr(n_elements(edges) + 2)
    printout[0] = 'Wavelength      Label'
    FOR n=0, n_elements(edges)-1 DO $
     printout[n+2] = string(FORMAT='(F10.3, 6X, A20)', $
                            edges[n], labels[n])
    return, printout
  ENDIF ELSE $
   return, edges
END
; -------- end ---------------------------- getEdges.pro ------------- ;

; -------- begin -------------------------- getLines.pro ------------- ;

FUNCTION getBackgrLines, TEXT=text, LABELSI=labelsi, LABELSJ=labelsj

@atmos.common

  IF (NOT H.active) THEN BEGIN
    index   = indgen(H.Nline)
    lines   = H.transition[index].lambda0
    labelsi = H.labels[H.transition[index].i]
    labelsj = H.labels[H.transition[index].j]
  ENDIF

  FOR n=0, n_elements(metals)-1 DO BEGIN
    index   = indgen((*metals[n]).Nline)
    IF (n_elements(lines) GT 0) THEN BEGIN
      lines   = [lines, $
                 (*metals[n]).transition[index].lambda0]
      labelsi = [labelsi, $
                 (*metals[n]).labels[(*metals[n]).transition[index].i]]
      labelsj = [labelsj, $
                 (*metals[n]).labels[(*metals[n]).transition[index].j]]
    ENDIF ELSE BEGIN
      lines   = (*metals[n]).transition[index].lambda0
      labelsi = (*metals[n]).labels[(*metals[n]).transition[index].i]
      labelsj = (*metals[n]).labels[(*metals[n]).transition[index].j]
    ENDELSE
  ENDFOR

  s = sort(lines)
  lines = lines[s]
  labelsi = labelsi[s]  &  labelsj = labelsj[s]

  IF (keyword_set(TEXT)) THEN BEGIN
    printout = strarr(n_elements(lines) + 2)
    printout[0] = 'Wavelength      Upper                    Lower'
    FOR n=0, n_elements(lines)-1 DO $
     printout[n+2] = string(FORMAT='(F10.3, 6X, A20, " --  ", A20)', $
                            lines[n], labelsj[n], labelsi[n])
    return, printout
  ENDIF ELSE $
   return, lines
END
; -------- end ---------------------------- getLines.pro ------------- ;

; -------- begin -------------------------- free_atom.pro ------------ ;

PRO free_atom, atom_ptr
  
  IF (ptr_valid((*atom_ptr).n_ptr))      THEN ptr_free, (*atom_ptr).n_ptr
  IF (ptr_valid((*atom_ptr).nstar_ptr))  THEN ptr_free, (*atom_ptr).nstar_ptr
  IF (ptr_valid((*atom_ptr).Cij_ptr))    THEN ptr_free, (*atom_ptr).Cij_ptr
  IF (ptr_valid((*atom_ptr).vbroad_ptr)) THEN ptr_free, (*atom_ptr).vbroad_ptr

  FOR kr=0, (*atom_ptr).Nline-1 DO BEGIN
    IF (ptr_valid((*atom_ptr).transition[kr].Adamp_ptr)) THEN $
     ptr_free, (*atom_ptr).transition[kr].Adamp_ptr
  ENDFOR

  Nrad = (*atom_ptr).Nline + (*atom_ptr).Ncont
  FOR kr=(*atom_ptr).Nline, Nrad - 1 DO BEGIN
    ptr_free, (*atom_ptr).transition[kr].lambda_ptr
    ptr_free, (*atom_ptr).transition[kr].alpha_ptr
  ENDFOR

  IF ((*atom_ptr).weight GT 0.0) THEN BEGIN
    FOR kr=0, Nrad - 1 DO BEGIN
      IF (ptr_valid((*atom_ptr).transition[kr].Rij_ptr)) THEN BEGIN
        ptr_free, (*atom_ptr).transition[kr].Rij_ptr
        ptr_free, (*atom_ptr).transition[kr].Rji_ptr
      ENDIF
    ENDFOR
  ENDIF

  ptr_free, atom_ptr
END
; -------- end ---------------------------- free_atom.pro ------------ ;
