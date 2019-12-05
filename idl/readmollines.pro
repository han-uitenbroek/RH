FUNCTION readKuruczMolecularLines, logicalUnit, Nlines

  molecular_rt = {KURUCZ_MOL_LINE_STR, $
                  wavelength: 0.0, log_gf: 0.0, Ji: 0.0, Ei: 0.0, $
                  Jj: 0.0, Ej: 0.0, code: 0, $
                  configi: 'X', vi: 0, parity_i: 'E1', $
                  configj: 'A', vj: 0, parity_j: 'E1'}

  rt_format = '(F10.4, F7.3, F5.1, F10.3, F5.1, F11.3, I4, ' + $
   'A1, I2, A2, 3X, A1, I2, A2)'

  lines = replicate(molecular_rt, Nlines)
  readf, logicalUnit, lines, FORMAT=rt_format

  return, lines
END

FUNCTION readKurucznewMolecularLines, logicalUnit, Nlines

  molecular_rt = {KURUCZ_MOL_LINE_STR2, $
                  wavelength: 0.0, log_gf: 0.0, Ji: 0.0, Ei: 0.0, $
                  Jj: 0.0, Ej: 0.0, code: 0, $
                  configi: '0X', vi: 0, parity_i: 'E1', $
                  configj: '0A', vj: 0, parity_j: 'E1', isotope: 0}

  rt_format = '(F10.4, F7.3, F5.1, F10.3, F5.1, F11.3, I4, ' + $
   'A1, I2, A2, 3X, A1, I2, A2, 3X, I2)'

  lines = replicate(molecular_rt, Nlines)
  readf, logicalUnit, lines, FORMAT=rt_format

  return, lines
END

FUNCTION readKuruczTiOMolecularLines, logicalUnit, Nlines

  molecular_rt = {KURUCZ_MOL_LINE_STR3, $
                  wavelength: 0.0, log_gf: 0.0, Ji: 0.0, Ei: 0.0, $
                  Jj: 0.0, Ej: 0.0, code: 0, $
                  configi: '0X', vi: 0, parity_i: 'E1', $
                  configj: '0A', vj: 0, parity_j: 'E1', isotope: 0}

  rt_format = '(F10.4, F7.3, F5.1, F10.3, F5.1, F10.3, I3, ' + $
   'A1, I2, A2, 3X, A1, I2, A2, 3X, I2)'

  lines = replicate(molecular_rt, Nlines)
  readf, logicalUnit, lines, FORMAT=rt_format

  return, lines
END

FUNCTION readKuruczCNMolecularLines, logicalUnit, Nlines, ISOTOPE=isotope

  molecular_rt = {KURUCZ_MOL_LINE_STR4, $
                  wavelength: 0.0, log_gf: 0.0, Ji: 0.0, Ei: 0.0, $
                  Jj: 0.0, Ej: 0.0, code: 0, $
                  configi: '0X', vi: 0, parity_i: 'E1', isotope: 0, $
                  configj: '0A', vj: 0, parity_j: 'E1'}

  rt_format = '(F10.4, F7.3, F5.1, F10.3, F5.1, F11.3, 1X, I3, ' + $
   'A1, I2, 1X, A1, 1X, I2, A1, I2, A2)'

  lines = replicate(molecular_rt, Nlines)
  readf, logicalUnit, lines, FORMAT=rt_format

  IF (keyword_set(ISOTOPE)) THEN $
     lines = lines[where(lines.isotope eq isotope)]

  Nnew = n_elements(lines)
  for n=0, Nnew-1 do begin
     if (lines[n].Ji eq lines[n].Jj) then begin
        case (lines[n].parity_j) of
           'E1':  lines[n].parity_i = 'F2'
           'E2':  lines[n].parity_i = 'F1'
           'F1':  lines[n].parity_i = 'E2'
           'F2':  lines[n].parity_i = 'E1'
        endcase
     endif else begin
        if (lines[n].parity_i eq $
            strmid(lines[n].parity_j, 1, 1)) then $
               lines[n].parity_i = lines[n].parity_j else $
                  lines[n].parity_i = strmid(lines[n].parity_j, 0, 1) + $
                                      lines[n].parity_i
     endelse
  endfor
  
  return, lines
END

FUNCTION readGoorvitchCOLines, logicalUnit, Nlines

  molecular_rt = {GOORVITCH_CO_LINE_STR, $
                  wavenumber: 0.0, dipole: 0.0, Aji: 0.0, Ei: 0.0, gf: 0.0, $
                  strength: 0.0, vj: 0, vi: 0, type: ' ', Ji: 0, $
                  isotope: '00'}

  rt_format = '(F9.4, 2E10.3, F11.4, 2E10.3, 2I3, 1X, A1, I4, 1X, A2)'

  lines = replicate(molecular_rt, Nlines)
  readf, logicalUnit, lines, FORMAT=rt_format

  return, lines
END

FUNCTION readmollines, filename, ISOTOPE=isotope

  openr, logicalUnit, filename, /GET_LUN

  Nline = 0L  &  type_and_format = ''
  readf, logicalUnit, Nline, type_and_format
  tokens = strsplit(type_and_format, ' ', /EXTRACT)
  type = tokens[0]
  format = tokens[1]

  Nlambda = 0  &  qwing = 0.0
  readf, logicalUnit, Nlambda, qwing 

  CASE (format) OF
    'GOORVITCH94': lines = readGoorvitchCOLines(logicalUnit, Nline)
    'KURUCZ_CD18': lines = readKuruczMolecularLines(logicalUnit, Nline)
    'KURUCZ_NEW':  lines = readKurucznewMolecularLines(logicalUnit, Nline)
    'KURUCZ_TIO':  lines = readKuruczTiOMolecularLines(logicalUnit, Nline)
    'KURUCZ_CN' :  lines = readKuruczCNMolecularLines(logicalUnit, $
                                                      Nline, ISOTOPE=isotope)

    ELSE: BEGIN
      print, "Unknown molecular line format: ", format
      return, 0
    END
  ENDCASE
  free_lun, logicalUnit

  print, n_elements(lines), type, FORMAT='("Found ", I8, " lines of type ", A)'
  return, lines
END

