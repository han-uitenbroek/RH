FUNCTION readnistlevels, inputfile, baselabel

  IF (n_params() LT 2) THEN BEGIN
    print, "Usage: levels = READNISTLEVELS(inputfile, baselabel)"
    return, 0
  END

  FALSE = 0
  TRUE  = 1

  HPLANCK = 6.6260755E-34
  CLIGHT  = 2.99792458E+08
  CM_TO_M = 1.0E-02

  openr, lun_in, /GET_LUN, inputfile
  inputline = ''

  Nlevel = 0
  WHILE (NOT eof(lun_in)) DO BEGIN
    readf, lun_in, inputline

    words = strsplit(inputline, "|", /EXTRACT)
    IF (strmid(words[0], 0, 6) NE "------"  AND $
        strmid(words[0], 0, 6) NE "Config"  ) THEN BEGIN

      IF (strmid(words[0], 0, 6) EQ "      "  AND $
          strmid(words[2], 0, 4) NE "    ") THEN $
       multiplet = TRUE $
      ELSE $
       multiplet = FALSE

      IF (strmid(words[0], 0, 6) NE "      "  OR multiplet) THEN BEGIN
        IF (Nlevel EQ 0) THEN BEGIN
          levels = [{label: string(FORMAT='(A20)', ""), $
                     g: 0.0D0, $
                     E: 0.0D0}]
        ENDIF ELSE BEGIN
          levels = [levels, levels[0]]
        ENDELSE

        IF (NOT multiplet) THEN BEGIN
          terms = strsplit(words[0], ".", COUNT=count, /EXTRACT)
          termlabel = strtrim(terms[0], 2)
          FOR n=1, count-1 DO termlabel += " " + strtrim(terms[n], 2)
          termlabel = strupcase(termlabel)

          orbit = strtrim(words[1], 2)
          IF (strlen(orbit) EQ 3) THEN $
           parity = strmid(orbit, 0, 2) + "O" $
          ELSE $
           parity = strmid(orbit, 0, 2) + "E"
        ENDIF

        fraction = strsplit(words[2], "/", /EXTRACT, COUNT=count)
        IF (count EQ 2) THEN $
         J = double(fraction[0]) / double(fraction[1]) $
        ELSE $
         J = double(words[2])

        label = baselabel + " " + termlabel + " " + parity + " " + $
                strtrim(fraction[0], 2)

        levels[Nlevel].label = string(label, FORMAT='(A-20)')
        levels[Nlevel].g     = 2.0*J + 1.0
        levels[Nlevel].E     = double(words[3]) * (HPLANCK * CLIGHT) / CM_TO_M

        Nlevel++
      ENDIF
    ENDIF

  ENDWHILE
  free_lun, lun_in

  print, "NREADNISTLEVELS: Read ", n_elements(levels), " levels"
  return, levels
END
