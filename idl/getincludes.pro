PRO getincludes, LIBRARY=library, OUTPUT_FILE=output_file, $
                 TEMPLATE=template, INCLUDE_DIR=include_dir

  ;;+
  ;;  List .h include files in c-source code files to define
  ;;  explicit dependencies in Makefile.
  ;;-


  IF (NOT keyword_set(TEMPLATE)) THEN template = '*.c'
  files = file_search(template, COUNT=Nfile)

  IF (keyword_set(OUTPUT_FILE)) THEN $
   openw, lun, output_file, /GET_LUN $
  ELSE $
   lun = -1

  FOR i=0, Nfile-1 DO BEGIN
    source_file = files[i]

    IF (strpos(source_file, ".c") NE strlen(source_file)-2) THEN BEGIN
      print, " Not a C source code file!"
      return
    ENDIF

    IF (keyword_set(INCLUDE_DIR)) THEN $
     spawn, ['cc', '-xM1', source_file, '-I', include_dir], result, /NOSHELL $
    ELSE $
     spawn, ['cc', '-xM1', source_file], result, /NOSHELL
    Nresult = n_elements(result)

    dot_h_files = strarr(Nresult)
    Ninclude = 0
    tokens = str_sep(result[0], ":", /TRIM)
    dot_o_file = tokens[0]

    match = 0
    FOR n=1, Nresult-1 DO BEGIN
      tokens = str_sep(result[n], ":", /TRIM)
      h_file = tokens[1]
      IF (strpos(h_file, ".h") EQ strlen(h_file)-2) THEN BEGIN
        FOR m=0, Ninclude-1 DO $
         IF ((match = (h_file EQ dot_h_files[m]))) THEN GOTO, BREAK
      ENDIF
      
BREAK:
      IF (NOT match) THEN BEGIN
        dot_h_files[Ninclude] = h_file
        Ninclude = Ninclude + 1         
      ENDIF
    ENDFOR

    IF (keyword_set(LIBRARY)) THEN $
      printf, lun, LIBRARY + "(" + dot_o_file + "):", FORMAT='(/A, T25, $)' $
    ELSE $
     printf, lun, dot_o_file + ":", FORMAT='(/A, T25, $)'
    length = 25

    FOR n=0, Ninclude-1 DO BEGIN
      length = length + strlen(dot_h_files[n])
      IF (length GT 65) THEN BEGIN
        fmt = (n EQ Ninclude - 1) ? $
         '("\", /T25, A)' : '("\", /T25, A, "  ", $)'
        printf, lun, dot_h_files[n], FORMAT=fmt
        length = 26
      ENDIF ELSE BEGIN
        printf, lun, dot_h_files[n], $
         FORMAT=(n EQ Ninclude - 1) ? '(A)' : '(A, "  ", $)'
        length = length + 2
      ENDELSE
    ENDFOR
  ENDFOR

  IF (keyword_set(OUTPUT_FILE)) THEN free_lun, lun
END


  
