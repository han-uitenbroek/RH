FUNCTION getLine, unit, COMMENT_CHAR=comment_char, EXIT_ON_EOF=exit_on_EOF

  IF (NOT keyword_set(COMMENT_CHAR)) THEN comment_char = '#'

  line = ''
  WHILE (NOT eof(unit)) DO BEGIN
    readf, unit, line
    line = strtrim(line, 2)
    IF ((strmid(line, 0, 1) NE '') AND $
        (strmid(line, 0, 1) NE comment_char)) THEN return, line
  ENDWHILE
  IF (keyword_set(EXIT_ON_EOF)) THEN BEGIN
    print, "Error: reached end of file"
    return, 0
  ENDIF ELSE return, -1

END
