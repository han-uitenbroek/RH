function replace_tag, s, tagname, tagvalue

  ;; --- Replace the value of tag tagname of structure s with tagvalue.
  ;;     Returns an anonymous version of the modified struture.
  ;;
  ;;     Note: This operation does not preserve the name of the structure

  Ntags = n_tags(s)
  IF (Ntags LE 0) THEN BEGIN
    print, "Expression must be structure: S"
    return, s
  ENDIF

  IF (NOT tag_present(s, tagname, POSITION=tag_position)) THEN BEGIN
    print, "Tag ", tagname, " is not present in structure S"
    return, s
  ENDIF

  names = tag_names(s)
  structname = tag_names(s, /STRUCTURE_NAME)
  IF (tag_position EQ 0) THEN $
   new_s = create_struct(names[0], tagvalue) $
  ELSE $
   new_s = create_struct(names[0], s.(0))

  FOR n=1, n_tags(s)-1 DO BEGIN
    IF (n EQ tag_position) THEN $
     new_s = create_struct(new_s, names[n], tagvalue) $
    ELSE $
     new_s = create_struct(new_s, names[n], s.(n))
  ENDFOR

  return, new_s
END
