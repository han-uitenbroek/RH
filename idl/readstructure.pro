PRO readstructure, s, unit, N_MAX_TAGS=N_max_tags, N_START=N_start

;+
; NAME:
;	READSTRUCTURE
;
; PURPOSE:
;	This procedure reads a structure from logical file unit.
;
; CATEGORY:
;	I/O
;
; CALLING SEQUENCE:
;	READSTRUCTURE, S, Unit, N_START=N_start, N_MAX_TAGS=N_max_tags
;
; INPUTS:
;	S:      The structure to be filled.
;       Unit:   The logical file unit.
;
; KEYWORD PARAMETERS:
;       N_START:     The first tag to be filled. Default is 0.
;	N_MAX_TAGS:  The maximum number of tags to be filled. Default is
;                    N_TAGS(S).
;
; MODIFICATION HISTORY:
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Fri Feb  4 10:19:37 2000 --
;-

  IF (NOT keyword_set(N_START)) THEN N_start = 0
  Nstart = MAX([0, N_start])
  IF (NOT keyword_set(N_MAX_TAGS)) THEN N_max_tags = n_tags(s)
  Ntags = MIN([N_max_tags, n_tags(s)])

  FOR n=Nstart, Nstart+Ntags-1 DO BEGIN
    tmp = s.(n)
    readu, unit, tmp
    s.(n) = tmp
  ENDFOR
END
