pro PSCLOSE, pname=pn, $
             transparency=trnsp, $
             print=print

@ps.common

  device, /CLOSE

  if (not keyword_set(pn)) then pn = 'lw'

  if (keyword_set(print)) then begin
    if (ColorPrint) then begin
      if (keyword_set(trnsp)) then $
        spawn, ['xprint', '-transparency', FileName], /NOSHELL $
      else $
        spawn, ['xprint', FileName], /NOSHELL
    endif else begin
      spawn, ['lpr', '-P'+pn, FileName], /NOSHELL
    endelse
  endif

  device, ENCAPSULATED=0
  set_plot, OldDevice
  !P = P_old
  return
end
