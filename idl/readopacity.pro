; -------- file: -------------------------- readopacity.pro ---------- ;

; -------- begin -------------------------- openOpacity.pro ---------- ;

FUNCTION openOpacity, fileName

@geometry.common
@atmos.common
@opacity.common
@files.common
@input.common

  ;; --- Open opacity file for active transitions --   --------------- ;

  WHILE (NOT existFile(fileName, UNIT=opacunit, /XDR)) DO BEGIN
    answer = dialog_message(/QUESTION, $
                             ["File " + fileName + " does not exist", $
                              "Continue?"])
    IF (answer EQ 'Yes') THEN $
     opacFile = dialog_pickfile(FILTER='*.out', TITLE='Find opacity file', $
                                /MUST_EXIST, /READ) $
    ELSE  return, 0
  ENDWHILE

  ;; --- Open opacity file for background transitions --  ------------ ;

  WHILE (NOT existFile(backgroundFile, UNIT=backgroundunit, $
                       SWAP_ENDIAN=(inputdata.big_endian NE $
                                    is_big_endian()))) DO BEGIN
    answer = dialog_message(/QUESTION, $
                             ["File " + backgroundFile + " does not exist", $
                              "Continue?"])
    IF (answer EQ 'Yes') THEN $
     backgroundFile = dialog_pickfile(FILTER='*.dat', $
                                      TITLE='Find background file', $
                                      /MUST_EXIST, /READ) $
    ELSE  return, 0
  ENDWHILE

  return, 1
END
; -------- end ---------------------------- openOpacity.pro ---------- ;

; -------- begin -------------------------- getTau.pro --------------- ;

FUNCTION getTau, path, chi

  N = n_elements(path)
  tau = dblarr(N)

  tau[0] = 0.0
  FOR k=1, N-1 DO BEGIN
    dtau = 0.5*(chi[k-1] + chi[k]) * (path[k-1] - path[k])
    tau[k] = tau[k-1] + dtau
  ENDFOR

  return, tau
END
; -------- end ---------------------------- getTau.pro --------------- ;

; -------- begin -------------------------- readOpacity.pro ---------- ;

PRO readOpacity, lambda_index, mu_index

@input.common
@geometry.common
@atmos.common
@spectrum.common
@opacity.common
@files.common

  IF (n_params() LT 2) THEN mu_index = 0

  lambdaDisplay = lambda_index
  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      as_rec_len = 2L * geometry.Ndep * 8L
      bg_reclen  = geometry.Ndep * 8L
    END
    "TWO_D_PLANE": BEGIN
      as_rec_len = 2L * (geometry.Nx*geometry.Nz) * 8L
      bg_reclen  = (geometry.Nx*geometry.Nz) * 8L
    END
    "THREE_D_PLANE": BEGIN
      as_rec_len = 2L * (geometry.Nx*geometry.Ny*geometry.Nz) * 8L
      bg_reclen  = (geometry.Nx*geometry.Ny*geometry.Nz) * 8L
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      as_rec_len = 2L * geometry.Nradius * 8L
      bg_reclen  = geometry.Nradius * 8L
    END
  ENDCASE

  IF (atmos.moving OR tag_present(atmos, 'STOKES') OR $
      inputData.PRD_angle_dep) THEN $
   as_index = lambdaDisplay*geometry.Nrays + mu_index $
  ELSE $
   as_index = lambdaDisplay

  IF (atmos.moving OR tag_present(atmos, 'STOKES')) THEN $
   backgr_index = 2*(lambdaDisplay * geometry.Nrays + mu_index) + 1 $
  ELSE $
   backgr_index = lambdaDisplay

  ;; --- Allocate space for one wavelength worth of data -- ---------- ;

  CASE geometryType OF
    "ONE_D_PLANE": BEGIN
      chi_as = dblarr(geometry.Ndep)
      eta_as = chi_as

      IF (tag_present(atmos, 'STOKES') AND $
          atmos.backgrflags[lambda_index].ispolarized) THEN BEGIN
        chi_c  = dblarr(geometry.Ndep, 4)
        IF (inputData.magneto_optical) THEN chip_c = dblarr(geometry.Ndep, 3)
        ENDIF ELSE $
       chi_c = dblarr(geometry.Ndep)

      eta_c = chi_c
      scatt = dblarr(geometry.Ndep)
    END
    "TWO_D_PLANE": BEGIN
      chi_as = dblarr(geometry.Nx, geometry.Nz)
      eta_as = chi_as
      IF (tag_present(atmos, 'STOKES') AND $
          atmos.backgrflags[lambda_index].ispolarized) THEN BEGIN
        chi_c  = dblarr(geometry.Nx, geometry.Nz, 4)
        chip_c = dblarr(geometry.Nx, geometry.Nz, 3)
      ENDIF ELSE $
       chi_c = dblarr(geometry.Nx, geometry.Nz)

      eta_c = chi_c
      scatt = dblarr(geometry.Nx, geometry.Nz)
    END
    "THREE_D_PLANE": BEGIN
      chi_as = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)
      eta_as = chi_as
      IF (tag_present(atmos, 'STOKES') AND $
          atmos.backgrflags[lambda_index].ispolarized) THEN BEGIN
        chi_c  = dblarr(geometry.Nx, geometry.Ny, geometry.Nz, 4)
        chip_c = dblarr(geometry.Nx, geometry.Ny, geometry.Nz, 3)
      ENDIF ELSE $
       chi_c = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)

      eta_c = chi_c
      scatt = dblarr(geometry.Nx, geometry.Ny, geometry.Nz)
    END
    "SPHERICAL_SYMMETRIC": BEGIN
      chi_as = dblarr(geometry.Nradius)  &  eta_as = chi_as
      chi_c  = dblarr(geometry.Nradius)  &  scatt  = chi_c  &  eta_c = chi_c
    END
  ENDCASE

  IF (spectrum.as_rn[as_index] GE 0) THEN BEGIN
    point_lun, opacunit, $
     long64(spectrum.as_rn[as_index]) * long64(as_rec_len)
    readu, opacunit, chi_as, eta_as
  ENDIF ELSE BEGIN
    chi_as = zero_array(chi_as)
    eta_as = chi_as
  ENDELSE

  point_lun, backgroundunit, $
   long64(atmos.backgrrecno[backgr_index]) * long64(bg_reclen)
  IF (inputData.magneto_optical AND $
      atmos.backgrflags[lambda_index].ispolarized) THEN $
   readu, backgroundunit, chi_c, chip_c, eta_c, scatt $
  ELSE $
   readu, backgroundunit, chi_c, eta_c, scatt

  IF ((tag_present(atmos, 'STOKES') AND $
       atmos.backgrflags[lambda_index].ispolarized)) THEN BEGIN
    CASE (geometryType) OF
      "ONE_D_PLANE": BEGIN
        chi_c = chi_c[*, 0]
        eta_c = eta_c[*, 0]
      END
      "TWO_D_PLANE": BEGIN
        chi_c = chi_c[*, *, 0]
        eta_c = eta_c[*, *, 0]
      END
     "THREE_D_PLANE": BEGIN
        chi_c = chi_c[*, *, *, 0]
        eta_c = eta_c[*, *, *, 0]
      END
      ELSE:
    ENDCASE
  ENDIF

  
END
; -------- end ---------------------------- readOpacity.pro ---------- ;
