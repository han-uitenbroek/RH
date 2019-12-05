!* ------- file: -------------------------- humlicek_.f90 -----------
!
!      Version:       rh1.0
!      Author:        Han Uitenbroek (huitenbroek@nso.edu)
!      Last modified: Fri Oct 13 14:21:14 2000 --
!
!      --------------------------                      ----------RH-- *!

!* --- Voigt function subroutines in different parts of parameter space.
!      Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
!      function needed in Stokes radiative transfer. Called from voigt.c
!
!      FORTRAN 90 version.
!
! See: Humlicek 1982, JQSRT 27, p. 437
!      --                                              -------------- *!

!* ------- begin -------------------------- humlicek.f90 ------------ *!

subroutine humlicek(a, v, W)

  implicit none

  real    (kind = 8), intent(in)  :: a, v
  complex (kind = 8), intent(out) :: W
  complex (kind = 8)              :: z, u
  real    (kind = 8)              :: s

  z = cmplx(a, -v)
  s = abs(v) + a

  if (s >= 15.0) then

     !* --- Approximation in region I --               -------------- *!

     W = (z * 0.5641896) / (0.5 + (z * z))
  else if (s >= 5.5) then

     !* --- Approximation in region II --              -------------- *!
     
     u = z * z
     W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
  else if (a >= 0.195*abs(v) - 0.176) then

     !* --- Approximation in region III --             -------------- *!

     W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
          0.5642236*z)))) / &
          (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
          z*(6.699398 + z)))))
  else
     !* --- Approximation in region IV --              -------------- *!

     u = z * z
     W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
          u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
          (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
          u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))
  endif

end subroutine humlicek
!* ------- end ---------------------------- humlicek.f90 ------------ *!
