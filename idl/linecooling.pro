FUNCTION linecooling, atom, kr

  CLIGHT  = 2.99792458D+08
  HPLANCK = 6.6260755D-34
  NM_TO_M = 1.0E-09

  Rij = *(atom.transition[kr].Rij_ptr)
  Rji = *(atom.transition[kr].Rji_ptr)

  i = atom.transition[kr].i
  j = atom.transition[kr].j

  sz  = size(Rij)
  CASE (sz[0]) OF
    1: BEGIN
      ni = (*(atom.n_ptr))[*, i]
      nj = (*(atom.n_ptr))[*, j]
    END
    2: BEGIN
      ni = (*(atom.n_ptr))[*, *, i]
      nj = (*(atom.n_ptr))[*, *, j]
    END
    3: BEGIN
      ni = (*(atom.n_ptr))[*, *, *, i]
      nj = (*(atom.n_ptr))[*, *, *, j]
    END
  ENDCASE

  hnu = (HPLANCK * CLIGHT) / (atom.transition[kr].lambda0 * NM_TO_M)
  return, hnu * (nj*Rji - ni*Rij)
END

