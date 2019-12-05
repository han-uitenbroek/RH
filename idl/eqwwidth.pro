FUNCTION eqwwidth, lambda, I

  Nlambda = n_elements(lambda)

  continuum = 0.5 * (I(0) + I(Nlambda-1))
  ref_width = continuum * (lambda(Nlambda-1) - lambda[0])

  return, (ref_width - int_tabulated(lambda, I, /SORT)) / continuum
END
