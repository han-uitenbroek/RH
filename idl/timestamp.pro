FUNCTION timestamp

  spawn, /NOSHELL, 'date', result
  return, strmid(result[0], strpos(result[0], ':') - 2, 8)
END
