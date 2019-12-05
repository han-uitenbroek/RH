FUNCTION is_big_endian

  return, (byte(1, 0, 2))[1]
END
