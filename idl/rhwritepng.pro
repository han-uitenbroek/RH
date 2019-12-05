PRO rhwritepng, WindowNo, baseName

  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  wset, WindowNo
  filename = baseName + timeStamp() + '.png'
  write_png, filename, tvrd(/ORDER), r_curr, g_curr, b_curr
  r = dialog_message(/INFORMATION, "Wrote graph to file: " + filename)
END
