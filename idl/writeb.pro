Ndep = 82
openw, 1, '~/src/rh/Atmos/FALC_82_1G_hor.B', /XDR
B = dblarr(Ndep) + 1.0E-4
gamma = dblarr(Ndep)-90.0/!radeg
chi = dblarr(Ndep) + 90.0/!radeg
writeu, 1, B, gamma, chi
close, 1

b= 0.1 + (0.3-0.1) * (alog(geometry.tau500)- alog(geometry.tau500[0])) / $
 (alog(geometry.tau500[81]) - alog(geometry.tau500[0]))

openw, 1, '~/src/rh/Atmos/FALC_82.B', /XDR
writeu, 1, b, dblarr(82)+0.0/!radeg, dblarr(82)
close, 1

