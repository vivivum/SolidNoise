pro generate_t_profile

dimx = 128
dimy = 128

init_milos,'6173',wl

Roughness_map = PERLINNOISE_2D(dimx,dimy,persistence=0.95,iorder=2)

step=10d0
Points=64.
STEP=STEP/1000d0
eje=wl(1)+Dindgen(Points)*step-0.310
eje_nm = (eje)/10.

t_filter = dblarr(dimx,dimy,points)
t_fwhm = dblarr(dimx,dimy)
t_lmax = dblarr(dimx,dimy)
t_tau = dblarr(dimx,dimy)
t_finesse = dblarr(dimx,dimy)
t_lfree = dblarr(dimx,dimy)

for i=0,dimx-1 do begin
for j=0,dimy-1 do begin

  etalon, eje_nm, g, rough = Roughness_map[i,j],theta=1.31d0

  t_filter[i,j,*] = g
  ;t_fwhm[i,j] = params.fwhm
  ;t_lmax[i,j] = params.lmax
  ;t_tau[i,j] = params.tau
  ;t_finesse[i,j] = params.finesse
  ;t_lfree[i,j] = params.lfree
  ;undefine,lam

endfor
print,i,j
endfor

save,filename='t_profile.save',t_filter,t_tau
end
