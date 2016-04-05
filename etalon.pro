;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   PRO ETALON
;
;   PURPOSE:
;            It calculates the transmission profile of an etalon
;            as a function of the wavelength. When the etalon is
;            in a collimated configuration, the angle of refraction
;            acts as a wavelength centering tool. When the etalon is
;            in a telecentric configuration, the resulting
;            transmission profile is the integral of all rays
;            belonging to the incidence cone of light.
;
;   INPUT PARAMETERS:
;           lam:  if lam is not define, set 1001 points with 1pm stepsize
;                 otherwise use it as input lambda vector (in nm)
;
;   KEYWORDS:
;
;            /COL: When this keyword is set, a collimated
;                  configuration is assumed.
;            /HRT: When this keyword is set, the HRT values are taken
;            /GET_PARAMS: When this keyword is set, params is NOT empty
;            ROUGH = ROUGH: error to add to thickness in nm
;            theta = tetha: if given, the tilt angle of the etalon (degree)
;   OUPUT:
;
;            params: Structure containing
;               params.tau: Transmission peak
;               params.fwhm: FWHM of the transmission profile (in pm)
;               params.lmax: Shift in mA with respect to the nominal position
;               params.lfree: Free spectral range in nm
;               params.finesse: Etalon finesse
;            lam: Array with wavelengths where g is calculated (input as well)
;            g: Transmission profile as a funciton of lambda
;            dg = dg: Derivative of g with respect to h (in nm^-1)
;
;   Created by: J.C. del Toro Iniesta
;               March, 2016
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;pro etalon, tau, fwhm, lmax, finesse, lam, g, lfree, dg, COLLIMATED = COLLIMATED, HRT = HRT, rough = rough,theta = theta

pro etalon, lam, g, dg = dg, COLLIMATED = COLLIMATED, $
    HRT = HRT, rough = rough,theta = theta,GET_PARAMS = GET_PARAMS

  if NOT(KEYWORD_SET(rough)) then ROUGH = 0d0
  if NOT(KEYWORD_SET(theta)) then theta = 0d0

  ;it is usefull in case I provide the wavelength axis becouse
  ;it needs more than one peak to calculate some of the parameters

  if KEYWORD_SET(GET_PARAMS) then begin
      fill_params = 1
      get_params = {tau:0.0, fwhm:0.0, lmax:0.0, finesse:0.0, lfree:0.0}
  endif else fill_params = 0

;VARIABLES

r_fine = 101 ;radio sampling ;NOTE line 101 change units!!!! (DO NOT MODIFY!!!)

lam0 = 617.33356d0   ; central wavelength in nm
n = 2.2908d0   ; refractive index (value for LiNbO_3 at 21 ºC)
;n = 2.292635d0   ; de Julián
h = 251.63d0+ROUGH/1d3   ; thickness in µm (for the SO/PHI etalon)
R = 0.92d0   ; internal reflectivity
Ab = 0.d0   ; absorptivity

;wavelegth axis (check whether it is defined or not)
if type(lam) eq 0 then begin
    l_fine = 1001 ;wavelength sampling
    lam = lam0 - 0.5d0 + findgen(l_fine)*1.d-3
    ; +/- 5 Angstroms from the central wavelength (in nm)
endif else begin
    lam = lam*1d0
    l_fine = n_elements(lam)
endelse

tau = (1.d0 - Ab/(1-R))^2    ; transmission peak
EFE = 4.d0*R/(1-R)^2

if KEYWORD_SET(HRT) then begin

   f = 7.92d3   ; focal at the etalon in mm for the HRT
;   f = 7.91d3   ; focal de Julián for the HRT
   rpup = 7.d1   ; pupil radius in mm for the HRT

endif else begin ;FDT CASE

   f = 1.1117d3   ; focal at the etalon in mm for the FDT
;   f = 9.982d2   ; focal de Julián for the FDT
   rpup = 8.75d0   ; pupil radius in mm for the FDT

endelse

; setting distances to m

h1 = h*1.d-6
f1 = f*1.d-3
rpup1 = rpup*1.d-3

a = 2.d0*!pi*n*h1/(lam*1.d-9)*cos(theta*!dpi/180d0)
lam_n = lam*1.d1 - 6173.3356d0 ; wavelength distance in Angstroms
;I just changed from lam to lam_n to not overide the input parameter

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; COLLIMATED CONFIGURATION ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if KEYWORD_SET(COLLIMATED) then begin

   theta0 = [0., 10., 20., 30.]*1.d-2   ; incidence angle in degree

; cosine of the refraction angle through Snell's law (external
; medium is vacuum)

   costheta =sqrt(1.d0 -(sin(theta0*!pi/180.d0)/n)^2)

   epsilon = 4/sqrt(EFE)
   fwhm = epsilon*lam0*lam0/4.d0/!pi/n/h ; in pm
   finesse = !pi*sqrt(EFE)/2.d0

   lfree = (lam0*1.d-9)^2/2.d0/n/h1/costheta * 1.d9 ; free spectral range in nm

   delta = 2.d0*costheta#a      ; phase difference per double reflection

   g = tau /(1.d0 + EFE*sin(delta/2.d0)^2d0)
   dg = -EFE*tau*delta*sin(delta)/(h*1.d3)/(1.d0+EFE*sin(delta/2.d0)^2d0)^2d0

   lmax = fltarr(n_elements(theta0))
   for i=0,n_elements(theta0)-1 do lmax(i) = 1.d3*lam_n(where(g(i,*) eq max(g(i,*))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; TELECENTRIC CONFIGURATION ;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endif else begin

   b = 1./(n*f1)   ; the external medium is assumed to be the vacuum

   rr = rpup1*findgen(r_fine)*1.d-3   ; distances across the pupil

   g = fltarr(l_fine)
   dg = g
   for i = 0,l_fine-1 do begin
      den = 1.d0 + EFE*sin(a(i)*sqrt(1-b^2*rr^2))^2
      g(i) = 2.d-6*tau * total(findgen(r_fine)/den)
      dg(i) = -2.d-6*tau*a(i)*EFE/(h*1.d3) * total(findgen(r_fine)*sqrt(1-b^2*rr^2)*sin(2.d0*a(i)*sqrt(1-b^2*rr^2))/den^2)
   endfor

;; YOUR ATTENTION!!!
;  the subsequent part of the procedure is not general. It must be
;  tuned according to the final transmission profile. Please plot
;  first and then tune the wavelength for a single peak

   if fill_params then begin
     donde = where(lam_n gt -3 and lam_n lt 0)   ; a single transmission peak
     gg = g(donde)
     ll = lam_n(donde)
     posmax = where(gg eq max(gg))
     get_params.lmax = ll(posmax)*1.d3   ; peak position in mA
     get_params.tau = (gg(posmax))[0]   ; re-evaluated transmission peak
     posmed = where(abs(gg - get_params.tau/2.d0) eq min(abs(gg - get_params.tau/2.d0)))
     get_params.fwhm = 2.d2*abs(ll(posmax)-ll(posmed))   ; in pm

     dd = deriv(g)
     wh = where(sign(deriv(g))-shift(sign(deriv(g)),1) lt 0)
     get_params.lfree = lam_n(wh(2)) - lam_n(wh(1))
     get_params.finesse = get_params.lfree/get_params.fwhm*1d2

 endif

endelse

end
