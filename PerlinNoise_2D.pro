Function Noise_2d,seed,x,y
   return,randomu(seed,x,y)
end

Function Noise_1d,seed,x
   return,randomu(seed,x)
end

Function PerlinNoise_2D,x,y,persistence=persistence,iorder = iorder

  if not(keyword_set(persistence)) then persistence = Noise_1d(seed,1)
  if not(keyword_set(iorder)) then iorder = 1d0 else iorder = iorder*1d0>1d0
       Map = dblarr(x,y)
       p = persistence[0]

       i = iorder
       frequency = 2d0^i
       amplitude = p^i
       repeat begin
           map = map + congrid(Noise_2d(seed,frequency,frequency),x,y,cubic=-0.5)*amplitude
           print,frequency,amplitude,p
           i = i + 1d0
           frequency = 2d0^i
           amplitude = p^i
       endrep until frequency gt x

       return,map/max(map)

   end
