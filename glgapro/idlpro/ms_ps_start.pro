pro ms_ps_start, filename=filename, color=color, bits=bits,$
 xsize=xsize,ysize=ysize,inches=inches, charsize=charsize,$
 xoffset=xoffset,yoffset=yoffset,landscape=landscape,true=true

 set_plot,'ps'
 device,/helv,/isolatin1,filename=filename,color=color,bits=bits,$
   xsize=xsize,ysize=ysize,inches=inches,xoffset=xoffset,yoffset=yoffset,$
   landscape=landscape;decomposed=true

 !p.font=0

 if keyword_set(charsize) then !p.charsize=charsize else !p.charsize=1.5

 !p.thick=4
 !p.symsize=4
 !x.thick=4
 !y.thick=4

return
end
