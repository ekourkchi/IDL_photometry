pro epsfile,file,color=color
;
set_plot,'ps'

device,file=file+'.eps'
device,/encap

if keyword_set(color) then $
	device,/color

return
end
