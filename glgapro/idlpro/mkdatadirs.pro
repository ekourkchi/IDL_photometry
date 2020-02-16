pro mkdatadirs
;subdirectory designation (00H to 23H)
deg = string(indgen(360), format='(i3.3)')+'D'

openw,ol,'mkdatadirs.csh',/get_lun

root=!GLGA_ROOT+'data/'
printf,ol,'mkdir '+root

for i=0,n_elements(deg)-1 do begin
	printf,ol,'mkdir '+root+deg[i]
	printf,ol,'mkdir '+root+deg[i]+'/aux'
	printf,ol,'mkdir '+root+deg[i]+'/plots'
	printf,ol,'mkdir '+root+deg[i]+'/photometry'
	printf,ol,'mkdir '+root+deg[i]+'/thumbs'
	printf,ol,'chmod 777 '+root+deg[i]+'/thumbs'
	printf,ol,'mkdir '+root+deg[i]+'/galex'
	printf,ol,'mkdir '+root+deg[i]+'/galex/fits'
	printf,ol,'mkdir '+root+deg[i]+'/galex/jpg'
	printf,ol,'mkdir '+root+deg[i]+'/dss'
	printf,ol,'mkdir '+root+deg[i]+'/dss/fits'
	printf,ol,'mkdir '+root+deg[i]+'/dss/jpg'
	printf,ol,'mkdir '+root+deg[i]+'/sdss'
	printf,ol,'mkdir '+root+deg[i]+'/sdss/fits'
	printf,ol,'mkdir '+root+deg[i]+'/sdss/jpg'
	printf,ol,'mkdir '+root+deg[i]+'/2mass'
	printf,ol,'mkdir '+root+deg[i]+'/2mass/fits'
	printf,ol,'mkdir '+root+deg[i]+'/2mass/jpg'
	printf,ol,'mkdir '+root+deg[i]+'/wise'
	printf,ol,'mkdir '+root+deg[i]+'/wise/fits'
	printf,ol,'mkdir '+root+deg[i]+'/wise/jpg'
	printf,ol,'mkdir '+root+deg[i]+'/irac'
	printf,ol,'mkdir '+root+deg[i]+'/irac/fits'
	printf,ol,'mkdir '+root+deg[i]+'/irac/jpg'
endfor

free_lun,ol

return
end
