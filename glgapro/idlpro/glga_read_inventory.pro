pro glga_read_inventory,glgainv, verbose=verbose
;+
; glga_read_inventory - read the glga inventory into a useful structure
;-
; test for inventory file
sep='/'
if strpos(!version.os,'Win32') ge 0 then sep='\'
invf = !GLGA_ROOT + 'data'+sep+'inventory'+sep+'inventory.txt'
if file_test(invf) then begin
	A = {glga_inven}
	A = struct_init(A)
	if keyword_set(verbose) then print,'GLGA Inventory ...',form='($,a)'
	readcol,invf,id,ra,dec,mjx,mnx,pa,ty,cat,elf,gxi,sdi,tmi,wii,iri,dsi, $
		qa0,qa1,qa2,qa3,qa4, $
		form='a,d,d,f,f,f,a,a,a,a,a,a,a,a,a,a,a,a,a,a',/silent
	if keyword_set(verbose) then print,' read',form='($,a)'
	ngal = n_elements(ra)
	glgainv = replicate(A,ngal)

	glgainv.id = id
	glgainv.catalog = cat
	glgainv.ra = ra
	glgainv.dec= dec
	glgainv.majax = mjx
	glgainv.minax = mnx
	glgainv.pa = pa
	glgainv.type = ty
	if keyword_set(verbose) then print,' props',form='($,a)'

	flags = intarr(ngal)
	t=where(strpos(elf,'elf') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.elfile = flags

	flags = intarr(ngal)
	t=where(strpos(dsi,'ok') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.dss_imgs = flags

	flags = intarr(ngal)
	t=where(strpos(gxi,'ok') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.galex_imgs = flags

	flags = intarr(ngal)
	t=where(strpos(qa0,'qa0') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	t=where(strpos(qa0,'qf0') ge 0, nt)
	if nt gt 0 then flags[t] = 2
	glgainv.galex_qa = flags

	flags = intarr(ngal)
	t=where(strpos(sdi,'ok') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.sdss_imgs = flags

	flags = intarr(ngal)
	t=where(strpos(qa1,'qa1') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	t=where(strpos(qa1,'qf1') ge 0, nt)
	if nt gt 0 then flags[t] = 2
	glgainv.sdss_qa = flags

	flags = intarr(ngal)
	t=where(strpos(tmi,'ok') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.twomass_imgs = flags

	flags = intarr(ngal)
	t=where(strpos(qa2,'qa2') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	t=where(strpos(qa2,'qf2') ge 0, nt)
	if nt gt 0 then flags[t] = 2
	glgainv.twomass_qa = flags

	flags = intarr(ngal)
	t=where(strpos(wii,'ok') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.wise_imgs = flags

	flags = intarr(ngal)
	t=where(strpos(qa3,'qa3') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	t=where(strpos(qa3,'qf3') ge 0, nt)
	if nt gt 0 then flags[t] = 2
	glgainv.wise_qa = flags

	flags = intarr(ngal)
	t=where(strpos(iri,'ok') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	glgainv.irac_imgs = flags

	flags = intarr(ngal)
	t=where(strpos(qa4,'qa4') ge 0, nt)
	if nt gt 0 then flags[t] = 1
	t=where(strpos(qa4,'qf4') ge 0, nt)
	if nt gt 0 then flags[t] = 2
	glgainv.irac_qa = flags

	if keyword_set(verbose) then print,' flags Done.'

	glgainv.mod_time = systime(1)

endif else print,'Inventory missing.'

return
end
