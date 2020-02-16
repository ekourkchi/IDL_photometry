pro glga_get_inventory
;+
; glga_get_inventory - make a *.glga file for every galaxy on the glga data disk
;-
;
	print,'Reading master list ... ',form='($,a)'
	readcol,!GLGA_ROOT+'data/inventory/master_list.glga', $
		mlid,mlra,mldec,mlmjx,mlmnx,mlpa,mlty, $
		format='a,d,d,f,f,f,a',/silent
	nml = n_elements(mlra)
	print,'Done.'

	surveys = ['galex','sdss','2mass','wise','irac']
	nsurveys = n_elements(surveys)
	srvlab = strarr(nsurveys)
	for i=0,nsurveys-1 do srvlab[i] = string('# ',i,surveys[i], $
		format='(a2,i3,2x,a-5)')

	invfil = !GLGA_ROOT+'data/inventory/inventory.txt'
	filestamp,invfil,/arch
	openw,li,invfil,/get_lun
	printf,li,'# GLGA_GET_INVENTORY: '+systime()
	for i=0,nsurveys-1 do printf,li,srvlab[i]

	gv1fil = !GLGA_ROOT+'data/inventory/glga_v1.txt'
	filestamp,gv1fil,/arch
	openw,l1,gv1fil,/get_lun
	printf,l1,'# GLGA_GET_INVENTORY: '+systime()
	for i=0,nsurveys-1 do printf,l1,srvlab[i]

	gv2fil = !GLGA_ROOT+'data/inventory/glga_v2.txt'
	filestamp,gv2fil,/arch
	openw,l2,gv2fil,/get_lun
	printf,l2,'# GLGA_GET_INVENTORY: '+systime()
	for i=0,nsurveys-1 do printf,l2,srvlab[i]

	xfil = !GLGA_ROOT+'data/inventory/extra.txt'
	filestamp,xfil,/arch
	openw,lx,xfil,/get_lun
	printf,lx,'# GLGA_GET_INVENTORY: '+systime()
	for i=0,nsurveys-1 do printf,lx,srvlab[i]

	for i=0L, nml-1L do begin

		id = mlid[i]
		ra = mlra[i]
		dec= mldec[i]
		mjx= mlmjx[i]
		mnx= mlmnx[i]
		pa = mlpa[i]
		ty = mlty[i]
		ddir = !GLGA_ROOT + 'data/' + glga_degdir(ra) + '/'
		;
		; get catalog status based on canonical diameters
		if mjx ge 1.5 then begin
			stat='glga_v1'
		endif else if mjx ge 0.8 then begin
			stat='glga_v2'
		endif else stat='extra'
		;
		; no anonymous galaxies in catalog
		if strmid(strtrim(id,2),0,1) eq 'A' then stat='extra'
		;
		; check ellipse file
		elf= ddir+'aux/'+id+'_ellipse.dat'
		if file_test(elf) then begin
			elfstat = 'elf'
			readcol,elf,majd_as,mind_as,ra,dec,pa, $
				form='f,f,d,d,f',/silent
			mjx = majd_as/60.
			mnx = mind_as/60.
		endif	else elfstat = 'elx'
		;
		; search for images
		flist = file_search(ddir+'*/fits/'+id+'[-_]*.fit*', count=nf)
		slist = extract_surveys(flist, nsurv)
		;
		; search for dss
		flist = file_search(ddir+'dss/fits/'+id+'_*.fit*', count=ndss)
		if ndss gt 0 then $
			dsstat = 'dss' $
		else	dsstat = 'dxx'
		;
		; search for jpegs
		flist = file_search(ddir+'*/jpg/'+id+'_*.jpg', count=nf)
		jlist = extract_surveys(flist, njsurv)
		;
		; get status of each image data set
		inven = (1 eq 0)	; are we in inventory?
		imstat = ['ix0','ix1','ix2','ix3','ix4']
		qastat = ['qx0','qx1','qx2','qx3','qx4']
		jpstat = ['jx0','jx1','jx2','jx3','jx4']
		for j=0,nsurveys-1 do begin
			;
			; check images
			idx = strcmp(slist,surveys[j])
			if total(idx) gt 0 then begin
				imstat[j] = 'im' + strn(j)
				inven = (1 eq 1)
			endif
			;
			; set qa status
			qaf = ddir+surveys[j]+'/fits/' + $
				id + '_qa.txt'
			qa = glga_read_qa_stat(qaf,stat=qas)
			qastat[j] = 'q' + strn(qas) + strn(j)
			;
			; set jpg status
			idx = strcmp(jlist,surveys[j])
			if total(idx) gt 0 then jpstat[j] = 'jp' + strn(j)
		endfor
		if inven then $
			printf,li,id,ra,dec,mjx,mnx,pa,ty,stat,elfstat,dsstat, $
				imstat,jpstat,qastat, $
		format='(a-25,2f13.8,3f9.3,2x,a-8,2x,a-8,2x,a-4,2x,a-4,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5)'
		if strcmp(stat,'glga_v1') eq 1 then $
			printf,l1,id,ra,dec,mjx,mnx,pa,ty,stat,elfstat,dsstat, $
				imstat,jpstat,qastat, $
		format='(a-25,2f13.8,3f9.3,2x,a-8,2x,a-8,2x,a-4,2x,a-4,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5)'
		if strcmp(stat,'glga_v2') eq 1 then $
			printf,l2,id,ra,dec,mjx,mnx,pa,ty,stat,elfstat,dsstat, $
				imstat,jpstat,qastat, $
		format='(a-25,2f13.8,3f9.3,2x,a-8,2x,a-8,2x,a-4,2x,a-4,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5)'
		if strcmp(stat,'extra') eq 1 then $
			printf,lx,id,ra,dec,mjx,mnx,pa,ty,stat,elfstat,dsstat, $
				imstat,jpstat,qastat, $
		format='(a-25,2f13.8,3f9.3,2x,a-8,2x,a-8,2x,a-4,2x,a-4,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5,2x,' + $
			strn(nsurveys) + 'a-5)'
		;
		;
		print,string(13B),i+1,'/',nml,id,stat, $
			format='($,a1,i7,a1,i7,2x,a-25,2x,a-7)'
	endfor
	free_lun,li,l1,l2,lx
	print,' '
	print,'Done.'
;
	return
end
