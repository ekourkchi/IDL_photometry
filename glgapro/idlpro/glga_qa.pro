pro glga_qa,lfile,sdss=sdss,galex=galex,twomass=twomass,irac=irac, $
	wise=wise, panstarrs=panstarrs, acs=acs, verbose=verbose,skip=skip,start=start,nobrowse=nobrowse, $
	pick_file=pick_file,annuli_size=annuli_size, wide=wide, $
	altsurv=altsurv, showall=showall, nofind=nofind, notv=notv, ra=ra
;+
; glga_qa - perform qa on image data
;
; lfile - list of objects with one line per object with these columns:
;	id
;	ra,dec	- degrees
;	majdiam,mindiam	- arcmin
;	pa	- degrees
;
; keywords:
;	sdss,galex,twomass,wise,irac- data to analyze
;	verbose - print out more details
;	skip	- skip all previously qa'd galaxies
;	start	- start at this galaxy (string name or sequence integer)
;	nobrowse- output qa results: if not set no output generated
;	pick_file - a filename to select objects of interest
;	annuli_size - size of aperture annuli in arcseconds
;	wide	- set to use wide jpg (urz, w1w2w4)
;	altsurv - alternate survey name (string: galex, sdss, 2mass, wise)
;	showall - set to display all alternate survey data
;	nofind - set to disable automatic star find
;
;-
; set defaults and check keywords
surveys = ['galex','sdss','2mass','wise','irac', 'panstarrs', 'acs']
;
;
srvn = 0
asrvn = [1,3,2]
ptail='_NUV.fits*'	; primary file tail
atail='_FUV.fits*'	; alternate file tail
jtail='_FUVNUV.jpg'	; jpeg file tail
if keyword_set(sdss) then begin
	srvn=1
	asrvn = [0,3,2]
	ptail='_r.fits*'
	atail='_g.fits*'
	bands=['u','g','r','i','z']
	if keyword_set(wide) then $
		jtail='_urz.jpg' $
	else	jtail='_gri.jpg'
endif
if keyword_set(twomass) then begin
	srvn=2
	asrvn = [0,3,1]
	ptail='_k.fits*'
	atail='_j.fits*'
	jtail='_jhk.jpg'
	bands = ['j','h','k']
endif
if keyword_set(wise) then begin
	srvn=3
	asrvn = [0,1,2]
	ptail='_w1.fits*'
	atail='_w2.fits*'
	bands = ['w1','w2','w3','w4']
	if keyword_set(wide) then $
		jtail='_w124.jpg' $
	else	jtail='_w123.jpg'
endif
if keyword_set(irac) then begin
	asrvn = [0,3,2]
	srvn=4
	bands = ['3p6um','4p5um']
	ptail='_4p5um.fits*'
	atail='_3p6um.fits*'
	jtail='_3p6um4p6um.jpg'
endif
if keyword_set(panstarrs) then begin
	srvn=5
	asrvn = [0,3,2]
	ptail='_r.fits*'
	atail='_g.fits*'
	jtail='_gri.jpg'
	bands=['g','r','i','z']
endif
if keyword_set(acs) then begin
	srvn=6
	asrvn = [0,3,2]
	ptail='_F814W.fits*'
	atail='_F606W.fits*'
	jtail='_color.jpg'
	bands=['F606W','F814W']
endif


srvy=surveys[srvn]
nbands = n_elements(bands)

if keyword_set(altsurv) then begin
	a = where(strcmp(surveys,strlowcase(altsurv)) eq 1, na)
	if na eq 1 then begin
		asrvn = [a[0]]
	endif else begin
		print,'GLGA_QA: Error - unknown alternate survey: ',altsurv
		return
	endelse
endif
defaltsrvy = surveys[asrvn[0]]
nalts = n_elements(asrvn)
;
; check input file
if not file_test(lfile) then begin
	print,'GLGA_QA: Error - file not found: ',lfile
	
	if keyword_set(ra) then begin
	   
	   id = [lfile]
	   ra = [ra]
	   dec = [0]
	   majdiam = [1.]
	   mindiam = [1.]
	   pa = [45.]
	   type = ['???']
	
	endif else return
endif
;
; read in sample data

if not keyword_set(ra) then begin
readcol,lfile, id, ra, dec, majdiam, mindiam, pa, type, $
	format='a,d,d,f,f,f,a', /silent
endif


if max(strlen(type)) le 0 then begin
	print,'GLGA_QA: Warning - no type field in input file (after pa)'
	q=''
	read,'Continue? (y/N): ',q
	if strupcase(strmid(q,0,1)) ne 'Y' then return
	readcol,lfile, id, ra, dec, majdiam, mindiam, pa, $
		format='a,d,d,f,f,f', /silent
	type = replicate('-',n_elements(id))
endif
;
; deal with negative sizes
a = where(majdiam le 0.0, count)
if count gt 0 then majdiam[a] = 0.5
if count gt 0 then mindiam[a] = 0.5
a = where(mindiam le 0.0,  count)
if count gt 0 then mindiam[a] = majdiam[a]
;
rat = majdiam/mindiam
;
;convert diam into arcsec and expand by 50%
d=majdiam*60.0*1.5
;
; define top level directory
deg = string(floor(ra), format='(i3.3)')+'D'
;
; used to generate file names
id = strcompress(id,/rem)

filebase=!GLGA_ROOT+'data/'+deg+'/'+srvy+'/fits/'+id
jpgbase=!GLGA_ROOT+'data/'+deg+'/'+srvy+'/jpg/'+id
;
; starting index
i0=0L
if keyword_set(start) then begin
	if size(start,/type) eq 7 and strtrim(start,2) ne '' then begin
		w=where(strcmp(id,strtrim(start,2)) eq 1, nw)
		if nw gt 1 then begin
			print,'GLGA_QA: Error - ambiguous id: ',start
			return
		endif else if nw lt 1 then begin
			print,'GLGA_QA: Error - not found: ',start
			return
		endif else i0 = w[0]
	endif else if (size(start,/type) ge 12 and size(start,/type) le 15) or $
		(size(start,/type) ge 1 and size(start,/type) le 3) then $
		i0 = start - 1L
endif
nloop = n_elements(id)
miscinfo=''
;
; set up pick file
if keyword_set(pick_file) then begin
	openw,pl,pick_file,/get_lun,/append
	printf,pl,'# GLGA_QA Selected Objects - '+systime(0)
	printf,pl,'# Selected from: '+lfile
	free_lun,pl
endif
;
; loop over object list


; this is the ulternate variable to do automatic photometry when 
; keyword 'notv' is chosen 
rerun = (1 eq 1)


for i=i0, nloop-1 do begin


;
; rerun target
	rerun:
;
; print host
	print,i+1,'/',nloop,id[i],deg[i],' start QA...', $
		format = '(i6,a1,i6,2x,a-25,a5,a)'
;
; directories
	outdir_plots = !GLGA_ROOT+'data/'+deg[i]+'/plots/'  
	outdir_aux = !GLGA_ROOT+'data/'+deg[i]+'/aux/'  
;
; log what happens
	qalogfile = filebase[i] + '_qa.txt'
	qa=glga_read_qa_stat(qalogfile)
	qastatus = ''
	
	
	if qa.ts gt 0. then begin
		if qa.complete then begin
			qastatus='QA complete'
			if keyword_set(skip) then goto, skipover
			if keyword_set(nobrowse) and not keyword_set(notv) then begin
				q=''
				read,'QA complete, <cr> to skip, Q to quit,'+$
					' p to go back, else re-do qa: ',q
				if q eq 'Q' then goto, done
				if q eq '' then begin 
				  miscinfo=''
				  goto, skipover
				endif
				if q eq 'p' then begin 
				   miscinfo='PREVIOUS'
				   goto, previous
				endif
			endif
			

			
		endif
	endif


	if keyword_set(notv) and not rerun then begin
	   rerun = (1 eq 1)
	   miscinfo=''
	   goto, skipover
	endif
	if keyword_set(notv) and rerun then begin
	   rerun = (1 eq 0)
	endif			
     
; generate filenames
	pfile=filebase[i]+ptail	; primary file
	afile=filebase[i]+atail	; alternate file
	jpgfile=jpgbase[i]+jtail	; jpeg file
	photplot = outdir_plots+id[i]+'_'+srvy+'_profile.jpg'
	imageplot = outdir_plots+id[i]+'_'+srvy+'_images.jpg'
	imageplot_masks = outdir_plots+id[i]+'_'+srvy+'_images_masks.jpg'
	ellipsefile = outdir_aux+id[i]+'_ellipse.dat'
	psrcfile = outdir_aux+id[i]+'_'+srvy+'_pointsrc.dat'
	roifile = outdir_aux+id[i]+'_'+srvy+'_roi.dat'
	mimgfile = outdir_aux+id[i]+'_'+srvy+'_mask.fits.gz'
	auxbase = outdir_aux+id[i]+'_'+srvy
	;
	
	; get refimplot
	for j = 0, nalts-1 do begin
		refimplot = outdir_plots+id[i]+'_'+ $
			surveys[asrvn[j]]+'_images.jpg'
		if file_test(refimplot) then begin
			refimn = asrvn[j]
			break
		endif else refimn = -1
	endfor
;
; check primary and alternate files
	nopfile = 0 &  noafile = 0
	finfo=file_info(pfile)
	if not finfo.exists then begin
		nopfile = 1
		qa.ts_primary_img = 0LL
	endif else qa.ts_primary_img = finfo.mtime
	finfo=file_info(afile)
	if not finfo.exists then begin
		noafile = 1
		qa.ts_secondary_img = 0LL
	endif else qa.ts_secondary_img = finfo.mtime
	usefile = pfile
	if nopfile then usefile = afile
;
; skip if neither is available
	if nopfile and noafile then begin
		qa.error = 1
		if keyword_set(verbose) then $
                        PRINT, pfile
                        print, afile
			print,'No images, skipping'
		goto,skiperr
	endif else qa.error = 0
;
; check jpg file
	finfo=file_info(jpgfile)
	if finfo.exists then begin
		qa.ts_jpg_img = finfo.mtime
		read_jpeg, jpgfile, jdata, /true
		nojfile = 0
	endif else begin
		qa.ts_jpg_img = 0LL
		jdata = -1.
		nojfile = 1
	endelse
;
; if it's not there, skip
	qa.require_jpg = 1
	if nojfile then begin
		qa.error = 1
		if keyword_set(verbose) then $
			print,'No jpgs, skipping'
		goto,skiperr
	endif else qa.error = 0
;
; set display variables
	device,get_screen_size=screen_size
	print, "Screen Size: ", screen_size
	if screen_size[0] gt 1900 then xmx = 1526 else xmx = 800 ; 2000 1526
	if screen_size[1] gt 1000 then ymx = 950 else ymx = 800 ; 1300 1526
;
; read in plots that are useful to see when doing QA
;
; radial plots
	finfo=file_info(photplot)  
	if finfo.exists and not keyword_set(notv) then begin
		read_jpeg, photplot, jplot, /true
		sz = size(jplot)
		yp = (screen_size[1] - 2.*sz[3] - 25) > 20
		window, 0, xsize = sz[2], ysize = sz[3], xpos = 1, ypos = yp, $
		    title=id[i] + ' Photometry ('+strn(i+1)+'/'+strn(nloop)+')'
		tvscl, jplot, /true
		qa.ts_phot_plots = finfo.mtime
	endif else qa.ts_phot_plots = 0LL
;
; image plots
	finfo=file_info(imageplot)
	if finfo.exists and not keyword_set(notv) then begin
		read_jpeg, imageplot, jplot2, /true
		sz = size(jplot2)
		window, 1, xsize = sz[2], ysize = sz[3], $
			xpos = 1, ypos = screen_size[1], $
			title=id[i] + ' Images ('+strn(i+1)+'/'+strn(nloop)+')'
		tvscl, jplot2, /true
		qa.ts_image_plots = finfo.mtime
	endif else qa.ts_image_plots = 0LL
;       

	finfo=file_info(imageplot_masks)
	if finfo.exists and not keyword_set(notv) then begin
	        
		read_jpeg, imageplot_masks, jplot3, /true
		sz = size(jplot3)
		window, 20, xsize = sz[2], ysize = sz[3], $
			xpos = 1.5*screen_size[0], ypos = screen_size[1]/3, $
			title=id[i] + ' MASKS ('+strn(i+1)+'/'+strn(nloop)+')'
		tvscl, jplot3, /true
	endif


	if keyword_set(sdss) or keyword_set(panstarrs) then begin
	   atail_filt = strmid(atail,0,strpos(atail,'.fit'))
	   single_phot = !GLGA_ROOT+'data/'+deg[i]+'/plots/'+id[i]+'_i.jpg'
	   finfo=file_info(single_phot)	
	endif else begin       
           ptail_filt = strmid(ptail,0,strpos(ptail,'.fit'))
           single_phot = !GLGA_ROOT+'data/'+deg[i]+'/plots/'+id[i]+ptail_filt+'.jpg'
	   finfo=file_info(single_phot)
	endelse
	
	if (keyword_set(sdss) or keyword_set(panstarrs)) and not finfo.exists then begin
	   ptail_filt = strmid(atail,0,strpos(atail,'.fit'))
	   single_phot = !GLGA_ROOT+'data/'+deg[i]+'/plots/'+id[i]+ptail_filt+'.jpg'
	   finfo=file_info(single_phot)	
	endif	
	
	if not finfo.exists then begin 
	   atail_filt = strmid(atail,0,strpos(atail,'.fit'))
	   single_phot = !GLGA_ROOT+'data/'+deg[i]+'/plots/'+id[i]+atail_filt+'.jpg'
	   finfo=file_info(single_phot)
	endif 
	

	
	if finfo.exists and not keyword_set(notv) then begin
	        
		read_jpeg, single_phot, jplot4, /true
		sz = size(jplot4)
		window, 30, xsize = sz[2], ysize = sz[3], $
			xpos = 1.0*screen_size[0], ypos = screen_size[1]/3, $
			title=id[i] + ' MASKS ('+strn(i+1)+'/'+strn(nloop)+')'
		tvscl, jplot4, /true
	endif
	
	
	
	
; read in image data to use
	data = mrdfits(usefile, 0, hdr, /fscale, /silent)
;MHS start addition - 82/2101
;merge irac images if possible for source detection
	if keyword_set(irac) and file_exist(afile) then begin 
		data2= mrdfits(afile, 0, hdr2, /fscale, /silent)
		m=where(finite(data2) and finite(data,/nan),count)
		if count gt 0 then data[m]=data2[m]
		m=where(finite(data2) and finite(data),count)
		data[m]=data[m]>data2[m]
	endif
;MHS end addition - 82/2101
	sz = size(data, /dim)
;
; image plot of a different wave-band
	finfo=file_info(refimplot)
	if finfo.exists  and not keyword_set(notv) then begin
		blab = ' GALEX'
		qa.altim_type='GALEX'
		if strpos(refimplot,'sdss') ge 0 then begin
			blab = ' SDSS'
			qa.altim_type='SDSS'
		endif
		if strpos(refimplot,'panstarrs') ge 0 then begin
			blab = ' PanSTARRS'
			qa.altim_type='PanSTARRS'
		endif	
		if strpos(refimplot,'acs') ge 0 then begin
			blab = ' ACS'
			qa.altim_type='ACS'
		endif			
		if strpos(refimplot,'2mass') ge 0 then begin
			blab= ' 2MASS'
			qa.altim_type='2MASS'
		endif
		if strpos(refimplot,'wise') ge 0 then begin
			blab= ' WISE'
			qa.altim_type='WISE'
		endif
		read_jpeg, refimplot, jplot3, /true
		asz = size(jplot3)
		if asz[0] ge 3 then begin
			xp = min([sz[0],xmx]) + 792
			yp = screen_size[1]
			window, 2, xsize = asz[2], ysize = asz[3], $
				xpos = xp, ypos = yp, $
				title=id[i] + ' Alt Waveband:' + blab
			tvscl, jplot3, /true
		endif
		qa.ts_altim_plots = finfo.mtime
	endif else qa.ts_altim_plots = 0LL
	;
	; check showall
	if keyword_set(showall) and nalts gt 1  and not keyword_set(notv) then begin
		;
		; get other refimplots
		for j = 1, nalts-1 do begin
			arefimplot = outdir_plots+id[i]+'_'+ $
				surveys[asrvn[j]]+'_images.jpg'
			if file_test(arefimplot) and $
			   asrvn[j] ne refimn then begin
				read_jpeg, arefimplot, jplot4, /true
				asz = size(jplot4)
				if asz[0] ge 3 then begin
					xp = min([sz[0],xmx]) + 792
					yp = j*20
					window, j+2, xsize = asz[2],  $
						ysize = asz[3], $
						xpos = xp, ypos = yp, $
						title=id[i]+' Alt Waveband: '+ $
						  strupcase(surveys[asrvn[j]])
					tvscl, jplot4, /true
				endif
			endif
		endfor
	endif	; keyword_set(showall) and nalts gt 1
;
; get astrometry
	extast, hdr, astr
	ad2xy, ra[i], dec[i], astr, xn, yn
	getrot, hdr, rot, cdelt
	as_pix = abs(cdelt[0])*3600.

	

; Taking a copy of the input file ellipse,  construct ellipse vector
ellipse_glga_file = [0.5*d[i]/as_pix, 0.5*(d[i]/rat[i]/as_pix), xn, yn, $
	pa[i],xn,yn]
; get ellipse defining photometric apertures
;
; is there a file modifying the default ellipse?
	if file_exist(ellipsefile) then begin
		qa.ellipse_update = 1
		readcol,ellipsefile,majdiam_as,mindiam_as,el_ra,el_dec, $
			pa_, majdiam_px, mindiam_px, x0_,y0_, astrom_bool, $
			el_as_pix, format='f,f,d,d,f,f,f,f,i,f', /silent
		if astrom_bool ne 1 and keyword_set(verbose) then $
			print,pre+'No astrometry, using default ellipse '+ $
				'parameters' $
		else begin
			ad2xy,el_ra,el_dec,astr,x0_,y0_
			x=x0_[0]
			y=y0_[0]
			pa[i]=pa_[0]
			d[i]=majdiam_as[0]
			rat[i]=majdiam_as[0]/mindiam_as[0]
			if keyword_set(verbose) then $
				print,'Using ellipse info: ',ellipsefile
		endelse
	endif else begin
		qa.ellipse_update = 0
		x=xn
		y=yn
		if keyword_set(verbose) then $
			print,'No ellipse file, using default ellipse '+ $
				'parameters'
	        
	        
	        if  strpos(srvy,'panstarrs')  ge 0 then begin
	        
	        id_rute = strmid(id[i],3,strlen(id[i])-1)
	        rute = strmid(auxbase,0,strpos(auxbase,'aux/'))
	        sdss_ellipse   = rute+'aux/'+id_rute+'_ellipse.dat'
	        
	        if file_exist(sdss_ellipse) then begin 
	          print, 'Using SDSS ellipse ...'
	          
		  readcol,sdss_ellipse,majdiam_as,mindiam_as,el_ra,el_dec, $
			pa_, majdiam_px, mindiam_px, x0_,y0_, astrom_bool, $
			el_as_pix, format='f,f,d,d,f,f,f,f,i,f', /silent	          
	          
	          ad2xy,el_ra,el_dec,astr,x0_,y0_
	          x=x0_[0]
	          y=y0_[0]
	          pa[i]=pa_[0]
	          d[i]=majdiam_as[0]
	          rat[i]=majdiam_as[0]/mindiam_as[0]
	          
	        endif ; if sdss ellipse exists
	        
	        endif  ; servey panstarrs
  
	endelse
        
        
; construct ellipse vector
        IF Finite(pa[i]) EQ 0  then begin 
            print, 'PA is not a defined numnber'
            print, 'PA = 45 deg by default'
            pa[i] = 45.
        endif
	ellipseinfo = [0.5*d[i]/as_pix, 0.5*(d[i]/rat[i]/as_pix), x, y, $
			pa[i],xn,yn]
	
	print,'ellipseinfo: ',ellipseinfo,format='(a,7f9.3)'
;
; set sky ellipses
	minskywid = 27.0 / as_pix	; minimum sky width  5.0 for ACS ; 
	factor = 1.5			; expand by 50%
	skyinner=factor*ellipseinfo[0]
	skyouter=ellipseinfo[0]*sqrt(1+factor^2) > (skyinner+minskywid)  ; ACS (skyinner+minskywid) ; ellipseinfo[0]*sqrt(1+factor^2) > (skyinner+minskywid) 
;
; get sky in sky ellipses
	dist_ellipse, im, sz, x, y, rat[i], pa[i]
	mask=bytarr(sz)*0b
	loc=where(im gt skyinner and im lt skyouter and $
		  finite(data,/nan) eq 0,nloc)
	if nloc gt 0 then mask[loc]=1b
	in=where(mask eq 1b, obj_npix)
	if obj_npix le 0 then in=lindgen(n_elements(data))
	meanclip, data[in], data_sky, data_sigma,  clipsig = 5
	if not keyword_set(galex) then begin
		mmm,data[in], mmm_sky, mmm_sigma,/silent
		if mmm_sigma gt 0. and finite(mmm_sky) eq 1 then begin
			data_sky=mmm_sky
			data_sigma=mmm_sigma
		endif
	endif
;
; set min and max display values
	min = (data_sky-data_sigma*0.8)
	max = data_sky+(7.*data_sigma)
	data[0]=min
	data[1]=max
	print,'sky,sig,min,max: ',data_sky,data_sigma,min,max, $
		format='(a,4g13.3)'
;
; set up color table
	loadct, 0,/silent
	rct
;
; clean up any previous runs
	spawn, 'rm '+auxbase+'_stv_cntrd_output.dat', res, err
	spawn, 'rm '+auxbase+'_stv_ellipse_output.dat', res, err
	spawn, 'rm '+auxbase+'_stv_roi_output.dat', res, err
;

        n_list = [i, nloop]
; display and mark
	pick_obj = (1 eq 0)	; default unselected
	
	if not keyword_set(notv) then begin
	stv, data, jdata, usefile, id[i], sz[0]<xmx , sz[1]<ymx , $
		hd = hdr, min = min, max = max,  $
		ellipseinfo = ellipseinfo, /block, miscinfo = miscinfo, n_list=n_list, $    
		asinh = 0, /true, smooth = 0, psrcfile = psrcfile, $
		roifile = roifile, nobrowse = nobrowse, qa = qa, $
		galex=galex,sdss=sdss,panstarrs=panstarrs,twomass=twomass,wise=wise,irac=irac, acs = acs, $
		pick_obj = pick_obj, nofind = nofind, jpgbase=jpgbase[i],auxbase=auxbase, ellipse_glga_file=ellipse_glga_file
	endif else begin
	  miscinfo = 'RUN'
	endelse
	   
 
;            
; clean up windows
	while !d.window ge 0 do wdelete
;
; check for stop QA-ing
	if miscinfo eq 'STOP' then goto, done
;
; picked?
	if pick_obj and keyword_set(pick_file) then begin
		openw,pl,pick_file,/get_lun,/append
		printf,pl,id[i],ra[i],dec[i],majdiam[i],mindiam[i], $
			pa[i],type[i], format='(a-25,2f13.8,3f9.3,2x,a)'
		free_lun,pl
		if keyword_set(verbose) then $
			print,'Wrote ',id[i],' to ',pick_file
	endif
;
; check for skip QA output
	if miscinfo eq 'SKIP' then goto, skipover
;
; what to do if mask files already exist?
	q='n'
;
; point source list
	if file_exist(auxbase+'_stv_cntrd_output.dat') then begin
		qa.psrc = 1
		if file_exist(psrcfile) then begin
			filestamp,psrcfile
			spawn,'mv '+auxbase+'_stv_cntrd_output.dat ' + psrcfile
		endif else $
			spawn,'mv '+auxbase+'_stv_cntrd_output.dat ' + psrcfile
	endif else qa.psrc = 0
;
; region of interest list
; 	if file_exist(auxbase+'_stv_roi_output.dat') then begin
; 		qa.roi = 1
; 		if file_exist(roifile) then begin
; 			read,'Append roi? <cr> - yes, else - no: ',q
; 			if strlen(strtrim(q,2)) eq 0 then begin
; 				spawn,'cat '+auxbase+'_stv_roi_output.dat >> ' + roifile
; 			endif else begin
; 				filestamp,roifile
; 				spawn,'mv '+auxbase+'_stv_roi_output.dat ' + roifile
; 			endelse
; 		endif else $
; 			spawn, 'mv '+auxbase+'_stv_roi_output.dat ' + roifile
; 	endif else qa.roi = 0
;
; replace previous ellipse info
	if file_exist(auxbase+'_stv_ellipse_output.dat') then begin
		qa.ellipse_update = 1
		filestamp,ellipsefile
		spawn, 'mv '+auxbase+'_stv_ellipse_output.dat ' + ellipsefile
		
		
		; making sure that Pan-STARRS ellipse is also available in the common file
		if  strpos(srvy,'panstarrs')  ge 0 then begin
		   
		   id_rute = strmid(id[i],3,strlen(id[i])-1)
		   rute = strmid(auxbase,0,strpos(auxbase,'aux/'))
		   common_ellipse   = rute+'aux/'+id_rute+'_ellipse.dat'
		   spawn, 'cp '+ellipsefile+' '+common_ellipse
		   
		
		endif
		
		
	endif
;
; re-generate mask list file
        
	if keyword_set(nobrowse) then begin
 		qa.mask = 1 
                if srvn eq 1 then $ ; sdss
		     glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose, filter=['u','g','r','i','z']
                if srvn eq 5 then $ ; panstarrs
		     glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose, filter=['g','r','i','z']
		if srvn eq 3 then $ ; wise
		     glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose, filter=['w1','w2','w3','w4']
		if srvn eq 2 or srvn eq 4 or srvn eq 6 then $
		     glga_genmask,sz,astr,as_pix,mimgfile,verbose=verbose 
		
	endif
;
; previously generated point source and roi lists
	if file_exist(psrcfile) then $
		qa.psrc = 1 $
	else	qa.psrc = 0
	if file_exist(roifile) then $
		qa.roi = 1 $
	else	qa.roi = 0
;
; go here if we can't do QA
	skiperr:
;
; update qa log
        
	if keyword_set(nobrowse) and qa.error ne 1 and miscinfo ne 'RUN' then $
		glga_write_qa_stat, qa, /hasNote
;
; go here if qa already done
	skipover:
;
; print status
	if strlen(miscinfo) gt 2 then $
		qastatus = miscinfo $
	else	  qastatus = 'Done.'
	print,i+1,'/',nloop,id[i],deg[i],qastatus, $
		format = '(i6,a1,i6,2x,a-25,a5,2x,a)'
	print,' '
	print,' '
;       
        
        previous:
	if miscinfo eq 'PREVIOUS' then begin 
	    i=i-2
	    if i eq -1 then i=nloop-2
	endif
	
	if miscinfo eq 'GOTO' then begin 
	    i=n_list[0]
	endif	
	
	if miscinfo eq 'SAVE' then i--
        
; re-run if requested
	if miscinfo eq 'RUN' then begin
		glga_measure,id[i],ra[i],dec[i],majdiam[i],mindiam[i], $
			pa[i], type[i],galex=galex,sdss=sdss,panstarrs=panstarrs,twomass=twomass, acs=acs, $
			wise=wise,irac=irac,verbose=verbose, $
			annuli_size=annuli_size,qa_run=i+1, background_radius=qa.backradius
		
		glga_write_qa_stat, qa, /hasNote
		
		
		;; Creating single-band light profiles
		auxpath=!GLGA_ROOT+'/data/'+deg[i]+'/aux/'
		dsspath=!GLGA_ROOT+'data/'+deg[i]+'/dss/fits/'
		photpath=!GLGA_ROOT+'data/'+deg[i]+'/photometry/'
		plotpath=!GLGA_ROOT+'data/'+deg[i]+'/plots/'
		fpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/fits/'
		jpath=!GLGA_ROOT+'data/'+deg[i]+'/'+srvy+'/jpg/'
		
		maskimgfile=auxpath+id[i]+'_'+srvy+'_mask.fits.gz'
		for j=0,nbands-1 do begin
		   if file_test(fpath+id[i]+'_'+bands[j]+'.fit*') eq 1 then begin
		        
		        bandmaskimgfile=auxpath+id[i]+'_'+srvy+'_mask_'+bands[j]+'.fits.gz'
		        
		         if srvn eq 6 then begin
                             
                            print, 'PLOTRADPROF: '+fpath+id[i]+'_'+bands[j]+'.fit*'
                            acs_plotradprof, id[i], bands[j], type=type[i], $
                               survey=strupcase(srvy), pathtoprofile=photpath, $
                               intfile=fpath+id[i]+'_'+bands[j]+'.fit*', $
                               maskimgfile=maskimgfile, bandmaskimgfile=bandmaskimgfile, outpath=plotpath, $
                               verbose=verbose, /yuan13, update=update, /rotatee
                         endif else begin
                            glga_plotradprof, id[i], bands[j], type=type[i], $
                               survey=strupcase(srvy), pathtoprofile=photpath, $
                               intfile=fpath+id[i]+'_'+bands[j]+'.fit*', $
                               maskimgfile=maskimgfile, bandmaskimgfile=bandmaskimgfile, outpath=plotpath, $
                               verbose=verbose, /yuan13, update=update, /rotatee 
                       
                        endelse  

		   endif
		endfor
		
		goto,rerun
		
		
	endif

endfor	; loop over object list
;
; go here if done QA-ing
done:
;
; clean up
while !d.window ge 0 do wdelete
print, 'Finis QA'
end
