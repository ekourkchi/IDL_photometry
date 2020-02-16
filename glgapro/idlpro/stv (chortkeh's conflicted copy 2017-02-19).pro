;STV  -  Widget slide TV display
;
; stv, image, width, height, log=log, min=min, max=max, same=same, $
;       cursor=cursor, auto=auto, losig=losig, hisig=hisig, $
;       hd=hd, true=true, cvfile=cvfile, block=block, ellipse=ellipse
;
; Optional Args
;    width, height : of display window (600,600 default)
;    log: logscale
;    min, max: scale vales
;    same: same values of log, min,max as last time    
;    auto: auto scale -3 to +7 sigma (ignores min,max,losig,hisig)
;    losig, hisig: set Nsigma scale values (ignores min,max)
;    true: true color image (buggy)
;    hd: fits header (useful for Curval coordinates etc.)
;    cursor: x cursor (24 default)
;    cvfile: name of file to send curval output to
;    block: disable widget xmanager no_block feature
;    ellipseinfo: array [semimajor,semiminor,x,y,pa] pixel units
;
; keyboard commands
;
;    'z': Zoom
;    'l': tvList
;    'q': Quit
;    'Q': Quit with miscinfo = 'STOP' 
;    'x': Quit with miscinfo = 'SKIP'
;    'v': Curval (optional output appended to cvfile)
;    'i': Invert colortable
;    'm': Measure photometry
;    'h': Help
;    'c': Centroid with mouse (appends to file aux_base+'stv_cntrd_output.dat')
;    'e': Ellipse tool
;
; by Mark Seibert
; last revised 10/2010
; requires astrolib, stv_curval.pro and cntrd.pro +...
;
;--------------------------------------------------------------------
;
; The features added by Ehsan Kourkchi 08/05/2016
;   
;  Keyboard Commands ...  
;
;   'b': adding regions to measure the background from
;        When a back-region is drawn, the anulus is not used for 
;        background estimation
;   '@': Remove all background regions and use the anulus for estimation 
;
;   'U': Displays SDSS  U-band image, if available
;   'G' 'R' 'I' 'Z' : displays g,r,i,z images respectively
;
;   The following actions need bash-scripts: (available in idlpro directory of the pipeline)
;   'a
;   remove_region.bash  
;   undo_region.bash
;   '~': Undo the last drawn mask region
;   '!': By clicking on an already available mask region, you may remove it
;       
;--------------------------------------------------------------------
pro ds9_reg_maker, reg_file, el, as_pix
        
        if file_exist(reg_file) then spawn, 'rm ' + reg_file
        
   	openw,  lun,reg_file,/get_lun
	printf, lun, 'image'
	
	xc =  strtrim(string(el[2]),2)
	yc =  strtrim(string(el[3]),2)
	a  =  strtrim(string(el[0]),2)
	b  =  strtrim(string(el[1]),2)
	pa =  strtrim(string(float(90.+el[4])),2)
	
	ellipse = 'ellipse('
	ellipse = ellipse + xc + ','
	ellipse = ellipse + yc + ','
	ellipse = ellipse + a + ','
	ellipse = ellipse + b + ','
	ellipse = ellipse + pa
	ellipse = ellipse + ') # color=red width=2'
	printf, lun, ellipse
	
	minskywid = 27.0 / as_pix ; 5.0 for ACS ; 27. otherwise
	factor=1.5
	skyinner=factor*el[0]
	skyouter=el[0]*sqrt(1+factor^2) > (skyinner+minskywid) ; ACS (skyinner+minskywid); el[0]*sqrt(1+factor^2) > (skyinner+minskywid)
    
	xc =  strtrim(string(el[2]),2)
	yc =  strtrim(string(el[3]),2)
	a  =  strtrim(string(skyinner),2)
	b  =  strtrim(string(skyinner/(el[0]/el[1])),2)
	pa =  strtrim(string(float(90.+el[4])),2)
	
	ellipse = 'ellipse('
	ellipse = ellipse + xc + ','
	ellipse = ellipse + yc + ','
	ellipse = ellipse + a + ','
	ellipse = ellipse + b + ','
	ellipse = ellipse + pa
	ellipse = ellipse + ') # color=cyan dash=1' 
	printf, lun, ellipse
	
	
	xc =  strtrim(string(el[2]),2)
	yc =  strtrim(string(el[3]),2)
	a  =  strtrim(string(skyouter),2)
	b  =  strtrim(string(skyouter/(el[0]/el[1])),2)
	pa =  strtrim(string(float(90.+el[4])),2)
	
	ellipse = 'ellipse('
	ellipse = ellipse + xc + ','
	ellipse = ellipse + yc + ','
	ellipse = ellipse + a + ','
	ellipse = ellipse + b + ','
	ellipse = ellipse + pa
	ellipse = ellipse + ') # color=cyan dash=1' 
	printf, lun, ellipse
	
	free_lun,lun  


 
end

;--------------------------------------------------------------------
function read_ds9_circle, reg_file

; read ds9 region file   
OPENR, lun, reg_file, /GET_LUN
array = ''
line = '' 
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array = [array, line] & $
ENDWHILE
FREE_LUN, lun


circles = []
for i=0, n_elements(array)-1 do begin
 
 
 if strcmp(array[i], 'circle', 6, /FOLD_CASE) EQ 1 and $ 
     strmatch(array[i], '*red*', /FOLD_CASE) NE 1 then begin
    line =  array[i]
    line =  STRSPLIT(line,',()',/EXTRACT)
    x = float(line[1])
    y = float(line[2])
    r = float(line[3])
    circles = [circles, [x,y,r]]
 endif
    
endfor


return, circles

end
;--------------------------------------------------------------------
function sign, x

if x gt 0 then return, 1
if x lt 0 then return, -1
if x eq 0 then return, 0

end 
;--------------------------------------------------------------------

function read_ds9_box, reg_file

; read ds9 region file   
OPENR, lun, reg_file, /GET_LUN
array = ''
line = '' 
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array = [array, line] & $
ENDWHILE
FREE_LUN, lun


box = []
for i=0, n_elements(array)-1 do begin
 
 
 if strcmp(array[i], 'box', 3, /FOLD_CASE) EQ 1 and $ 
     strmatch(array[i], '*red*', /FOLD_CASE) NE 1 then begin
    line =  array[i]
    line =  STRSPLIT(line,',()',/EXTRACT)
    x     = float(line[1])
    y     = float(line[2])
    Lx    = float(line[3])
    Ly    = float(line[4])
    theta = float(line[5])
    box = [box, [x,y,Lx,Ly,theta]]
 endif
    
endfor


return, box

end
;--------------------------------------------------------------------
function read_ds9_ellipse, reg_file



; read ds9 region file   
OPENR, lun, reg_file, /GET_LUN
array = ''
line = '' 
WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array = [array, line] & $
ENDWHILE
FREE_LUN, lun


for i=0, n_elements(array)-1 do begin
 
 if strcmp(array[i], 'ellipse', 7, /FOLD_CASE) EQ 1 and $ 
    strmatch(array[i], '*red*', /FOLD_CASE) EQ 1 then begin
    line =  array[i]
    break
    endif
    
endfor

line =  STRSPLIT(line,',()',/EXTRACT)

el = fltarr(5)

el[0] = float(line[3])
el[1] = float(line[4])
el[2] = float(line[1])
el[3] = float(line[2])
el[4] = float(line[5])-90.

; make sure pa is positive
 while el[4] lt 0. do el[4] = el[4] + 180.


return, el


END
;--------------------------------------------------------------------

PRO stv_display_plots, nty

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
	     
     print, 'It is still under developmnet ....'
     print, nty

END	; stv_display_plots
;--------------------------------------------------------------------

PRO stv_find,galaxy=galaxy,outside=outside,noprompt=noprompt

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
;
p = ''
p = DIALOG_MESSAGE('Are you sure, finding new point sources?', /QUESTION,  /CENTER, /DEFAULT_NO, TITLE='New point sources?')
if  strcmp(p, 'No') then return else begin


; data type
if strpos(dtype,'galex') ge 0 then galex = 1 else galex = 0
if strpos(dtype,'sdss') ge 0 then sdss = 1 else sdss = 0
if strpos(dtype,'panstarrs') ge 0 then panstarrs = 1 else panstarrs = 0
if strpos(dtype,'acs') ge 0 then acs = 1 else acs = 0
if strpos(dtype,'2mass') ge 0 then twomass = 1 else twomass = 0
if strpos(dtype,'wise') ge 0 then wise = 1 else wise = 0
if strpos(dtype,'irac') ge 0 then irac = 1 else irac = 0
;
; get inputs
zp=0.
defmult = 4.6
defmagl = 22.0
if irac then begin
	defmult = 2.5
endif
if sdss then begin
	defmult = 2.5
	defmagl = 25.0
endif
if panstarrs then begin
	defmult = 2.5
	defmagl = 25.0
endif
if acs then begin
	defmult = 2.5
	defmagl = 25.0
endif
if galex then begin
	defmult = 2.5
	defmagl = 25.0
endif
if twomass then begin
	defmagl = 18.0
endif
if wise then begin
	zp = sxpar(hdr,'MAGZP')
	defmult = 4.6
	defmagl = 21.0
endif
;
; sigma multiplier
if keyword_set(noprompt) then begin
	sigmult=defmult
	maglim=defmagl
endif else begin
	prompt='Sigma multiplier for thresshold (default = '+string(defmult,form='(f4.1)')+'): '
	imult=''
	print,'Enter negative multiplier to cancel.'
	read,imult,prompt=prompt
	if strtrim(imult,2) eq '' then $
		sigmult=defmult $
	else	sigmult=float(imult)
	;
	if sigmult gt 0. then begin
		;
		; magnitude faint limit
  		prompt='Magnitude faint limit (default = '+string(defmagl,form='(f5.1)')+'): '
  		imagl=''
  		read,imagl,prompt=prompt
  		if strtrim(imagl,2) eq '' then $
			maglim=defmagl $
  		else	maglim=float(imagl)
	endif
endelse

if sigmult gt 0. then begin

  glga_find_ptsrc,im,el,astr,as_pix,asrc,nobj,sigmult=sigmult,maglim=maglim, $
	galex=galex,sdss=sdss,panstarrs=panstarrs,twomass=twomass,wise=wise,irac=irac, acs=acs, $
	galaxy=galaxy, outside=outside, zpmag=zp, /verbose

  print,'STV_FIND - found this many: ',nobj

  vec = reform(psrc[6,*])
  test = 2.
  if keyword_set(outside) then test = 3.
  if keyword_set(galaxy) then test = 4.
;
; erase previous find results
  w=where(vec eq test, nw)
  if nw gt 0 then $
	psrc[6,w] = 0.
;
; insert after last good source
  vec = reform(psrc[6,*])
  w=where(vec gt 0., nw)
  if nw le 0 then $
	i0 = 0L $
  else	i0 = max(w) + 1L
;
; insert
  if nobj gt 0 then begin
	p=0L
	for i=i0,maxp do begin
		psrc[*,i] = asrc[*,p]
		p = p + 1L
		if p ge nobj then break
	endfor
  endif

  wsrc = (1 eq 1)	; write out new point source file

  stv_display
;
endif	; sigmult gt 0.

endelse  ; dialoge yes/no box

END	; stv_find
;--------------------------------------------------------------------

PRO read_psrc, pfile, npread

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
;
; read list
npread = 0l
if file_test(pfile) then begin
	readcol,pfile, a, d, xin, yin, rad, psc_astr, $
		format='d,d,f,f,f,i',/silent

	rpx = rad / as_pix	; convert radius from arcsec into pixels
	npread = n_elements(a)<maxp
	for i=0l,npread-1 do begin
		if have_astr and psc_astr[i] then begin
			ad2xy,a[i],d[i],astr,x,y ; convert ra,dec to x,y px
		endif else begin
			x = xin[i]
			y = yin[i]
		endelse
		psrc[*,i] = [a[i],d[i],x,y,rpx[i],psc_astr[i],1.]
	endfor
endif

END	; read_psrc
;--------------------------------------------------------------------

PRO stv_event, ev

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

Widget_Control, ev.top, Get_UValue=info, /No_Copy
whichButton = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
eventTypes = ['DOWN', 'UP', 'MOTION']

rute = strmid(imfil,0,strpos(imfil,'.fit'))
rute =  strmid(rute, 0, strlen(rute)-2)
rute = strmid(rute,0,strpos(rute,id))
b = strpos(imfil,'.fit')
a = strpos(imfil,'_', /REVERSE_SEARCH)
imfil_band = strmid(imfil,a+1,b-a-1)
filter = tv_filter
if strcmp(strtrim(tv_filter,2), 'w_com') or strcmp(strtrim(tv_filter,2), 'n_com') then $
     filter = imfil_band




if not tag_exist(ev,'type') then begin 
     
     WIDGET_CONTROL, ev.id, get_UVALUE = UVALUE
     
     if strcmp(UVALUE, 'r_sex') then Command_handler, 'sex'
       
     if strcmp(UVALUE, 'seg_sex') then begin

         segment_file =  rute+id+'_'+filter+'.segment.fits'
         if file_exist(segment_file) then begin
             spawn, 'ds9 '+segment_file+' &' 
         endif else $
            p = DIALOG_MESSAGE("Segmentation: not found for '"+ filter +"' band ... ! (RUN SExtractor first)",  /INFORMATION,  /CENTER)
     
     endif
     
     
     if strcmp(UVALUE, 'reg_sex_on') then begin

         reg_segment_file =  rute+id+'_'+filter+'.sex.reg'
         spawn, "xpaset -p ds9 regions load " + reg_segment_file
     endif
     
     if strcmp(UVALUE, 'reg_sex_off') then $
         spawn, "xpaset -p ds9 regions delete all"

     if strcmp(UVALUE, 'tmp_sex') then begin
         
         
         if file_exist(rute+id+'*tmp.fits') then $
             spawn, 'ds9 '+rute+id+'*tmp.fits &' $
         else $
             p = DIALOG_MESSAGE("Not found for '"+ filter +"' band ! (RUN SExtractor first)",  /INFORMATION,  /CENTER)

     endif 
     if strcmp(UVALUE, 'use_sex') then stv_reg_creator, /use_sex

          
     if strcmp(UVALUE, 'force_sex') then stv_reg_creator, /force_sex

     if strcmp(UVALUE, 'clear_sex') then begin
        q = DIALOG_MESSAGE('Remove all SExtractgor masks?', /QUESTION,  /CENTER, /DEFAULT_NO)
        if  strcmp(q, 'Yes') then begin
          
          if file_exist(aux_base+'_stv_usesex_vertices.dat') then spawn, 'rm ' + aux_base+'_stv_usesex_vertices.dat'
          if file_exist(aux_base+'_stv_forcesex_vertices.dat') then spawn, 'rm ' + aux_base+'_stv_forcesex_vertices.dat'
          if file_exist(aux_base+'_stv_forcesex_output.dat') then spawn, 'rm ' + aux_base+'_stv_forcesex_output.dat'
          if file_exist(aux_base+'_stv_usesex_output.dat') then spawn, 'rm ' + aux_base+'_stv_usesex_output.dat'
          stv_display
        endif
     endif
     
     if strcmp(UVALUE, 'gen_back')    then Command_handler, 'b'
     if strcmp(UVALUE, 'remove_back') then Command_handler, '@'
     
     
     
     if strcmp(UVALUE, 'open_pyl') then begin
       
         findpro, 'glga_qa', DirList=DirList, ProList=ProList, /Noprint
         py_root = GETENV('IDL_GLGAPRO') + '/pypro/'
         py_scr = py_root+'pylipse.py'
         rute = strmid(imfil,0,strpos(imfil,'.fit'))
         fits_base =  strmid(rute, 0, strlen(rute)-2)
         a = el[0]          ; major axis  
         b = el[1]          ; minor axis  
         xcenter = el[2]    ; xc
         ycenter = el[3]    ; yc
         pa = el[4]         ; pa
         

         command = 'python ' + py_scr 
         command += ' -p ' + ID
         command += ' -a ' + strtrim(string(a),2)
         command += ' -b ' + strtrim(string(b),2)
         command += ' -x ' + strtrim(string(xcenter),2)
         command += ' -y ' + strtrim(string(ycenter),2)
         command += ' -g ' + strtrim(string(pa),2)
         command += ' -s ' + dtype
         command += ' -j ' + jpg_base
         command += ' -f ' + fits_base
         command += ' -r ' + py_root
         command += ' -B ' + aux_base
         command += ' &'
         
         print, ''
         print, 'Opening Pylipse .... :)'
         print, 'Please waite ...'
         print, ''
         spawn, command
         
     endif     
     
     
     if strcmp(UVALUE, 'load_pyl') then begin
          
          ellipsefile = aux_base+'_stv_ellipse_output.dat'
          
          if file_exist(ellipsefile) then begin

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
			pa=pa_[0]
			d=majdiam_as[0]
			rat=majdiam_as[0]/mindiam_as[0]
			if keyword_set(verbose) then $
				print,'Using ellipse info: ',ellipsefile
		endelse				
                ellipseinfo = [0.5*d/as_pix, 0.5*(d/rat/as_pix), x, y, $
			pa,el[5],el[6]]
	        print,'ellipseinfo: ',ellipseinfo,format='(a,7f9.3)'
	        
	        el = ellipseinfo
	        stv_display
			
			
          endif else q = DIALOG_MESSAGE('NO temporary ellipse file !',  $ 
                   /INFORMATION,  /CENTER, TITLE='Missing ellipse file !')
          
     
     
     endif 
     
     if strcmp(UVALUE, 'measure') then  Command_handler, 'm'
     if strcmp(UVALUE, 'refresh') then  Command_handler, '^'
     if strcmp(UVALUE, 'quit') then  Command_handler, 'q'
     if strcmp(UVALUE, 'save') then  Command_handler, 'save'
     if strcmp(UVALUE, 'skip') then  Command_handler, 'x'
     if strcmp(UVALUE, 'previous') then  Command_handler, 'previous'
     if strcmp(UVALUE, 'exit') then  Command_handler, 'Q'

     if strcmp(UVALUE, 'reg_undo') then  Command_handler, '~'
     if strcmp(UVALUE, 'reg_add') then  Command_handler, 'r'
     if strcmp(UVALUE, 'reg_remove') then  Command_handler, '!'  
     if strcmp(UVALUE, 'reg_clean') then  Command_handler, 'C'     
     if strcmp(UVALUE, 'reg_save') then  Command_handler, 'save_reg'
     if strcmp(UVALUE, 'reg_change') then  Command_handler, '$'
     
     if strcmp(UVALUE, 'ell_edit') then  Command_handler, 'e'
     if strcmp(UVALUE, 'ell_orig') then  Command_handler, '#'
     if strcmp(UVALUE, 'ell_sdss') then  Command_handler, 'ell_sdss'
     if strcmp(UVALUE, 'ell_save') then  Command_handler, 'ell_save'
     if strcmp(UVALUE, 'ell_catal') then  Command_handler, 'ell_catal'
     
     if strcmp(UVALUE, 'star_find') then  Command_handler, 'f'
     if strcmp(UVALUE, 'star_del') then  Command_handler, 'd'
     if strcmp(UVALUE, 'star_edit') then  Command_handler, 'c' 

     if strcmp(UVALUE, 'ds9_open') then  Command_handler, '9'  
     if strcmp(UVALUE, 'ds9_contour') then  Command_handler, '\'
     if strcmp(UVALUE, 'ds9_imp') then  Command_handler, '*'
     if strcmp(UVALUE, 'ds9_exp') then  Command_handler, '='
     if strcmp(UVALUE, 'ds9_close') then  Command_handler, 'c_ds9' 
     
     if strcmp(UVALUE, 'Help_Program') then begin
         
         print, ''
         print, 'Open this URL in your browser: http://svn.pan-starrs.ifa.hawaii.edu/trac/ipp/wiki/glga_manual'
         print, ''
         p = DIALOG_MESSAGE("Open this URL in your web browser: http://svn.pan-starrs.ifa.hawaii.edu/trac/ipp/wiki/glga_manual",  /INFORMATION,  /CENTER)

     
     endif     
   
goto, finish
endif


IF ev.type GT 2 THEN RETURN

if strcmp(eventTypes[ev.type],'DOWN') and strcmp(whichButton(ev.press),'RIGHT') then begin
   
   if lock_reg ne 0 then begin
      lock_button =  WIDGET_INFO(sbase, FIND_BY_UNAME='lockreg')
      lock_reg = 0
      WIDGET_CONTROL, lock_button, SET_VALUE='lock reg'
      goto, finish
   endif

endif
   
if strcmp(eventTypes[ev.type],'DOWN') and strcmp(whichButton(ev.press),'RIGHT') and not strcmp(last_command,'') then begin
   
   
   ; This line is essential, don't pass by reference, i.e last_command is a global variable
    new_command = last_command 
    Command_handler, new_command
    if (xregistered('stv', /noshow)) then $
       widget_control, keyboard_text_id, /clear_events
       goto, finish
endif     

WIDGET_CONTROL, ev.top, GET_UVALUE=slide_window

ename=tag_names(ev,/structure_name)

;just a mouse movement event
IF ename EQ 'WIDGET_DRAW' THEN BEGIN 
        
        if strcmp(eventTypes[ev.type],'DOWN') and strcmp(whichButton(ev.press),'LEFT') then begin
            
          if ev.clicks eq 2 then begin  ; double_click
;               XYOUTS, ev.x, ev.y, 'O', /DEVICE
              tvcirc, ev.x, ev.y, 3, color=colordex('Y')
              p = DIALOG_MESSAGE('Do you really want to re-center the ellispe?', /QUESTION,  /CENTER, /DEFAULT_NO)
              if strcmp(p, 'Yes') then begin
                  el[2] = ev.x
                  el[3] = ev.y 
                  stv_display
              endif
              goto, finish
          endif
          
          
           if ev.clicks eq 1 and lock_reg eq 1 then begin
           

              tvcirc, ev.x, ev.y, reg_size, color=colordex('G')
              xcen = ev.x
              ycen = ev.y
              sz = size(im,/dim)
              sx = sz[0]
              sy = sz[1] 
		xc = ev.x
		yc = ev.y
		rc = reg_size
		x_min = round(xc-rc) 
		x_max = round(xc+rc) 
		y_min = round(yc-rc) 
		y_max = round(yc+rc) 		
		Roi = []
		for xx=x_min, x_max do begin
		  for yy=y_min, y_max do begin
		      
		      if (xx-xc)^2+(yy-yc)^2 le rc*rc then begin
			index = yy*sx + xx
			Roi = [Roi, index]
		      endif
		      
		  endfor
		endfor
		
		if n_elements(Roi) gt 0 then begin
		s=size(im)
		openw, lun, edit_region, /append, /get_lun
		printf,lun,'# STV ROI output (image index): '+systime(0)
		printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

		if Roi[0] ge -1 then $
			for i=0L,n_elements(Roi)-1L do $
				printf,lun,Roi[i]

		endif ; Roi               
                free_lun,lun              
              
           
           endif ; signle left click - locked-reg
          

            stv_update_info_widgets, ev.x, ev.y
        
        endif
        
        
        if strcmp(eventTypes[ev.type],'DOWN') and strcmp(whichButton(ev.press),'MIDDLE') then begin
            
            if ev.clicks eq 2 then begin ; double midle click

                tvcirc, ev.x, ev.y, reg_size, color=colordex('G')
                xcen = ev.x
                ycen = ev.y
                sz = size(im,/dim)
                sx = sz[0]
                sy = sz[1]   
                
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                
; 		if have_astr then begin
; 			XY2AD,xcen,ycen,astr,a,d
; 			psrc_ast = 1
; 		endif else begin
; 			a = xcen
; 			d = ycen
; 			psrc_ast = 0
; 		endelse                
; 
; 		radius = reg_size
;    	
;                 ; install in new array
;                 
; 		w=where(nsrc[0,*] eq 0.)
; 		w=w[0]
; 		nsrc[0,w] = a
; 		nsrc[1,w] = d
; 		nsrc[2,w] = xcen
; 		nsrc[3,w] = ycen
; 		nsrc[4,w] = radius
; 		nsrc[5,w] = psrc_ast
; 		nsrc[6,w] = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		xc = ev.x
		yc = ev.y
		rc = reg_size
		x_min = round(xc-rc) 
		x_max = round(xc+rc) 
		y_min = round(yc-rc) 
		y_max = round(yc+rc) 		
		Roi = []
		for xx=x_min, x_max do begin
		  for yy=y_min, y_max do begin
		      
		      if (xx-xc)^2+(yy-yc)^2 le rc*rc then begin
			index = yy*sx + xx
			Roi = [Roi, index]
		      endif
		      
		  endfor
		endfor
		
		if n_elements(Roi) gt 0 then begin
		s=size(im)
		openw, lun, edit_region, /append, /get_lun
		printf,lun,'# STV ROI output (image index): '+systime(0)
		printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

		if Roi[0] ge -1 then $
			for i=0L,n_elements(Roi)-1L do $
				printf,lun,Roi[i]

		endif ; Roi               
                free_lun,lun
                           
            endif else begin ; middle double click
                 tvcirc, ev.x, ev.y, reg_size, color=colordex('R')
            endelse   ; middle click
         
        endif
        
        
	if ev.type lt 2 then begin
 	    WIDGET_CONTROL, keyboard_text_id, /input_focus
	endif

ENDIF 

finish:

END

;--------------------------------------------------------------------

Pro Command_handler, Command

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, $
	     has_counter, last_command, reg_size, lock_reg, ellipse_glga
	     

eventchar = Command
last_command = ''

browsemess = 'In browse mode, set /nobrowse for output'
common_reg = aux_base + '_stv_roi_output.dat'

rute = strmid(imfil,0,strpos(imfil,'.fit'))
rute =  strmid(rute, 0, strlen(rute)-2)

change_region = (0 eq 1)
if not strcmp(strtrim(edit_region,2), strtrim(common_reg,2)) then change_region = (1 eq 1)

case eventchar of

    'a': begin
	    choi = -1
	    read,'0 - galex, 1 - sdss, 2 - 2mass, 3 - wise, 4 - panstarrs, 5 - acs: ',choi
	    stv_display_plots,choi
         end
    'U':if strpos(dtype,'sdss') ge 0 then begin
          if not file_exist(jpg_base+'_u.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_u.jpg'
	        if file_exist(rute+'_u.fits') then $ 
	           fspec  = rute+'_u.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then $
	               rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/uonly          
          endif 
          if file_exist(jpg_base+'_u.jpg') then begin 
	    jpg_file = jpg_base+'_u.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'u'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
        endif  
    'G': if strpos(dtype,'sdss') ge 0 or strpos(dtype,'panstarrs') ge 0 then begin
          if not file_exist(jpg_base+'_g.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_g.jpg'
	        if file_exist(rute+'_g.fits') then $ 
	           fspec  = rute+'_g.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
	               if strpos(dtype,'sdss') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/gonly  
	               if strpos(dtype,'panstarrs') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/gonly  
	               endif
          endif 
          if file_exist(jpg_base+'_g.jpg') then begin 
	    jpg_file = jpg_base+'_g.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'g'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif           
    'R': if strpos(dtype,'sdss') ge 0 or strpos(dtype,'panstarrs') ge 0  then begin
          if not file_exist(jpg_base+'_r.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_r.jpg'
	        if file_exist(rute+'_r.fits') then $ 
	           fspec  = rute+'_r.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
	               if strpos(dtype,'sdss') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ronly  
	               if strpos(dtype,'panstarrs') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ronly  
	               endif
          endif 
          if file_exist(jpg_base+'_r.jpg') then begin 
	    jpg_file = jpg_base+'_r.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'r'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif           
    'I': if strpos(dtype,'sdss') ge 0 or strpos(dtype,'panstarrs') ge 0  then begin
          if not file_exist(jpg_base+'_i.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_i.jpg'
	        if file_exist(rute+'_i.fits') then $ 
	           fspec  = rute+'_i.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
	               if strpos(dtype,'sdss') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ionly  
	               if strpos(dtype,'panstarrs') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/ionly  
	               endif
          endif 
          if file_exist(jpg_base+'_i.jpg') then begin 
	    jpg_file = jpg_base+'_i.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'i'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif           
    'Z': if strpos(dtype,'sdss') ge 0 or strpos(dtype,'panstarrs') ge 0  then begin
          if not file_exist(jpg_base+'_z.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_z.jpg'
	        if file_exist(rute+'_z.fits') then $ 
	           fspec  = rute+'_z.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
	               if strpos(dtype,'sdss') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/zonly  
	               if strpos(dtype,'panstarrs') ge 0 then rgb_sdss_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/zonly  
	               endif
          endif 
          if file_exist(jpg_base+'_z.jpg') then begin 
	    jpg_file = jpg_base+'_z.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'z'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif          
    'W': if strpos(dtype,'sdss') ge 0 then begin
          if file_exist(jpg_base+'_urz.jpg') then begin 
	  jpg_file = jpg_base+'_urz.jpg'
	  read_jpeg, jpg_file, tim, /true
	  tv_filter = 'w_com'
	  edit_region = common_reg
 	  stv_display
 	  stv_printhelpline
 	  endif else print, 'Warning: Could NOT find: ' + jpg_base+'_urz.jpg'
         endif              
    'N': if strpos(dtype,'sdss') ge 0 or strpos(dtype,'panstarrs') ge 0 then begin
          if file_exist(jpg_base+'_gri.jpg') then begin 
	  jpg_file = jpg_base+'_gri.jpg'
	  read_jpeg, jpg_file, tim, /true
	  tv_filter = 'n_com'
	  edit_region = common_reg
 	  stv_display 
 	  stv_printhelpline
 	  endif else print, 'Warning: Could NOT find: ' + jpg_base+'_gri.jpg'
         endif  
         
    'W1': if strpos(dtype,'wise') ge 0 then begin
          if not file_exist(jpg_base+'_w1.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_w1.jpg'
	        if file_exist(rute+'w1.fits') then begin
	           fspec  = rute+'w1.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
                       print, outdir
                       print, fspec
	               rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w1only,scales=[0.044,0.24,0.12],brite=5
	           endif
	        endif
          endif 
          if file_exist(jpg_base+'_w1.jpg') then begin 
	    jpg_file = jpg_base+'_w1.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'w1'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif  
         

    'W2': if strpos(dtype,'wise') ge 0 then begin
          if not file_exist(jpg_base+'_w2.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_w2.jpg'
	        if file_exist(rute+'w2.fits') then begin
	           fspec  = rute+'w2.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
                       print, outdir
                       print, fspec
	               rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w2only,scales=[0.044,0.24,0.12],brite=5
	           endif
	        endif
          endif 
          if file_exist(jpg_base+'_w2.jpg') then begin 
	    jpg_file = jpg_base+'_w2.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'w2'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif           
         
 
 
    'W3': if strpos(dtype,'wise') ge 0 then begin
          if not file_exist(jpg_base+'_w3.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_w3.jpg'
	        if file_exist(rute+'w3.fits') then begin
	           fspec  = rute+'w3.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
                       print, outdir
                       print, fspec
	               rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w3only,scales=[0.044,0.24,0.12],brite=5
	           endif
	        endif
          endif 
          if file_exist(jpg_base+'_w3.jpg') then begin 
	    jpg_file = jpg_base+'_w3.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'w3'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif   
 

    'W4': if strpos(dtype,'wise') ge 0 then begin
          if not file_exist(jpg_base+'_w4.jpg') then begin 
	        print, 'Warning: Could NOT find: ' + jpg_base+'_w4.jpg'
	        if file_exist(rute+'w4.fits') then begin
	           fspec  = rute+'w4.fits'
	           outdir = strmid(jpg_base,0,strpos(jpg_base,id))
	           q = DIALOG_MESSAGE('Could not find the jpg file, but the fits file is available !',  /INFORMATION,  /CENTER, TITLE='Missing JPG file !')
                   p = DIALOG_MESSAGE('Do you want to generate this jpg file?', /QUESTION,  /CENTER, TITLE='Regenerate JPG file?')
                   if  strcmp(p, 'Yes') then begin
                       print, outdir
                       print, fspec
	               rgb_wise_asinh,fspec=fspec,outdir=outdir,update=update,verbose=verbose,/w4only,scales=[0.03,0.24,0.12],brite=3
	           endif
	        endif
          endif 
          if file_exist(jpg_base+'_w4.jpg') then begin 
	    jpg_file = jpg_base+'_w4.jpg'
	    read_jpeg, jpg_file, tim, /true
	    tv_filter = 'w4'
	    filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
	    if change_region then edit_region=filter_reg
	    stv_display 
	    stv_printhelpline
	   endif 
         endif    
 
    'W123': if strpos(dtype,'wise') ge 0 then begin
          if file_exist(jpg_base+'_w123.jpg') then begin 
	  jpg_file = jpg_base+'_w123.jpg'
	  read_jpeg, jpg_file, tim, /true
	  tv_filter = 'n_com'
	  edit_region = common_reg
 	  stv_display 
 	  stv_printhelpline
 	  endif else print, 'Warning: Could NOT find: ' + jpg_base+'_w123.jpg'
         endif  

    'W124': if strpos(dtype,'wise') ge 0 then begin
          if file_exist(jpg_base+'_w124.jpg') then begin 
	  jpg_file = jpg_base+'_w124.jpg'
	  read_jpeg, jpg_file, tim, /true
	  tv_filter = 'w_com'
	  edit_region = common_reg
 	  stv_display 
 	  stv_printhelpline
 	  endif else print, 'Warning: Could NOT find: ' + jpg_base+'_w124.jpg'
         endif          
         
         
    '9': begin
         has_counter =  (1 eq 1)
         contour_id =  WIDGET_INFO(sbase, FIND_BY_UNAME='contour')
         WIDGET_CONTROL, contour_id, SET_VALUE='Contours On'
         ds9_image
         end 
    'c': if browse eq 0 then begin
                last_command = Command
         	stvcntrd
         	stv_printhelpline
	 endif else print,browsemess

    'd': if browse eq 0 then begin
            last_command = Command
	    stv_delete
	    stv_printhelpline
	 endif else print,browsemess
    'e': if browse eq 0 then $
		    stv_ellipse $
	 else    print,browsemess
    'e+': if browse eq 0 then $
		    stv_ellipse_size, 'e+' $
	 else    print,browsemess
    'e-': if browse eq 0 then $
		   stv_ellipse_size, 'e-' $
	 else    print,browsemess	 
    'f': if browse eq 0 then $
		    stv_find $
	 else    print,browsemess
    'g': if browse eq 0 then $
		    stv_find,/galaxy $
	 else    print,browsemess
    'h': stv_printhelpline
    '?': stv_printhelpline
    'i': invert_colortable
    'k': begin
	 for i=0l,maxp-1l do $
		 if psrc[6,i] ge 2. then $
		 	psrc[6,i] = 0.
	 stv_display
	 end
    'l': begin
         tvlist, im
         stv_printhelpline
         end
    'm': if browse eq 0 then begin
	 	misc = 'RUN'
	 	stv_shutdown
	 endif else print,browsemess
    'o': if browse eq 0 then $
		    stv_find,/outside $
	 else    print,browsemess
    'p': stv_print_qa_flags,/help,/values,/note,/details
    'q': stv_shutdown
    'Q': begin
         misc = 'STOP'
         stv_shutdown
         end
    '$': if browse eq 0 then begin
                   stv_change_edit_region
                   stv_printhelpline
         endif else print,browsemess
    'r': if browse eq 0 then begin
                    last_command = Command
		    stv_roi
         	    stv_printhelpline
	 endif else print,browsemess
    'C': if browse eq 0 then begin
		    stv_clear_roi
         	    stv_printhelpline
	 endif else print,browsemess	
	 
    'save_reg': if browse eq 0 then begin
		    stv_save_roi
         	    stv_printhelpline
	 endif else print,browsemess		 
    '!': if browse eq 0 then begin
                    last_command = Command
		    rmvrgn
         	    stv_printhelpline
	 endif else print,browsemess
    '~': if browse eq 0 then begin
		    stv_undo_roi
         	    stv_printhelpline
	 endif else print,browsemess
    '^': if browse eq 0 then begin
                    print, ''
                    print, 'Refreshing STV ....'
                    print, ''
		    stv_display
         	    stv_printhelpline         	    
	 endif else print,browsemess	 
    'b': if browse eq 0 then begin
                    last_command = Command
		    stv_reg_creator
         	    stv_printhelpline
         endif else print,browsemess	 	    
    '@': if browse eq 0 then begin
		    stv_remove_background
         	    stv_printhelpline
	 endif else print,browsemess	
	 
    '*': if browse eq 0 then begin
		    stv_grab_region, action='update'
         	    stv_printhelpline
	 endif else print,browsemess	
    '#': if browse eq 0 then begin
		    stv_grab_region, action='undo'
         	    stv_printhelpline
         endif else print,browsemess	    
    'ell_sdss': if browse eq 0 then begin
		    stv_grab_region, action='sdss'
         	    stv_printhelpline         	    
         endif else print,browsemess	
    'ell_catal': if browse eq 0 then begin
		    stv_grab_region, action='ell_catal'
         	    stv_printhelpline         	    
         endif else print,browsemess	         
    'ell_save': if browse eq 0 then begin
		    stv_grab_region, action='save_ellipse'
         	    stv_printhelpline         	    
         endif else print,browsemess	         
    '=': if browse eq 0 then begin
		    stv_grab_region, action='to_ds9'
         	    stv_printhelpline         	    
	 endif else print,browsemess	
    'sex': if browse eq 0 then begin
		    run_sextract
         	    stv_printhelpline         	    
	 endif else print,browsemess	
    's': if plog then begin
		    plog = (1 eq 0)
		    print,'De-selecting : ',id
	 endif else begin
		    plog = (1 eq 1)
		    print,'Selecting    : ',id
	 endelse
    'v': begin
         IF n_elements(hdr) NE 1 THEN stv_curval,hdr,im,filename=cfile ELSE $
           stv_curval,im, filename=cfile
         stv_printhelpline
         end
    'x': begin
         
         p = DIALOG_MESSAGE('Do you really want to SKIP?', /QUESTION,  /CENTER, /DEFAULT_NO)
         if strcmp(p, 'Yes') then begin
	   misc = 'SKIP'
	   stv_shutdown
	 endif
	 end
    'z': begin
         zoom, /cont;, /keep
         stv_printhelpline
         end
    'previous': begin
         
         p = DIALOG_MESSAGE('Do you really want to go back?', /QUESTION,  /CENTER, /DEFAULT_NO)
         if strcmp(p, 'Yes') then begin
	   misc = 'PREVIOUS'
	   stv_shutdown
	 endif
	 end  
    'save': begin
	   misc = 'SAVE'
	   stv_shutdown
	 end   	 
    else:  ;any other key press does nothing
ENDCASE


if strcmp(eventchar, 'lock') then begin
   
   bol = 1
   lock_info =  WIDGET_INFO(sbase, FIND_BY_UNAME='lock_info')
   WIDGET_CONTROL, lock_info, SET_VALUE='unlock'
   print, 'Press right mouse mouse button to unlock ...'

   while bol eq 1 do begin
       !MOUSE.BUTTON = 0
       cursor,x,y,2,/DEVICE,/CHANGE 
       cr_err = !MOUSE.BUTTON
       
       stv_update_info_widgets, x, y
       
       
       if cr_err EQ 4 then begin 
          bol = 0
          WIDGET_CONTROL, lock_info, SET_VALUE='lock'
       endif  
   endwhile  
    
endif


if strcmp(eventchar, 'clr') then begin
    
  spawn, 'xpaaccess ds9', results
  if strcmp(results[0], 'yes') then begin
     spawn, "xpaset -p ds9 regions delete all"
  endif
  

endif


if strcmp(eventchar, 'import_reg') then begin
    
   spawn, 'xpaaccess ds9', results
   if strcmp(results[0], 'yes') then begin 
       
       spawn, "xpaset -p ds9 regions system  image"
       spawn, "xpaset -p ds9 regions format ds9"
       spawn, "xpaset -p ds9 regions save idl_foo.reg"
   
     
      circles = read_ds9_circle('idl_foo.reg')
      boxes   = read_ds9_box('idl_foo.reg')
      sz = size(im,/dim)
      sx = sz[0]
      sy = sz[1]
      
      
      
    n_circles = n_elements(circles)/3
    n_boxes   = n_elements(boxes)/5
    
    s=size(im)
    openw, lun, edit_region, /append, /get_lun
    
    if n_circles ge 1 then begin 
    
	for n=0, n_circles-1 do begin
	      xc = circles[0+3*n]
	      yc = circles[1+3*n]
	      rc = circles[2+3*n]
	      
	      x_min = round(xc-rc) 
	      x_max = round(xc+rc) 
	      y_min = round(yc-rc) 
	      y_max = round(yc+rc) 
	      Roi = []
	      for xx=x_min, x_max do begin
		for yy=y_min, y_max do begin
		    
		    if (xx-xc)^2+(yy-yc)^2 le rc*rc then begin
		      index = yy*sx + xx
		      Roi = [Roi, index]
		    endif
		    
		endfor
	      endfor
	  if n_elements(Roi) gt 0 then begin

	  printf,lun,'# STV ROI output (image index): '+systime(0)
	  printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

	  if Roi[0] ge -1 then $
		  for i=0L,n_elements(Roi)-1L do $
			  printf,lun,Roi[i]

          endif ; Roi
	endfor
    endif  ; n_circles ge 1
    
    if n_boxes ge 1 then begin 

	for n=0, n_boxes-1 do begin
	      xc = boxes[0+5*n]
	      yc = boxes[1+5*n]
	      Lx = boxes[2+5*n]
	      Ly = boxes[3+5*n]
	      theta = boxes[4+5*n]*!PI/180.
   
	      x_min = floor(xc-Lx*0.5) 
	      x_max = ceil(xc+Lx*0.5) 
	      y_min = floor(yc-Ly*0.5) 
	      y_max = ceil(yc+Ly*0.5) 
	      print, x_min, x_max, y_min, y_max
	      print, xc, yc, Lx, Ly
	      print, '----------------'
	      x1 = (x_min-xc)*cos(theta)-(y_min-yc)*sin(theta)+xc
	      y1 = (x_min-xc)*sin(theta)+(y_min-yc)*cos(theta)+yc
	      
	      x2 = (x_max-xc)*cos(theta)-(y_min-yc)*sin(theta)+xc
	      y2 = (x_max-xc)*sin(theta)+(y_min-yc)*cos(theta)+yc
	      
	      x3 = (x_max-xc)*cos(theta)-(y_max-yc)*sin(theta)+xc
	      y3 = (x_max-xc)*sin(theta)+(y_max-yc)*cos(theta)+yc
	      
	      x4 = (x_min-xc)*cos(theta)-(y_max-yc)*sin(theta)+xc
	      y4 = (x_min-xc)*sin(theta)+(y_max-yc)*cos(theta)+yc
	      
	      m1 = (y2-y1)/(x2-x1)
	      m2 = (y3-y2)/(x3-x2)
	      m3 = (y3-y4)/(x3-x4)
	      m4 = (y4-y1)/(x4-x1)
	      if (x2-x1)*(x3-x2)*(x3-x4)*(x4-x1) ne 0 then begin 

	      x_min = min([x1,x2,x3,x4])
	      x_max = max([x1,x2,x3,x4])
	      y_min = min([y1,y2,y3,y4])
	      y_max = max([y1,y2,y3,y4])

	      
	      Roi = []
	      for xx=x_min, x_max do begin
		for yy=y_min, y_max do begin
		    
		    c = 0
		    c+= sign(yy - (m1*(xx-x1)+y1))
		    c+= sign(yy - (m2*(xx-x2)+y2))
		    c+= sign(yy - (m3*(xx-x4)+y4))
		    c+= sign(yy - (m4*(xx-x1)+y1))
		    
		    if c eq 0 then begin 
         		    index = round(yy)*sx + round(xx)
	         	    Roi = [Roi, index]
		    endif 
		    
		endfor
	      endfor
	    if n_elements(Roi) gt 0 then begin

	  printf,lun,'# STV ROI output (image index): '+systime(0)
	  printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

	if Roi[0] ge -1 then $
		  for i=0L,n_elements(Roi)-1L do $
			  printf,lun,Roi[i]

        endif ; Roi
	endif else begin
	      Roi = [] 
	      for xx=x_min, x_max do begin
		for yy=y_min, y_max do begin
         	  index = round(yy)*sx + round(xx)
	          Roi = [Roi, index]
		endfor
	      endfor  	  
	      
	  if n_elements(Roi) gt 0 then begin

	    printf,lun,'# STV ROI output (image index): '+systime(0)
	    printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

	    if Roi[0] ge -1 then $
		for i=0L,n_elements(Roi)-1L do $
		  printf,lun,Roi[i]

          endif ; Roi
	
	endelse 


  endfor ; for n_boxes
  endif  ; n_boxes ge 1    
  free_lun,lun
  stv_display 
  endif ; if there is an open ds9

endif ; import_reg


if strcmp(eventchar, 'add_centroid') then begin
   
   spawn, 'xpaaccess ds9', results
   if strcmp(results[0], 'yes') then begin 
   reg_file = 'idl_foo.reg'
   openw,  lun,reg_file,/get_lun
   printf, lun, 'image'
   
   sz = size(psrc,/dim)
   sx = sz[0]
   sy = sz[1]
   print, sx,sy 
   
   for w=0, sy-1 do begin
      circle = 'circle('
      circle = circle + strtrim(string(psrc[2,w]),2)  + ','
      circle = circle + strtrim(string(psrc[3,w]),2)  + ','
      circle = circle + strtrim(string(psrc[4,w]),2)
      circle = circle + ') # color=red'
      
      ; if the centroid is currentyl turned-on
      if psrc[6,w] eq 1 then printf, lun, circle
      
   
   endfor
   free_lun,lun  
   spawn, "xpaset -p ds9 regions load idl_foo.reg"
   
   endif 
   

endif



if strcmp(eventchar, '\') or strcmp(eventchar, '9') then begin
      contour_id =  WIDGET_INFO(sbase, FIND_BY_UNAME='contour')
      if not has_counter then begin; contour is NOT activated

          if strcmp(dtype,'sdss') then begin
            spawn, "xpaset -p ds9 contour nlevels 30"
            spawn, "xpaset -p ds9 contour smooth 10"
          endif else begin
            spawn, "xpaset -p ds9 contour nlevels 50"
            spawn, "xpaset -p ds9 contour smooth 5"
          endelse
           spawn, "xpaset -p ds9 contour scale log"
           spawn, "xpaset -p ds9 contour generate"
           spawn, "xpaset -p ds9 contour color green"
           spawn, "xpaset -p ds9 contour"
           WIDGET_CONTROL, contour_id, SET_VALUE='Contours Off'
           has_counter =  (1 eq 1) ; True
      endif else begin  ; contour is activated
           spawn, "xpaset -p ds9 contour clear"
           spawn, "xpaset -p ds9 contour close"
	   if strcmp(eventchar, '\') then WIDGET_CONTROL, contour_id, SET_VALUE='Contours On'
	   has_counter =  (1 eq 0) ; False
      endelse
endif

if strcmp(eventchar, '9\') then begin
  spawn, 'xpaaccess ds9', results
  if strcmp(results[0], 'yes') then spawn, 'xpaset -p ds9 exit'
  Command_handler, '9'
  wait, 1.5
  Command_handler, '\'

endif


spawn, 'xpaaccess ds9', results ; check if ds9 is open
if strcmp(eventchar, 'c_ds9') and strcmp(results[0], 'yes') then begin
  p = DIALOG_MESSAGE('Do you really want to close all ds9 windows?', /QUESTION,  /CENTER, /DEFAULT_NO)  
  if strcmp(p, 'Yes') then spawn, 'xpaset -p ds9 exit'
endif




END

;--------------------------------------------------------------------
PRO stv_update_info_widgets, x, y

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
	     
	     
	     WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='cor_x'), set_value='X: '+ strtrim(string(x),2)
             WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='cor_y'), set_value='Y: '+ strtrim(string(y),2)
	     
	     if have_astr  then begin
			XY2AD, x, y, astr, ra, dec ; convert ra,dec to x,y px
			radec, ra,dec,h,m,s,dd,mm,ss
			ra_hex = strtrim(string(h,format='(I02)'),2)+':'+strtrim(string(m,format='(I02)'),2)+':'+strtrim(string(s,format='(f06.3)'),2)
			dec_hex = strtrim(string(dd,format='(I03)'),2)+':'+strtrim(string(mm,format='(I02)'),2)+':'+strtrim(string(ss,format='(f06.3)'),2)

			WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='cor_ra'), set_value='RA : '+ strtrim(string(ra,format='(f9.5)'),2)
			WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='cor_ra0'), set_value=ra_hex
			
			WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='cor_dec'), set_value='DEC: '+ strtrim(string(dec,format='(f9.5)'),2)
			WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='cor_dec0'), set_value=dec_hex
			
			x = round(x)  & y = round(y)
			value = im[x,y]
			WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='im_value'), set_value='Value: '+ strtrim(string(value),2)
			
			bscale = sxpar(hdr,'BSCALE')

			if (bscale ne 0)  then begin 
			    bzero = sxpar(hdr,'BZERO')
			    flux = bscale*value + bzero  
			    WIDGET_CONTROL, WIDGET_INFO(sbase, FIND_BY_UNAME='im_flux'), set_value='Flux: '+ strtrim(string(value),2)
		       endif	     
	     
; 			print , 'RA : '+ strtrim(string(ra,format='(f9.5)'),2)
; 			print , ra_hex
; 			print , 'DEC: '+ strtrim(string(dec,format='(f9.5)'),2)
; 			print , dec_hex
	     endif


END
;--------------------------------------------------------------------

PRO panel_event, event

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
Common indices, last_band_ind, last_tv_filter


WIDGET_CONTROL, keyboard_text_id, event_pro = 'stv_keyboard_event'
stv_printhelpline
if (el[0] ne elorig[0] or el[1] ne elorig[1] or $
    el[2] ne elorig[2] or el[3] ne elorig[3] or $
    el[4] ne elorig[4]) then begin

 if have_astr then begin
	 XY2AD,el[2],el[3],astr,el_a,el_d
	 el_ast = 1
 endif else begin
	 as_pix = 1.
	 el_a = el[2]
	 el_d = el[3]
	 el_ast = 0
 endelse
 ;
 ; make sure pa is positive
 while el[4] lt 0. do el[4] = el[4] + 180.

 openw,lun,aux_base+'_stv_ellipse_output.dat',/get_lun
 printf,lun,'# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix ): '+systime(0)
 printf, lun, el[0]*2.*as_pix, el[1]*2.*as_pix, el_a, el_d, el[4], $
	 el[0]*2, el[1]*2, el[2], el[3], el_ast, as_pix, $
   format = '(2f9.1,2f13.8,f7.1,4f9.2,i5,f9.5)'
 free_lun,lun
 
endif


values_A = ['none','none'] ; for later development for other surveys


if strcmp(strtrim(dtype), 'sdss') then $
  values_A = ['U', 'G', 'R', 'I', 'Z', 'N', 'W'] 
if strcmp(strtrim(dtype), 'panstarrs') then $
  values_A = ['G', 'R', 'I', 'Z', 'N'] 
if strcmp(strtrim(dtype), 'wise') then $
  values_A = ['W1', 'W2', 'W3', 'W4', 'W123', 'W124']   
  
Bgroup_A = 'filter'


WIDGET_CONTROL, event.id, get_UVALUE = UVALUE

if strcmp(UVALUE, Bgroup_A) and not Event.select then begin 

last_band_ind = Event.value
last_tv_filter = tv_filter


endif


if strcmp(UVALUE, Bgroup_A) and Event.select then begin 
  
  i = Event.value
  Command_handler, values_A[i]
  
  if strcmp(last_tv_filter, tv_filter) then $
       WIDGET_CONTROL, event.id, SET_VALUE=last_band_ind

       
  Return
  endif

 
if strcmp(UVALUE, 'reg_circle')  and Event.select then begin
   
  sizes = [3,5,7,10,15,20,25,30,35] 
  i = Event.value
  reg_size = sizes[i]

return 
endif



if strcmp(UVALUE, 'lock_reg') then begin
  
  lock_button =  WIDGET_INFO(sbase, FIND_BY_UNAME='lockreg')
  if lock_reg eq 0 then begin
      WIDGET_CONTROL, lock_button, SET_VALUE='unlock'
      lock_reg = 1
  endif else begin
      lock_reg = 0
      WIDGET_CONTROL, lock_button, SET_VALUE='lock reg'
  endelse
  

return 
endif

; run the command
Command_handler, UVALUE   
  

  
RETURN
END

;--------------------------------------------------------------------

PRO stv_keyboard_event, event

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

eventchar = string(event.ch)


Command_handler, eventchar

if (xregistered('stv', /noshow)) then $
  widget_control, keyboard_text_id, /clear_events

END

;--------------------------------------------------------------------

PRO stv_printhelpline

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

if browse eq 0 then begin
  print, '--------------------------------------------------------------------'
  print, 'STV: (z)oom|(l)ist|cur(v)al|(d)elete|(c)entroid|(e)llipse|(r)egion|'
  print, '     (s)elect host|ds(9)|(p)rint qa flags|e(x)it-no output|(m)easure|'
  print, '     (f)ind|(g)alaxy|(o)utside|(k)ill|(h)elp|(q)uit-next|(Q)uit-exit'
  print, '--------------------------------------------------------------------'
  print, '     (b)ackground region|(@) remove background region |'
  print, '     (~) Undo the last region | (!) Remove a region | (^) Refresh STV |'
  print, '     (C)lear all masks | ($) Change region of interest'
  print, '     (*) import ellipse from ds9 | (#) back to the original ellipse |'
  print, '     (=) Export ellipse to ds9 | '
  print, '--------------------------------------------------------------------'
  
  
  
  if strcmp(strtrim(tv_filter,2), 'w_com') or strcmp(strtrim(tv_filter,2), 'n_com') then begin 
     foi = 'Composite Image' 
     editable_reg = 'All Bands'
  endif else begin 
     foi = tv_filter
     common_reg = aux_base + '_stv_roi_output.dat'
     filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
     if strcmp(strtrim(edit_region,2), strtrim(common_reg,2)) then begin
         editable_reg = 'All Bands'
     endif else editable_reg = 'Filter ' + tv_filter
  endelse   
  
  print, ''
  print, 'STV Image       : ', foi
  print, 'Editable Region : ', editable_reg
  label_id =  WIDGET_INFO(sbase, FIND_BY_UNAME='label')
  WIDGET_CONTROL, label_id, SET_VALUE='Current: '+editable_reg
  
  
  if strpos(dtype,'sdss') ge 0 then begin
      print, ''
      print, '--------------------------------------------------------------------'
      print, 'To change the display image (SDSS only)'
      print, '     (U,G,R,I,Z) sdss bands | (W) urz | (N) gri'
      print, '--------------------------------------------------------------------'
      print, ''
  endif
  
  if strpos(dtype,'panstarrs') ge 0 then begin
      print, ''
      print, '--------------------------------------------------------------------'
      print, 'To change the display image (Pan-STARRS only)'
      print, '     (G,R,I,Z) sdss bands | (N) gri'
      print, '--------------------------------------------------------------------'
      print, ''
  endif  

  endif else begin
  stv_print_qa_flags,/help,/values,/note,/details
  print, 'STV: (z)oom|(l)ist|cur(v)al|(s)elect host|ds(9)|(p)rint qa flags|(h)elp|(q)uit-next|(Q)uit-exit'
endelse

END

;--------------------------------------------------------------------
;--------------------------------------------------------------------

function get_roi_sep,im,rfile

out_lst = list()

array = ''
line = ''
; file = 'stv_background_vertices.dat'

OPENR, lun, rfile, /GET_LUN
WHILE NOT EOF(lun) DO BEGIN
  
  indices = []
  READF, lun, line
  first = strmid(line, 0, 1)
  while strcmp(first , '#') and NOT EOF(lun) do begin
      if EOF(lun) then break
      READF, lun, line
      first = strmid(line, 0, 1)
      
  endwhile
  while first ne '#' and NOT EOF(lun) do begin
       
       S = STRSPLIT(line,' ', /extract)
       indices = [indices, long(S[0])]
       if EOF(lun) then break
       READF, lun, line
       first = strmid(line, 0, 1)
  endwhile  

  
  if n_elements(indices) gt 1 then out_lst.add, indices
ENDWHILE 


free_lun,lun
return, out_lst

end
;--------------------------------------------------------------------

function undo_stv_region, filename

array = ''
line = ''

sharps = [0]


findpro, 'glga_qa', DirList=DirList, ProList=ProList, /Noprint
Undo_Script = GETENV('IDL_GLGAPRO') + '/scripts/undo_region.bash'

if file_test(filename) then begin
   if file_test(Undo_Script) then begin
       spawn, 'bash '+ Undo_Script + ' ' + filename
   endif else begin
        
        print, ' '
        print, 'Warning ... '
        print, 'Could not find '+ Undo_Script
        print, 'Undoing may takes a while depending on how big is the region file ...'
        print, ' '
        
	OPENR, lun, filename, /GET_LUN

	WHILE NOT EOF(lun) DO BEGIN & $
	  READF, lun, line & $
	  array = [array, line] & $
	ENDWHILE
	free_lun,lun


	i = 0
	while i lt n_elements(array) do begin
	  first = strmid(array[i], 0, 1)
	  if strcmp(first , '#') then sharps=[sharps,i]
	  i++
	endwhile

	N = n_elements(sharps)

	if N eq 3 then begin
	  spawn, 'rm '+filename+' '
	  print, 'removed completely ...'
	endif else begin
	  spawn, 'rm '+filename+' '
	  openw, lun, filename,/append,/get_lun
	  print, sharps[N-2], n_elements(array)
	  for p=0, sharps[N-2]-1 do begin
	      print, p
	      printf, lun, array[p]
	  endfor
	  free_lun,lun
	  print, 'removed the last region ...'
	endelse
	  
   endelse

endif else begin
    print, ' '
    print, ' '
    print, '*****************'
    print, 'No new region file is available .... '
    print, '*****************'
    print, ' '
    print, ' '
endelse

return, sharps

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO invert_colortable

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


TVLCT, R, G, B, /get
tvlct,reverse(r),reverse(g),reverse(b)
IF !d.n_colors Gt 256 THEN stv_display ;stv, im, hd=hdr ,/same

END

;--------------------------------------------------------------------

PRO stv_setflags_keyboard_event, event

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

common stvc, flag_w, pnale_w

eventchar = strlowcase(string(event.ch))

cancel = 0
exit = 0

case eventchar of
    'f' : if qa.complete eq 0 then qa.complete = 1 else qa.complete = 0
    'm' : if qa.multiple eq 0 then qa.multiple = 1 else qa.multiple = 0
    'b' : if qa.bright_star eq 0 then qa.bright_star = 1 else qa.bright_star = 0
    'u' : if qa.uncertain_mask eq 0 then qa.uncertain_mask = 1 else qa.uncertain_mask = 0
    'v' : if qa.fov_expand eq 0 then qa.fov_expand = 1 else qa.fov_expand = 0
    '1' : if qa.band1_edge eq 0 then qa.band1_edge = 1 else qa.band1_edge = 0
    '2' : if qa.band1_sn_grad eq 0 then qa.band1_sn_grad = 1 else qa.band1_sn_grad = 0
    '3' : if qa.band1_artifact eq 0 then qa.band1_artifact = 1 else qa.band1_artifact = 0
    '4' : if qa.band1_missing eq 0 then qa.band1_missing = 1 else qa.band1_missing = 0
    '5' : if qa.band1_other eq 0 then qa.band1_other = 1 else qa.band1_other = 0
    '6' : if qa.band2_edge eq 0 then qa.band2_edge = 1 else qa.band2_edge = 0
    '7' : if qa.band2_sn_grad eq 0 then qa.band2_sn_grad = 1 else qa.band2_sn_grad = 0
    '8' : if qa.band2_artifact eq 0 then qa.band2_artifact = 1 else qa.band2_artifact = 0
    '9' : if qa.band2_missing eq 0 then qa.band2_missing = 1 else qa.band2_missing = 0
    '0' : if qa.band2_other eq 0 then qa.band2_other = 1 else qa.band2_other = 0
    '+' : if qa.band3_edge eq 0 then qa.band3_edge = 1 else qa.band3_edge = 0
    '-' : if qa.band3_sn_grad eq 0 then qa.band3_sn_grad = 1 else qa.band3_sn_grad = 0
    '*' : if qa.band3_artifact eq 0 then qa.band3_artifact = 1 else qa.band3_artifact = 0
    '=' : if qa.band3_other eq 0 then qa.band3_other = 1 else qa.band3_other = 0
    '/' : if qa.band3_missing eq 0 then qa.band3_missing = 1 else qa.band3_missing = 0
    '?' : stv_print_qa_flags,/help
    'q' : begin
          if (xregistered ('stv')) then widget_control, sbase, /destroy
          if (xregistered ('flag_base')) then begin
             
             widget_control, flag_w, /destroy
             endif
	      cm=''
	      read,'Enter QA note: ',cm
	      if strlen(cm) gt 0 then $
	      qa.note = cm
          end
    'c' : begin
             
             WIDGET_CONTROL, keyboard_text_id, event_pro = $
		'stv_keyboard_event'
             WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col_widg'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col2_widg'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='import'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='export'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='sextract'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col4_widg'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='wLRow1_widg'), sensitive=1
	     if (xregistered ('flag_base')) then  widget_control, flag_w, /destroy
	     stv_printhelpline
	     goto, finish
    
          end
    else:  ;any other key press does nothing
ENDCASE

stv_print_qa_flags,/values

if (xregistered('stv', /noshow)) then $
  widget_control, keyboard_text_id, /clear_events

  

if (xregistered ('flag_base')) then begin  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='flag_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.complete,qa.uncertain_mask,qa.fov_expand,qa.multiple,qa.bright_star]
  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='band1_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.band1_edge, qa.band1_sn_grad, qa.band1_artifact, qa.band1_missing, qa.band1_other]
  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='band2_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.band2_edge, qa.band2_sn_grad, qa.band2_artifact, qa.band2_missing, qa.band2_other]  
  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='band3_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.band3_edge, qa.band3_sn_grad, qa.band3_artifact, qa.band3_missing, qa.band3_other]  
endif  

finish:
  
END

;--------------------------------------------------------------------

PRO stv_print_qa_flags,help=help,values=values,note=note,details=details

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

if keyword_set(details) then begin
	print,' '
	print,'User     : ',qa.user_name
	print,'Machine  : ',qa.machine_name
	print,'# of runs: ',qa.nqa
endif
if keyword_set(note) and strlen(qa.note) gt 0 then begin
	print,' '
	print,'QA NOTE: ', qa.note
endif
if keyword_set(help) then begin
	print,' '
	print,'FLAGS: (F)inished | (U)ncertain | FO(V)expand | (M)ultiple| (B)rightstar'
	print,'BAND1[FUV,g,j]: (1)edge | (2)sn grad | (3)artifact | (4)missing | (5)other'
	print,'BAND2[NUV,r h]: (6)edge | (7)sn grad | (8)artifact | (9)missing | (0)other'
	print,'BAND3    [i,k]: (+)edge | (-)sn grad | (*)artifact | (/)missing | (=)other'
	print,'(c) to cancel, (q) to exit, (?) for this help printout'
	print,'  F  U  V  M  B  1  2  3  4  5  6  7  8  9  0  +  -  *  /  ='
endif
if keyword_set(values) then begin
  print,qa.complete,qa.uncertain_mask,qa.fov_expand,qa.multiple,qa.bright_star,$
	qa.band1_edge,qa.band1_sn_grad,qa.band1_artifact,qa.band1_missing, $  
	qa.band1_other, $
	qa.band2_edge,qa.band2_sn_grad,qa.band2_artifact,qa.band2_missing, $  
	qa.band2_other, $
	qa.band3_edge,qa.band3_sn_grad,qa.band3_artifact,qa.band3_missing, $  
	qa.band3_other, $
	format='(20i3)'
endif

END

;--------------------------------------------------------------------
PRO stv_ellipse_size, command


COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
    
    ratio = el[0]/el[1]
    case command of
    'e+' : begin 
         el[0] *= 1.1
         el[1] = el[0]/ratio      
         end
    'e-' : begin ; 
         el[0] /= 1.1
         el[1] = el[0]/ratio  
         end 

ENDCASE


stv_display
tvellipse,el[0],el[1],el[2],el[3],float(90.+el[4]),$
   color=colordex('R'),linestyle=1
plots,el[5],el[6],psym=4,symsi=3.,color=colordex('G'),/device
plots,el[2],el[3],psym=1,symsi=2.,color=colordex('R'),/device

minskywid = 27. / as_pix   ; 5 for ACS ; 27. otherwise
factor=1.5
skyinner=factor*el[0]
skyouter= el[0]*sqrt(1+factor^2) > (skyinner+minskywid) ; ACS (skyinner+minskywid) ; el[0]*sqrt(1+factor^2) > (skyinner+minskywid)

tvellipse,skyinner,skyinner/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
   color=colordex(anulus_col),linestyle=1
tvellipse,skyouter,skyouter/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
   color=colordex(anulus_col),linestyle=1


if  (el[0] ne elorig[0] or el[1] ne elorig[1] or $
    el[2] ne elorig[2] or el[3] ne elorig[3] or $
    el[4] ne elorig[4]) then begin

 if have_astr then begin
	 XY2AD,el[2],el[3],astr,el_a,el_d
	 el_ast = 1
 endif else begin
	 as_pix = 1.
	 el_a = el[2]
	 el_d = el[3]
	 el_ast = 0
 endelse
 ;
 ; make sure pa is positive
 while el[4] lt 0. do el[4] = el[4] + 180.

 openw,lun,aux_base+'_stv_ellipse_output.dat',/get_lun
 printf,lun,'# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix ): '+systime(0)
 printf, lun, el[0]*2.*as_pix, el[1]*2.*as_pix, el_a, el_d, el[4], $
	 el[0]*2, el[1]*2, el[2], el[3], el_ast, as_pix, $
   format = '(2f9.1,2f13.8,f7.1,4f9.2,i5,f9.5)'
 free_lun,lun
 
endif



END


;--------------------------------------------------------------------


PRO stv_ellipse_keyboard_event, event

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


eventchar = string(event.ch)

cancel = 0
exit = 0

case eventchar of
    'z' : begin ;dec ellipse
         ratio = el[0]/el[1]
         el[0] = el[0]-1 > (el[0]-1)/ratio > 0 
         el[1] = el[0]/ratio        
         end
    'x' : begin ;inc ellipse
         ratio = el[0]/el[1]
         el[0] = el[0]+1 > el[1] > 0 
         el[1] = el[0]/ratio        
         end
    '%' : begin ; multiply ellipse size by percentage
         ratio = el[0]/el[1]
	 pc = 100.
	 read,'Enter % ( > 100 to expand, < 100 to contract ): ',pc
	 if abs(pc) le 5. then pc=pc*100.  ; if user make a mistake and forgets about the % sign
	 if pc lt 0. then $
		 el[0] = el[0] / ( abs(pc) / 100. ) $
	 else	 el[0] = el[0] * ( abs(pc) / 100. )
         el[1] = el[0]/ratio        
         end
    '/' : begin ; set a/b axial ratio
         print,'  Current b/a: ',el[1]/el[0],form='(a,f7.3)'
	 read,'Enter new b/a: ',ba
	 if ba gt 0. and ba le 1. then $
		 el[1] = el[0] * ba $
	 else	 print,'b/a must > 0. and <= 1.'
         end
    ;inc or dec major diam
    "'" : el[0] = el[0]+1 > el[1] > 0;'
    ";" : el[0] = el[0]-1 > el[1] > 0;;
    ;inc or dec min diam
    ']' : el[1] = el[1]+1 < el[0] > 0;[
    '[' : el[1] = el[1]-1 < el[0] > 0;]
    ;move center
    'w': el[3] = el[3]+1 ;up arrow
    's': el[3] = el[3]-1 ;down arrow
    'a': el[2] = el[2]-1 ;left arrow
    'd': el[2] = el[2]+1 ;right arrow
    ;change PA
    '<' : el[4] =  (el[4] +1) < 360;rotate pa ccw <
    '>' : el[4] =  (el[4] -1) >   0;rotate pa cw  >
    ',' : el[4] =  (el[4] +1) < 360;rotate pa ccw ,
    '.' : el[4] =  (el[4] -1) >   0;rotate pa cw  .
    ;print, exit, cancel, help
    'p': begin
         print, 'semimajor[pix], semiminor[pix], xcen, ycen, pa, xnom, ynom' 
         print, el
         print, "Ellipse: z|x|%|/|;|'|[|]|w|s|a|d|<|>|(r)e-center|(p)rint|(h)elp|(c)ancel|(q)uit"
         end
    'q': begin
         exit = 1      
         WIDGET_CONTROL, keyboard_text_id, event_pro = 'stv_keyboard_event'
         stv_printhelpline
         end
    'c': begin
         cancel = 1     
         WIDGET_CONTROL, keyboard_text_id, event_pro = 'stv_keyboard_event'
         stv_printhelpline
         end
    'r': begin
         print,'Mark the center of galaxy.'
         cursor,xx,yy,3,/dev 
         ;IF !mouse.button EQ 4 THEN GOTO, finish
         delpos = 7.0	; pixels
         if strpos(dtype,'sdss') ge 0 then delpos = 6.0
         if strpos(dtype,'panstarrs') ge 0 then delpos = 12.0
         if strpos(dtype,'acs') ge 0 then delpos = 60.0
         sz = size(im,/dim)
         cntrd, im, xx, yy, xcen, ycen, delpos, /extendbox, /silent
         delx = abs(xcen - float(xx))
         dely = abs(ycen - float(yy))
         IF xcen lt 0. OR ycen lt 0. or xcen ge sz[0] or ycen ge sz[1] or $
	   delx gt delpos/2. or dely gt delpos/2. THEN BEGIN
            print,"Can't find centroid! - using x,y"
	    xcen=xx
	    ycen=yy
         ENDIF
	 print,"Center X,Y: "+string(xx,yy,"(2f9.2)")
         el[3]=ycen
         el[2]=xcen          
         end
    'h': begin
          print, 'z|x : dec, inc ellipse'
          print, ";|' : dec, inc majordiam"
          print, "[|] : dec, inc minordiam"
          print, "w|s|a|d : move center up,down,left,right"
          print, "<|> or ,|. : rotate ccw/cw"
	  print, "% : multiply size by percentage"
	  print, "/ : enter b/a axial ratio"
	  print, "r : centroid on mouse position"
          print, "p : print ellipse parameters"
	  print, "h : print this help"
          print, "c : cancel all changes"
          print, "q : quit and return"
         end
    else:  ;any other key press does nothing
ENDCASE

if cancel then begin
  p = DIALOG_MESSAGE('Are you sure, cancel changes? (Yes: Cancel / NO: Save Changes)', /QUESTION,  /CENTER, /DEFAULT_NO, TITLE='cancel?')
  if  strcmp(p, 'Yes') then el = elorig
endif

stv_display
tvellipse,el[0],el[1],el[2],el[3],float(90.+el[4]),$
   color=colordex('R'),linestyle=1
plots,el[5],el[6],psym=4,symsi=3.,color=colordex('G'),/device
plots,el[2],el[3],psym=1,symsi=2.,color=colordex('R'),/device

minskywid = 27. / as_pix   ; 5 for ACS ; 27. otherwise
factor=1.5
skyinner=factor*el[0]
skyouter= el[0]*sqrt(1+factor^2) > (skyinner+minskywid) ; ACS (skyinner+minskywid) ; el[0]*sqrt(1+factor^2) > (skyinner+minskywid)

tvellipse,skyinner,skyinner/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
   color=colordex(anulus_col),linestyle=1
tvellipse,skyouter,skyouter/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
   color=colordex(anulus_col),linestyle=1


if exit and (el[0] ne elorig[0] or el[1] ne elorig[1] or $
    el[2] ne elorig[2] or el[3] ne elorig[3] or $
    el[4] ne elorig[4]) then begin

 if have_astr then begin
	 XY2AD,el[2],el[3],astr,el_a,el_d
	 el_ast = 1
 endif else begin
	 as_pix = 1.
	 el_a = el[2]
	 el_d = el[3]
	 el_ast = 0
 endelse
 ;
 ; make sure pa is positive
 while el[4] lt 0. do el[4] = el[4] + 180.

 openw,lun,aux_base+'_stv_ellipse_output.dat',/get_lun
 printf,lun,'# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix ): '+systime(0)
 printf, lun, el[0]*2.*as_pix, el[1]*2.*as_pix, el_a, el_d, el[4], $
	 el[0]*2, el[1]*2, el[2], el[3], el_ast, as_pix, $
   format = '(2f9.1,2f13.8,f7.1,4f9.2,i5,f9.5)'
 free_lun,lun
 
endif

if (xregistered('stv', /noshow)) then $
  widget_control, keyboard_text_id, /clear_events

END

;--------------------------------------------------------------------

PRO stv_ellipse

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

if total(el) eq 0 then begin
 print, 'Ellipse: no initial parameters'
 goto, finish
endif


if total(el) eq 0 then goto,  finish

tvellipse,el[0],el[1],el[2],el[3],float(90.+el[4]),$
   color=colordex('R'),linestyle=1
plots,el[5],el[6],psym=4,symsi=2.,color=colordex('G'),/device
plots,el[2],el[3],psym=1,symsi=1.,color=colordex('R'),/device

minskywid = 27.0 / as_pix ; 5. for ACS ; 27. otherwise
factor=1.5
skyinner=factor*el[0]
skyouter=el[0]*sqrt(1+factor^2) > (skyinner+minskywid)  ; ACS (skyinner+minskywid) ; el[0]*sqrt(1+factor^2) > (skyinner+minskywid)

tvellipse,skyinner,skyinner/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
   color=colordex(anulus_col),linestyle=1
tvellipse,skyouter,skyouter/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
   color=colordex(anulus_col),linestyle=1

print, "Ellipse: z|x|%|/|;|'|[|]|w|s|a|d|<|>|(r)e-center|(p)rint|(h)elp|(c)ancel|(q)uit"
WIDGET_CONTROL, keyboard_text_id, event_pro = 'stv_ellipse_keyboard_event'

finish:

END 

;--------------------------------------------------------------------

PRO stv_delete

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

s=size(im,/dimensions)
print, ''
print, 'Deleting point sources in a defined region.'
print, 'Choosing the region ... '
print, 'Right click when you are done ...'
print, ''
Roi = defroi(s[0],s[1])

q=''
q = DIALOG_MESSAGE('Delete point sources in region?', /QUESTION, /DEFAULT_NO,  /CENTER)

if strcmp(q, 'Yes') then begin

	xy  = array_indices(im,Roi)
	dx  = reform(xy[0,*])
	dy  = reform(xy[1,*])
	;
	; psrc
	if total(psrc) ne 0. then begin
 		t=where(psrc[6,*] gt 0., nt)
 		if nt gt 0 then $
   			for i=0l,nt-1 do begin
				rsrc=sqrt((dx-psrc[2,t[i]])^2 + $
					  (dy-psrc[3,t[i]])^2)
				if min(rsrc) le 1. then $
					psrc[6,t[i]] = 0.
   			endfor
	endif

	;
	; nsrc
	if total(nsrc) ne 0. then begin
 		t=where(nsrc[6,*] gt 0., nt)
 		if nt gt 0 then $
   			for i=0l,nt-1 do begin
				rsrc=sqrt((dx-nsrc[2,t[i]])^2 + $
					  (dy-nsrc[3,t[i]])^2)
				if min(rsrc) le 1. then $
					nsrc[6,t[i]] = 0.
   			endfor
	endif

	wsrc = (1 eq 1)	; write out new point source list

endif 

stv_display

END 
;--------------------------------------------------------------------

PRO stv_undo_roi

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

print, ' '
print, 'The last created region is being removed ....'
print, 'Please wait ...'
print, ' '
print, ' '
	    
sharps = undo_stv_region(aux_base+'_stv_roi_output.dat')
stv_display

END 

;--------------------------------------------------------------------

PRO stv_change_edit_region

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
	     
	     

if strcmp(strtrim(tv_filter,2), 'w_com') or strcmp(strtrim(tv_filter,2), 'n_com') then begin
   
   edit_region = aux_base + '_stv_roi_output.dat'
   print, 'Editing the common mask region: ', edit_region
   
endif else begin
   
   common_reg = aux_base + '_stv_roi_output.dat'
   filter_reg = aux_base + '_roi_' + tv_filter + '.dat'
   
   if strcmp(strtrim(edit_region,2), strtrim(common_reg,2)) then begin
      edit_region = filter_reg
   endif else begin
      edit_region = common_reg
   endelse
   
   
   print, "STV fliter: ", tv_filter
   print, 'Editing the mask region: ', edit_region

endelse


end
;--------------------------------------------------------------------
PRO stv_roi

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


s=size(im)

print, ''
print, 'Choosing a mask region. '
print, 'Right click when you are done ...'
print, ''

Roi = defroi(s[1],s[2], Xverts, Yverts)

q=''
q = DIALOG_MESSAGE('Keep region?', /QUESTION,  /CENTER)

if strcmp(q, 'Yes') then begin

openw, lun, edit_region, /append, /get_lun
printf,lun,'# STV ROI output (image index): '+systime(0)
printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

if Roi[0] ge -1 then $
	for i=0L,n_elements(Roi)-1L do $
		printf,lun,Roi[i]
free_lun,lun





endif else stv_display

END 
;--------------------------------------------------------------------

PRO stv_clear_roi

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
q=''

common_reg = aux_base + '_stv_roi_output.dat'
filter_reg = aux_base + '_roi_' + tv_filter + '.dat'

if strcmp(strtrim(edit_region,2), strtrim(common_reg,2)) then b=1 else b=0
if b eq 1 then foi = 'Common Filter' else foi=tv_filter

print, 'Warning: removing all regions ... '
print, 'I am about to remove all regions (filter): ' + foi
print, "ehsan:", edit_region


q = DIALOG_MESSAGE('Keep regions?', /QUESTION,  /CENTER)

if  strcmp(q, 'No') then begin
    p = ''
    p = DIALOG_MESSAGE('Are you sure, you want to remove all masks?', /QUESTION,  /CENTER, /DEFAULT_NO)
    if  strcmp(p, 'Yes') then begin
         if file_test(edit_region) then spawn, 'rm '+ edit_region
         if file_test(rfile) and b eq 1 then spawn, 'rm '+ rfile
         stv_display
    endif

endif


END
;--------------------------------------------------------------------

Pro stv_save_roi


COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
        
        roifile = aux_base+'_roi.dat'
	if file_exist(aux_base+'_stv_roi_output.dat') then begin

		if file_exist(roifile) then begin
			q = DIALOG_MESSAGE('Append new regions?', /QUESTION,  /CENTER)
			if strcmp(q, 'Yes') then begin
				spawn,'cat '+aux_base+'_stv_roi_output.dat >> ' + roifile
				spawn,'rm '+aux_base+'_stv_roi_output.dat'
			endif 
		endif else begin $
			spawn, 'mv '+aux_base+'_stv_roi_output.dat ' + roifile
			endelse
			
        endif




END
;--------------------------------------------------------------------

PRO stv_remove_background

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
	     
	     
q=''
q = DIALOG_MESSAGE('Remove background region?', /QUESTION,  /CENTER, /DEFAULT_NO)

vertices = aux_base+'_stv_background_vertices.dat'
background = aux_base+'_stv_background_output.dat'


if  strcmp(q, 'Yes') then begin
   
   if file_test(vertices) then begin
       spawn, 'rm '+vertices+' '
       endif
   if file_test(background) then begin
       spawn, 'rm '+background+' '
       endif
   stv_display
   
endif 

END
;--------------------------------------------------------------------


;--------------------------------------------------------------------
 ; default is used for creating background region
PRO stv_reg_creator, use_sex=use_sex, force_sex=force_sex

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

	     
	     
vertices = aux_base+'_stv_background_vertices.dat'
background = aux_base+'_stv_background_output.dat'
region_name = 'background'

if keyword_set(use_sex) then begin
vertices = aux_base+'_stv_usesex_vertices.dat'
background = aux_base+'_stv_usesex_output.dat'
region_name = 'Use_SEx'
endif

if keyword_set(force_sex) then begin
vertices = aux_base+'_stv_forcesex_vertices.dat'
background = aux_base+'_stv_forcesex_output.dat'
region_name = 'Force_SEx'
endif     


print, ''
print, 'Choosing '+region_name+' region. '
print, 'Right click when you are done ...'
print, ''

s=size(im)
; reading the region, by choosing its vertices
Back = defroi(s[1],s[2], Xverts, Yverts)







q=''
q = DIALOG_MESSAGE('Keep '+region_name+' region?', /QUESTION,  /CENTER)

if strcmp(q, 'Yes') then begin

openw,lun,background,/append,/get_lun
printf,lun,'# STV '+region_name+' output (image index): '+systime(0)
printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

if Back[0] ge -1 then $
	for i=0L,n_elements(Back)-1L do $
		printf, lun, Back[i]
free_lun,lun

openw,lun,vertices,/append,/get_lun
printf,lun,'# STV '+region_name+' vertices (image index): '+systime(0)
printf,lun,'# IMAGE NX,NY: ',s[1],s[2]

if Xverts[0] ge -1 then $
	for i=0L,n_elements(Xverts)-1L do $
		printf, lun, Xverts[i],  Yverts[i]
free_lun,lun
stv_display

endif else stv_display

END 
;--------------------------------------------------------------------
Pro run_sextract


COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga



findpro, 'glga_qa', DirList=DirList, ProList=ProList, /Noprint
rute = strmid(imfil,0,strpos(imfil,'.fit'))
rute =  strmid(rute, 0, strlen(rute)-2)
rute = strmid(rute,0,strpos(rute,id))

b = strpos(imfil,'.fit')
a = strpos(imfil,'_', /REVERSE_SEARCH)
imfil_band = strmid(imfil,a+1,b-a-1)
filter = tv_filter



if strcmp(strtrim(tv_filter,2), 'w_com') or strcmp(strtrim(tv_filter,2), 'n_com') then begin 

  filter = imfil_band
  
  q = DIALOG_MESSAGE("Note: Using '"+ imfil_band +"' band ... !",  /INFORMATION,  /CENTER)

endif  




sextract_file =  rute+id+'_'+filter+'_sex.txt'

if file_exist(sextract_file) then begin
   
   p = DIALOG_MESSAGE('Use existing SExtract info?', /QUESTION,  /CENTER, TITLE='Use exisiting file?')
   if  strcmp(p, 'Yes') then begin
   
   spawn, 'awk '+"'" +'(NR>2){print($6" "$8" "$16" "$18" "$22)}'+"' "+ sextract_file, results
   results = STRSPLIT(results[n_elements(results)-1], ' ', /extract)
   print, 'flag: ', results
   
      ; testing the SEx-info
      if n_elements(results) ne 5 then begin
	    p = DIALOG_MESSAGE('Warning, something went wrong, Want to re-run SExtractor?', /QUESTION,  /CENTER)
	    if  strcmp(p, 'Yes') then goto, run_again else goto, finish
      endif
   
   
   sc = 60./as_pix    ; major and minor axis are in minutes (convert to pixel)
   goto, apply_results
   endif
   
endif


run_again:


Sex_Script = GETENV('IDL_GLGAPRO') + '/scripts/ellipse_sextract.csh'
if not file_test(Sex_Script) then begin
        print, ' '
        print, 'Warning ... '
        print, 'Could not find '+ Sex_Script
        print, 'Could Run Sextractor ...'
        print, ' '
        goto, finish
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Creating a temporary filtered images
;;; to be used by SExtractor


badmap = aux_base+'_'+filter+'_badmap.fits.gz'
maskimgfile = rute+'/../../aux/'+id+'_'+strtrim(dtype)+'_mask'+'.fits.gz'
bandmask = rute+'/../../aux/'+id+'_'+strtrim(dtype)+'_mask_'+filter+'.fits.gz'


intfile = rute + '/'+id+'_'+filter+'.fits'
int_fi = file_info(intfile)
if not int_fi.exists then begin 
    print,pre+'Cannot find the input file.'
    goto, finish
endif

print, ''
print, 'Generating the temeporary masked fits file to used by SExtractor.'
print, 'This may take a while based on the image size.'
print, 'Please hold on ... :)'
print, ''


int=mrdfits(int_fi.name,0,inthdr,/fscale,/silent)
sz=size(int,/dim)
; mask -- 1:good pixels, 0:bad pixels
mask=bytarr(sz)+1b


use = (1 eq 0)
use_file = aux_base+'_stv_usesex_output.dat'
if file_test(use_file) then begin
	readcol,use_file,useindx,form='l',/silent    
	use = (1 eq 1)
        endif
; only use this region; start from this region
if use then begin
  mask *= 0.
  mask[useindx] = 1b  
  endif

maskidx=-1
maskidx_band=-1

if file_exist(maskimgfile) then begin
  mskimg = glga_getmask(maskimgfile,sz,astr,as_pix,/update)
  maskidx = where(mskimg ge 1)
  delvarx,mskimg,/free
endif

if keyword_set(bandmask) and file_test(bandmask) then begin 
  mskimg2 = glga_getmask(bandmask,sz,astr,as_pix,/update)
  maskidx_band = where(mskimg2 ge 1)
  delvarx,mskimg2,/free
endif  


force = (1 eq 0)
force_file = aux_base+'_stv_forcesex_output.dat'
if file_test(force_file) then begin
	readcol,force_file,forceindx,form='l',/silent    
	force = (1 eq 1)
        endif


if maskidx[0] ge 0 then $
 	mask[maskidx]=0b
 
 ; taking care of masking of each individual mask
if maskidx_band[0] ge 0 then $   
 	mask[maskidx_band]=0b 	

if force then $ 
  mask[forceindx] = 1b  ; These pixel would be forced to be used by SExtractor


  
in=where(mask eq 0, bad_npix)


if bad_npix gt 0 then begin
   spawn, 'rm '+ rute+'*tmp.fits'
   int[in] = !VALUES.F_NAN
   optional_file = rute + '/'+id+'_'+filter+'_tmp.fits'
   writefits, optional_file, int
endif else optional_file='optional_file'



;;;; END Test ZONE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


scr_command  = 'csh '
scr_command += Sex_Script+ ' '
scr_command += id + ' '
scr_command += filter + ' '
scr_command += rute + ' '
scr_command += GETENV('IDL_GLGAPRO') + '/sextract/ '
scr_command += GETENV('IDL_GLGAPRO') + '/bin/ '


print, ''
print, "Click on the galaxy center, to run SExtractor on '"+ imfil_band +"' band"
print, ''
q = DIALOG_MESSAGE("Click on the galaxy center, to run SExtractor on '"+ imfil_band +"' band",  /INFORMATION,  /CENTER, TITLE='choose the center')
cursor,xx,yy,3,/dev 
XYOUTS, xx,yy, 'X', /DEVICE
scr_command += strtrim(string(xx),2) + ' '
scr_command += strtrim(string(xx),2) + ' '
scr_command += strtrim(string(as_pix),2) + ' '
scr_command += optional_file
spawn, scr_command, results

; if not strcmp(optional_file, 'optional_file') then spawn, 'rm ' + optional_file

results = STRSPLIT(results[n_elements(results)-1], ' ', /extract)
sc = 1   ; major and minor axis are already in pixel

if n_elements(results) ne 5 then begin

  print, ' '
  q = DIALOG_MESSAGE('Error: Something went wrong ...',  /INFORMATION,  /CENTER, TITLE='choose the center')
  print, 'Error: Could not find a good ellipse ...'
  print, ' '
  goto, finish

endif

apply_results:

if n_elements(results) eq 5 then begin
 el_2 = float(results[0]) ; xc
 el_3 = float(results[1]) ; yc
 el_0 = float(results[2])*sc ; major axis
 el_1 = float(results[3])*sc ; minor axis
 el_4 = float(results[4]) ; pa

 if el_2 eq 0 or el_3 eq 0 or el_0 eq 0 or el_1 eq 0 then begin 
  print, ' '
  q = DIALOG_MESSAGE('Error: Could not find an ellipse ...',  /INFORMATION,  /CENTER, TITLE='choose the center')
  print, 'Error: Could not find a good ellipse ...'
  print, ' '
  goto, finish
 endif  
    
    
    
 el[0] = 1.3*el_0 ; major axis  - adding 30% (SExtrctor ellipses look small)
 el[1] = 1.3*el_1 ; minor axis  - adding 30% (which is choses arbitrarily)
 el[2] = el_2 ; xc
 el[3] = el_3 ; yc
 el[4] = el_4 ; pa
 ; make sure pa is positive
 while el[4] lt 0. do el[4] = el[4] + 180. 
 

	if have_astr then begin
		XY2AD,el[2],el[3],astr,el_a,el_d
		el_ast = 1
	endif else begin
		as_pix = 1.
		el_a = el[2]
		el_d = el[3]
		el_ast = 0
	endelse  
	

	  openw,lun,aux_base+'_stv_ellipse_output.dat',/get_lun
	  printf,lun,'# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix ): '+systime(0)
	  printf, lun, el[0]*2.*as_pix, el[1]*2.*as_pix, el_a, el_d, el[4], $
	  el[0]*2, el[1]*2, el[2], el[3], el_ast, as_pix, $
	  format = '(2f9.1,2f13.8,f7.1,4f9.2,i5,f9.5)'
	  free_lun,lun
	  
	  ds9_reg_maker, 'idl_foo.reg', el, as_pix
	  spawn, "xpaset -p ds9 regions delete all"
	  spawn, "xpaset -p ds9 regions load idl_foo.reg"
	  stv_display
	  
	  
	  print, ' '
	  print, 'STV updated with new SExtractor generated ellipse ...'
	  print, ''
	  
	  
endif


finish:

END
;--------------------------------------------------------------------

PRO rmvrgn


COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

print,"Click inside the region you'd like to remove ..."

region_file = edit_region
common_reg = aux_base + '_stv_roi_output.dat'


findpro, 'glga_qa', DirList=DirList, ProList=ProList, /Noprint

Remove_Script = GETENV('IDL_GLGAPRO') + '/scripts/remove_region.bash'
if not file_test(Remove_Script) then begin
        print, ' '
        print, 'Warning ... '
        print, 'Could not find '+ Remove_Script
        print, 'Could not remove the region ...'
        print, ' '
        goto, fini
endif

if not file_test(region_file) then begin
        print, ' '
        print, 'Warning ... '
        print, 'Could not find '+ region_file
        print, 'Could not remove the region ...'
        print, ' '
        
        if file_test(rfile) and (strcmp(strtrim(tv_filter,2), 'w_com') or strcmp(strtrim(tv_filter,2), 'n_com')) then begin
        
           region_file = rfile
        
        endif else begin
           if strcmp(strtrim(tv_filter,2), 'w_com') or strcmp(strtrim(tv_filter,2), 'n_com') then begin
 	     print, ' '
	     print, 'Warning ... '
	     print, 'Could not find old region file:'+ rfile
	     print, 'Could not remove the region ...'
	     print, ' '       
             goto, fini
           endif
        endelse
        
endif


cursor,xx,yy,3,/dev 


XYOUTS, xx,yy, 'X', /DEVICE
sz = size(im,/dim)
sx = sz[0]
sy = sz[1]
index = yy*sx + xx

;  print, sz
;  print, xx,yy,index

spawn, 'bash '+ Remove_Script + ' ' + region_file + ' ' + strtrim(string(index),2)

if not strcmp(region_file , rfile) and strcmp(strtrim(edit_region,2), strtrim(common_reg,2))   then begin
    if file_test(rfile) then begin
    spawn, 'bash '+ Remove_Script + ' ' + rfile + ' ' + strtrim(string(index),2)
    endif
    
endif


stv_display

fini:

END 

;--------------------------------------------------------------------

PRO stvcntrd 

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


print,'Mark the center of star (right button to exit).'
print,'Mark an old star to modify.'
expmd = (1 eq 0)

xloop:

cursor,xx,yy,3,/dev 

IF !mouse.button EQ 4 THEN GOTO, finish

wsrc = (1 eq 1)	; write out point source list

delpos = 7.0	; pixels
if strpos(dtype,'sdss') ge 0 then delpos = 6.0
if strpos(dtype,'panstarrs') ge 0 then delpos = 12.0
if strpos(dtype,'acs') ge 0 then delpos = 60.0

sz = size(im,/dim)
cntrd, im, xx, yy, xcen, ycen, delpos, /extendbox, /silent
delx = abs(xcen - float(xx))
dely = abs(ycen - float(yy))
IF xcen lt 0. OR ycen lt 0. or xcen ge sz[0] or ycen ge sz[1] or $
	delx gt delpos/2. or dely gt delpos/2. THEN BEGIN
 print,"Can't find centroid! - using x,y"
	 xcen=xx
	 ycen=yy
ENDIF
;
; are we old?
;
; check read in point sources
roff = sqrt( (psrc[2,*] - xcen)^2 + (psrc[3,*] - ycen)^2 )
w = where(roff le 2., nw)	; within two pixels
; if we are old, modify
if nw ge 1 then begin
	w=w[0]
	q = ' '
	tvcirc,psrc[2,w],psrc[3,w],psrc[4,w],color=colordex('Q')
	while strupcase(strmid(q,0,1)) ne 'Q' do begin
		read,'(x)-expand, (z)-shrink, (o)-on/off, (q)-quit: ',q
		case strupcase(strmid(q,0,1)) of
		'X': psrc[4,w] = psrc[4,w]+1.
		'Z': psrc[4,w] = psrc[4,w]-1.
		'O': if psrc[6,w] eq 0. then $
			  psrc[6,w] = 1.  $ ; turn it on
		     else psrc[6,w] = 0.    ; turn it off
		else:
		endcase
		stv_display
		if psrc[6,w] gt 0. and strupcase(strmid(q,0,1)) ne 'Q' then $
			tvcirc,psrc[2,w],psrc[3,w],psrc[4,w],color=colordex('Q')
	endwhile
	print,'Mark the center of star (right button to exit).'
;
; check new point sources
endif else begin
	roff = sqrt( (nsrc[2,*] - xcen)^2 + (nsrc[3,*] - ycen)^2 )
	w = where(roff le 2., nw)	; within two pixels
	if nw ge 1 then begin
		w=w[0]
		if nsrc[6,w] eq 0. then begin
			nsrc[6,w] = 1.	; turn it on
			tvcirc,nsrc[2,w],nsrc[3,w],nsrc[4,w],color=colordex('G')
		endif else begin
			nsrc[6,w] = 0.	; turn it off
			tvcirc,nsrc[2,w],nsrc[3,w],nsrc[4,w],color=colordex('R')
		endelse
	endif else begin
;
; must be new
		if have_astr then begin
			XY2AD,xcen,ycen,astr,a,d
			psrc_ast = 1
		endif else begin
			a = xcen
			d = ycen
			psrc_ast = 0
		endelse
;
; use subim to speed things up
		x0=fix(xcen-63.)>0
		x1=x0+127<(sz[0]-1)
		y0=fix(ycen-63.)>0
		y1=y0+127<(sz[1]-1)
		radius = get_mask_radius(im[x0:x1,y0:y1],xcen-float(x0),ycen-float(y0),dtype)
;
; install in new array
		w=where(nsrc[0,*] eq 0.)
		w=w[0]
		nsrc[0,w] = a
		nsrc[1,w] = d
		nsrc[2,w] = xcen
		nsrc[3,w] = ycen
		nsrc[4,w] = radius
		nsrc[5,w] = psrc_ast
		nsrc[6,w] = 1
;
; display new one
		tvcirc,xcen,ycen,radius, color=colordex('G')

		radius = radius * as_pix	; convert to arcsec

		print, 'Centroid is (a,d,x,y,r): ',a,d,xcen,ycen,radius, $
			format='(a,2f13.8,2f9.2,f9.3)'
	endelse	; are we from new point sources?
endelse	; are we from read in point sources?

GOTO, xloop

finish:

END 

;--------------------------------------------------------------------

function get_roi

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

; is roi file name set?
if strlen(rfile) gt 1 then begin
	; does roi file exist?
	if file_test(rfile) then begin
		; get image dimensions
		spawn,'grep IMAGE '+rfile,sout
		; do we have the requisite record?
		if strlen(sout[0]) gt 0 then begin
			rx = 0L & ry = 0L
			; use last (most recent) image size
			for i=0l,n_elements(sout)-1 do begin
				cmnt = sout[i]
				for j=0,2 do jnk = gettok(cmnt,' ')
				rx = long(gettok(cmnt,' '))
				ry = long(gettok(cmnt,' '))
			endfor
			; read in roi data
			readcol,rfile,indx,form='l',/silent
			; get current image size
			s = size(im)
			; do we need to adjust roi?
			; we expanded
			if s[1] gt rx then begin
				x0 = (s[1] - rx) / 2
				y0 = (s[2] - ry) / 2
				x1 = x0 + rx - 1
				y1 = y0 + ry - 1
				large = bytarr(s[1],s[2])
				small = bytarr(rx,ry)
				small[indx] = 1b
				large[x0:x1,y0:y1] = small[*,*]
				indx = where(large eq 1)
				; re-write file with new img dims
				print,'Re-writing ROI file: '+rfile
				openw,lun,rfile,/get_lun
				printf,lun,'# STV ROI output (image index): '+systime(0)
				printf,lun,'# IMAGE NX,NY: ',s[1],s[2]
				for i=0L,n_elements(indx)-1L do $
					printf,lun,indx[i]
				free_lun,lun
			endif
			; we shrank
			if s[1] lt rx then begin
				x0 = (rx - s[1]) / 2
				y0 = (ry - s[2]) / 2
				x1 = x0 + s[1] - 1
				y1 = y0 + s[2] - 1
				large = bytarr(rx,ry)
				large[indx] = 1b
				small = large[x0:x1,y0:y1]
				indx = where(small eq 1)
				; re-write file with new img dims
				print,'Re-writing ROI file: '+rfile
				openw,lun,rfile,/get_lun
				printf,lun,'# STV ROI output (image index): '+systime(0)
				printf,lun,'# IMAGE NX,NY: ',s[1],s[2]
				for i=0L,n_elements(indx)-1L do $
					printf,lun,indx[i]
				free_lun,lun
			endif
		; else blindly read it in
		endif else $
			readcol,rfile,indx,form='l',/silent
	endif else	indx = -1
endif else indx = -1

return,indx

END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function read_back_roi, filename

array = ''
line = ''
; file = 'stv_background_vertices.dat'

OPENR, lun, filename, /GET_LUN

WHILE NOT EOF(lun) DO BEGIN & $
  READF, lun, line & $
  array = [array, line] & $
ENDWHILE


xvertices = [-1]
yvertices = [-1]
i = 0 
while i lt n_elements(array) do begin
  

  first = strmid(array[i], 0, 1)
  while strcmp(first , '#') and i lt n_elements(array) do begin
      i++
      if i eq n_elements(array) then break
      first = strmid(array[i], 0, 1)
      
  endwhile
  
  S = STRSPLIT(array[i],' ', /extract)
  x0 = long(S[0])
  y0 = long(S[n_elements(S)-1])

  while first ne '#' and i lt n_elements(array) do begin
       S = STRSPLIT(array[i],' ', /extract)
       xvertices = [xvertices, long(S[0])]
       yvertices = [yvertices, long(S[n_elements(S)-1])]
       i++
       if i eq n_elements(array) then break
       first = strmid(array[i], 0, 1)
  endwhile
  
  ; to have a closed polygon
  xvertices = [xvertices, x0, -1]
  yvertices = [yvertices, y0, -1]
  

endwhile 

i = 0
while xvertices[i] le 0 do i++

xvertices =  xvertices[i-1:n_elements(xvertices)-1]
yvertices =  yvertices[i-1:n_elements(yvertices)-1]

free_lun,lun


return, [[xvertices], [yvertices]]

end ; function ends here


;;;;;;;;;;;;;;;;;;;;;;;;;;;;



PRO stv_display

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga

	     

vertices = aux_base+'_stv_background_vertices.dat'
background = aux_base+'_stv_background_output.dat'	     

; Background
backindx = -1
if file_test(background) then begin
	readcol,background,rndx,form='l',/silent    ; additional mask
; 	if backindx[0] ge 0 then $
; 		indx = [indx,rndx] $
; 	else	indx = rndx
	backindx = rndx
endif



; taking SDSS masks 
if strpos(dtype,'panstarrs') ge 0 and not file_test(aux_base+'_stv_roi_output.dat') and not file_test(aux_base+'_roi.dat') then begin
   
   id_rute = strmid(id,3,strlen(id)-1)
   rute = strmid(aux_base,0,strpos(aux_base,'aux/'))
   
   sdss_rfile   = rute+'aux/'+id_rute+'_sdss_roi.dat'
   ps_rfile = aux_base+'_stv_roi_output.dat'
   
   im_ps = im
   hdr_ps = hdr
   extast,hdr_ps,astr_ps
   sdss_file = rute+'sdss/fits/'+id_rute+'_r.fits'
   

   
   if file_test(sdss_file) and file_test(sdss_rfile) then begin
   g = DIALOG_MESSAGE('Do you want to use SDSS masks (speed-up)?', /QUESTION,  /CENTER)
   if strcmp(g, 'Yes') then begin
      
      print, ''
      print, ' OK ... :)'
      print, 'I am now translating SDSS masks for this current image ...'
      print, 'You can remove each region seprately.'
      print, '' 
       
      im_sdss  = mrdfits(sdss_file, 0, hdr_sdss, /fscale, /silent)
      extast,hdr_sdss,astr_sdss
   
      indx_lst = get_roi_sep(im_sdss, sdss_rfile)
      N = indx_lst.count()
      for j=0, N-1 do begin

            indx = indx_lst[j]


            s = size(im_sdss,/dim)
            sx = s[0]
            sy = s[1]


            y_sdss = indx/sx
            x_sdss = indx mod sx



            XY2AD, x_sdss, y_sdss, astr_sdss, ra, dec 
            AD2XY, ra, dec, astr_ps, x_ps, y_ps

            s = size(im_ps,/dim)
            sx = s[0]
            sy = s[1]

            y_ps = round(y_ps)
            x_ps = round(x_ps)
            indx_ps = y_ps*sx+x_ps

            M = n_elements(indx_ps) 
            indx_ps_new = lonarr(9L*M)


            p = 0L
            for i=0L,M-1L do begin

                  QX = x_ps[i]
                  QY = y_ps[i]

                  indx_ps_new[p] = (QY)*sx+(QX)
                  p++
                  indx_ps_new[p] = (QY+1)*sx+(QX+1)
                  p++
                  indx_ps_new[p] = (QY)*sx+(QX+1)
                  p++
                  indx_ps_new[p] = (QY+1)*sx+(QX)
                  p++
                  indx_ps_new[p] = (QY)*sx+(QX-1)
                  p++
                  indx_ps_new[p] = (QY-1)*sx+(QX-1)
                  p++
                  indx_ps_new[p] = (QY+1)*sx+(QX-1)
                  p++
                  indx_ps_new[p] = (QY-1)*sx+(QX)
                  p++
                  indx_ps_new[p] = (QY-1)*sx+(QX+1)
                  p++

            endfor


            ; ; Getting rid of redundancy
            ind = sort(indx_ps_new)
            indx_ps_new =  indx_ps_new[ind]

            tmp = -10000
            p = 0L
            for i=0L,n_elements(indx_ps_new)-1L do begin
              if indx_ps_new[i] ne tmp then begin 
                 indx_ps_new[p] = indx_ps_new[i]
                 p++
                 tmp = indx_ps_new[i]
              endif 
            endfor



            s = size(im_ps)
            if j eq 0 then openw,lun,ps_rfile,/get_lun else openw,lun,ps_rfile,/append,/get_lun
            printf,lun,'# STV ROI output (image index): '+systime(0)
            printf,lun,'# IMAGE NX,NY: ',s[1],s[2]
            for i=0L,p-1L do $
            printf,lun,indx_ps_new[i]
            free_lun,lun


      endfor   ; for all regions
   endif ; question box
   endif ; if sdss info exists

endif   ; end taking SDSS masks


; Masking
indx = get_roi()   ; Initial mask

if file_test(aux_base+'_stv_roi_output.dat') then begin
	readcol,aux_base+'_stv_roi_output.dat',rndx,form='l',/silent    ; additional mask
	if indx[0] ge 0 then $
		indx = [indx,rndx] $
	else	indx = rndx
endif



if not strcmp(strtrim(tv_filter,2), 'w_com') and not strcmp(strtrim(tv_filter,2), 'n_com') then begin

   filter_reg = aux_base + '_roi_' + tv_filter + '.dat'

   if file_test(filter_reg) then begin
        readcol, filter_reg,f_indx,form='l',/silent  
   endif
   
   
        

endif



if file_test(vertices) then begin
	readcol,vertices,xvertices, yvertices,form='l, l',/silent    ; additional mask
endif


if strue then begin 
 dim = tim
 if indx[0] ge 0 then begin
	 dimr = dim[0, *, *]
	 dimg = dim[1, *, *]
	 dimb = dim[2, *, *]
	 
	 dimr[indx] = 100
	 dimg[indx] = 160
	 dimb[indx] = 100
	 
         dim[0, *, *] = dimr
	 dim[1, *, *] = dimg
	 dim[2, *, *] = dimb
 endif
 if n_elements(f_indx) gt 0 then begin

 	 dimr = dim[0, *, *]
	 dimg = dim[1, *, *]
	 dimb = dim[2, *, *]
	 
	 dimr[f_indx] = 255
	 dimg[f_indx] = 0
	 dimb[f_indx] = 127
	 
         dim[0, *, *] = dimr
	 dim[1, *, *] = dimg
	 dim[2, *, *] = dimb 
 
 endif
  
 if n_elements(xvertices) gt 0 then begin
	 dimr = dim[0, *, *]
	 dimg = dim[1, *, *]
	 dimb = dim[2, *, *]

	 timr = tim[0, *, *]
	 timg = tim[1, *, *]
	 timb = tim[2, *, *]
	 
	 dimr[backindx] = timr[backindx]
	 dimg[backindx] = timg[backindx]
	 dimb[backindx] = timb[backindx]
	 
	 dim[0, *, *] = dimr
	 dim[1, *, *] = dimg
	 dim[2, *, *] = dimb
	 
 endif 

 
 
tv, dim,true=1 
endif else begin
 dim = im
 if indx[0] ge 0 then $
	 dim[indx] = smin
 if sasinh then tvscl,asinh(smooth(dim > smin < smax, ssmooth,/nan))
 if slog then tvscl,alog10(smooth(dim > smin < smax, ssmooth,/nan) +abs(smin)+1)
 if not slog and not sasinh then tvscl,smooth(dim > smin < smax, ssmooth,/nan)
endelse

; Background regions 
if file_test(vertices) then begin
  polygons = read_back_roi(vertices)
  xvertices = polygons[*,0]
  yvertices = polygons[*,1]

  split = where(xvertices eq -1)

    for i = 0, n_elements(split)-2 do begin
    plots,xvertices[split[i]+1: split[i+1]-1], yvertices[split[i]+1: split[i+1]-1],color=colordex('S'),linestyle=2, thick=2, /device
  endfor
  
  anulus_col = 'Y'
endif else anulus_col = 'S'


; usesex regions 
vertices = aux_base+'_stv_usesex_vertices.dat'
if file_test(vertices) then begin
  polygons = read_back_roi(vertices)
  xvertices = polygons[*,0]
  yvertices = polygons[*,1]

  split = where(xvertices eq -1)

    for i = 0, n_elements(split)-2 do begin
    plots,xvertices[split[i]+1: split[i+1]-1], yvertices[split[i]+1: split[i+1]-1],color=colordex('H'),linestyle=3, thick=1, /device
  endfor
  
endif 


; forcesex regions 
vertices = aux_base+'_stv_forcesex_vertices.dat'
if file_test(vertices) then begin
  polygons = read_back_roi(vertices)
  xvertices = polygons[*,0]
  yvertices = polygons[*,1]

  split = where(xvertices eq -1)

    for i = 0, n_elements(split)-2 do begin
    plots,xvertices[split[i]+1: split[i+1]-1], yvertices[split[i]+1: split[i+1]-1],color=colordex('F'),linestyle=3, thick=1, /device
  endfor
  
endif 




; Ellipse
if total(el) ne 0 then begin
 tvellipse,el[0],el[1],el[2],el[3],float(90.+el[4]),$
   color=colordex('R'),linestyle=1
 plots,el[5],el[6],psym=4,symsi=3.,color=colordex('G'),/device
 plots,el[2],el[3],psym=1,symsi=1.,color=colordex('R'),/device

 minskywid = 27.0 / as_pix ; 5.0 for ACS ; 27. otherwise
 factor=1.5
 skyinner=factor*el[0]
 skyouter=el[0]*sqrt(1+factor^2) > (skyinner+minskywid) ; ACS (skyinner+minskywid) ;   el[0]*sqrt(1+factor^2) > (skyinner+minskywid)

 tvellipse,skyinner,skyinner/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
    color=colordex(anulus_col),linestyle=1
 tvellipse,skyouter,skyouter/(el[0]/el[1]),el[2],el[3],float(90.+el[4]),$
    color=colordex(anulus_col),linestyle=1
endif

; point Sources
if total(psrc) ne 0. then begin
 t=where(psrc[6,*] gt 0., nt)
 if nt gt 0 then $
   for i=0l,nt-1 do begin
     if psrc[6,t[i]] gt 1. then $
	     col = colordex('B') $
     else    col = colordex('O')
     tvcirc,psrc[2,t[i]], psrc[3,t[i]], psrc[4,t[i]], color=col
   endfor
endif

if total(nsrc) ne 0. then begin
 t=where(nsrc[6,*] gt 0., nt)
 if nt gt 0 then $
   for i=0l,nt-1 do $
     tvcirc,nsrc[2,t[i]], nsrc[3,t[i]], nsrc[4,t[i]], color=colordex('G')
endif

END

;--------------------------------------------------------------------
pro flag_event, event
   
 COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth, $
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
 common stvc, flag_w, pnale_w
  
WIDGET_CONTROL, event.id, get_UVALUE = UVALUE


if strcmp(UVALUE, 'flag_gr') then begin 
    if Event.value eq 0 then if Event.select then qa.complete = 1 else qa.complete = 0             ; Finished
    if Event.value eq 1 then if Event.select then qa.uncertain_mask = 1 else qa.uncertain_mask = 0 ; Uncertain
    if Event.value eq 2 then if Event.select then qa.fov_expand = 1 else qa.fov_expand = 0         ; FOV
    if Event.value eq 3 then if Event.select then qa.multiple = 1 else qa.multiple = 0             ; Multiple
    if Event.value eq 4 then if Event.select then qa.bright_star = 1 else qa.bright_star = 0       ; Bright-star
    stv_print_qa_flags,/values
endif  

if strcmp(UVALUE, 'band1_gr') then begin 
    if Event.value eq 0 then if Event.select then qa.band1_edge = 1 else qa.band1_edge = 0           ; edge
    if Event.value eq 1 then if Event.select then qa.band1_sn_grad = 1 else qa.band1_sn_grad = 0     ; sn grad
    if Event.value eq 2 then if Event.select then qa.band1_artifact = 1 else qa.band1_artifact = 0   ; artifact
    if Event.value eq 3 then if Event.select then qa.band1_missing = 1 else qa.band1_missing = 0     ; missing
    if Event.value eq 4 then if Event.select then qa.band1_other = 1 else qa.band1_other = 0         ; other
    stv_print_qa_flags,/values
endif  

if strcmp(UVALUE, 'band2_gr') then begin 
    if Event.value eq 0 then if Event.select then qa.band2_edge = 1 else qa.band2_edge = 0           ; edge
    if Event.value eq 1 then if Event.select then qa.band2_sn_grad = 1 else qa.band2_sn_grad = 0     ; sn grad
    if Event.value eq 2 then if Event.select then qa.band2_artifact = 1 else qa.band2_artifact = 0   ; artifact
    if Event.value eq 3 then if Event.select then qa.band2_missing = 1 else qa.band2_missing = 0     ; missing
    if Event.value eq 4 then if Event.select then qa.band2_other = 1 else qa.band2_other = 0         ; other
    stv_print_qa_flags,/values
endif  

if strcmp(UVALUE, 'band3_gr') then begin 
    if Event.value eq 0 then if Event.select then qa.band3_edge = 1 else qa.band3_edge = 0           ; edge
    if Event.value eq 1 then if Event.select then qa.band3_sn_grad = 1 else qa.band3_sn_grad = 0     ; sn grad
    if Event.value eq 2 then if Event.select then qa.band3_artifact = 1 else qa.band3_artifact = 0   ; artifact
    if Event.value eq 3 then if Event.select then qa.band3_missing = 1 else qa.band3_missing = 0     ; missing
    if Event.value eq 4 then if Event.select then qa.band3_other = 1 else qa.band3_other = 0         ; other
    stv_print_qa_flags,/values
endif  

;;; If enter-key is pressed when tryping the note
; if strcmp(UVALUE, 'flag_note') then begin 
;     print, event.ch
; endif

if strcmp(UVALUE, 'save') then begin 
   if xregistered ('flag_base') then begin
        WIDGET_CONTROL, WIDGET_INFO(flag_w, FIND_BY_UNAME='flag_qa_note'), get_value=text
        qa.note = text
        if (xregistered ('stv')) then widget_control, flag_w, /destroy
        if (xregistered ('stv')) then widget_control, sbase, /destroy
   endif
endif


if strcmp(UVALUE, 'cancel') then begin 
             WIDGET_CONTROL, keyboard_text_id, event_pro = $
		'stv_keyboard_event'
             WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col_widg'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col2_widg'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='import'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='export'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='sextract'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col4_widg'), sensitive=1
	     WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='wLRow1_widg'), sensitive=1
	     if (xregistered ('flag_base')) then  widget_control, flag_w, /destroy
	     stv_printhelpline
endif



end
;--------------------------------------------------------------------


PRO flag_widget

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


common stvc, flag_w, pnale_w
   
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col_widg'), sensitive=0
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col2_widg'), sensitive=0
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='col4_widg'), sensitive=0
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='wLRow1_widg'), sensitive=0
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='import'), sensitive=0
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='export'), sensitive=0
   WIDGET_CONTROL, WIDGET_INFO(pnale_w, FIND_BY_UNAME='sextract'), sensitive=0
   
   device,get_screen_size=screen_size
   xoff = screen_size[0]/2
   yoff = screen_size[1]/3
   
   flag_base =  WIDGET_BASE(TITLE='Flags & Note', uname='flag_base_uname', Row=6,  space=10, /BASE_ALIGN_CENTER, xoffset = xoff, yoffset = yoff)
      
      rowBase0 = widget_base(flag_base, COLUMN=2, /frame, space=20) 
        rowBase01 = widget_base(rowBase0, /col, /frame) 
         flag_label1  = WIDGET_Label(rowBase01, value = 'User       : ' + qa.user_name, /ALIGN_LEFT) 
         flag_label2  = WIDGET_Label(rowBase01, value = 'Machine    : ' + qa.machine_name, /ALIGN_LEFT) 
         flag_label3  = WIDGET_Label(rowBase01, value = 'No of runs : ' + strtrim(string(qa.nqa),2), /ALIGN_LEFT)
         flag_label4  = WIDGET_Label(rowBase01, value = '', /ALIGN_LEFT) 
        rowBase02 = widget_base(rowBase0, /col, /frame) 
         flag_label5  = WIDGET_Label(rowBase02, value = 'Obj:  ' + id, /ALIGN_LEFT)
         flag_label6  = WIDGET_Label(rowBase02, value = 'Survey: ' + dtype, /ALIGN_LEFT)
         flag_label7  = WIDGET_Label(rowBase02, value = '', /ALIGN_LEFT) 
         
      
      rowBase1 = widget_base(flag_base, /row, /frame) 
         flag_row1  = WIDGET_Label(rowBase1, value = 'FLAGS:') 
         values = ['Finished', 'Uncertain', 'FOV', 'Multiple', 'Bright-star'] 
         flag_gr = CW_BGROUP(rowBase1, values, UVALUE='flag_gr', /Row, /NONEXCLUSIVE, $ 
               /FRAME, SET_VALUE=set_val, uname='flag_gr_uname') 
      rowBase2 = widget_base(flag_base, /row, /frame) 
         flag_row2  = WIDGET_Label(rowBase2, value = 'BAND1 [FUV,g,j]:')  
         values = ['edge', 'sn grad', 'artifact', 'missing', 'other'] 
         band1_gr = CW_BGROUP(rowBase2, values, UVALUE='band1_gr', /Row, /NONEXCLUSIVE, $ 
               /FRAME, SET_VALUE=set_val, uname='band1_gr_uname') 
      rowBase3 = widget_base(flag_base, /row, /frame) 
         flag_row3  = WIDGET_Label(rowBase3, value = 'BAND2 [NUV,r h]:') 
         band2_gr = CW_BGROUP(rowBase3, values, UVALUE='band2_gr', /Row, /NONEXCLUSIVE, $ 
               /FRAME, SET_VALUE=set_val, uname='band2_gr_uname') 
      rowBase4 = widget_base(flag_base, /row, /frame) 
         flag_row4  = WIDGET_Label(rowBase4, value = 'BAND3 [i,k]    :') 
         band3_gr = CW_BGROUP(rowBase4, values, UVALUE='band3_gr', /Row, /NONEXCLUSIVE, $ 
               /FRAME, SET_VALUE=set_val, uname='band3_gr_uname')
      rowBase5 = widget_base(flag_base, /row, /frame, UVALUE='flag_base5')    
         flag_qa_lab  = WIDGET_Label(rowBase5, value = 'QA Note: ') 
         flag_note    = WIDGET_TEXT( rowBase5, value = qa.note, /editable, uname='flag_qa_note', UVALUE='flag_note')
         flag_row5  = WIDGET_BUTTON(rowBase5, value = 'Close and Continue ...', UVALUE='save') 
         flag_row6  = WIDGET_BUTTON(rowBase5, value = 'Cancel', UVALUE='cancel') 
  
  WIDGET_CONTROL, flag_base, /realize, TAB_MODE=0
  XMANAGER, 'flag_base', flag_base, event_handler='flag_event'
  
  flag_w = flag_base
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='flag_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.complete,qa.uncertain_mask,qa.fov_expand,qa.multiple,qa.bright_star]
  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='band1_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.band1_edge, qa.band1_sn_grad, qa.band1_artifact, qa.band1_missing, qa.band1_other]
  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='band2_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.band2_edge, qa.band2_sn_grad, qa.band2_artifact, qa.band2_missing, qa.band2_other]  
  
  group = WIDGET_INFO(flag_w, FIND_BY_UNAME='band3_gr_uname')
  WIDGET_CONTROL, group, SET_VALUE=[qa.band3_edge, qa.band3_sn_grad, qa.band3_artifact, qa.band3_missing, qa.band3_other]  
  

END 





;--------------------------------------------------------------------

PRO stv_shutdown

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


if strcmp(misc, 'STOP') then begin
 
 p = DIALOG_MESSAGE('Are you sure, you want to QUIT?', /QUESTION,  /CENTER, /DEFAULT_NO, TITLE='Quite?')
 if strcmp(p, 'No') then begin
    misc=''
    goto, finish
 endif

endif
	     
stv_save_roi


 spawn, 'xpaaccess ds9', results ; check if ds9 is open
 if strcmp(results[0], 'yes') then begin
   spawn, 'xpaset -p ds9 exit'
 endif

	     
if strpos(misc,'STOP') lt 0  and $
	browse eq 0 then begin
   if (total(psrc[6,*]) gt 0. or total(nsrc[6,*]) gt 0.) and wsrc then begin
	openw,lun,aux_base+'_stv_cntrd_output.dat',/get_lun
	printf,lun,'# STV centroid output (a,d,x,y,r_asec,astr?): '+systime(0)
	for i=0L,maxp-1 do $
		if psrc[6,i] gt 0. then $
			printf,lun,psrc[0,i],psrc[1,i],psrc[2,i],psrc[3,i], $
				psrc[4,i]*as_pix,fix(psrc[5,i]), $
					format='(2f13.8,3f9.3,i5)'
	for i=0L,maxp-1 do $
		if nsrc[6,i] gt 0. then $
			printf,lun,nsrc[0,i],nsrc[1,i],nsrc[2,i],nsrc[3,i], $
				nsrc[4,i]*as_pix,fix(nsrc[5,i]), $
					format='(2f13.8,3f9.3,i5)'
	free_lun,lun
   endif
   if strpos(misc,'RUN') lt 0 and strpos(misc,'PREVIOUS') and strpos(misc,'SAVE') lt 0 and strpos(misc,'SKIP') lt 0 then begin
	stv_print_qa_flags,/help,/values,/note,/details
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	flag_widget   ; This opens the pop-up closing window
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	
	WIDGET_CONTROL, keyboard_text_id, event_pro = $
		'stv_setflags_keyboard_event'

   endif else begin 
       if (xregistered ('stv')) then begin 
           widget_control, sbase, /destroy
       endif
   endelse
endif else begin
	if (xregistered ('stv')) then begin 
           widget_control, sbase, /destroy
       endif
endelse


finish:

END

;-------------------------------------------------------------------

PRO stv_cleanup,windowid

COMMON stvb

IF !d.name EQ 'X' AND !d.window NE -1 THEN curs,34
 delvarx, stvb

END
; 


;--------------------------------------------------------------------
pro stv_grab_region, action=action

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga
	     

   
   
   if strcmp(action, 'to_ds9') then begin 
      
        ds9_reg_maker, 'idl_foo.reg', el, as_pix
        
        spawn, "xpaset -p ds9 regions delete all"
        spawn, "xpaset -p ds9 regions load idl_foo.reg"
   
   endif
   
   
   if strcmp(action, 'update') then begin 
       
       spawn, "xpaset -p ds9 regions system  image"
       spawn, "xpaset -p ds9 regions format ds9"
       spawn, "xpaset -p ds9 regions save idl_foo.reg"
   
   
   endif 

   

   
   
   
   if strcmp(action, 'update') and file_exist('idl_foo.reg') then begin
      
	el[0:4] = read_ds9_ellipse('idl_foo.reg')

	if have_astr then begin
		XY2AD,el[2],el[3],astr,el_a,el_d
		el_ast = 1
	endif else begin
		as_pix = 1.
		el_a = el[2]
		el_d = el[3]
		el_ast = 0
	endelse      
	  
	  openw,lun,aux_base+'_stv_ellipse_output.dat',/get_lun
	  printf,lun,'# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix ): '+systime(0)
	  printf, lun, el[0]*2.*as_pix, el[1]*2.*as_pix, el_a, el_d, el[4], $
	  el[0]*2, el[1]*2, el[2], el[3], el_ast, as_pix, $
	  format = '(2f9.1,2f13.8,f7.1,4f9.2,i5,f9.5)'
	  free_lun,lun
	  
	  ds9_reg_maker, 'idl_foo.reg', el, as_pix
	  spawn, "xpaset -p ds9 regions delete all"
	  spawn, "xpaset -p ds9 regions load idl_foo.reg"
	  stv_display
      
   endif
   
   
   if strcmp(action, 'undo') then begin
      
      q = ''
      q = DIALOG_MESSAGE('Use original ellipse info?', /QUESTION,  /CENTER, /DEFAULT_NO, TITLE='recover the original ellipse?')
      
      if  strcmp(q, 'Yes') then begin
          el = elorig
          if file_exist(aux_base+'_stv_ellipse_output.dat') then spawn, 'rm  '+aux_base+'_stv_ellipse_output.dat'
          
          ds9_reg_maker, 'idl_foo.reg', el, as_pix
          spawn, "xpaset -p ds9 regions delete all"
          spawn, "xpaset -p ds9 regions load idl_foo.reg"
          stv_display
      endif
   
   endif
   
   
   if strcmp(action, 'sdss') then begin   
   
        id_rute = strmid(id,3,strlen(id)-1)
        rute = strmid(aux_base,0,strpos(aux_base,'aux/'))
        sdss_ellipse   = rute+'aux/'+id_rute+'_ellipse.dat'
        
	if file_exist(sdss_ellipse) then begin 
	    
	    readcol,sdss_ellipse,majdiam_as,mindiam_as,el_ra,el_dec, $
			pa_, majdiam_px, mindiam_px, x0_,y0_, astrom_bool, $
			el_as_pix, format='f,f,d,d,f,f,f,f,i,f', /silent
			
	    ad2xy,el_ra,el_dec,astr,x0_,y0_
            x=x0_[0]
            y=y0_[0]
            pa=pa_[0]
            d=majdiam_as[0]
            rat=majdiam_as[0]/mindiam_as[0]
	    el = [0.5*d/as_pix, 0.5*(d/rat/as_pix), x, y, $
			pa,el[5],el[6]]
	    stv_display
	    
	endif
   endif
   
   if strcmp(action, 'ell_catal') then begin   
       
      q = DIALOG_MESSAGE('Use the input catalog ellipse info?', /QUESTION,  /CENTER, /DEFAULT_NO, TITLE='recover the catalog ellipse?')
      
      if  strcmp(q, 'Yes') then begin
        el = ellipse_glga
        stv_display
      endif
   endif ; end ell_catal
   
   if strcmp(action, 'save_ellipse') then begin   
   
	; make sure pa is positive
	while el[4] lt 0. do el[4] = el[4] + 180.

		if have_astr then begin
			XY2AD,el[2],el[3],astr,el_a,el_d
			el_ast = 1
		endif else begin
			as_pix = 1.
			el_a = el[2]
			el_d = el[3]
			el_ast = 0
		endelse   
	openw,lun,aux_base+'_stv_ellipse_output.dat',/get_lun
	printf,lun,'# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix ): '+systime(0)
	printf, lun, el[0]*2.*as_pix, el[1]*2.*as_pix, el_a, el_d, el[4], $
		el[0]*2, el[1]*2, el[2], el[3], el_ast, as_pix, $
	  format = '(2f9.1,2f13.8,f7.1,4f9.2,i5,f9.5)'
	free_lun,lun   
	  
	      if file_exist(aux_base+'_stv_ellipse_output.dat') then begin
		  
		  rute = strmid(aux_base,0,strpos(aux_base,'aux/'))
		  ellipsefile = rute+'aux/'+id+'_ellipse.dat'
		  filestamp,ellipsefile
		  print, ellipsefile
		  spawn, 'cp '+aux_base+'_stv_ellipse_output.dat ' + ellipsefile
		  
	      endif
   endif   ; end save_ellipse

   
   
end

;--------------------------------------------------------------------

pro ds9_image

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base, aux_base, anulus_col, tv_filter, edit_region, has_counter, last_command, reg_size, lock_reg, ellipse_glga


  ds9_image = imfil
  
  if not strcmp(strtrim(tv_filter,2), 'w_com') and not strcmp(strtrim(tv_filter,2), 'n_com') and $ 
     ((strpos(dtype,'sdss') ge 0 or strpos(dtype,'panstarrs') ge 0))  then begin
     rute = strmid(imfil,0,strpos(imfil,'.fit'))
     rute =  strmid(rute, 0, strlen(rute)-2)
     
     if strcmp(strtrim(tv_filter,2), 'u') then tmp_image = rute+'_u.fits'
     if strcmp(strtrim(tv_filter,2), 'g') then tmp_image = rute+'_g.fits'
     if strcmp(strtrim(tv_filter,2), 'r') then tmp_image = rute+'_r.fits'
     if strcmp(strtrim(tv_filter,2), 'i') then tmp_image = rute+'_i.fits'
     if strcmp(strtrim(tv_filter,2), 'z') then tmp_image = rute+'_z.fits'

     if file_exist(tmp_image) then ds9_image=tmp_image
     
  endif

  
  if not strcmp(strtrim(tv_filter,2), 'w_com') and not strcmp(strtrim(tv_filter,2), 'n_com') and $ 
     strpos(dtype,'wise') ge 0 then begin
     rute = strmid(imfil,0,strpos(imfil,'.fit'))
     rute =  strmid(rute, 0, strlen(rute)-2)
     
     if strcmp(strtrim(tv_filter,2), 'w1') then tmp_image = rute+'w1.fits'
     if strcmp(strtrim(tv_filter,2), 'w2') then tmp_image = rute+'w2.fits'
     if strcmp(strtrim(tv_filter,2), 'w3') then tmp_image = rute+'w3.fits'
     if strcmp(strtrim(tv_filter,2), 'w4') then tmp_image = rute+'w4.fits'

     if file_exist(tmp_image) then ds9_image=tmp_image
     
  endif  
  
  
  if file_exist(ds9_image) then begin 
  
    ds9_reg_maker, 'idl_foo.reg', el, as_pix
    xc =  strtrim(string(el[2]),2)
    yc =  strtrim(string(el[3]),2)
   
    ds9_command = "ds9 " +  ds9_image
    ds9_command += " -height 800"
    ds9_command += " -width 800"
    ds9_command += " -scale log"
    ds9_command += " -scale minmax"
    ds9_command += " -zoom to fit"
    ds9_command += " -pan to "+xc+" "+yc
    ds9_command += ' -region load idl_foo.reg & '
    
    spawn, ds9_command
    
  endif

end

;--------------------------------------------------------------------



PRO stv, image, jimg, usefile, object, width, height,  $
	 log=log, min=min, max=max, same=same, $
         cursor=cursor, auto=auto, losig=losig, hisig=hisig, $
         hd=hd, true=true, cvfile=cvfile, block=block,  $
         ellipseinfo = ellipseinfo,  miscinfo = miscinfo, $
         asinh = asinh,  smooth = smooth, psrcfile = psrcfile, $
	 roifile = roifile, nobrowse = nobrowse, qa_stat = qa_stat, $
	 galex=galex,sdss=sdss,twomass=twomass,wise=wise,irac=irac, acs=acs, $
	 panstarrs=panstarrs, pick_obj = pick_obj, nofind = nofind, $
	 jpgbase=jpgbase,auxbase=auxbase,ellipse_glga_file=ellipse_glga_file

COMMON stvb, im, tim, sbase, draw, smin, smax, slog, sasinh, strue, ssmooth,$
	     have_astr, as_pix, astr, dtype, rfile, browse, qa, wsrc, plog, $
             keyboard_text_id, hdr, cfile, nsrc, psrc, maxp, el, elorig, misc, $
	     id, imfil, jpg_base,aux_base, anulus_col, tv_filter, edit_region, $ 
	     has_counter, last_command, reg_size, lock_reg, ellipse_glga

common stvc, flag_w, pnale_w
Common indices, last_band_ind, last_tv_filter

IF n_params() LT 1  THEN BEGIN 
 print,"stv, image, jpg_image, imfilename, object_id [, width, height, $ "
 print,"     log=log, min=min, max=max, same=same, $"
 print,"     cursor=cursor, auto=auto, losig=losig, hisig=hisig, $"
 print,"     hd=hd, true=true, cvfile=cvfile, psrcfile=psrcfile, block=block]"
 print," "
 print,"     width, height : of display window (600,600 default)"
 print,"     log: use logscale"
 print,"     min, max: scale vales"
 print,"     same: same values of log, min,max as last time"    
 print,"     auto: auto scale -3 to +7 sigma (ignores min,max,losig,hisig)"
 print,"     losig, hisig: set Nsigma scale values (ignores min,max)"
 print,"     true: true color image (buggy)"
 print,"     hd: fits header (useful for Curval coordinates etc.)"
 print,"     cursor: x cursor (24 default)"
 print,"     cvfile: name of file to send curval output to"
 print,"     psrcfile: name of existing list of point sources to mask"
 print,"     block: disable widget xmanager no_block feature"

 print,"keyboard:  'z': Zoom"
 print,"           'l': tvList"
 print,"           'q': Quit"
 print,"           'c': Curval (optional output appended to cvfile)"
 print,"           'i': Invert colortable"
 print,"           'h': Help"
 print,"           'n': ceNtroid with mouse (appends to file aux_base+'stv_cntrd_output.dat')"

 GOTO, finish
ENDIF 

IF NOT keyword_set(cursor) THEN cursor = 24
IF NOT keyword_set(width) THEN width = 600
IF NOT keyword_set(height) THEN height = 600

cfile = '/dev/tty'
IF keyword_set(cvfile) THEN cfile = cvfile

maxp = 90000L
psrc = fltarr(7,maxp)	; point sources: a,d,x,y,r_asec,astr?,use?
nsrc = fltarr(7,maxp)

imfil = usefile
id = object

plog = pick_obj

hdr = ''
IF keyword_set(hd) THEN hdr = hd

if n_elements(hdr) ge 1 then begin
	extast,hdr,astr
	getrot,hdr,rot,cdelt
	as_pix = abs(cdelt[0])*3600.
	have_astr = 1
endif else begin
	astr = -1
	as_pix = 1.
	have_astr = 0
endelse

dtype='galex'
if keyword_set(sdss) then dtype='sdss'
if keyword_set(panstarrs) then dtype='panstarrs'
if keyword_set(acs) then dtype='acs '
if keyword_set(twomass) then dtype='2mass'
if keyword_set(wise) then dtype='wise'
if keyword_set(irac) then dtype='irac'

if keyword_set(psrcfile) then read_psrc, psrcfile, nps

wsrc = (1 eq 0)	; default to not write out point source file

if keyword_set(roifile) then begin
	if file_test(roifile) then $
		rfile = roifile $
	else	rfile = ''
endif else rfile = ''

if keyword_set(nobrowse) then $
	browse = 0 $
else	browse = 1

if keyword_set(qa_stat) then $
	qa = qa_stat $
else	qa = {qa_stat}

el = [0,0,0,0,0,0,0]
elorig = [0,0,0,0,0,0,0]
IF keyword_set(ellipseinfo) THEN begin
 el = ellipseinfo
 elorig = ellipseinfo
endif

misc = ''

image_info = size(image)

IF image_info[0] NE 2 AND NOT keyword_set(true) THEN BEGIN
   print,'ERROR: image must be 2 dimensional'
   GOTO ,finish
ENDIF

imx=image_info[1]
imy=image_info[2]

image1=image
rmin=min(image1)
rmax=max(image1)

ssmooth = 1
IF keyword_set(smooth) then ssmooth = smooth

IF n_elements(max) EQ 0 THEN max = rmax
IF n_elements(min) EQ 0 THEN min = rmin

IF keyword_set(losig) AND keyword_set(hisig) THEN BEGIN
 sky,image1, sky, sigma,/silent
 max = sky + hisig*sigma
 min = sky + losig*sigma
endif

IF keyword_set(auto) THEN BEGIN
 losig=-3
 hisig=7
 sky,image1, sky, sigma,/silent
 max = sky + hisig*sigma
 min = sky + losig*sigma
endif

IF  (xregistered ('stv')) AND keyword_set(same) THEN BEGIN
 min=smin
 max=smax
 log=slog
ENDIF

smin=min
smax=max
slog=keyword_set(log)
strue=keyword_set(true)
sasinh = keyword_set(asinh)
im=image
tim=jimg
jpg_base = jpgbase
aux_base = auxbase

tv_filter = 'n_com'
last_tv_filter = 'n_com'
last_band_ind = 0
edit_region = aux_base + '_stv_roi_output.dat'
ellipse_glga = ellipse_glga_file

if file_exist(aux_base+'_stv_usesex_vertices.dat') then spawn, 'rm ' + aux_base+'_stv_usesex_vertices.dat'
if file_exist(aux_base+'_stv_forcesex_vertices.dat') then spawn, 'rm ' + aux_base+'_stv_forcesex_vertices.dat'
if file_exist(aux_base+'_stv_forcesex_output.dat') then spawn, 'rm ' + aux_base+'_stv_forcesex_output.dat'
if file_exist(aux_base+'_stv_usesex_output.dat') then spawn, 'rm ' + aux_base+'_stv_usesex_output.dat'
rute = strmid(imfil,0,strpos(imfil,'.fit'))
rute =  strmid(rute, 0, strlen(rute)-2)
spawn, 'rm '+ rute+'*tmp.fits'
     
     
anulus_col = 'Y'   ; default anulus color
if (xregistered ('stv')) THEN BEGIN 

 if nps le 0 and browse eq 0 then $
	stv_find,/noprompt $
 else	stv_display
 curs, cursor

ENDIF ELSE BEGIN ; otherwise create widget
 
 
 ;---------------------------------------------------------------
; Create a base

 ; This generates the main window
 baseID = WIDGET_BASE( COLUMN=2, $
       MBAR=wMBarBase, xoffset = 760, yoffset = 0, title='Obj: ' + object)

 wFileMenu = Widget_Button( wMBarBase, VALUE='File', /MENU)
 wMeasure = Widget_Button( wFileMenu, VALUE='Measure >', UVALUE='measure')
 wQuitButton = Widget_Button( wFileMenu, VALUE='Refresh TV', UVALUE='refresh')
 wquit = Widget_Button( wFileMenu, VALUE='quit-next', UVALUE='quit') 
 wsave = Widget_Button( wFileMenu, VALUE='w-save', UVALUE='save') 
 wSkip = Widget_Button( wFileMenu, VALUE='Skip ...', UVALUE='skip') 
 wPrevious = Widget_Button( wFileMenu, VALUE='Previous', UVALUE='previous') 
 wQuitButton = Widget_Button( wFileMenu, VALUE='Quit', UVALUE='exit') 
 
 wReg = Widget_Button( wMBarBase, VALUE='Region', /MENU)  ; 
 wRegUndo  = Widget_Button( wReg, VALUE='Undo', UVALUE='reg_undo') 
 wRegAdd  = Widget_Button( wReg, VALUE='Add ...', UVALUE='reg_add')    
 wRegRemove  = Widget_Button( wReg, VALUE='Remove one', UVALUE='reg_remove')    
 wRegClean = Widget_Button( wReg, VALUE='Clean all', UVALUE='reg_clean')    
 wRegSave  = Widget_Button( wReg, VALUE='Save', UVALUE='reg_save') 
 wRegChange  = Widget_Button( wReg, VALUE='Change Mask', UVALUE='reg_change')    
   

 wELLIPSE     = Widget_Button( wMBarBase, VALUE='Ellipse', /MENU) ;   
 wELLedit     = Widget_Button( wELLIPSE, VALUE='Edit ...', UVALUE='ell_edit')  
 wELLoriginal = Widget_Button( wELLIPSE, VALUE='Original', UVALUE='ell_orig')  
 wELLcatal = Widget_Button( wELLIPSE, VALUE='Catalog', UVALUE='ell_catal')  
 wELLsdss     = Widget_Button( wELLIPSE, VALUE='SDSS', UVALUE='ell_sdss')  
 wELLsave     = Widget_Button( wELLIPSE, VALUE='SAVE', UVALUE='ell_save')  
   
 wStar = Widget_Button( wMBarBase, VALUE='Stars', /MENU)  ;  
 wStarFind  = Widget_Button( wStar, VALUE='Find ...', UVALUE='star_find')   
 wStarDel  = Widget_Button( wStar, VALUE='Delete ...', UVALUE='star_del')  
 wStarEdit  = Widget_Button( wStar, VALUE='Edit/Add ...', UVALUE='star_edit')  
   

   
 wBACK = Widget_Button( wMBarBase, VALUE='Background', /MENU)
 genBACK = Widget_Button( wBACK, VALUE='Back. Region', UVALUE='gen_back') 
 genBACK = Widget_Button( wBACK, VALUE='Clear Region', UVALUE='remove_back')    

 
 wDS9 = Widget_Button( wMBarBase, VALUE='ds9', /MENU)  ; 
 wDS9open  = Widget_Button( wDS9, VALUE='Open', UVALUE='ds9_open') 
 wDS9Cont  = Widget_Button( wDS9, VALUE='Contours on/off', UVALUE='ds9_contour')  
 wDS9Imp  = Widget_Button( wDS9, VALUE='Import ellipse', UVALUE='ds9_imp')  
 wDS9Exp  = Widget_Button( wDS9, VALUE='Export ellipse', UVALUE='ds9_exp') 
 wDs9Close  = Widget_Button( wDS9, VALUE='Close all', UVALUE='ds9_close')    
 


 
 wSEX = Widget_Button( wMBarBase, VALUE='SExtractor', /MENU)
 runSEX = Widget_Button( wSEX, VALUE='Run', UVALUE='r_sex') 
 segSEX = Widget_Button( wSEX, VALUE='Segmentation', UVALUE='seg_sex')  
 segSEX = Widget_Button( wSEX, VALUE='Masked Image', UVALUE='tmp_sex') 
 regSEXon = Widget_Button( wSEX, VALUE='Region on', UVALUE='reg_sex_on') 
 regSEXoff = Widget_Button( wSEX, VALUE='Region off', UVALUE='reg_sex_off')
 useSEXoff = Widget_Button( wSEX, VALUE='Use-Mask', UVALUE='use_sex')
 forceSEXoff = Widget_Button( wSEX, VALUE='Force-Mask', UVALUE='force_sex')
 resetSEXoff = Widget_Button( wSEX, VALUE='Clear-Mask', UVALUE='clear_sex')
 
 
 wElliptcon = Widget_Button( wMBarBase, VALUE='Pylipse', /MENU)
 openEllipticon = Widget_Button( wElliptcon, VALUE='Open', UVALUE='open_pyl') 
 reload_ell = Widget_Button( wElliptcon, VALUE='Load Ellipse', UVALUE='load_pyl')  
 
 

 
 wHelpMenu = Widget_Button( wMBarBase, VALUE='Help', /HELP)
 wHelpButton = Widget_Button( wHelpMenu, VALUE='Help', UVALUE='Help_Program')

 


wLRow = WIDGET_BASE(baseID,  /COLUMN, /BASE_ALIGN_LEF, space=5) 


 values = ['not available']
 if strcmp(strtrim(dtype), 'sdss') then begin 
    values = ['u', 'g', 'r', 'i', 'z', 'gri', 'urz'] 
    set_val = 5 ; initialize to 'gri'
    endif
 if strcmp(strtrim(dtype), 'panstarrs') then begin
    values = ['g', 'r', 'i', 'z', 'gri'] 
    set_val = 4 ; initialize to 'gri'
    endif
 if strcmp(strtrim(dtype), 'wise') then begin
    values = ['w1', 'w2', 'w3', 'w4', 'w123', 'w124'] 
    set_val = 4 ; initialize to 'w123'
    endif

infowidget = widget_base(wLRow, COLUMN=2, /frame, uname='info_widg') 

  infowidget_ = widget_base(infowidget,Row=2, /frame, uname='info_widg_')  

  bgroup1 = CW_BGROUP(infowidget_, values, UVALUE='filter', /Column, /EXCLUSIVE, $ 
    LABEL_TOP='Filter', /FRAME, SET_VALUE=set_val) 
    
  infowidget__ = widget_base(infowidget_,/col, /frame, uname='info_widg__')  
  button16 = WIDGET_BUTTON(infowidget__, VALUE='Refresh', UVALUE='^')
  
  info_label = widget_base(infowidget,/col, /frame) 
;   info_label1  = WIDGET_Label(info_label, value = 'Obj/Survey:') 
  info_label2  = WIDGET_Label(info_label, value = object) 
  info_label4  = WIDGET_Label(info_label, value = dtype)
  info_label3  = WIDGET_Label(info_label, value = '---------')
  info_label5  = WIDGET_Label(info_label, value = 'X:        ', /ALIGN_LEFT, uname='cor_x')
  info_label6  = WIDGET_Label(info_label, value = 'Y:         ', /ALIGN_LEFT, uname='cor_y')
  info_label7  = WIDGET_Label(info_label, value = 'RA :               ', /ALIGN_LEFT, uname='cor_ra')
  info_label70 = WIDGET_Label(info_label, value = '              ', /ALIGN_LEFT, uname='cor_ra0')
  info_label8  = WIDGET_Label(info_label, value = 'DEC:               ', /ALIGN_LEFT, uname='cor_dec')
  info_label80 = WIDGET_Label(info_label, value = '               ', /ALIGN_LEFT, uname='cor_dec0')
  info_label9  = WIDGET_Label(info_label, value = 'Vlue:               ', /ALIGN_LEFT, uname='im_value')
  info_label10 = WIDGET_Label(info_label, value = 'Flux:               ', /ALIGN_LEFT, uname='im_flux')
  info_but0 = WIDGET_BUTTON(info_label, VALUE='lock', UVALUE='lock', uname='lock_info')
  info_lock = (1 eq 0)
  
colBase = widget_base(wLRow,COLUMN=2, /frame, uname='col_widg')  

  colBase00 = widget_base(colBase,/column, /frame) 
  button0 = WIDGET_BUTTON(colBase00, VALUE='Change Mask ($)', UVALUE='$')
  label0  = WIDGET_Label(colBase00, value = 'Current: All Bands', uname='label')
  button1 = WIDGET_BUTTON(colBase00, VALUE='Add Region (r)', UVALUE='r')
  button2 = WIDGET_BUTTON(colBase00, VALUE='Undo Region (~)', UVALUE='~')
  button3 = WIDGET_BUTTON(colBase00, VALUE='Remove a Region (!)', UVALUE='!')
  button4 = WIDGET_BUTTON(colBase00, VALUE='Clean all regions (C)', UVALUE='C')
  button45 = WIDGET_BUTTON(colBase00, VALUE='Save regions', UVALUE='save_reg')
  
  
  colBase_ = widget_base(colBase,ROW=2, /frame, uname='col_widg_')   
  bgroup2 = CW_BGROUP(colBase_, ['5','7','10','15','20','25','30','35'], UVALUE='reg_circle', /Column, /EXCLUSIVE, $ 
    LABEL_TOP='reg size', /FRAME, SET_VALUE=0) 
  reg_size = 5
  colBase__ = widget_base(colBase_,/row, /frame, uname='col_widg__')   
  button_lock_reg = WIDGET_BUTTON(colBase__, VALUE='lock reg', UVALUE='lock_reg', UNAME='lockreg')
  lock_reg = 0
  
colBase2 = widget_base(wLRow,/col, /frame, uname='col2_widg') 

  button5 = WIDGET_BUTTON(colBase00, VALUE='Find Stars (f)', UVALUE='f')
  button6 = WIDGET_BUTTON(colBase00, VALUE='Delete Centeroid (d)', UVALUE='d')
  button7 = WIDGET_BUTTON(colBase00, VALUE='Edit/Add Centeroid (c)', UVALUE='c')

wLRow3 = widget_base(wLRow,COLUMN=2, /frame, uname='col3_widg') 
  
  colBase3 = widget_base(wLRow3,/column, /frame) 
  button8 = WIDGET_BUTTON(colBase3, VALUE='ds(9)', UVALUE='9')
  button85 = WIDGET_BUTTON(colBase3, VALUE='Contours On', UVALUE='\', uname='contour')
  button86 = WIDGET_BUTTON(colBase3, VALUE='Close all ds9', UVALUE='c_ds9')
  has_counter = (1 eq 0)
  button9 = WIDGET_BUTTON(colBase3, VALUE='Import Ellipse (*)', UVALUE='*', uname='import')
  button10 = WIDGET_BUTTON(colBase3, VALUE='Export Ellipse (=)', UVALUE='=', uname='export')
  button101 = WIDGET_BUTTON(colBase3, VALUE='Run SExtractor', UVALUE='sex', uname='sextract')

  
  colBase31 = widget_base(wLRow3, /column,/frame)
  button150 = WIDGET_BUTTON(colBase31, VALUE='+', UVALUE='9\', xsize=30, ysize=55)
  button160 = WIDGET_BUTTON(colBase31, VALUE='CL-reg', UVALUE='clr')
  button161 = WIDGET_BUTTON(colBase31, VALUE='Cntrd', UVALUE='add_centroid')
  button161 = WIDGET_BUTTON(colBase31, VALUE='-->', UVALUE='import_reg')
  
  
wLRow4 = widget_base(wLRow,COLUMN=3, /frame, uname='col4_widg') 
  
  colBase4 = widget_base(wLRow4,/column, /frame) 
  button12 = WIDGET_BUTTON(colBase4, VALUE='Edit Ellipse (e)', UVALUE='e')	
;   button11 = WIDGET_BUTTON(colBase4, VALUE='Original Ellipse (#)', UVALUE='#')

  colBase41 = widget_base(wLRow4, /column,/frame)
  button121 = WIDGET_BUTTON(colBase41, VALUE='-', UVALUE='e-')
  
  colBase42 = widget_base(wLRow4, /column,/frame)
  button122 = WIDGET_BUTTON(colBase42, VALUE='+', UVALUE='e+')  


wLRow1 = widget_base(wLRow, COLUMN=2, /frame, uname='wLRow1_widg') 
  
  colBase6 = widget_base(wLRow1,/column, /frame) 
  
  button17 = WIDGET_BUTTON(colBase6, VALUE='quit-next (q)', UVALUE='q')
  button17 = WIDGET_BUTTON(colBase6, VALUE='SAVE', UVALUE='save')
  button19 = WIDGET_BUTTON(colBase6, VALUE='SKIP (x)', UVALUE='x')
;   button20 = WIDGET_BUTTON(colBase6, VALUE='PREVIOUS', UVALUE='previous')
  button18 = WIDGET_BUTTON(colBase6, VALUE='Quit-Exit (Q)', UVALUE='Q')
  
  colBase7 = widget_base(wLRow1, /column,/frame)
  button15 = WIDGET_BUTTON(colBase7, VALUE='>', UVALUE='m', xsize=30, ysize=100)

;---------------------------------------------------------------
 ; Create the Draw widget
 baseRow1 = WIDGET_BASE(baseID, /COLUMN , /BASE_ALIGN_RIGHT) 
 
 base = WIDGET_BASE(baseRow1, RESOURCE_NAME = 'stv',title='stv '+id,uname='stv',$
                    xoffset = 760, yoffset = 0)
 draw = widget_draw(base, xsize=imx, ysize=imy, /scroll,/no_copy, $
          x_scroll_size=width < imx, y_scroll_size= height < imy, $
	  /button_events) 

 keyboard_text_id = widget_text(base, $
                                     /all_events, $
                                     scr_xsize = 1, $
                                     scr_ysize = 1, $
                                     units = 0, $
                                     uvalue = 'keyboard_text', $
                                     value = '')
   
; ;  WIDGET_CONTROL, base, /realize
 WIDGET_CONTROL, baseID, /realize
 WIDGET_CONTROL, draw, get_value=SLIDE_WINDOW
 WIDGET_CONTROL, base, SET_UVALUE=slide_window
 WIDGET_CONTROL, keyboard_text_id, event_pro = 'stv_keyboard_event'
 WIDGET_CONTROL, wLRow, event_pro = 'panel_event'
 WSET, SLIDE_WINDOW

 sbase  = baseID
 flag_w  = 0 
 pnale_w = wLRow
 last_command = ''
 
 if not keyword_set(nofind) and nps le 0 and browse eq 0 then $
	stv_find,/noprompt $
 else	stv_display
 curs, cursor

 stv_printhelpline
 
 

 
 
  
 IF keyword_set(block) THEN  BEGIN 
  XMANAGER, 'stv', baseID, event_handler='stv_event', $
          cleanup = 'stv_cleanup'
 ENDIF ELSE BEGIN 
  XMANAGER, 'stv', baseID, event_handler='stv_event',/no_block, $
           cleanup = 'stv_cleanup'
 ENDELSE 

 ENDELSE 

miscinfo = misc
qa_stat = qa
pick_obj = plog

finish:
END

