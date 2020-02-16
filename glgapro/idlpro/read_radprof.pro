pro read_radprof, id, band, $
	ntot_a, ntot_mag, ntot_mag_e, nann_mu, nann_mu_e, $
	nra_cen, ndec_cen, nsemimajor, nsemiminor, npa, nscale, $
	ntf_mag, ntf_mag_e, naf_mag, naf_mag_e, nmu_bg, nmu_bg_e, $
	nskyradius_in, nskyradius_out, ntf_a, asymag, asymag_e, asyma, $
	pathtoprofile=pathtoprofile, silent=silent, $
	nbg_flx=nbg_flx, nbg_e=nbg_e, tot_int=ntot_int, phot_ts=phot_ts, $
	annuli_size=annuli_size, int_to_mjy=int_to_mjy

;
;
if not keyword_set(pathtoprofile) then pathtoprofile='./'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; files names from radprof

ntotfile =pathtoprofile+'/'+id+'_'+band+'_totprofile.dat'
nannfile =pathtoprofile+'/'+id+'_'+band+'_annprofile.dat'
nbgfile  =pathtoprofile+'/'+id+'_'+band+'_background.dat'
nellfile =pathtoprofile+'/'+id+'_'+band+'_ellipsepar.dat'
ntotflxfile =pathtoprofile+'/'+id+'_'+band+'_total.dat'
asymflxfile =pathtoprofile+'/'+id+'_'+band+'_asymptotic.dat'
naperflxfile =pathtoprofile+'/'+id+'_'+band+'_aperture.dat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get photometry file timestamps
totfi = file_info(ntotfile)
annfi = file_info(nannfile)
phot_ts = max([totfi.mtime,annfi.mtime])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; initialize
ntot_a=!values.f_nan
ntot_mag=!values.f_nan
ntot_mag_e=!values.f_nan
nann_mu=!values.f_nan
nann_mu_e=!values.f_nan
nra_cen=!values.f_nan
ndec_cen=!values.f_nan
nsemimajor=!values.f_nan
nsemiminor=!values.f_nan
npa=!values.f_nan
nscale=!values.f_nan
ntf_mag=!values.f_nan
ntf_mag_e=!values.f_nan
naf_mag=!values.f_nan
naf_mag_e=!values.f_nan
nmu_bg=!values.f_nan
nmu_bg_e=!values.f_nan
nskyradius_in=!values.f_nan
nskyradius_out=!values.f_nan
ntf_a=!values.f_nan
asymag=!values.f_nan
asymag_e=!values.f_nan
asyma=!values.f_nan

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; annuli size and mJy conversion?

annuli_size=-9.99
int_to_mjy=-9.99
if file_exist(nannfile) then begin
	openr,ilun,nannfile,/get_lun
	rec=''
	while not eof(ilun) do begin
		readf,ilun,rec
		if strpos(rec,'AnnSz:') ge 0 then begin
			for i=0,2 do $
				val = gettok(rec,' ')
			annuli_size = float(val)
		endif
		if strpos(rec,'mJy:') ge 0 then begin
			for i=0,2 do $
				val = gettok(rec,' ')
			int_to_mjy = float(val)
		endif
	endwhile
	free_lun,ilun
endif

;;;;;;;;;;;;;;;;;;;;;;;
; read in profile data

if file_exist(ntotfile) and file_exist(nannfile) and $
   file_exist(nbgfile) and file_exist(nellfile) then begin

 readcol, ntotfile, ntot_a, ntot_mag, ntot_mag_e, ntot_int, ntot_int_e, $
          ntot_int_e_cnt, ntot_int_e_bg, ntot_bg, ntot_area,/silent

 readcol, ntotflxfile, ntf_a, ntf_mag, ntf_mag_e, ntf_int, ntf_int_e, $
          ntf_int_e_cnt, ntf_int_e_bg, ntf_bg, ntf_area,/silent

 readcol, asymflxfile, asyma, asymag, asymag_e, nsf_int, nsf_int_e, $
          nsf_int_e_cnt, nsf_int_e_bg, nsf_bg, nsf_area,/silent

 readcol, naperflxfile, naf_a, naf_mag, naf_mag_e, naf_int, naf_int_e, $
          naf_int_e_cnt, naf_int_e_bg, naf_bg, naf_area,/silent

 readcol, nannfile, nann_a, nann_mu, nann_mu_e, nann_int, nann_int_e, $
          nann_int_e_cnt, nann_int_e_bg, nann_bg, nann_area,/silent

 readcol, nellfile, nra_cen, ndec_cen, nsemimajor, nsemiminor, npa,/silent

 readcol, nbgfile, nbg_flx, nbg_e, nmu_bg, nmu_bg_e, nscale, $
          nskyradius_in, nskyradius_out,/silent

endif else begin
	if not keyword_set(silent) then $
		print,'No '+band+' profiles found for: ',id
	ntot_a = -1.0
	return
endelse

;;;;;;;;;;;;
; all done

return
end
