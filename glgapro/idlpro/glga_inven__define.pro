pro glga_inven__define
;+
; glga_inven__define - define glga inventory structure
;
; NOTE: remember to keep structure name on a line by itself 
;	(see struct_init.pro)
;
; Catalog data are first followed by analysis data
; Analysis is tracked as follows:
;	0 - initial state
;	1 - image data assembled
;	2 - integrated magnitude and profile derived
;
; see structure snhphsrc for host photometry sources:
;	higher source numbers are more reliable
;
;-
tmp = { glga_inven, $
	id:'', $		; ID
	catalog:'', $		; catalog: glga_v1, glga_v2, or extra
	ra:-9.d0, dec:-99.d0, $	; coords (J2000)
	majax:-9., minax:-9., $	; major/minor axes (arcsec)
	pa:-99.0, $		; PA (degrees)
	type:'', $		; type
	elfile: 0, $		; ellipse.dat file? 0, 1 - no, yes
	dss_imgs: 0, $		; dss images? 0,1 - no, yes
	galex_imgs: 0, $	; galex images? 0,1 - no, yes
	galex_qa: 0, $		; galex qa status? 0,1,2 - no, yes, done
	sdss_imgs: 0, $		; sdss images? 0,1 - no, yes
	sdss_qa: 0, $		; sdss qa status? 0,1,2 - no, yes, done
	twomass_imgs: 0, $	; 2mass images? 0,1 - no, yes
	twomass_qa: 0, $	; twomass qa status? 0,1,2 - no, yes, done
	wise_imgs: 0, $		; wise images? 0,1 - no, yes
	wise_qa: 0, $		; wise qa status? 0,1,2 - no, yes, done
	irac_imgs: 0, $		; irac images? 0,1 - no, yes
	irac_qa: 0, $		; irac qa status? 0,1,2 - no, yes, done
	mod_time:0.d0 $		; time stamp in seconds systime(1)
	}
end
