function glga_getextin,ebv,band,yuan13=yuan13
;+
;	glga_getextin - return the extinction for the given band
;-
	extin=0.
	;
	; did we select Yuan et al. 2013 Table 2 values?
	; use column 3 except for Johnson and GALEX
	if keyword_set(yuan13) then begin
		case strupcase(band) of
			'FUV': extin = ebv*6.892
			'NUV': extin = ebv*6.738	; column 5 (CCM)
			'UJ': extin = ebv*5.434
			'BJ': extin = ebv*4.315
			'VJ': extin = ebv*3.315	; Johnson SFD98 values
			'U': extin = ebv*4.35
			'G': extin = ebv*3.31
			'R': extin = ebv*2.32
			'I': extin = ebv*1.72
			'Z': extin = ebv*1.28
			'J': extin = ebv*0.72
			'H': extin = ebv*0.46
			'K': extin = ebv*0.306
			'W1': extin = ebv*0.19
			'W2': extin = ebv*0.15
			'W3': extin = ebv*0.
			'W4': extin = ebv*0.	; unknown
			'3P6UM': extin = ebv*0.
			'4P5UM': extin = ebv*0.	; unknown
			else: print,'Unknown band: ',band
		endcase
	;
	; just use canonical values
	endif else begin
		case strupcase(band) of
			'FUV': extin = ebv*8.24	; Wyder et al. 2007 GALEX
			'NUV': extin = ebv*8.24 - 0.67*ebv^2
			'UJ': extin = ebv*5.434
			'BJ': extin = ebv*4.315
			'VJ': extin = ebv*3.315	; Johnson SFD98 values
			'U': extin = ebv*5.155
			'G': extin = ebv*3.793
			'R': extin = ebv*2.751
			'I': extin = ebv*2.086
			'Z': extin = ebv*1.479	; original SFD98 values
			'J': extin = ebv*0.900
			'H': extin = ebv*0.576
			'K': extin = ebv*0.365	; original SFD98 values
			'W1': extin = ebv*0.171
			'W2': extin = ebv*0.096	; CCM values
			'W3': extin = ebv*0.
			'W4': extin = ebv*0.	; unknown
			'3P6UM': extin = ebv*0.
			'4P5UM': extin = ebv*0.	; unknown
			else: print,'Unknown band: ',band
		endcase
	endelse
	;
	return,extin
end
