function glga_getimplotscl,band
;+
;	glga_getimplotscl - return the scale factor for plotting images
;-
	scale=1.
	;
	case strupcase(band) of
		'FUV': scale = 390.
		'NUV': scale = 100.
		'UJ': scale = 1.
		'BJ': scale = 1.
		'VJ': scale = 1.
		'U': scale = 590.
		'G': scale = 1200.
		'R': scale = 590.
		'I': scale = 590.
		'Z': scale = 500.
		'J': scale = 30.
		'H': scale = 20.
		'K': scale = 20.
		'W1': scale = 0.6
		'W2': scale = 1.20
		'W3': scale = 0.22
		'W4': scale = 0.09
		'3P6UM': scale = 300.
		'4P5UM': scale = 400.
		else: print,'Unknown band: ',band
	endcase
	;
	return,scale
end
