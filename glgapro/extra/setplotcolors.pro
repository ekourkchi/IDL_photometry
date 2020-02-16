PRO setplotcolors, background, foreground, basetable, list=list, $
                   test=test, gray=gray, colors=colors, low=low

;syntax: setplotcolors [background, foreground, basetable, /list, $
;                      /test, gray=, colors=, /low]
;
;Parameters:
; background = color to set background to
; foreground = color to set foreground to
; basetable = colortable to load and modify (default=0)
; 
;Keywords:
; list   - print list of color names available
; test   - plot colors available in a window
; gray   - number of grayscale colors to establish
; colors - name of long array to store index values of colors
; low    - put colors at low end of color table (24 bit contour plots)
;          do not use low of using 256 or less colors.
;
;Description:
; Loads a base color table (linear B&W color table (0) is default)
; and then modifies the color table starting with !d.n_colors-2
; (unless low keyword is set) from the values specified in the 'd'
; structure. To add or remove basic colors, modify 'd' structure.
;
; Sets system variables to the index of each color specified.
; i.e. !red = color index for red etc. If grayscale option is selected
; the color names are gray1, gray2 etc. (use test or list to view names)
;
;Example:
;  setplotcolors <- setup the colors from setplotcolors_rgb.dat
;  setplotcolors, /test <- setup the colors and plots them
;  setplotcolors, !gray, /test <- plots colors on gray background
;
;  setplotcolors, !gray,!purple  
;  plot,[0,0],[1,1] ; gray and purple plots...why not
;  oplot,[0.5],[0.5],psym=4, symsize=3, thick=3,color=!orange
;
; set up 10 shades of gray for a contour plot:
;  pcolors=lonarr(10)
;  setplotcolors,gray=10,colors=pcolors,/low
;  (contour does not like high end of color table in 24 bit mode)
;  contour,...,c_colors=pcolors  
;
; to remove color table modifications reload a color table:
;  loadct,0
;
;
; Mark Seibert (July,2000)
;
; Now supports 24 bit true color (MS 3/01)
;
; Added grayscale, low & colors keyword options (MS 6/01)
;
; converted external data table to internal structure. add background
; and foreground keywords (MS 7/03)


;define color names & RGB values (0-255)

c_R =[255,  0,112,211,190,169,173,  0,  0,144,  0,  0,144,  0,  0,255,255,255,210,165,255,255,255,139,255,255,255,255,139,160]

c_G=[255,  0,128,211,190,169,216,  0,  0,255,255,139,238,255,190,255,255,215,180, 42,165,140,  0,  0,182,192, 20,  0,  0, 32]

c_B=[255,  0,144,211,190,169,230,255,139,255,255,139,144,  0,  0,100,  0,  0,140, 42,  0,  0,  0,  0,193,203,147,255,139,240]

c_NAME=['white','black','slate','lgray','gray','dgray','lblue','blue','dblue','lcyan','cyan','dcyan','lgreen','green','dgreen','lyellow','yellow','dyellow','tan','brown','orange','dorange','red','dred','lpink','pink','dpink','magenta','dmagenta','purple']

IF keyword_set(gray) THEN BEGIN
 c_R = indgen(gray+2)
 c_name=strarr(gray+2)
 c_R[0]=255
 c_R[1]=0
 c_name[0]='white'
 c_name[1]='black'
 FOR i=2,n_elements(c_R)-1 DO BEGIN
    c_R[i]=200 - fix(200./gray)*(i-2)
    c_name[i]='gray'+strn(i-1)
 endfor  
 c_G = c_R
 c_B = c_R
ENDIF

IF keyword_set(colors) THEN colors=lonarr(n_elements(c_name))

;load color table & get RGB values
table = 0 ;default B&W linear
IF n_params() gt 2 THEN table = basetable
loadct,table,/silent
TVLCT, R, G, B, /get

n_colors = !d.n_colors < 256
IF keyword_set(low) THEN n_colors = n_elements(c_R)+2

;setup the colors & system variables specified above
FOR i = 0, n_elements(c_name)-1 DO begin
 cvar = '!'+c_name[i]
; defsysv,cvar,n_colors-(2+i)  
 R[n_colors-(2+i)]=c_R[i] 
 G[n_colors-(2+i)]=c_G[i]
 B[n_colors-(2+i)]=c_B[i]

 IF !d.n_colors LE 256 THEN BEGIN
   defsysv,cvar,n_colors-(2+i) ; 8 bit mode
   IF keyword_set(colors) THEN colors[i]=n_colors-(2+i)
 ENDIF ELSE BEGIN ; 24 bit mode
  cval=n_colors-(2+i) + 256L * (n_colors-(2+i) + 256L * n_colors-(2+i))
  defsysv,cvar,cval 
  device,decomposed=0
  IF keyword_set(colors) THEN colors[i]=cval
 endelse

ENDFOR

;modify the color table
TVLCT, R, G, B

IF keyword_set(list) THEN BEGIN
print
print,'The following colors may be accesed via system variables:'
print
print,c_name
print
print,'!colorname = colortable index for colorname'
print,'l = light; d = dark'
print,'Ex: plot,x,y,color=!red'
print
endif

IF n_params() GT 0 THEN !p.background=background
IF n_params() GT 1 THEN !p.color=foreground

IF keyword_set(test) THEN BEGIN
window,6,xsize=100.*ceil(n_elements(c_name)/10.),ysize=170,$
       title='Plot Colors Added'
FOR i = 0,n_elements(c_name)-1 DO BEGIN
 x = 10 + 10*(i - (i MOD 10))
 y = 10+((i MOD 10)*15.)
 IF !d.n_colors LE 256 THEN BEGIN 
  xyouts,x,y,c_name[i],color=n_colors-(2+i),$
        charsize=1.25,charthick=1.25,/device
 ENDIF ELSE BEGIN
  cval=n_colors-(2+i) + 256L * (n_colors-(2+i) + 256L * n_colors-(2+i))  
  xyouts,x,y,c_name[i],color=cval,$
        charsize=1.25,charthick=1.25,/device
 ENDELSE

ENDFOR
ENDIF
END


