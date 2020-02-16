!quiet=1
device,retain=2
udir = 'Users'	; default for Mac OS X
if strpos(!version.os,'linux') ge 0 then udir = 'users'
loginfo = get_login_info()
unam = loginfo.user_name
;
astrolib
;
defsysv,'!top_color',0
defsysv,'!top_colorx',0
defsysv,'!top_colorps',255
defsysv,'!colorstr',''
defsysv,'!colorstr_ps',''
defsysv,'!colorstr_x',''
;
defsysv,'!PHYS_C',2.99792458d5  ; speed o' light
defsysv,'!COSMO_OM',0.27	; Omega matter
defsysv,'!COSMO_OL',0.73	; Omega Lambda
defsysv,'!COSMO_H0',73.0	; Hubble constant

defsysv,'!GLGA_ROOT','/home/ehsan/PanStarrs/'



!quiet=0
