PRO rct

TVLCT, R, G, B, /get

;r_new=[r[0],reverse(r[1:n_elements(r)-2]),r[n_elements(r)-1]]
;g_new=[g[0],reverse(g[1:n_elements(g)-2]),g[n_elements(g)-1]]
;b_new=[b[0],reverse(b[1:n_elements(b)-2]),b[n_elements(b)-1]]

;iF !d.n_colors LE 256 THEN 
;tvlct,r_new,g_new,b_new

tvlct,reverse(r),reverse(g),reverse(b) 

end
