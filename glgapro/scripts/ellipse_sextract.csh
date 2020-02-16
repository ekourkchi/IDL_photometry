#!/usr/bin/csh

setenv  pgc_no $argv[1]   # galaxy ID, PGC number
setenv  filter $argv[2]   # filter 
setenv  location $argv[3] # where is the fits file

setenv sextract $argv[4]      # where to find necessary executables
setenv bin $argv[5] # where is SExtractor files

setenv xc $argv[6]        # center of the galaxy (x,y)
setenv yc $argv[7]

setenv pixel_scael $argv[8]   # arcsec / pix
setenv optional_fname $argv[9]   # optional file name


setenv root ${location}'/'${pgc_no}'_'${filter}

if ( -f ${optional_fname} ) then
   setenv filename ${optional_fname}
else 
   
   setenv filename ${root}'.fits'
   
endif


setenv cataloge ${root}'.cat'
setenv segment ${root}'.segment.fits'



if ( -f ${filename} ) then
    
    
    setenv backsize 512
    


    sex ${filename} -c ${sextract}'/default.sex' \
    -DETECT_THRESH   0.9  \
    -ANALYSIS_THRESH  0.9  \
    -PARAMETERS_NAME ${sextract}'/default_aper.param'	\
    -CATALOG_NAME ${cataloge} \
    -CATALOG_TYPE	ASCII_HEAD \
    -STARNNW_NAME	${sextract}'/default.nnw' \
    -CHECKIMAGE_TYPE  SEGMENTATION  \
    -CHECKIMAGE_NAME   $segment \
    -BACKPHOTO_THICK 64  \
    -BACKPHOTO_TYPE LOCAL \
    -BACK_FILTERSIZE 3  \
    -BACK_SIZE ${backsize}  \
    -CLEAN Y  \
    -CLEAN_PARAM 1.0  \
    -DEBLEND_MINCONT 0.03   \
    -DEBLEND_NTHRESH 32  \
    -DETECT_MINAREA 150  \
    -FILTER Y \
    -FILTER_NAME ${sextract}'/gauss_2.5_5x5.conv' \
    -FLAG_TYPE OR \
    -MAG_ZEROPOINT 16.40006562 \
    -MASK_TYPE CORRECT \
    -PHOT_APERTURES 2.4,8 \
    -PHOT_AUTOPARAMS 2.5,3.5 \
    -PHOT_FLUXFRAC 0.2,0.3,0.5,0.8,0.9  \
    -PHOT_PETROPARAMS 2.0,3.5 \
    -WEIGHT_GAIN Y 
    
    # Making the segmentation region file,
    # These are basically the ellipses found by SExtractor
    awk '(NR>31) {print("image;ellipse("$2","$3","$20*$33","$21*$33","$19")  #text={"$1"}")}' ${cataloge} > ${root}'.sex.reg'

    setenv n `${bin}'/pixval' $segment $xc $yc`
    
 
    set n_line=`echo "$n+31" | bc -l`
    setenv xy `awk -F: -v n="$n_line" '(NR==n) {print($0)}' ${cataloge}`

    setenv RA `echo $xy | awk '{print ($7)}'`
    setenv DEC `echo $xy | awk '{print ($8)}'`

    setenv x0 `echo $xy | awk '{print ($2)}'`
    setenv y0 `echo $xy | awk '{print ($3)}'`
    setenv A `echo $xy | awk '{print ($20)}'`
    setenv B `echo $xy | awk '{print ($21)}'`
    setenv THETA `echo $xy | awk '{print ($19)}'`
    setenv KRON_R `echo $xy | awk '{print ($33)}'`

    setenv mag `echo $xy | awk '{print ($14)}'`

    #re --> pixel
    setenv re  `echo $xy | awk '{print ($25)}'`
    
    # ellipticity
    setenv elips `echo $xy | awk '{print ($32)}'`
    set ratio=`echo "1-$elips" | bc -l`
    # position angle (PA)  [Degrees: Up=0, Left=90]"
    set pa=`echo "${THETA}+90" | bc -l`
    set a=`echo "$A*$KRON_R" | bc -l`
    set b=`echo "$B*$KRON_R" | bc -l`
    
    set a_min=`echo "$A*$KRON_R*$pixel_scael/60" | bc -l`
    set b_min=`echo "$B*$KRON_R*$pixel_scael/60" | bc -l`
    set Re_min=`echo "$re*$pixel_scael/60" | bc -l`
    
    
    printf "||filter || Mag || X0 || Y0 || RA || DEC || Re || Ra  || Rb  || b/a || PA  || \n"  > ${root}_sex.txt
    printf "|| glga || mag || pix || pix || deg || deg || arcmin || arcmin || arcmin ||  ||  deg || \n"  >> ${root}_sex.txt
    printf "|| %s || %.2f || %.1f || %.1f || %.4f || %.4f || %.4f || %.4f || %.4f || %.2f || %.2f || \n" $root $mag $x0 $y0 $RA $DEC $Re_min $a_min $b_min $ratio $pa   >> ${root}_sex.txt

    
    
    # this is used by IDL
    echo $x0 $y0 $a $b $pa    
    
#     rm ${segment}
#     rm ${cataloge}
    
endif



