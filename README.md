[[PageOutline]]


![stv_preview2](https://user-images.githubusercontent.com/13570487/74598526-42f18c00-5030-11ea-84a8-50dcf856ec23.jpg =250x)


## SETUP 

Before launching the program in IDL:

 1. Set the following environment variable to point at the address where '''GLGA''' codes are stored.
{{{
export IDL_GLGAPRO='$HOME/PanStarrs/glgapro/'
}}}

This is what you may have in "glgapro" folder:

{{{
 * astrolib
 * bin          > binary files
 * coyote 
 * cpp          > C programs (compiled binary files in bin folder)
 * dust         > necessary files to calculate extinction
 * extra
 * documents    > dcumentations and manual
 * idlpro       > all glga IDL clodes
 * plot  
 * pypro        > python codes 
 * scripts      > shell scripts
 * sextract     > SExtractor config. files
 * startup.pro
}}}

2. In IDL, define where the data-base is:
{{{
IDL> !GLGA_ROOT = '~/data_base/' 
}}}

3. Under ''~/data_base/data/'', you have your data organised by their coordinates (RA):
{{{
000D  015D  030D  045D  060D  075D  090D  105D  120D  135D  150D  165D  180D  195D  210D  225D  240D  255D  270D  285D  300D  315D  330D  345D
001D  016D  031D  046D  061D  076D  091D  106D  121D  136D  151D  166D  181D  196D  211D  226D  241D  256D  271D  286D  301D  316D  331D  346D
002D  017D  032D  047D  062D  077D  092D  107D  122D  137D  152D  167D  182D  197D  212D  227D  242D  257D  272D  287D  302D  317D  332D  347D
003D  018D  033D  048D  063D  078D  093D  108D  123D  138D  153D  168D  183D  198D  213D  228D  243D  258D  273D  288D  303D  318D  333D  348D
004D  019D  034D  049D  064D  079D  094D  109D  124D  139D  154D  169D  184D  199D  214D  229D  244D  259D  274D  289D  304D  319D  334D  349D
005D  020D  035D  050D  065D  080D  095D  110D  125D  140D  155D  170D  185D  200D  215D  230D  245D  260D  275D  290D  305D  320D  335D  350D
006D  021D  036D  051D  066D  081D  096D  111D  126D  141D  156D  171D  186D  201D  216D  231D  246D  261D  276D  291D  306D  321D  336D  351D
007D  022D  037D  052D  067D  082D  097D  112D  127D  142D  157D  172D  187D  202D  217D  232D  247D  262D  277D  292D  307D  322D  337D  352D
008D  023D  038D  053D  068D  083D  098D  113D  128D  143D  158D  173D  188D  203D  218D  233D  248D  263D  278D  293D  308D  323D  338D  353D
009D  024D  039D  054D  069D  084D  099D  114D  129D  144D  159D  174D  189D  204D  219D  234D  249D  264D  279D  294D  309D  324D  339D  354D
010D  025D  040D  055D  070D  085D  100D  115D  130D  145D  160D  175D  190D  205D  220D  235D  250D  265D  280D  295D  310D  325D  340D  355D
011D  026D  041D  056D  071D  086D  101D  116D  131D  146D  161D  176D  191D  206D  221D  236D  251D  266D  281D  296D  311D  326D  341D  356D
012D  027D  042D  057D  072D  087D  102D  117D  132D  147D  162D  177D  192D  207D  222D  237D  252D  267D  282D  297D  312D  327D  342D  357D
013D  028D  043D  058D  073D  088D  103D  118D  133D  148D  163D  178D  193D  208D  223D  238D  253D  268D  283D  298D  313D  328D  343D  358D
014D  029D  044D  059D  074D  089D  104D  119D  134D  149D  164D  179D  194D  209D  224D  239D  254D  269D  284D  299D  314D  329D  344D  359D
}}}

4. Make sure you have '''awk''', '''ds9''' and '''SExtractor''' already installed to be able to use all the features of the program.

5. In order to use "Python" features, you need to have the following python libraries. Other libraries may have been also used.

{{{
 - numpy
 - pylab
 - matplotlib
 - astropy
 - mpl_toolkits
 - optparse
 - PIL
}}}

6. Most often the provided executable files in ''bin'' folder should work on Linux systems. However you can recompile them in other systems. In '''glgapro''' folder, where ''Makefile'' exists, run the following command to regenerate executable binary files.
{{{
$ make all
}}}

You need to have '''cfitsio''' installed before running making binary files. `cfitsio` is a C libaray to interact with fits files in C/C++ programs.

To remove all binary files in the ''bin'' folder, you can use the following command. Be advised, this removes everything in ''bin'' folder. Please use it wisely.
{{{
$ make clean
}}}



[[Image(stv_masking.jpg)]][[br]]
Green: common mask, Magenta: band mask

= Manual =

== Features ==
 * using masks for individual bands
 * using masks for all bands
 * getting bad pixels from the original image. Bad pixels must have NaN value
 * using ds9 to open all fits images of the same object inside the TV
 * using SExtractor to fit the ellipse to the galaxy
 * using Pylipse, a python program to manually fit an ellipse over the galaxy
 * Pylipse also interacts with ds9. ds9 can be used as a bridge to transfer ellipse parameters between GLGA TV and Pylipse



== TV-Image ==

 * Mouse 
  * Double-click: re-centering the ellipse
  * Left-Click: choosing the vertices of a Region Of Interest (ROI)
  * Right-Click: 
    a. enclosing a region of interest, when choosing a Region Of Interest (ROI)
    b. repeating the highly demanding actions, which are
      1. 'c' adding/editing centeroid
      2. 'd' deleting centeroids
      3. 'r' adding regions (masks)
      4. '!' removing regions (masks)
      5. 'b' adding background regions
 * Keyboard:
  * Note: Short keys may be the same, but having different actions, based on process you are in

 * Line-Styles and Colors:
  * dotted red ellipse: Photometry Aperture 
  * dotted sky blue ellipses: Photometry Annulus (used to estimate the background level/statistics)
  * dotted yellow ellipses: The annulus that is not used for background estimations (when background regions are in use)
  * orange circles: Already Masked point sources
  * green circles: new selected point sources
  * [[span(style=color: #43BFC7	, '''Thick dashed sky blue open regions''' )]]: These are the optional regions to measure the background value. When they are defined, the elliptical annulus turns into yellow. 
  * [[span(style=color: #E2A76F	, '''Dash-dotted peach-color open region''' )]]: If defined '(menu > SExtractor > Use-Mask)', SExtractor only uses this area for measurements. It is useful when half of the image is missing and only the area around the galaxy is useful.
  * [[span(style=color: #4AA02C	, '''Dash-dotted green open region''' )]]: If defined '(menu > SExtractor > Force-Mask)', SExtractor does not filter this area for measurements. It is useful for instance when a masked spike divides a galaxy into two halves. We do not want SExtractor fits a separate ellipse for each half.
  * [[span(style=color: #008000, '''Green filled regions''' )]]: these are the masks used for all bands (i.e. common-masks)  
  * [[span(style=color: #FF00FF	, '''Magenta  filled regions''' )]]: these are the masks used for each individual mask (i.e. band-masks)  

 * Notes:
  * Background estimation regions are always blue
  * Ellipse aperture/annulus are displayed with dotted lines
  * SExtractor related regions are displayed with dash-dotted
  * Background related regions are displayed with dashed lines

== Control Panel ==

=== Filter ===

  1. Filter: You can change the image band. If the JPG file does not exist but the fits file exists, GUI asks the user to generate the JPG file and displays it.
  2. If the JPG file cannot be found anyway, TV does not change to the switched band.

=== Buttons ===

For all the buttons, you have a counterpart on the menu bar, but the opposite is not true. The letter in the parenthesis shows the short key to the same action, once you click on the image. A full set of all actions is also available on the terminal. 

 1. Change Mask: Generally, when drawing a new region to mask the image, it would be applied on all band, and its easier to look at the composite image (however you can still choose to look at other band and edit the common mask file). If a band is selected (and not a composite image), then you will have an opportunity to generate masks for the chosen band. The file that contains each band mask is stored separately. Under this button, you will see the status of the masking. 
{{{
 Current: All Bands: the generated mask would be used for all bands and it is displayed with Green color.
 Current: '_band_' : the mask is generated for the specified band, and it is displayed with Magenta color.
}}}

 2. Add Region: Adding a new mask, based on the status of the label above. To choose a region, left click on the images to specify the region vertices, (or simply drag your mouse while the left mouse button is pressed). Once you are done press right mouse button to enclose the region.

 3. Undo Region: Removing the last added region. Currently this only works for common regions, i.e. green masks.
 
 4. Remove a region: After clicking on this button you need to lick inside a region to remove it. This works for both band-masks and common-mask (i.e. Green and Magenta regions). If the regions does not exist (in case of any mistake), the action would be ignored. 

 5. Clean all regions: After user verification, it removes all the region. This works for both band-masks and common-mask (i.e. Green and Magenta regions).
  
 6. Save region: Normally all common-regions (Green reg.) would be saved after doing the measurement. However, once can save any change before doing the measurement. This is helpful when working with big galaxies. One should note that all band-masks are stored instantly and users should not be worried about saving their pink masks. This is because, we don not expect to have a lot of masking for each individual mask.

 7. Find Star: This helps to find stars and point sources by defining the proper PSF and limiting magnitude. the inputs must be entered in the terminal.

 8. Delete Centeroid: Removing all the previously masked point sources inside a specified region. 

 9: Edit/Add Centeroid: Changing the size of the centeroids, turning them off or adding new ones. All inputs are taken from the terminal


[[Image(stv_pylipse.jpeg, 60%)]]


 10: ds9: Opens the fits file of the displayed image in the TV. If not band has been specified (in case of using the composite image), the use-file ('r' band for /sdss and /panstarrs) would be opened. In ds9, one can play around with the [[span(style=color: #FF0000,  red )]] elliptical aperture, (e.g. re-center it, or resize it). The parameter of this [[span(style=color: #FF0000,  red )]] aperture can be later on captured by GLGA TV and update the ellipse 

 11: Contours On/Off: If ds9 is open, then isophots can be displayed. This helps to locate a proper border around the galaxy. To rotate the [[span(style=color: #FF0000,  red )]] aperture, first select it by clicking inside it on ds9. You see four selection points around it. Press Shift-key on keyboard and choose one of the four point to rotate the ellipse by dragging the mouse. To move the ellipse, click and drag ...

[[Image(ds9_demo.jpeg, 30%)]]

 12: Import Ellipse: To import the [[span(style=color: #FF0000,  red )]] aperture from ds9 to the TV (when ds9 is open). All ellipse information would be updated and be used for the measurement. Note: ds9 can also be used a bridge between Pylipes (i.e. the python code to modify the ellipse) and GLGA TV. However, Saving (in Pylipse) and loading (in TV: Menu > Pyliupse > Load Ellipse) is a way to import from Pylipse to GLGA TV.


 13: Export Ellipse: To update the ellipse in ds9 when it's open. If you have already opened ds9, but you have edited the ellipse in the TV, use this action to export the TV ellipse to ds9. Note: ds9 can also be used a bridge between Pylipes (i.e. the python code to modify the ellipse) and GLGA TV. To export an ellipse from TV to Pylipse, you can re-open Pylipse (in TV: Menu > Pylipse > Open)

 14: Run SExtractor: Running SExtractor to update the ellipse information. If SExtractor has been already run and the information exists in the data-base, the GUI asks whether to use the previous results, or run SExtractor again. Note: one can run SExtractor automatically for the entire data-base, and if there is any adjustment needed, it is possible to re-run it in the TV. In this case, TV first creates a temporary masked image (applying both band-mask and common-mask files) and then runs SExtractor. This way one can take care of any foreground/background problems that could have been already resolved. Also there are other masking options just to run SExtractor (see: menu > SExtractor)

 15: Edit Ellipse: Editing the ellipse parameter in the TV. To see all the options, look at the terminal and follow the right commands. Note: if in the middle of the editing, user clicks on any button on the control-panel, the ellipse editing action would be ignored and the new command would be run.

 16: Original Ellipse: Use the original ellipse paramter and update the TV

 17: Refresh TV: To reload the image and all masks. For larger galaxies, reloading takes time, but it helps to see the color of the masks (in case of having both Magenta and Green masks)

[[Image(stv_note_flag.jpeg, 40%)]]

 18: quit-next: quit the GLGA TV and save the flags and note about the image. You can simply use the pop-up window which is synchronized with terminal. You can also use the short keys suggested in terminal. If you happen to close the pop-up window, don't worry, you can still follow-up with terminal in old-fashion. To prevent users to accidentally mess around with the important TV parameters (e.g. masks, ellipse etc.), all unnecessary actions on the control-panel would remain deactivated unless this process is cancelled. Please don't use the menu-bar when you are in the middle of quitting process. Updated ellipse parameter would be saved. 

 19: Quit-Exit: Exiting the TV. Updated ellipse parameter would not be saved.

 20: Skip: Skipping the analysis of the current galaxy and move on to the next galaxy on the list. Updated ellipse parameter would not be saved.

 21: ''>'' button: The long key to the left of quit/skip, is used to trigger a measurement event. This would shutdown the entire TV to do the optometry.



== Menu Bar ==

=== File ===

 a. Measure >
   * run photometry
 b. Refresh TV
   * reload the JPG image
   * reload ellipse info
   * reload all mask
 c. quit-next
   * saving ellipse info
   * entering flags and note 
   * go to the next galaxy in the list
 d. Skip ...
   * skipping the current galaxy
   * updated ellipse info would not be saved
   * ask user if he/she wants to save the newly added regions
   * jumping to the next galaxy in the list
 e. Quit
   * exiting the program
   * updated ellipse info would not be saved
   * ask user if he/she wants to save the newly added regions


=== Region ===

 a. Undo
   * Remove the last added region  
   * Only works on common regions (all-bands, green)
 b. Add ...
   * Add new masking regions (by left clicking to define vertices, right click to close the region)
   * Works on both common and band regions (green and magenta)
   * To switch between common and band regions, go to ''change mask'' option, under the same menu
 c. Remove one
   * Remove masking regions one by one (click inside the region of interest)
   * Works on both common and band regions (green and magenta)
   * To switch between common and band regions, go to ''change mask'' option, under the same menu
 d. Clean all
   * Remove all masking regions based on the current mask status
   * Works on both common and band regions (green and magenta)
   * To switch between common and band regions, go to ''change mask'' option, under the same menu
 e. Save
   * Saves newly added '''common''' regions
 f. Change Mask
   * Switch between common and band regions
   * It is effective when a particular band is chosen


=== Ellipse ===
 a. Edit ...
   * Edit ellipse information
   * Look at terminal for the effective short key 
   * You need to switch back and force between TV and terminal 
 b. Original
   * Updates the ellipse info. and switched back to the original ellipse
   * Refresh the TV

=== Stars === 
 a. Find ...
   * Find new point sources based on the entered PSF and limiting magnitudes
   * Use terminal to enter the input variables
 b. Delete ...
   * Delete the masked point sources
   * Define a region (left/right mouse clicks)
   * All point sources inside the region would be removed from the point source mask
   * Useful to deselect the star forming regions of spiral galaxy arms
 c. Edit/Add ...
   * Add more point sources manually
   * Edit the existing point sources, e.g. increase/decrease their sizes


=== Background === 
 a. Back. Region
   * Defining regions to measure the background
   * When using these regions, the elliptical annulus would no longer be used
 b. Clear Region
   * Clean the previously defined background regions and go back to the annulus background estimation


=== ds9 ===
 a. Open
   * Opens the fits imaged of the current TV band on ds9
   * If TV is in the composite image mode, the main fits file (use-file) would be used
 b. Contours on/off
   * turning on/off the contours on the ds9 image
 c. Import ellipse
   * Imports the red elliptical aperture from ds9 to TV
 d. Export ellipse
   * Exports the current ellipse from TV to ds9
 e. Close all
   * Closes all open ds9 windows
   * be careful this would close all ds9 windows, not necessarily those opened by GLGA TV

 
=== SExtractor === 
 a. Run
   * Run SEWxtractor on the fits images of the current band
   * If TV is in the composite image mode, the main fits file (use-file) would be used
   * If SExtractor has been previously run, this asks the user to either use the existing results, or re-run SExtractor
   * Sometimes, masking the image helps SExtractor to fit a better ellipse
   * If users chooses to re-run SExtractor, a temporary fits image would be created, which includes all the common-/band- masks and the SExtractor uses that file to extract the ellipse paramters
   * Note: if the image is large, this process may take time. A large temporary file needs to be prepared for SExtractor, and SExtracting also takes time
 b. Segmentation
   * This would display the Segmentation image that is produced by SExtractor. If this is not available, you need to run SExtractor first.
 c. Masked Image
   * This would display the temporary masked region that has been already used in this session
 d. Region on
   * This would turn on all the fitted ellipses by SExtractor. To see them, you need to open ds9
   * Segmentation, Maksed Image, TV Image can be used in ds9 to see these regions
 e. Region off
   * This would turn off all regions previously openned in all ds9 windows
 f. Use-Mask
   * User can define a regions to be specifically used by SExtractor in this session
   * When TV starts, all previous SExtractor regions would be deleted, so this region is only valid in the current session
   * The rest of the image would be totally masked for SExtractor in the process of making the temporary file
   * The un-maksed defined region, then later would be masked by common/band/centroid masks
 g. Force-Mask
   * User can define regions to be specifically used by SExtractor in this session
   * This mask would not be affected by any other common/band/centeroid masks
   * For example, this is useful when a spike cuts a galaxy in to two halves and we do not want SExtractor fits two separates ellipses for one single object
 h. Clear-Mask
   * Clear both Use-Mask and Force-Mask in the current session
   * If you don't clear these masks, they would be removed on the next run, because users are expected to run SExtractor only a few times to find the proper ellipse
   * Once a good ellipse is found, no more messing up with its parameters

=== Pylipse === 
 a. Open
   * Opens Pylipse, a python program
   * There is a complete help on how to use this program
   * This can also be used stand-alone
   * ds9 can be used as a bridge to import/export ellipse to/from this program
   * To export an ellipse from TV to Pylipse, you can re-open Pylipse (Menu > Pylipse > Open)


 b. Load Ellipse
   * Saving (in Pylipse) and loading (Menu > Pyliupse > Load Ellipse) is a way to transfer ellipse from Pylipse to GLGA TV


=== Help === 
 * The URL to this page ... :)




= OUTPUTS =

[[Image(pgc055_PS_panstarrs_images.jpg)]]

 * Images have this format: <ID>_<survey>_images.<pdf/jpg>
 * Bad pixels would be displayed in red in output images. 
 * GLGA does not get any separate bad-pixel map. If there is any badpixel in the image, it should have NaN value ...
 * Pan-STARRS images already include bad pixels (NaN value)
 * The common-mask that is used for all bands are displayed over the composite image, in orange.

[[Image(pgc055_PS_panstarrs_images_masks.jpg)]]

 * The band mask image has this format: <ID>_<survey>_images_mask.<pdf/jpg>
 * If there is any band-mask designed, this would be generated. Otherwise, there would not be such file.
 * Bad pixels are in red
 * Band-masks are displayed in green

[[Image(pgc055_PS_panstarrs_profile.jpg)]]

 * The profile image has this format: <ID>_<survey>_profile.<pdf/jpg>
