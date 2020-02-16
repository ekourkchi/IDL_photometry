#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2016"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"


import sys
import os
import subprocess
import math
import numpy as np
from datetime import *
from pylab import *
import matplotlib
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column 

from mpl_toolkits.axes_grid1 import make_axes_locatable
from optparse import OptionParser
from PIL import Image, ImageTk
from subprocess import Popen, PIPE

  
#################################################
def help():
  return '''
Version : v1.0 (March 2016)
email   : ehsan@ifa.hawaii.edu
Copyright 2016 Ehsan Kourkchi

********** Mouse Actions ***********
Use mouse middle wheel for scrolling
************************************

A) When mouse pointer is on the image

 a1) scroll-down  = zoom-in
 a2) scroll-up    = zoom-outer
 a3) middle-click = re-center the image

B) When anywhere on the GUI window
 
 b1) Ctrl+Scroll_up   = increase the semi-major axis
 b2) Ctrl+Scroll_down = decrease the semi-major axis
 
 b3) Shift+Scroll_up   = increase the semi-minor axis
 b4) Shift+Scroll-down = decrease the semi-minor axis

 b5) Alt+Scroll_up   = increase the position angle (PA)
 b6) Alt+Scroll_down = decrease the position angle (PA)
 
 b7) Ctrl+a          = Increase the ellipse size by 10% (*1.1)
 b8) Ctrl+d          = Decrease the ellipse size by 10% (/1.1)
 
 i  ) PA increases in clock-wise direction. 
 ii ) PA = 0 if semi-major axis is horizontal
 iii) Step-size is displayed in the box next to the a-/b- and PA-control bars.
     - "Left  click" on the box: increases the step size by one pixel
     - "Right click" on the box: decreases the step size by one pixel
   
 b7) "Left/Middle double click" = choose a new center for the ellipse
 b8) "Right click" or "enter-key" = re-center the ellipse
 b9) "Right double click" = choose the new center and re-draw the ellipse at once (b7+b8)
 b10) z+Left_click = choose a new center for the ellipse (b7)
 b11) q/Esc = Ignoring the new chosen center 
'''



#################################################

def arg_parser():
    parser = OptionParser(usage="""\

 - A GUI for manual ellipse fitting for Elliptical and Spiral galaxies ...

 - How to run: 
     %prog [options]

 - Use the -h option to see all options ...


""")
    
    parser.add_option('-p', '--pgc',
                      type='string', action='store',
                      help="""PGC ID""")
    parser.add_option('-a', '--major',
                      type='float', action='store',
                      help="""major-axis""")
    parser.add_option('-b', '--minor',
                      type='float', action='store',
                      help="""minor-axis""")
    parser.add_option('-x', '--xc',
                      type='float', action='store',
                      help="""X-center""")
    parser.add_option('-y', '--yc',
                      type='float', action='store',
                      help="""Y-center""")
    parser.add_option('-g', '--pa',
                      type='float', action='store',
                      help="""position angle""")
    parser.add_option('-s', '--survey',
                      type='string', action='store',
                      help="""survey""")         
    parser.add_option('-j', '--jpg_base',
                      type='string', action='store',
                      help="""jpg base""")     
    parser.add_option('-f', '--fits_base',
                      type='string', action='store',
                      help="""fits base""")  
         
    parser.add_option('-r', '--py_root',
                      type='string', action='store',
                      help="""The root of this code""")  
     
    parser.add_option('-B', '--auxbase',
                      type='string', action='store',
                      help="""The root of this code""")     
         
    (opts, args) = parser.parse_args()
    return opts, args
########

def split_func(sex_file):
  
  output = Popen(["awk", "(NR>2){print($6\" \"$8\" \"$16\" \"$18\" \"$22)}", sex_file], stdout=PIPE)
  s = output.stdout.read()
  spl = s.split()
  xcenter = np.float(spl[0])
  ycenter = np.float(spl[1])  
  ra = np.float(spl[2])*60/.3960
  rb = np.float(spl[3])*60/.3960
  pa = np.float(spl[4])
  
  Status = True
  if ra*rb*pa == 0: Status = False
  
  pa+=90
  if pa > 180: pa-=180
  return Status, xcenter, ycenter, ra, rb, pa
#################################


def ellipse_param(line):
  
  separator = []
  for i in range(len(line)):
    if line[i] == '(' or line[i] == ')':
      separator.append(i)
    elif line[i] == ',':
      separator.append(i)
  
  if len(separator) != 6: 
    return None
  
  xcenter = line[separator[0]+1:separator[1]]
  ycenter = line[separator[1]+1:separator[2]]
  a       = line[separator[2]+1:separator[3]]
  b       = line[separator[3]+1:separator[4]]
  pa      = line[separator[4]+1:separator[5]]
  
  xcenter = np.float(xcenter)
  ycenter = np.float(ycenter)
  a = np.float(a)
  b = np.float(b)
  pa = np.float(pa)
  
  state = [xcenter, ycenter, a, b, pa]
  
  return state 
      
    
  


def load_ellipse_ds9(reg_file):
  
   try: 
      fo = open(reg_file, "rw+")
   except: 
      print "\""+reg_file+"\": No such file ..."
      return None
   
   for i in range(20):
     
     line = fo.readline()
     if line[0:7] == 'ellipse':
       return ellipse_param(line)
     
   return None  

#################################################
   #Circlx, Circly = myCircle(936, 1861, r)
   #circle, = ax.plot(Circlx, Circly, color='red', lw=1)
   #circle.set_dashes([2,3])
def myCircle(xcenter, ycenter, r):
   
   theta = np.arange(0.0, 360.0, 1.0)*np.pi/180.0
   Circlx = r*np.cos(theta) + xcenter
   Circly = r*np.sin(theta) + ycenter
   
   return Circlx, Circly
   
#################################################
def myEllipse(xcenter, ycenter, a, b, angle):


   theta = np.arange(0.0, 360.0, 1.0)*np.pi/180.0
   
   x = a * np.cos(theta)
   y = b * np.sin(theta)

   rtheta = np.radians(angle)
   R = np.array([
    [np.cos(rtheta), -np.sin(rtheta)],
    [np.sin(rtheta),  np.cos(rtheta)],
    ])
   
   x, y = np.dot(R, np.array([x, y]))
   x += xcenter
   y += ycenter

   return x, y

#################################################
class ImDisp:
  
  def __init__(self, Xmin, Xmax, Ymin, Ymax):
    self.Xmin = Xmin
    self.Xmax = Xmax
    self.Ymin = Ymin
    self.Ymax = Ymax
    
    self.x1 = Xmin
    self.x2 = Xmax
    self.y1 = Ymin
    self.y2 = Ymax
  
  def zoom(self, xc, yc, ratio=1):
    
    delta_x = self.x2 - self.x1
    delta_y = self.y2 - self.y1
    
    dx = 0.5 * ratio * delta_x
    dy = 0.5 * ratio * delta_y
    
    ###
    if xc + dx > self.Xmax:
      self.x2 = self.Xmax
      self.x1 = self.Xmax - 2. * dx
    elif xc - dx < self.Xmin:
      self.x1 = self.Xmin
      self.x2 = self.Xmin + 2. * dx
    else: 
      self.x1 = xc - dx
      self.x2 = xc + dx
    
    if self.x1 < self.Xmin:
      self.x1 = self.Xmin
    if self.x2 > self.Xmax:
      self.x2 = self.Xmax
    ###
    if yc + dy > self.Ymax:
      self.y2 = self.Ymax
      self.y1 = self.Ymax - 2. * dy
    elif yc - dy < self.Ymin:
      self.y1 = self.Ymin
      self.y2 = self.Ymin + 2. * dy
    else: 
      self.y1 = yc - dy
      self.y2 = yc + dy  
      
    if self.y1 < self.Ymin:
      self.y1 = self.Ymin
    if self.y2 > self.Ymax:
      self.y2 = self.Ymax
    ###
    
  def zoom_IN(self, xc, yc, ratio=0.75):
    self.zoom(xc, yc, ratio=ratio)
    return self.x1, self.x2, self.y1, self.y2
      
  def zoom_OUT(self, xc, yc, ratio=4./3):
    self.zoom(xc, yc, ratio=ratio)
    return self.x1, self.x2, self.y1, self.y2  

  def reset(self):
    self.x1 = self.Xmin
    self.x2 = self.Xmax
    self.y1 = self.Ymin
    self.y2 = self.Ymax
    return self.x1, self.x2, self.y1, self.y2
  
  def pan(self, xc, yc):
    self.zoom(xc, yc, ratio=1)
    return self.x1, self.x2, self.y1, self.y2 
    
    
    


#################################################
class Undo_Redo:
  
  
  def __init__(self):
    
    self.data = []
    self.n = 0
    self.iter = 0
    self.size = 100
    
  
  def push(self, item):
    self.data = self.data[0:self.iter]
    self.iter += 1
    self.data.append(item)
    self.n = self.iter
    
    
    if self.n > self.size:
      self.data = self.data[self.n-self.size:self.n]
      self.n = self.size
      self.iter = self.size
    
  
  def undo(self):
    if self.iter > 1:
       self.iter -= 1
       return self.data[self.iter-1]
       
    return None
  
  
  def redo(self):
    if self.iter < self.n:
       self.iter += 1
       return self.data[self.iter-1]
    return None


  def current(self):
    return self.data[self.iter-1]
#################################################



class objEllipse():
  
  def __init__(self, ax, xcenter, ycenter, a, b, pa, color='red'):
    
    self.xcenter = xcenter
    self.ycenter = ycenter
    self.a = a
    self.b = b
    self.pa = pa
    
    ellX, ellY = myEllipse(xcenter, ycenter, a, b, pa)
    self.ellipse, = ax.plot(ellX, ellY, color=color, lw=1)
    self.ellipse_c, = ax.plot(xcenter, ycenter, 'r+')
    self.ellipse.set_dashes([2,3])
    
    
    annotate("Center ", (0.016,0.42), xycoords='figure fraction', size=14, color='maroon')
    
    self.uXc = annotate("Xc: "+'{:.1f}'.format(xcenter), (0.015,0.39), xycoords='figure fraction', size=10)
    self.uYc = annotate("Yc: "+'{:.1f}'.format(ycenter), (0.015,0.36), xycoords='figure fraction', size=10)
    
    state = [xcenter, ycenter, a, b, pa]
    self.state_iter = Undo_Redo()
    self.state_iter.push(state)
    
    b_a = float(b)/float(a)
    self.b_a = annotate("b/a: "+'{:.2f}'.format(b_a), (0.20,0.20), xycoords='figure fraction', size=12)
    
    
    
  def update(self, xcenter, ycenter, a, b, pa):
    
    self.xcenter = np.float(xcenter)
    self.ycenter = np.float(ycenter)
    self.a = np.float(a)
    self.b = np.float(b)
    self.pa = np.float(pa)
    
    ellX, ellY = myEllipse(xcenter, ycenter, a, b, pa)
    self.ellipse.set_xdata(ellX)
    self.ellipse.set_ydata(ellY)
    self.ellipse_c.set_xdata([xcenter])
    self.ellipse_c.set_ydata([ycenter])
    
    self.uXc.set_text("Xc: "+'{:.1f}'.format(xcenter))
    self.uYc.set_text("Yc: "+'{:.1f}'.format(ycenter))
    
    state = [xcenter, ycenter, a, b, pa]
    self.state_iter.push(state)
    b_a = float(b)/float(a)
    self.b_a.set_text("b/a: "+'{:.2f}'.format(b_a))
    
    draw()
  
  def undo_redo(self, undo=True):
    
    if undo == True:
      state = self.state_iter.undo()
    else:
      state = self.state_iter.redo()
    
    if state != None:
      [xcenter, ycenter, a, b, pa] = state

      self.xcenter = xcenter
      self.ycenter = ycenter
      self.a = a
      self.b = b
      self.pa = pa
    
      ellX, ellY = myEllipse(xcenter, ycenter, a, b, pa)
      self.ellipse.set_xdata(ellX)
      self.ellipse.set_ydata(ellY)
      self.ellipse_c.set_xdata([xcenter])
      self.ellipse_c.set_ydata([ycenter])
      self.uXc.set_text("Xc: "+'{:.1f}'.format(xcenter))
      self.uYc.set_text("Yc: "+'{:.1f}'.format(ycenter))
      b_a = float(b)/float(a)
      self.b_a.set_text("b/a: "+'{:.2f}'.format(b_a))
      draw()  
    return state

  def update_color(self, color):
    self.ellipse.set_color(color)
    self.ellipse_c.set_color(color)
  
    
    
############################################
############################################
############################################  
    
has_counter = False

a = 443
b = 208
pa = 23
xcenter = 1361
ycenter = 1367

delta_axis = 10
delta_angle = 5
delta_position = 10

x_new = 0
y_new = 0
center_change = False

plus_minus = r'$\pm$'


py_root = ''
jpg_base  = ''
fits_base  = ''
auxbase = ''
dtype = ''
ID = ''



############################################
############################################
############################################

def main():
   global a, b, pa, xcenter, ycenter
   
   #print a
   #print b
   #print xcenter
   #print ycenter
   #print pa
  
   #print dtype      # sdss
   #print ID         # pgc55
   #print jpg_base   # +  '_u.jpg'
   #print fits_base  # + '_u.fits'  
   
   u_file = fits_base+'_u.fits'
   g_file = fits_base+'_g.fits'
   r_file = fits_base+'_r.fits'
   i_file = fits_base+'_i.fits'
   z_file = fits_base+'_z.fits'
   HaveFitsFile = True
   
   gri_image = None
   urz_image = None
   
   has_u_fits = os.path.isfile(u_file)
   has_g_fits = os.path.isfile(g_file)
   has_r_fits = os.path.isfile(r_file)
   has_i_fits = os.path.isfile(i_file)
   has_z_fits = os.path.isfile(z_file)
   
   if not has_u_fits:
       has_u_fits = os.path.isfile(u_file+'.gz')
       u_file+='.gz'
   if not has_g_fits:
       has_g_fits = os.path.isfile(g_file+'.gz')  
       g_file+='.gz'
   if not has_r_fits:
       has_r_fits = os.path.isfile(r_file+'.gz')
       r_file+='.gz'
   if not has_i_fits:
       has_i_fits = os.path.isfile(i_file+'.gz')
       i_file+='.gz'
   if not has_z_fits:
       has_z_fits = os.path.isfile(z_file+'.gz')  
       z_file+='.gz' 
   
   jpeg_gri = jpg_base + '_gri.jpg'
   jpeg_urz = jpg_base +'_urz.jpg'
   
   fig = plt.figure(figsize=(10, 8), dpi=100)
   ax = fig.add_axes([0.2, 0.25, 0.7,  0.70]) 
   subplots_adjust(left=0.25, bottom=0.25)
   fig.patch.set_facecolor('lightgray')
   #print fig.patch.get_facecolor()
   
   if has_g_fits == False:
     print g_file + '  does not exist ... !!! \n'
     if has_r_fits:
        fits_file = r_file
     elif has_i_fits:
        fits_file = i_file
     elif has_z_fits:
        fits_file = z_file
     elif has_u_fits:
        fits_file = u_file       
     else:
        HaveFitsFile = False
   else:
     fits_file = g_file
   
   
   if HaveFitsFile:
      hdu_list = fits.open(fits_file)
      image_data = hdu_list[0].data
      w = wcs.WCS(hdu_list[0].header) 
      fits_x_max, fits_y_max = image_data.shape
      

      
      
   
   if os.path.isfile(jpeg_gri):
      img = Image.open(jpeg_gri)
      rsize = img.resize((img.size[0],img.size[1])) # Use PIL to resize 
      rsizeArr = np.asarray(rsize)  # Get array back
      #rsizeArr = np.fliplr(rsizeArr) # left/right flip
      gri_image = np.flipud(rsizeArr) # up/down flip
      gri_x_max, gri_y_max, dimension = gri_image.shape
   else: 
      gri_image = None
      
      
   if os.path.isfile(jpeg_urz):
      img = Image.open(jpeg_urz)
      rsize = img.resize((img.size[0],img.size[1])) # Use PIL to resize 
      rsizeArr = np.asarray(rsize)  # Get array back
      urz_image = np.flipud(rsizeArr) # up/down flip
      urz_x_max, urz_y_max, dimension = urz_image.shape
   else: 
      urz_image = None
         
   
   if gri_image is not None:
      imgplot = plt.imshow(gri_image)
      x_max = gri_x_max
      y_max = gri_y_max
      x_min = 0
      y_min = 0
      axis([x_min, x_max, y_min, y_max])
   elif urz_image is not None:
      imgplot = plt.imshow(urz_image)
      x_max = urz_x_max
      y_max = urz_y_max
      x_min = 0
      y_min = 0
      axis([x_min, x_max, y_min, y_max])
   elif HaveFitsFile: 
     scale_min = np.min(image_data)
     scale_max = np.max(image_data)  
     x_max = fits_x_max
     y_max = fits_y_max
     x_min = 0
     y_min = 0     
     imgplot = plt.imshow(image_data, cmap='gray', vmin=scale_min, vmax=scale_max, norm=LogNorm())
     axis([x_min, x_max, y_min, y_max])
   else:
     print "Nither JPG file nor FITS file ... !!!!"
     sys.exit(0)
   
   
   
   
   
   disp = ImDisp(x_min, x_max, y_min, y_max)

   root_file = fits_base
   save_name = auxbase+'_stv_ellipse_output.dat'
   
   
   is_done = False   
   
   

   pa += 90
   if pa > 180: pa-=180





   semi_a_max = np.round(0.75 * x_max)
   semi_b_max = np.round(0.75 * y_max)


   resetax = axes([0.01, 0.34, 0.18, 0.43])
   info_bts = Button(resetax, '', color='lightgray', hovercolor='lightgray') 
   
   resetax = axes([0.86, 0.65, 0.13, 0.31])
   ds9_bts = Button(resetax, '', color='lightgray', hovercolor='lightgray')  
   
   annotate("object: " + ID, (0.45,0.02), xycoords='figure fraction', size=12, color='maroon')
   
   

        
   user_center, = ax.plot([], [],  'ro')
        

   ellipse = objEllipse(ax, np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))

   axcolor = 'lightgoldenrodyellow'

   


    
        
        
   #######################
   scale = 1
   
   ax_a   = axes([0.22, 0.15, 0.60, 0.02], axisbg=axcolor)
   ax_b   = axes([0.22, 0.1, 0.60, 0.02], axisbg=axcolor)
   ax_pa  = axes([0.22, 0.05, 0.60, 0.02], axisbg='tan')
   

   slider_a = Slider(ax_a, 'a', 0, semi_a_max, valinit=a, dragging=True)
   slider_b = Slider(ax_b, 'b', 0, semi_b_max, valinit=b, dragging=True)
   slider_pa = Slider(ax_pa, 'PA', 0, 180, valinit=pa, dragging=True)
   
   def update(val):
      global a, b, pa, xcenter, ycenter
      a = slider_a.val
      b = slider_b.val
      pa = slider_pa.val
      ellipse.update(np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))

      
   slider_a.on_changed(update)
   slider_b.on_changed(update)
   slider_pa.on_changed(update)
   #######################   
   #######################   
   rax = axes([0.020, 0.045, 0.15, 0.15], axisbg='tan')
   radio = RadioButtons(rax, ('red', 'cyan', 'lawngreen', 'yellow'), active=0)
   annotate('ellipse color', (0.05,0.025), xycoords='figure fraction', size=8, color='black')
 
   
   def colorfunc(label):
      ellipse.update_color(label)
      draw()
      
   radio.on_clicked(colorfunc)
   #######################   
   #######################   
   #######################   
   
   if gri_image is not None and urz_image is not None: 
      rax_image = axes([0.02, 0.22, 0.15, 0.08], axisbg='tan')
      radio_image = RadioButtons(rax_image, ('gri', 'urz'), active=0)
   
   def change_image(color):
      if color == 'gri':
	imgplot = ax.imshow(gri_image)
	draw()
      elif color == 'urz':
        imgplot = ax.imshow(urz_image)
        draw()
     
   if gri_image is not None and urz_image is not None: 
      radio_image.on_clicked(change_image)
   ####################### 
   
   resetax = axes([0.93, 0.125, 0.05, 0.05])
   delt_axis_button = Button(resetax, plus_minus+str(delta_axis), color=axcolor, hovercolor='navajowhite')  
   
   def addOne_axis(event):
      global delta_axis
      if event.button == 1:
	 delta_axis +=1
	 delt_axis_button.label.set_text(plus_minus+str(delta_axis))
	 draw()
      elif event.button == 3:
         delta_axis -=1
         delt_axis_button.label.set_text(plus_minus+str(delta_axis))            
         draw()
        
   delt_axis_button.on_clicked(addOne_axis)
   #######################  
   tmp = Image.open(py_root+'/icons/Help_icon.png')
   help_label = annotate("HELP", (0.035,0.87), xycoords='figure fraction', size=12)
   tmp_rsize = tmp.resize((tmp.size[0],tmp.size[1]))
   help_icon = np.asarray(tmp_rsize)  
   help_icon = np.flipud(help_icon)
   
   resetax = axes([0.02, 0.9, 0.07, 0.07])
   help_button = Button(resetax, '?', color='white', hovercolor='yellow', image=help_icon) 
   
   def help_me(event):
     help_fi = auxbase+'_pylipse_help.txt'
     
     if not os.path.isfile(help_fi):
      fout=open(help_fi,'w')
      fout.write(help()+'\n')
      fout.close()
     
     if os.path.isfile(help_fi):
       subprocess.Popen('python ' + py_root + 'help_widget.py  -x -f '+help_fi+' &', shell=True)
   
   help_button.on_clicked(help_me) 
   #######################   
   #######################  
   tmp = Image.open(py_root+'/icons/home_icon.png')
   home_label = annotate("Home", (0.115,0.87), xycoords='figure fraction', size=12)
   tmp_rsize = tmp.resize((tmp.size[0],tmp.size[1]))
   home_icon = np.asarray(tmp_rsize)  
   home_icon = np.flipud(home_icon)
   
   resetax = axes([0.1, 0.9, 0.07, 0.07])
   home_button = Button(resetax, ' ', color='white', hovercolor='yellow', image=home_icon) 
   
   def home_disp(event):
     i1, i2, j1, j2 = disp.reset()
     ax.set_xlim(i1,i2)
     ax.set_ylim(j1,j2)
     draw()
   
   home_button.on_clicked(home_disp) 
   #######################   
   resetax = axes([0.86, 0.25, 0.055, 0.055])
   
   undo_button = Button(resetax, 'Undo', color='tan', hovercolor='navajowhite')  
   
   def undo(event):
      global a, b, pa, xcenter, ycenter
      state = ellipse.undo_redo(undo=True)
      if state != None: 
	 [xcenter, ycenter, a, b, pa] = state
         
         slider_a.disconnect(0)
         slider_a.set_val(a)
         slider_a.observers[0] = update
         
         slider_b.disconnect(0)
         slider_b.set_val(b)
         slider_b.observers[0] = update
         
         slider_pa.disconnect(0)
         slider_pa.set_val(pa)
         slider_pa.observers[0] = update
        
   undo_button.on_clicked(undo) 
   #######################    
   #######################   

   resetax = axes([0.93, 0.25, 0.055, 0.055])
   redo_button = Button(resetax, 'Redo', color='tan', hovercolor='navajowhite')  
   
   def redo(event):
      global a, b, pa, xcenter, ycenter
      state = ellipse.undo_redo(undo=False)
      if state != None: 
	 [xcenter, ycenter, a, b, pa] = state
         
         slider_a.disconnect(0)
         slider_a.set_val(a)
         slider_a.observers[0] = update
         
         slider_b.disconnect(0)
         slider_b.set_val(b)
         slider_b.observers[0] = update
         
         slider_pa.disconnect(0)
         slider_pa.set_val(pa)
         slider_pa.observers[0] = update
        
   redo_button.on_clicked(redo) 
   #######################   
   tmp = Image.open(py_root+'/icons/Save-icon.png')
   tmp_rsize = tmp.resize((tmp.size[0],tmp.size[1]))
   save_icon = np.asarray(tmp_rsize)  
   save_icon = np.flipud(save_icon)
   resetax = axes([0.92, 0.35, 0.07, 0.07])
   annotate("Save", (0.865,0.39), xycoords='figure fraction', size=10)
   annotate("Ellipse", (0.86,0.365), xycoords='figure fraction', size=10)
  
   
   save_button = Button(resetax, '', color='white', hovercolor='navajowhite', image=save_icon)  
   
   def save(event):
      global a, b, pa, xcenter, ycenter
      
      save_name = auxbase+'_stv_ellipse_output.dat'
      
      world = w.wcs_pix2world([[xcenter, ycenter]], 1)
      RA = world[0][0]
      DEC = world[0][1]
      
      #ellipse_string = '{:.5f}'.format(RA)+ ', ' + '{:.5f}'.format(DEC)+ ', ' + '{:.2f}'.format(xcenter)+', '+'{:.2f}'.format(ycenter)+', '+'{:.2f}'.format(a)+', '+'{:.2f}'.format(b)+', '+'{:.2f}'.format(pa)+', '+'{:.2f}'.format(a*0.396/60)+', '+'{:.2f}'.format(b*0.396/60)+', '+'{:.2f}'.format(pa-90)
      ellipse_string = '{:.2f}'.format(2*a*0.396)+ '  ' + '{:.2f}'.format(2*b*0.396) + '  ' + '{:.5f}'.format(RA)+ '  ' + '{:.5f}'.format(DEC)+ '  ' + '{:.2f}'.format(pa-90)+'  ' + '{:.2f}'.format(a)+'  '+'{:.2f}'.format(b) + '  ' + '{:.2f}'.format(xcenter)+'  '+'{:.2f}'.format(ycenter)+'  1  0.396'
      
      ellipse_string += datetime.datetime.now().strftime('  %b-%d-%Y  %H:%M:%S')

      #header = "ellipse, RA_deg, DEC_deg, xcenter, ycenter, a_pix, b_pix, pa_deg_x, a_arcmin, b_arcmin, pa_deg_y, date, time \n"
      header = "# STV ellipse output (majordiam_as, minordiam_as, ra_deg, dec_deg, PA_deg, majordiam_px, minordiam_px, x_px, y_px, astrom_bool, as_pix, date, time ) \n"

      
      if os.path.isfile(save_name):
	with open(save_name, "r+") as f:
	  old = f.read() # read everything in the file
	  l_old = len(old)
	  i = 0 
	  while i < l_old:
	    if old[i] == '\n': break
	    i += 1
	  
	  if i+1 < l_old:
	    old = old[i+1:l_old]
	    
	  new_str = ''
	  n_old = len(old)
	  if n_old> 1 and old[0]!='#':
	    new_str += '# '
	    
	  for i in range(n_old):
	    new_str += old[i]
	    if old[i] == '\n' and i<n_old-1 and old[i+1] != '#':
	      new_str += "# "
	  
	  new_str = new_str[0:len(new_str)-2]
	  f.seek(0) # rewind
	  f.write(header + new_str + '\n' + ellipse_string + '\n') # write the new line after
      
      else:
        with open(save_name, "w+") as f:
          f.write(header + ellipse_string+"\n")
      #sys.exit(0)
        
   save_button.on_clicked(save)    
   
   #######################   

   
   
   
   ds9_label = annotate("DS9", (0.94,0.92), xycoords='figure fraction', size=12)
   
   
   resetax = axes([0.87, 0.90, 0.05, 0.05])
   if has_u_fits:
      ds9_u_button = Button(resetax, 'u', color='lightblue', hovercolor='navajowhite')  
   else: 
      ds9_u_button = Button(resetax, 'u', color='grey', hovercolor='grey')  
   
   def ds9_u_open(event):
      global a, b, pa, xcenter, ycenter
      
      ds9_command = "ds9 " +  u_file
      #ds9_command += " -height 500"ds9
      ds9_command += " -width 500"
      ds9_command += " -scale log"
      ds9_command += " -scale minmax"
      ds9_command += " -zoom to fit"
      ds9_command += " -pan to "+str(xcenter)+" "+str(ycenter)
      #ds9_command += " -regions command \"ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2\""
      os.system(ds9_command+" &")
      
      os.system("sleep 3")
      ##os.system("xpaset -p ds9 pan to "+str(xcenter)+" "+str(ycenter))
      #os.system("xpaset -p ds9 scale minmax")
      #os.system("xpaset -p ds9 scale log")
      #os.system("xpaset -p ds9 height 500")
      #os.system("xpaset -p ds9 width 500")
      os.system("xpaset -p ds9 regions command \"{ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2}\"")
      #os.system("xpaset -p ds9 zoom to fit")
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")

      
   if has_u_fits:    
      ds9_u_button.on_clicked(ds9_u_open) 
   ####################### 
   #######################   
   #######################   

   resetax = axes([0.87, 0.83, 0.05, 0.05])
   if has_g_fits:
      ds9_g_button = Button(resetax, 'g', color='lightgreen', hovercolor='navajowhite')  
   else: 
      ds9_g_button = Button(resetax, 'g', color='grey', hovercolor='grey')  
   
   def ds9_g_open(event):
      global a, b, pa, xcenter, ycenter
      
      ds9_command = "ds9 " +  g_file
      ds9_command += " -height 500"
      ds9_command += " -width 500"
      ds9_command += " -scale log"
      ds9_command += " -scale minmax"
      ds9_command += " -zoom to fit"
      ds9_command += " -pan to "+str(xcenter)+" "+str(ycenter)
      #ds9_command += " -regions command \"ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2\""
      os.system(ds9_command+" &")
      
      os.system("sleep 3")
      ##os.system("xpaset -p ds9 pan to "+str(xcenter)+" "+str(ycenter))
      #os.system("xpaset -p ds9 scale minmax")
      #os.system("xpaset -p ds9 scale log")
      #os.system("xpaset -p ds9 height 500")
      #os.system("xpaset -p ds9 width 500")
      os.system("xpaset -p ds9 regions command \"{ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2}\"")
      #os.system("xpaset -p ds9 zoom to fit")
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")

      
   if has_g_fits:        
     ds9_g_button.on_clicked(ds9_g_open) 
   ####################### 
   #######################   

   resetax = axes([0.93, 0.83, 0.05, 0.05])
   if has_r_fits:
      ds9_r_button = Button(resetax, 'r', color='red', hovercolor='navajowhite')  
   else: 
      ds9_r_button = Button(resetax, 'r', color='grey', hovercolor='grey')  
   
   def ds9_r_open(event):
      global a, b, pa, xcenter, ycenter
      
      ds9_command = "ds9 " +  r_file
      ds9_command += " -height 500"
      ds9_command += " -width 500"
      ds9_command += " -scale log"
      ds9_command += " -scale minmax"
      ds9_command += " -zoom to fit"
      ds9_command += " -pan to "+str(xcenter)+" "+str(ycenter)
      #ds9_command += " -regions command \"ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2\""
      os.system(ds9_command+" &")
      
      os.system("sleep 3")
      ##os.system("xpaset -p ds9 pan to "+str(xcenter)+" "+str(ycenter))
      #os.system("xpaset -p ds9 scale minmax")
      #os.system("xpaset -p ds9 scale log")
      #os.system("xpaset -p ds9 height 500")
      #os.system("xpaset -p ds9 width 500")
      os.system("xpaset -p ds9 regions command \"{ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2}\"")
      #os.system("xpaset -p ds9 zoom to fit")
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")

      
   if has_r_fits:     
      ds9_r_button.on_clicked(ds9_r_open) 
   #######################    
   ####################### 
   #######################   

   resetax = axes([0.87, 0.76, 0.05, 0.05])
   if has_i_fits:
      ds9_i_button = Button(resetax, 'i', color='orange', hovercolor='navajowhite')  
   else: 
      ds9_i_button = Button(resetax, 'i', color='grey', hovercolor='grey')  
   
   def ds9_i_open(event):
      global a, b, pa, xcenter, ycenter
      
      ds9_command = "ds9 " +  i_file
      ds9_command += " -height 500"
      ds9_command += " -width 500"
      ds9_command += " -scale log"
      ds9_command += " -scale minmax"
      ds9_command += " -zoom to fit"
      ds9_command += " -pan to "+str(xcenter)+" "+str(ycenter)
      #ds9_command += " -regions command \"ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2\""
      os.system(ds9_command+" &")
      
      os.system("sleep 3")
      ##os.system("xpaset -p ds9 pan to "+str(xcenter)+" "+str(ycenter))
      #os.system("xpaset -p ds9 scale minmax")
      #os.system("xpaset -p ds9 scale log")
      #os.system("xpaset -p ds9 height 500")
      #os.system("xpaset -p ds9 width 500")
      os.system("xpaset -p ds9 regions command \"{ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2}\"")
      #os.system("xpaset -p ds9 zoom to fit")
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")

      
   if has_i_fits:
      ds9_i_button.on_clicked(ds9_i_open) 
   #######################   
   #######################   

   resetax = axes([0.93, 0.76, 0.05, 0.05])
   if has_z_fits:
      ds9_z_button = Button(resetax, 'z', color='gold', hovercolor='navajowhite')  
   else: 
      ds9_z_button = Button(resetax, 'z', color='grey', hovercolor='grey')  
   
   def ds9_z_open(event):
      global a, b, pa, xcenter, ycenter
      
      ds9_command = "ds9 " +  z_file
      ds9_command += " -height 500"
      ds9_command += " -width 500"
      ds9_command += " -scale log"
      ds9_command += " -scale minmax"
      ds9_command += " -zoom to fit"
      ds9_command += " -pan to "+str(xcenter)+" "+str(ycenter)
      #ds9_command += " -regions command \"ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2\""
      os.system(ds9_command+" &")
      
      os.system("sleep 3")
      ##os.system("xpaset -p ds9 pan to "+str(xcenter)+" "+str(ycenter))
      #os.system("xpaset -p ds9 scale minmax")
      #os.system("xpaset -p ds9 scale log")
      #os.system("xpaset -p ds9 height 500")
      #os.system("xpaset -p ds9 width 500")
      os.system("xpaset -p ds9 regions command \"{ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2}\"")
      #os.system("xpaset -p ds9 zoom to fit")
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")

      
   if has_z_fits:     
      ds9_z_button.on_clicked(ds9_z_open) 
   #######################
   
   resetax = axes([0.92, 0.55, 0.07, 0.07])
   annotate("Export", (0.86,0.59), xycoords='figure fraction', size=10)
   annotate("to ds9", (0.865,0.565), xycoords='figure fraction', size=8)

     
   tmp = Image.open(py_root+'/icons/right_arrow.png')
   tmp_rsize = tmp.resize((tmp.size[0],tmp.size[1]))
   right_arrow = np.asarray(tmp_rsize)   
   ds9_sync_button = Button(resetax, '', color='cornsilk', hovercolor='navajowhite', image=right_arrow)
   ds9_sync_button.label.set_fontsize(10)
   
   def ds9_sync(event):
      global a, b, pa, xcenter, ycenter
    
      
      os.system("xpaset -p ds9 regions delete all")
      os.system("xpaset -p ds9 regions command \"{ellipse "+str(xcenter)+" "+str(ycenter)+" "+str(a)+" "+str(b)+" "+str(pa)+" # color=red width=2}\"")
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")
      
        
   ds9_sync_button.on_clicked(ds9_sync) 
   #######################   
   resetax = axes([0.92, 0.45, 0.07, 0.07])
   annotate("Import", (0.86,0.49), xycoords='figure fraction', size=10)
   annotate("from ds9", (0.86,0.465), xycoords='figure fraction', size=8)
   
   
   tmp = Image.open(py_root+'/icons/left_arrow.png')
   tmp_rsize = tmp.resize((tmp.size[0],tmp.size[1]))
   left_arrow = np.asarray(tmp_rsize)
   esn_sync_button = Button(resetax, '', color='cornsilk', hovercolor='navajowhite', image=left_arrow)  
   esn_sync_button.label.set_fontsize(10)
   
   def ds9_sync(event):
      global a, b, pa, xcenter, ycenter
    
      os.system("xpaset -p ds9 regions system  image")
      os.system("xpaset -p ds9 regions format ds9")
      os.system("xpaset -p ds9 regions save python_foo.reg")
      state = load_ellipse_ds9("python_foo.reg")
      if state != None:
	#print state
	[xcenter, ycenter, a, b, pa] = state
        
        
        ellipse.update(np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))
        
        slider_a.disconnect(0)
        slider_a.set_val(a)
        slider_a.observers[0] = update
         
        slider_b.disconnect(0)
        slider_b.set_val(b)
        slider_b.observers[0] = update
         
        slider_pa.disconnect(0)
        slider_pa.set_val(pa)
        slider_pa.observers[0] = update
        
   esn_sync_button.on_clicked(ds9_sync) 
   #######################
   resetax = axes([0.89, 0.67, 0.07, 0.07])
   ds9_contour_button = Button(resetax, 'Contours\nOn', color='lightgrey', hovercolor='navajowhite')  
   ds9_contour_button.label.set_fontsize(10)
   
   def ds9_counter(event):
      global has_counter
      
      if has_counter == False:
         os.system("xpaset -p ds9 contour nlevels 30")
         os.system("xpaset -p ds9 contour smooth 10")
         os.system("xpaset -p ds9 contour scale log")
         os.system("xpaset -p ds9 contour generate")
         os.system("xpaset -p ds9 contour color green")
         os.system("xpaset -p ds9 contour")
         
         ds9_contour_button.label.set_text('Contours\nOff')
         ds9_contour_button.color ='cornsilk'
         has_counter =  True
      else:
	 os.system("xpaset -p ds9 contour clear")
	 os.system("xpaset -p ds9 contour close")
	 ds9_contour_button.label.set_text('Contours\nOn')
	 ds9_contour_button.color = 'lightgrey'
	 has_counter =  False



        
   ds9_contour_button.on_clicked(ds9_counter) 
   #######################   
   
   #######################   

   resetax = axes([0.93, 0.045, 0.05, 0.05])
   delt_angle_button = Button(resetax, plus_minus+str(delta_angle), color='tan', hovercolor='navajowhite')  
   
   def addOne_angle(event):
      global delta_angle
      if event.button == 1:
	 delta_angle +=1
	 delt_angle_button.label.set_text(plus_minus+str(delta_angle))
	 draw()
      elif event.button == 3:
         delta_angle -=1
         delt_angle_button.label.set_text(plus_minus+str(delta_angle))            
         draw()
        
   delt_angle_button.on_clicked(addOne_angle) 
   ####################### 
   #######################   

   resetax = axes([0.12, 0.355, 0.05, 0.05])
   delt_position_button = Button(resetax, plus_minus+str(delta_position), color='gainsboro', hovercolor='navajowhite')  
   
   def addOne_position(event):
      global delta_position
      if event.button == 1:
	 delta_position +=1
	 delt_position_button.label.set_text(plus_minus+str(delta_position))
	 draw()
      elif event.button == 3:
         delta_position -=1
         delt_position_button.label.set_text(plus_minus+str(delta_position))            
         draw()
        
   delt_position_button.on_clicked(addOne_position) 
   #######################    
   
   ####################### 
    
   def scroll_event(event):
         #print 'you pressed', event.key, event.button, event.xdata, event.ydata, event.key
            
         global a, b, pa, xcenter, ycenter, delta_axis, delta_angle
         
         
         if event.inaxes == ax: 
	   if event.key is None and event.button == 'up':
	     #print "Zoom out"
	     i1, i2, j1, j2 = disp.zoom_OUT(event.xdata, event.ydata)
	     ax.set_xlim(i1,i2)
             ax.set_ylim(j1,j2)
             draw()
           elif event.key is None and event.button == 'down':
	     #print "Zoom in"
	     i1, i2, j1, j2 = disp.zoom_IN(event.xdata, event.ydata)
             ax.set_xlim(i1,i2)
             ax.set_ylim(j1,j2)
             draw()
         
         
         change = False
         
         if event.button == 'up' and event.key == 'control':
	    a += delta_axis
	    slider_a.set_val(a)
	    draw()
	    
	 elif event.button == 'down' and event.key == 'control':
	    a -= delta_axis
	    slider_a.set_val(a)
	    draw()

	 elif event.button == 'up' and event.key == 'shift':
	    b += delta_axis
	    slider_b.set_val(b)
	    draw()
	    
	 elif event.button == 'down' and event.key == 'shift':
	    b -= delta_axis
	    slider_b.set_val(b)
	    draw()
	    
	 elif event.button == 'up' and event.key == 'alt':
	    pa += delta_angle
	    slider_pa.set_val(pa)
	    draw()
	    
	 elif event.button == 'down' and event.key == 'alt':
	    pa -= delta_angle
	    slider_pa.set_val(pa)
	    draw()
	    
	 elif event.button == 'up' and (event.key == 'alt+control' or  event.key == 'ctrl+alt'):
	    pa += 0.5
	    slider_pa.set_val(pa)
	    draw()
	    
	 elif event.button == 'down' and (event.key == 'alt+control' or event.key == 'ctrl+alt'):
	    pa -= 0.5
	    slider_pa.set_val(pa)
	    draw()

	    	    
	    
 
	    
	 else:
	    change = False
	    
	 if change:
            ellipse.update(np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))
            
         
         
         
   fig.canvas.mpl_connect('scroll_event', scroll_event)
   ####################### 
   
   
   
   def press_key(event):
         global a, b, pa, xcenter, ycenter, delta_axis, delta_angle, center_change, x_new, y_new
         #print 'you pressed', event.key
         
         if event.key == 'up':
            ycenter += delta_position
            change = True
         elif event.key == 'down':
	    ycenter -= delta_position
	    change = True
	 elif event.key == 'right':
	    xcenter += delta_position
	    change = True
	 elif event.key == 'left':
	    xcenter -= delta_position
	    change = True	    
  
	 elif event.key == 'ctrl+z':
	    undo(event)
	    change = False
	 elif event.key == 'ctrl+y':
	    redo(event)
	    change = False
	    
	 elif event.key == 'ctrl+a':
	    a = 1.1*a
	    b = 1.1*b
	    slider_a.disconnect(0)
            slider_a.set_val(a)
            slider_a.observers[0] = update
         
            slider_b.disconnect(0)
            slider_b.set_val(b)
            slider_b.observers[0] = update
	    change = True
	 elif event.key == 'ctrl+d':
	    a = 10.*a/11.
	    b = 10.*b/11.
	    slider_a.disconnect(0)
            slider_a.set_val(a)
            slider_a.observers[0] = update
         
            slider_b.disconnect(0)
            slider_b.set_val(b)
            slider_b.observers[0] = update
	    change = True	    
	    
	    
	 else: 
	    change = False
	 
	 if change:
            ellipse.update(np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))
            
         if event.key == 'q' or event.key == 'escape':
	    center_change = False
	    user_center.set_xdata([])
            user_center.set_ydata([])
            draw()

         
	 if center_change == True and event.key == 'enter':
	   xcenter = x_new
	   ycenter = y_new
           ellipse.update(np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))
           center_change = False
           user_center.set_xdata([])
           user_center.set_ydata([])
           draw()	        
        
        
   def on_click(event):    
         global a, b, pa, xcenter, ycenter, delta_axis, delta_angle, center_change, x_new, y_new
         #print 'you pressed', event.key, event.button, event.xdata, event.ydata, event.key
         
         if event.key == 'control' or event.key == 'shift':
	   addOne_axis(event)
	 elif event.key == 'alt':
	   addOne_angle(event)
	 	 
	 if event.dblclick or (event.key == 'z' and event.button == 1):
	   x_new = event.xdata
	   y_new = event.ydata
	   center_change = True
	   user_center.set_xdata([x_new])
	   user_center.set_ydata([y_new])
	   draw()
	   
	 if center_change == True and event.button == 3:
	   xcenter = x_new
	   ycenter = y_new
           ellipse.update(np.float(xcenter), np.float(ycenter), np.float(a), np.float(b), np.float(pa))
           center_change = False
           user_center.set_xdata([])
           user_center.set_ydata([])
           draw()	   
         
         if event.inaxes == ax:
	   if event.key is None and event.button == 2:
	     #print "Pan"
	     i1, i2, j1, j2 = disp.pan(event.xdata, event.ydata)
             ax.set_xlim(i1,i2)
             ax.set_ylim(j1,j2)
             draw()
        
        
   fig.canvas.mpl_connect('button_press_event', on_click)
   
   fig.canvas.mpl_connect('key_press_event', press_key)
   
   
   annotate("Coordinates ...", (0.016,0.73), xycoords='figure fraction', size=14, color='maroon')
   
   Ux = annotate(" ", (0.016,0.69), xycoords='figure fraction', size=12)
   Uy = annotate(" ", (0.016,0.65), xycoords='figure fraction', size=12)
   
   annotate("WCS ...", (0.016,0.60), xycoords='figure fraction', size=14, color='maroon')
   
   Ura = annotate(" ", (0.016,0.56), xycoords='figure fraction', size=11)
   Ualf = annotate(" ", (0.016,0.54), xycoords='figure fraction', size=10)
   Udec = annotate(" ", (0.016,0.50), xycoords='figure fraction', size=11)
   Udelt = annotate(" ", (0.016,0.47), xycoords='figure fraction', size=10)
   
   
   
   def in_motion(event):
     #print 'you pressed', event.key, event.button, event.xdata, event.ydata, event.key, event.inaxes
     x = event.xdata
     y = event.ydata
     if event.inaxes == ax: 
       
        Ux.set_text("X: "+'{:.2f}'.format(x))
        Uy.set_text("Y: "+'{:.2f}'.format(y))
        
        if HaveFitsFile:
	  world = w.wcs_pix2world([[x, y]], 1)
	  RA = world[0][0]
	  DEC = world[0][1]
	  Ura.set_text("RA: "+'{:.4f}'.format(RA))
	  Udec.set_text("DEC: "+'{:.4f}'.format(DEC))
	  c = SkyCoord(ra=RA, dec=DEC, unit=(u.degree, u.degree))
	  wcs_hmsdms =  c.to_string('hmsdms', precision=4, sep=':')
	  wcs_hmsdms = wcs_hmsdms.split(" ")
	  Ualf.set_text(r"$\alpha: $"+wcs_hmsdms[0])
	  Udelt.set_text(r"$\delta: $"+wcs_hmsdms[1])        
        
        draw()
     else:
        Ux.set_text(" ")
        Uy.set_text(" ")
        Ura.set_text(" ")
        Udec.set_text(" ")
        Ualf.set_text(" ")
	Udelt.set_text(" ")
        draw()
   
   
   fig.canvas.mpl_connect('motion_notify_event', in_motion)
   

   
   
   

   plt.show()
        



#################################################################
def db_id(ra):       
  
  ra_id = str(int(np.floor(ra)))
  if ra < 10:
    ra_id = '00'+ra_id+'D'
  elif ra < 100:
    ra_id = '0'+ra_id+'D'
  else:
    ra_id = ra_id+'D'        
  
  return ra_id

#################################################################

if __name__ == '__main__':
  
  if (len(sys.argv) < 2): 
     print >> sys.stderr, "\n Use \""+sys.argv[0]+" -h\" for help ...\n"
     exit(1)
     
     
  opts, args =  arg_parser()
  
  ID =  opts.pgc
  a = opts.major
  b = opts.minor
  xcenter = opts.xc
  ycenter = opts.yc
  pa = opts.pa
  dtype = opts.survey
  jpg_base = opts.jpg_base
  fits_base = opts.fits_base
  py_root = opts.py_root
  auxbase = opts.auxbase
  
  if ID == None: exit(1)
  if a == None: exit(1)
  if b == None: exit(1)
  if xcenter == None: exit(1)
  if ycenter == None: exit(1)
  if pa == None: exit(1)
  if dtype == None: exit(1)
  if jpg_base == None: exit(1)
  if fits_base == None: exit(1)
  if auxbase == None: exit(1)
  
  if py_root == None: 
    p = repr(os.getcwd())
    py_root = p[1:len(p)-1]    # removing quotation signs
  
  main()





