import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from astropy.visualization import make_lupton_rgb
import sep
from matplotlib.widgets import Slider
plt.close('all')

bias_list = np.loadtxt("bias.csv", dtype = str)
h_alpha_flats = np.loadtxt("halp_flat.csv", dtype = str)
r_flats = np.loadtxt("r_flat.csv", dtype = str)
u_flats = np.loadtxt("u_flat.csv", dtype = str)
g_flats = np.loadtxt("g_flat.csv", dtype = str)

def bias(list):

   first_frame = fits.open(list[0])
   data = np.array(first_frame[1].data)
   ysize = data.shape[0]
   xsize = data.shape[1]

   nframes = len(list)

   bigbias = np.zeros((nframes,ysize,xsize),float)
   for i,file in enumerate(list):
      rawbias = fits.open(file)
      data = np.array(rawbias[1].data)
      bigbias[i-1,0:ysize-1,0:xsize-1] = data[0:ysize-1,0:xsize-1]

   medianbias = np.median(bigbias,axis=0)

   return  medianbias

BIAS = bias(bias_list)

def flats(list):
    first_frame = fits.open(list[0])
    data = np.array(first_frame[1].data)
    ysize = data.shape[0]
    xsize = data.shape[1]

    nframes = len(list)

    bigflat = np.zeros((nframes,ysize,xsize),float)
    for i,file in enumerate(list):
        rawflat = fits.open(file)
        data = np.array(rawflat[1].data)
        data = data-BIAS
        bigflat[i-1,0:ysize-1,0:xsize-1] = data[0:ysize-1,0:xsize-1]

    medianflat = np.median(bigflat,axis=0)

    medianflat = medianflat/np.mean(medianflat[500:1750, 500:1750])

    return  medianflat

h_alpha_flat = flats(h_alpha_flats)
r_flat = flats(r_flats)
u_flat = flats(u_flats)
g_flat = flats(g_flats)

H_ALPH = fits.open("ALIa050122.fits")[1].data
R = fits.open("ALIa050123.fits")[1].data
U = fits.open("ALIa050125.fits")[1].data
G = fits.open("ALIa050125.fits")[1].data

# Clip the images to remove extreme values
H_ALPH[h_alpha_flat > 1.5] = 0
H_ALPH[h_alpha_flat < 0.5] = 0

R[r_flat > 1.5] = 0
R[r_flat < 0.5] = 0

U[u_flat > 1.5] = 0
U[u_flat < 0.5] = 0

G[g_flat > 1.5] = 0
G[g_flat < 0.5] = 0

# this is not necessary for flats, but it can remove some warnings,
# and code runs much smoother without NaNs and Infs
h_alpha_flat[h_alpha_flat > 1.5] = 1
h_alpha_flat[h_alpha_flat < 0.5] = 1

r_flat[r_flat > 1.5] = 1
r_flat[r_flat < 0.5] = 1

u_flat[u_flat > 1.5] = 1
u_flat[u_flat < 0.5] = 1

g_flat[g_flat > 1.5] = 1
g_flat[g_flat < 0.5] = 1

# subtract bias and divide by flats
H_ALPH = (H_ALPH-BIAS)/h_alpha_flat
R = (R-BIAS)/r_flat
U = (U-BIAS)/u_flat
G = (G-BIAS)/g_flat


window_start_x = 200
window_start_y = 125
window_end_x = 1950
window_end_y = 1850

# crop the images
H_ALPH = np.ascontiguousarray(H_ALPH[window_start_y:window_end_y, window_start_x:window_end_x]) # H-alpha not used at the moment
B = np.ascontiguousarray(U[window_start_y:window_end_y, window_start_x:window_end_x]) # ultraviolet as blue
G = np.ascontiguousarray(G[window_start_y:window_end_y, window_start_x:window_end_x]) 
R = np.ascontiguousarray(R[window_start_y:window_end_y, window_start_x:window_end_x]) 

bkgB = sep.Background(B)
bkgG = sep.Background(G)
bkgR = sep.Background(R)
B = (B-bkgB.globalback)/bkgB.globalrms
G = (G-bkgG.globalback)/bkgG.globalrms
R = (R-bkgR.globalback)/bkgR.globalrms/1.5

# these will be used to scale intensity
R_factor = 1
G_factor = 1
B_factor = 1

# Interactive plot
def update(val=None):
    Q = 5.
    stretch = 25.
    minimum = [-1, -1, -1]
    image = make_lupton_rgb(R_factor * R, G_factor * G, B_factor * B, Q=Q, stretch=stretch, minimum=minimum)
    # Plot image
    ax.imshow(image, origin='lower')
    fig.canvas.draw_idle()

# Create the figure and axis
fig, ax = plt.subplots(figsize=(12, 12), dpi=80)
plt.subplots_adjust(left=0.1, bottom=0.25)

# Initial plot
update()

# Add sliders for intensity adjustment
axcolor = 'lightgoldenrodyellow'
ax_r_factor = plt.axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_g_factor = plt.axes([0.1, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_b_factor = plt.axes([0.1, 0.2, 0.65, 0.03], facecolor=axcolor)

s_r_factor = Slider(ax_r_factor, 'Red', 0, 10.0, valinit=1)
s_g_factor = Slider(ax_g_factor, 'Green', 0, 10.0, valinit=1)
s_b_factor = Slider(ax_b_factor, 'Blue', 0, 10.0, valinit=1)

def update_factors(val):
    global R_factor, G_factor, B_factor
    R_factor = s_r_factor.val
    G_factor = s_g_factor.val
    B_factor = s_b_factor.val
    update()

s_r_factor.on_changed(update_factors)
s_g_factor.on_changed(update_factors)
s_b_factor.on_changed(update_factors)

plt.show()

Q = 5.
stretch = 25.
minimum = [-1, -1, -1]
image = make_lupton_rgb(R_factor * R, G_factor * G, B_factor * B, Q=Q, stretch=stretch, minimum=minimum)
plt.imshow(image, origin='lower')
plt.axis('off')
plt.savefig('rgb.png', bbox_inches='tight')
