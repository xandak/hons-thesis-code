"""
Set pyplot parameters for nice plots
"""

import matplotlib as mpl
from scipy import sqrt

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]

# Explicitly set fontsizes:
font_size = 20
tick_size = 18
params = {
          #'backend': 'ps',
          'axes.labelsize': font_size,
          'font.size': font_size,
          'legend.fontsize': font_size,
          'xtick.labelsize': tick_size,
          'ytick.labelsize': tick_size,
          'text.usetex': True,
          'figure.figsize': fig_size,
          }
mpl.rcParams.update(params)
adjustprops = dict(right=0.95, left=0.14, bottom=0.17, top=0.93, 
                   wspace=0.0, hspace=0.0)

figprops = dict(figsize=(1.0*fig_width,1.0*fig_height))

def subplots_adjust(fig, props):
    fig.subplots_adjust(left=props['left'])
    fig.subplots_adjust(right=props['right'])
    fig.subplots_adjust(bottom=props['bottom'])
    fig.subplots_adjust(top=props['top'])
    fig.subplots_adjust(wspace=props['wspace'])
    fig.subplots_adjust(hspace=props['hspace'])
