
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from gaslight.grid import Grid
from gaslight.line import (
    get_diagram_labels,
    get_ratio_label,)
from gaslight import line_ratios
import streamlit as st

# set style
plt.style.use('matplotlibrc.txt')


st.image('elements/banner.png')

# set up streamlit
col1, col2 = st.columns(spec=[0.4, 0.6])

col2.markdown("## Gaslight BPT explorer")
col2.markdown("Use the sliders on the right to explore the effect of various parameters on the BPT diagram.")


grid_dir = 'data/'
grid_name = 'bpass-2.2.1-bin_chabrier03-0.1,300.0-ages:6.,7.,8.-c23.01-full-BPT-NII'

grid = Grid(grid_dir=grid_dir, grid_name=grid_name)


default_point = {
    'log10age': 0,
    'alpha': 1,
    'abundance_scalings.nitrogen_to_oxygen': 2,
    'abundance_scalings.carbon_to_oxygen': 2,
    'depletion_scale': 2,
    'ionisation_parameter': 2,
    'hydrogen_density': 2,
}

diagram_limits = {
    'BPT-NII': [[-4.,1.],[-3.,2.0]]
}

diagram_id = 'BPT-NII'

axes = copy.deepcopy(grid.axes)
axes.remove('metallicity')

slider_values = {}
for axis in axes:
    default_value = grid.axes_values[axis][default_point[axis]]
    slider_values[axis] = col1.select_slider(
        axis,
        value=default_value,
        options=grid.axes_values[axis])


fig = plt.figure(figsize=(3.5, 3.5))

bottom = 0.15
height = 0.8
left = 0.15
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


if diagram_id == 'BPT-NII':

    # plot Kewley and Kauffmann lines 
    for f, ls, limit, label in zip([line_ratios.get_bpt_kewley01], # , line_ratios.get_bpt_kauffman03
                                    ['-', '--'],
                                    [0.47, 0.05],
                                    ['Kewley+2001', 'Kauffman+2003']):
        log10x = np.arange(-5., limit, 0.01)
        ax.plot(log10x, f(log10x), ls=ls, lw=2, c='k', alpha=0.2, label=label)



x = []
y = []

colours = cm.magma(np.linspace(0, 0.7, len(grid.metallicity)))

for metallicity, colour in zip(grid.metallicity, colours):

    grid_value_dict = {'metallicity': metallicity} | slider_values

    grid_point = grid.get_nearest_grid_point(grid_value_dict)

    lines = grid.get_line_collection(grid_point)

    x_, y_ = lines.get_diagram(diagram_id)

    x.append(x_)
    y.append(y_)
    ax.scatter(np.log10(x_), np.log10(y_), s=15, alpha=1, zorder=2, c=colour, label=f'Z={metallicity}')

ax.plot(np.log10(x), np.log10(y), c='k', lw=2, alpha=0.5, zorder=1)

if diagram_id in diagram_limits:
    xlim, ylim = diagram_limits[diagram_id]
else:
    xlim = [-5., 1.5]
    ylim = [-3., 1.5]

ax.set_xlim(xlim)
ax.set_ylim(ylim)

x_label, y_label = get_diagram_labels(diagram_id)

# add axes labels
ax.set_xlabel(rf'${x_label}$')
ax.set_ylabel(rf'${y_label}$')

ax.legend(fontsize=6, labelspacing=0.1, loc='upper left')

col2.pyplot(fig)
