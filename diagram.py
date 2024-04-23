
import copy
import numpy as np
import matplotlib.pyplot as plt
from gaslight.grid import Grid
from synthesizer.line import (
    get_diagram_labels,
    get_ratio_label,)
from synthesizer import line_ratios
import streamlit as st

# set style
plt.style.use('matplotlibrc.txt')

# set up streamlit
col1, col2 = st.columns(spec=[0.4, 0.6])




grid_dir = 'data/'
grid_name = 'bpass-2.2.1-bin_chabrier03-0.1,300.0-ages:6.,7.,8.-c23.01-full-BPT-NII'

grid = Grid(grid_dir=grid_dir, grid_name=grid_name)


diagram_limits = {
    'BPT-NII': [[-4.,1.],[-4.,1]]
}

diagram_id = 'BPT-NII'

axes = copy.deepcopy(grid.axes)
axes.remove('metallicity')

slider_values = {}
for axis in axes:
    slider_values[axis] = col1.select_slider(axis, options=grid.axes_values[axis])


fig = plt.figure(figsize=(3.5, 3.5))

bottom = 0.15
height = 0.8
left = 0.15
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

x = []
y = []

for metallicity in grid.metallicity:

    grid_value_dict = {'metallicity': metallicity} | slider_values

    grid_point = grid.get_nearest_grid_point(grid_value_dict)

    lines = grid.get_line_collection(grid_point)

    x_, y_ = lines.get_diagram(diagram_id)

    x.append(x_)
    y.append(y_)

ax.plot(np.log10(x), np.log10(y))

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

col2.pyplot(fig)
