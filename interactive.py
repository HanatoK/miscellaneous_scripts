#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 15:29:12 2018

@author: hanatok
"""

import numpy as np
import plotly
import plotly.graph_objs as go

plotly.offline.init_notebook_mode(connected=True)

xi, yi, zi, err = np.genfromtxt('/home/hanatok/HDD/documents/calixarene/pmf/error.pmf', unpack = True)
px, py = np.genfromtxt('/home/hanatok/HDD/documents/calixarene/pmf/mepsa_res1.txt', unpack = True)
epx, epy, epz = np.genfromtxt('/home/hanatok/HDD/documents/calixarene/pmf/path_experiment.txt', unpack = True)
data = [go.Contour(z = zi, x = xi, y = yi, colorscale='Jet', autocontour=False, contours=dict(start=0, end=20, size=1)), go.Scatter(x = px, y = py), go.Scatter(x = epx, y = epy)]
layout = go.Layout(autosize=False, width=1000, height=1000)

fig = go.Figure(data=data, layout=layout)
plotly.offline.iplot(fig)
