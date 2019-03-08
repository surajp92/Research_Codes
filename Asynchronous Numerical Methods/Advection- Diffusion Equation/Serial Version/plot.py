 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt

with open('field_final.dat') as f:
    lines = (line for line in f if not line.startswith('#'))
    FH = np.loadtxt(lines, delimiter=',', skiprows=2)

plt.plot(FH[:,0],FH[:,1], FH[:,0], FH[:,2])