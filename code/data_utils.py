import os
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.constants import k_B, h, c, m_p, m_e
k_B = k_B.to(u.erg / u.K)
h = h.to(u.erg * u.s)
c = c.to(u.cm / u.s)
m_p = m_p.to(u.g)
m_e = m_e.to(u.g)


# Loading the data from models
fold = '../data/'
files = ['t4000.dat', 't6000.dat', 't8000.dat']
data_dic = {file.split('.')[0]: file for file in files}
for file in files:
    file_name = file.split('.')[0]
    data_dic[file_name] = {}

    path = fold + file
    with open(path, 'r') as file:
        lines = file.readlines()

    d1 = lines[:11]
    d2 = lines[11:23]

    d1 = ''.join(d1)
    d2 = ''.join(d2)

    data_dic[file_name]['d1'] = d1
    data_dic[file_name]['d2'] = d2

    d3 = lines[24:]
    d3 = ''.join(d3)
    data_io = StringIO(d3)
    df = pd.read_csv(data_io, delim_whitespace=True)
    df = df.iloc[:, 1:]

    data_dic[file_name]['data'] = df


atomic_data = {'HI': {}, 'HII': {}, 'H-': {}}
# Values of partition function
Us = [2, 1, 1]
gs = [[2, 8, 18], [], []]

wave_numbers_HI = [0,
                   82259.158 ,
                   97492.304 ,]
wave_number = [wave_numbers_HI, 
               [], 
               [],
               [],[],[]]

ion_en = [13.59844 * u.eV, 0 * u.eV, 0.755 * u.eV]

for i, el in enumerate(atomic_data.keys()):
    atomic_data[el]['U'] = Us[i]
    atomic_data[el]['g'] = gs[i]
    atomic_data[el]['wave_number'] = wave_number[i]
    atomic_data[el]['ex_en'] = [(h*c * wave_number[i][k]*u.cm**(-1)).to(u.eV) for k in range(len(wave_number[i]))]
    atomic_data[el]['ion_en'] = ion_en[i]


atomic_data_ext = {'HI': {}, 'HII': {}, 'H-': {}, 'HeI': {}, 'HeII': {}, 'HeIII': {}}
# Values of partition function
Us = [2, 1, 1, 1, 2, 1]
gs = [[2, 2, 4], [], [], [], [], []]

wave_numbers_HI = [0* u.cm**(-1),
                   82258.9543992821* u.cm**(-1),
                   1200000* u.cm**(-1)]

#82258.9191133

wave_number = [wave_numbers_HI, 
               [], 
               [],
               [],[],[]]

ion_en = [13.59844 * u.eV, 0 * u.eV, 0.755 * u.eV,
          24.58741 * u.eV, 54.41778 * u.eV, 0 * u.eV]

for i, el in enumerate(atomic_data_ext.keys()):
    atomic_data_ext[el]['U'] = Us[i]
    atomic_data_ext[el]['g'] = gs[i]
    atomic_data_ext[el]['wave_number'] = wave_number[i]
    atomic_data_ext[el]['ion_en'] = ion_en[i]