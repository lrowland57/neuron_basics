#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:51:25 2021

@author: lrow692
"""


#LIBRARIES
from neuron import h
try: #
    from neuron.units import ms, mV
except ModuleNotFoundError:
    ms = 1
    mV = 1   
import matplotlib.pyplot as plt
import pandas as pd
import csv
import pprint

soma = h.Section(name="soma") #create a neuron with one soma
h.topology()
pprint.pprint(soma.psection()) #properties
print(soma.psection()['morphology']['L'])

soma.L = 20
soma.diam = 20

soma.insert('hh')

print(type(soma(0.5)))

mech = soma(0.5).hh
print(type(mech))
print(soma(0.5).hh.gkbar) #variables

iclamp = h.IClamp(soma(0.5))
print([item for item in dir(iclamp) if not item.startswith('__')])
iclamp.delay =2
iclamp.dur = 0.1
iclamp.amp = 0.9
pprint.pprint(soma.psection())

v = h.Vector().record(soma(0.5)._ref_v)
t = h.Vector().record(h._ref_t)

#run
h.load_file('stdrun.hoc')
h.finitialize(-65*mV)
h.continuerun(40*ms)

#plotting

plt.figure()
plt.plot(t,v)
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.show()

#save and load
with open('data.csv', 'w') as f:
    csv.writer(f).writerow(zip(t,v))
      
data = pd.read_csv('data.csv', header = None, names = ['t','v'])
 
plt.figure()
plt.plot(data['t'],data["v"])
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.show()
    
    