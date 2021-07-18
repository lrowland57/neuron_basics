#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 15:36:37 2021

@author: lrow692
"""
import matplotlib.pyplot as plt

try: #
    from neuron.units import ms, mV
except ModuleNotFoundError:
    ms = 1
    mV = 1
    
h.load_file("stdrun.hoc")

from neuron import h

class BallAndStick:
    def __init__(self, gid):
        self._gid = gid #Make the relationship
        self._setup_morphology()
        self._setup_biophysics()
    def _setup_morphology(self):
        self.soma = h.Section(name = "soma", cell = self) #soma
        self.dend = h.Section(name="dend", cell=self) #dendrite
        self.all = [self.soma,self.dend]
        self.dend.connect(self.soma) #connect soma
        self.soma.L = self.soma.diam = 12.6157 #diameter
        self.dend.L = 200
        self.dend.diam = 1
        
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100 #axial radius
            sec.cm = 1 #membrance capitance
        self.soma.insert("hh")
        for seg in self.soma:
            seg.hh.gnabar = 0.12 #sodium
            seg.hh.gkbar = 0.036 #potassium
            seg.hh.gl = 0.0003 #leak conductance
            seg.hh.el = -54.3 #reversal potential
            
            #passive current in dendrite
            self.dend.insert("pas")
            for seg in self.dend:
                seg.pas.g = 0.001 #passive conductance
                seg.pas.e = -65 # leak reversal potential
        
    def __repr__(self):
        return "BallAndStick[{}]".format(self._gid)
    
my_cell = BallAndStick(0)
#my_other_cell = BallAndStick(1)


try:
    h.PlotShape(False).plot(plt)
    ps = h.PlotShape(True)
except AttributeError:
    ps = h.PlotShape(True)
ps.show(0)


stim = h.IClamp(my_cell.dend(1))
stim.delay = 5
stim.dur = 1
stim.amp = 0.1

soma_v = h.Vector().record(my_cell.soma(0.5)._ref_v)
t = h.Vector().record(h._ref_t)
dend_v = h.Vector().record(my_cell.dend(0.5)._ref_v)

f1 = plt.figure()
plt.xlabel("t (ms)")
plt.ylabel("v (mV)")
plt.plot(t, soma_v, linewidth=2)

f = plt.figure()
plt.xlabel("t (ms)")
plt.ylabel("v (mV)")
amps = [0.075*i for i in range(1,5)]
colors = ["green","blue","red","black"]
for amp, color in zip(amps, colors):
    stim.amp = amp
    for my_cell.dend.nseg, width in [(1,2),(101,1)]:
        h.finitialize(-65)
        h.continuerun(25)
        plt.plot(t,list(soma_v),
                 linewidth=width,
                 label = "amp=%g" % amp if my_cell.dend.nseg == 1 else None,
                 color = color)
        plt.plot(t, list(dend_v), "--",
                 linewidth = width,
                 color = color)
plt.legend()
plt.show(f)




for sec in h.allsec():
    print("%s: %s" % (sec, ", ".join(sec.psection()["density_mechs"].keys())))
    
    
    
















