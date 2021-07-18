#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:18:32 2021

@author: lrow692
"""

#LIBRARIES
from neuron import h, gui
try: #
    from neuron.units import ms, mV
except ModuleNotFoundError:
    ms = 1
    mV = 1
    
import matplotlib.pyplot as plt
    
h.load_file("stdrun.hoc")


class Cell:
    def __init__(self, gid, x, y, z, theta):
        """initializes cell object"""
        self._gid = gid
        self._setup_morphology()
        try:
            self.all = self.soma.wholetree()
        except AttributeError:
            self.all = []
            for key in self.__dict__:
                if str(type(self.__dict__[key])) == "<class 'nrn.section'>":
                    self.all.append(self.__dict__[key])
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        self._rotate_z(theta)
        self._set_position(x, y, z)
        
        self._spike_detector = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        
        self._ncs = []
        
        self.soma_v = h.Vector().record(self.soma(0.5)._ref_v)
        
    def __repr__(self):
        """represents cell as a string"""
        return "{}[{}]".format(self.name, self._gid)
    
    def _set_position(self, x, y, z):
        """assigns/sets position"""
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3change(i,
                              x - self.x + sec.x3d (i),
                              y - self.y + sec.y3d (i),
                              z - self.z + sec.z3d (i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
            
    def _rotate_z(self, theta):
        """Rotate the cell about the z axis."""
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = h.cos(theta)
                s = h.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))
                
class BallAndStick(Cell):
    name = 'BallAndStick'
    def _setup_morphology(self):
        """sets up cell morphology"""
        self.soma = h.Section(name = "soma", cell = self) #soma
        self.dend = h.Section(name="dend", cell=self) #dendrite
        #self.all = [self.soma,self.dend]
        self.dend.connect(self.soma) #connect soma
        self.soma.L = self.soma.diam = 12.6157 #diameter
        self.dend.L = 200
        self.dend.diam = 1
        
    def _setup_biophysics(self):
        """sets up cell biophysics"""
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
            self.syn = h.ExpSyn(self.dend(0.5))
            self.syn.tau = 2 * ms
            
            
class Ring:
    """A network of N ball-and-stick cells wjere cell n makes an excitatory synapse
    onto cell n+ 1 and the ;ast Nth cell in the network projects to teh first cell"""
    def __init__(self, N = 5, stim_w = 0.04, stim_t = 9, 
                 stim_delay = 1, syn_w = 0.01, syn_delay = 5, r = 50):
        """
        :Param N: Number of cells
        :param stim_w: weight of the stimulus
        :param stime_t: time of the stimulus
        :param stim_delay: delay of the stimulus
        :[aram syn_w: synaptc weight
        :param syn_delay: Delay of the syna[se
        :param r: radius of the netwrok
        """
        
        self._syn_w = syn_w
        self._syn_delay = syn_delay
        self._create_cells(N, r)
        self._connect_cells()
        #add stimulus
        self._netstim = h.NetStim()
        self._netstim.number = 1
        self._netstim.start = stim_t
        self._nc = h.NetCon(self._netstim, self.cells[0].syn)
        self._nc.delay = stim_delay
        self._nc.weight[0] = stim_w
        
    def _create_cells(self, N, r):
        self.cells = []
        for i in range(N):
            theta = i * 2 * h.PI / N
            self.cells.append(BallAndStick(i, h.cos(theta) * r, h.sin(theta) * r, 0, theta))
            
    def _connect_cells(self):
        for source, target in zip(self.cells, self.cells[1:] + [self.cells[0]]):
            nc = h.NetCon(source.soma(0.5)._ref_v, target.syn, sec = source.soma)
            nc.weight[0] = self._syn_w
            nc.delay = self._syn_delay
            source._ncs.append(nc)
            
            
ring = Ring(N=5)
shape_window = h.PlotShape(True)
shape_window.show(0)

t = h.Vector().record(h._ref_t)
h.finitialize(-65*mV)
h.continuerun(100)
            
plt.plot(t, ring.cells[0].soma_v)
plt.show()

plt.figure()
for syn_w, color in [(0.01, "black"), (0.005, 'red')]:
    ring = Ring(N=5, syn_w=syn_w)
    h.finitialize(-65*mV)
    h.continuerun(100*ms)
    for i, cell in enumerate(ring.cells):
        plt.vlines(cell.spike_times, i + 0.5, i + 1.5, color = color)

plt.show()
                
                
                
                
                
                
                