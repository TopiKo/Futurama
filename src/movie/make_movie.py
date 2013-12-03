'''
Created on Oct 17, 2013

@author: tohekorh
'''
import numpy as np

from main.bender_rw import read_h_consts, read_energies, read_umat
from main.deform_energies import deform_energies
from main.tests import test_prime_maps
from main.strain import u
from main.surface import surf
from main.plots import plot_all, plot_energies
    
pi          =   np.pi
opts        =   [[2, 3, 1], [1,0,0]]
ir          =   [8, 16]
ih          =   5
n           =   None

all_data    =   []
heights     =   []
energies    =   []

for opt in opts:
    data    =   read_h_consts(ir, opt)
    all_data.append(data)
    heights.append(data[('height')])
    energies.append(data[('energies')])

heights_hb, energies_hb =   read_energies(ir)[:2]

H           =   heights[0][ih]
hangle      =   H/2/pi

#energies        =   deform_energies(ue)

plot_energies(ir, heights, energies, heights_hb, energies_hb)

for all_data_sys in all_data:
    for ida, data in enumerate(all_data_sys):
        
        if ida == n or n == None:
        
            H                       =   data[0]
            consts                  =   data[1]
            [nr, rmin, rmax]        =   data[3]
            [nphi, phimin, phimax]  =   data[4]
            
            hangle                  =   H/2/pi
            
            rs                      =   np.linspace(rmin, rmax, nr)
            phis                    =   np.linspace(phimin, phimax, nphi)  
            asurf                   =   surf(rs, phis)
            
            print consts, H
            ue                      =   u(opt, hangle, asurf.get_all_surf(), \
                                          phimax - phimin, 'roto_trans', consts)
            
            print ir
            umat                    =   read_umat(ir, opt, H)
    
            ue.set_u(umat) 
            energies                =   deform_energies(ue)
        
            if opt[1] == 3:
                print 'testing phi map..' 
                test_prime_maps(ue)
            
            print 'E_s'
            E_s_surf, E_s           =   energies.F_s()
        
            print 'E_b'
            E_b_surf, E_b, normals  =   energies.F_b()
            
            plot_all(ue, E_s_surf, E_b_surf, normals)
    
