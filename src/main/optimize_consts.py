'''
Created on Oct 17, 2013

@author: tohekorh
'''

from help_classes import  valminmax #rs, phis,
#from strain import u
#from tests import test_prime_maps
from surface import surf #, e_surf
from datetime import datetime
#from plots import plot_all
from help_classes import get_quess
import numpy as np
from deform_energies import deform_energies

class optimize():
    
    def __init__(self, ue, optimizer): #, energies
        self.opt                    =   ue.opt
        self.hangle                 =   ue.hangle
        [self.rmin, self.rmax],     \
        [self.phimin, self.phimax]  =   valminmax(ue.phys_surf)
        self.nr, self.nphi          =   ue.phys_surf.shape
        self.H                      =   self.hangle*2*np.pi
        #self.energies               =   energies
        self.energies               =   deform_energies()
        self.optimizer              =   optimizer
        self.ue                     =   ue
        
        
    def optimize_consts(self):

        
    
        from scipy.optimize import fmin_l_bfgs_b  # fmin, fmin_tnc,  bfgs, lbfgs optimointiin
        print 'optimizing.. consts'
        
        def energy(consts):

            print str(consts) + ',',
            
            E_b, E_s                =   self.energies.calc_energies(consts)[:2]
            grad                    =   self.energies.calc_grad(consts, E_b, E_s)
            
            #print 'grad = ' + str(grad)        
            print 'tot E = ' + str(E_b + E_s) + ', E_b = ' +  str(E_b) + ', E_s = ' +  \
                    str(E_s) + ', time: ' + str(datetime.now()) 
            #print ''
            
            return E_s + E_b, grad
        
        
        #n                          =   [self.nr, self.nphi]
        
        #for iin, n in enumerate(ns):
            
        #print 'grid n    = '    +   str([self.nr, self.nphi]) + ', optimizer = ' + self.optimizer
        #print 'Height    = '    +   str(self.H)
        
        #nr_o, nphi_o            =   n    
        asurf                   =   surf(self.rmin, self.rmax, self.nr, \
                                         self.phimin, self.phimax, self.nphi)
        
        self.ue.set_ainit_surf(asurf.get_all_surf())  
        self.energies.set_u(self.ue)
                
        #if iin == 0:
        
        #o_consts=   get_quess(self.opt, asurf.get_ext_surf(), \
        #                     self.hangle, self.phimax- self.phimin, 1)[0] 
        
            
            
        #self.ue.set_const(o_consts)
        #self.energies.set_u(self.ue)
                  
        #o_consts                =   fmin_l_bfgs_b(energy, o_consts)[0]                      
        i = 0
        E_wave      =   [0, -1]
        data_arr    =   []
        out         =   False
        
        while i < 2 or not out:
            
            o_consts, bounds    =   get_quess(self.opt, asurf.get_ext_surf(), \
                                              self.hangle, self.phimax- self.phimin, i) 
                                
            #bounds              =   get_quess(self.opt, asurf.get_ext_surf(), \
            #                                  self.hangle, self.phimax- self.phimin, i)[1] 
            
            
            #o_consts[4]         =   i #number of waves in unit cell
            self.ue.set_const(o_consts)
            
            print 'num of waves = ' + str(i)
            #print 'o_consts     = ' + str(o_consts)
            #print 'bounds amplit= ' + str(bounds[1])
            
            o_consts            =   fmin_l_bfgs_b(energy, o_consts[:-1], \
                                                  bounds = bounds)[0]
            o_consts            =   np.append(o_consts, [int(i)])
            
            E_b, E_s, E_b_surf, E_s_surf, normals \
                                =   self.energies.calc_energies(o_consts)
            
            if i in [0, 1]:
                E_wave[i]           =   E_b + E_s 
                
                if i == 0:
                    data_arr.append([E_b, E_s, E_b_surf, E_s_surf, normals, o_consts])
                    print data_arr[i][5]
            
            else:
                E_wave.append(E_b + E_s)

            if 0 < i:
                if E_wave[i] < E_wave[i - 1] and 1e-3 < o_consts[1]:
                    data_arr.append([E_b, E_s, E_b_surf, E_s_surf, normals, o_consts])
                    
                    #print 'write consts = ' + str(data_arr[i][5])
                    
                else:
                    out         =   True
                    data        =   data_arr[-1][:5]
                    set_consts    =   data_arr[-1][5] 
                    #o_consts[4]-=   1
                    
            i                  +=   1
                
                
        
        
        #print 'set_consts' + str(set_consts)
        
        self.ue.set_const(set_consts)

        E_b, E_s, E_b_surf, E_s_surf, normals  = data
        #E_b, E_s, E_b_surf, E_s_surf, normals =  \
        #                                    self.energies.calc_energies(o_consts)
        
        
        return [E_b_surf, E_s_surf], [E_b, E_s], self.ue, normals
    
    
    
    def optimize_moldy(self):

        from scipy.optimize import fmin, fmin_tnc, fmin_l_bfgs_b  # fmin, fmin_tnc,  bfgs, lbfgs optimointiin
        print 'optimizing.. moldy'
        self.count                  =   0
        
        def energy(flat_u):
            
            E_b, E_s        =   self.energies.calc_energies_moldy(flat_u, self.count)[:2]
                      
            if self.count % 1000 == 0:
                print   'count '    +   str(self.count) + ', tot E = '   +   str(E_b + E_s) +  \
                        ', E_b  = ' +   str(E_b) + ' , E_s = '  +   str(E_s)        
            
            self.count     +=   1
            return E_s + E_b
        
        
        
        if self.energies.u == None:
            self.energies.set_u(self.ue)
        
        if not self.ue.moldy:
            self.ue.set_moldy(True)  
        
        if self.optimizer == 'l_bfgs':
            o_flat_umat     =   fmin_l_bfgs_b(energy, self.ue.flat_umat, bounds = self.ue.bounds, \
                                                  approx_grad = True, epsilon = 1e-5)[0]
        
        elif self.optimizer == 'tnc':
            o_flat_umat     =   fmin_tnc(energy, self.ue.flat_umat, bounds = self.ue.bounds, approx_grad = True)[0]
        
        else:
            o_flat_umat     =   fmin(energy, self.ue.flat_umat)
        

        E_b, E_s, E_b_surf, E_s_surf, normals =  \
                                    self.energies.calc_energies_moldy(o_flat_umat)
        
        self.ue.parse_u(o_flat_umat)
        
        return [E_b_surf, E_s_surf], [E_b, E_s], self.ue, normals
    
    