'''
Created on Oct 17, 2013

@author: tohekorh
'''

from help_classes import  valminmax, reduce_bounds, get_quess
from datetime import datetime
from deform_energies import deform_energies

import numpy as np

class optimize():
    
    def __init__(self, ue, optimizer): #, energies
        self.opt                    =   ue.opt
        self.hangle                 =   ue.hangle
        [self.rmin, self.rmax],     \
        [self.phimin, self.phimax]  =   valminmax(ue.phys_surf)
        self.nr, self.nphi          =   ue.phys_surf.shape
        self.H                      =   self.hangle*2*np.pi
        self.energies               =   deform_energies()
        self.optimizer              =   optimizer
        self.ue                     =   ue
        
        
    def optimize_consts(self):
    
        from scipy.optimize import fmin_l_bfgs_b #,  fmin_tnc #,  bfgs, lbfgs optimointiin
        print 'optimizing.. consts'
        
        def energy(consts):

            print str(consts) + ',',
            
            E_b, E_s                =   self.energies.calc_energies(consts)[:2]
            grad                    =   self.energies.calc_grad(consts, E_b, E_s)
            
            print 'tot E = ' + str(E_b + E_s) + ', E_b = ' +  str(E_b) + ', E_s = ' +  \
                    str(E_s) + ', time: ' + str(datetime.now()) 
            
            return E_s + E_b, grad
        
        self.energies.set_u(self.ue)
        
        if self.ue.system != 'spiral':
                
            i           =   0
            E_wave      =   [0, -1]
            data_arr    =   []
            out         =   False
            
            while i < 2 or not out:
                print 'num of waves = ' + str(i)
                
                o_consts_set, bounds=   get_quess(self.opt, self.ue.ext_surf, \
                                                  self.hangle, self.phimax- self.phimin, i) 
                
                E_min               =   1e8
                for consts in o_consts_set:
                    
                    print 'Amplitude, A = %.2f, bounds A = ' %consts[1] + str(bounds[1])
                    
                    consts_set          =   False
                    n                   =   0
                    while not consts_set:
                        if n    == 30:
                            consts[0]   =   0.
                            consts[1]   =   0.
                        elif n  == 31:
                            raise ValueError
                        
                        try:
                            self.ue.set_const(consts)
                            consts_set  =   True
                        except ValueError:
                            n          +=   1
                            consts      =   reduce_bounds(consts, bounds)[0]
                        
                    
                    converged           =   False
                    n                   =   0
                    while not converged:
                        if n    == 30:
                            consts[0]   =   0.
                            consts[1]   =   0.
                        elif n  == 31:
                            raise ValueError
                        
                        try:
                            ol_consts, E=   fmin_l_bfgs_b(energy, consts[:-1], \
                                                          bounds = bounds)[:2]
                            
                            converged   =   True
                        except ValueError:
                            #print 'red bounds'
                            n          +=   1
                            bounds      =   reduce_bounds(consts, bounds)[1]
                
                    if E < E_min:
                        E_min           =   E
                        min_consts      =   ol_consts    
                        
                        
                o_consts            =   np.append(min_consts, [i])
                
                E_b, E_s, E_b_surf, E_s_surf, normals \
                                    =   self.energies.calc_energies(o_consts)
                
                if i in [0, 1]:
                    E_wave[i]       =   E_b + E_s 
                    if i == 0:
                        data_arr.append([E_b, E_s, E_b_surf, E_s_surf, normals, o_consts])
                else:
                    E_wave.append(E_b + E_s)
    
                if 0 < i:
                    if E_wave[i] < E_wave[i - 1] and 1e-3 < o_consts[1]:
                        data_arr.append([E_b, E_s, E_b_surf, E_s_surf, normals, o_consts])
                    else:
                        out         =   True
                        data        =   data_arr[-1][:5]
                        set_consts  =   data_arr[-1][5] 
                i                  +=   1
                    
            #print 'set_consts' + str(set_consts)
    
            self.ue.set_const(set_consts)
    
            E_b, E_s, E_b_surf, E_s_surf, normals  = data
            
            return [E_b_surf, E_s_surf], [E_b, E_s], self.ue, normals
        
        elif self.ue.system == 'spiral':
            
            quess_consts, bounds=   get_quess(self.opt, self.ue.ext_surf, \
                                              self.hangle, self.phimax- self.phimin, 0) 
                
            o_consts, E=   fmin_l_bfgs_b(energy, quess_consts, \
                                              bounds = bounds)[:2]
                
            E_b, E_s, E_b_surf, E_s_surf, normals \
                               =   self.energies.calc_energies(o_consts)
                
            self.ue.set_const(o_consts)
            
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
        
            
        o_flat_umat     =   fmin_l_bfgs_b(energy, self.ue.flat_umat, bounds = self.ue.bounds, approx_grad = True)[0]
        
        E_b, E_s, E_b_surf, E_s_surf, normals =  \
                                    self.energies.calc_energies_moldy(o_flat_umat)
        
        self.ue.parse_u(o_flat_umat)
        
        return [E_b_surf, E_s_surf], [E_b, E_s], self.ue, normals
    
    