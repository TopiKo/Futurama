'''
Created on Oct 3, 2013

@author: tohekorh
'''

import numpy as np
from help_classes import make_periodic

class surf():
    
    
    
    def __init__(self, rmin, rmax, nr, phimin, delta_phi, nphi):
        
        phimax      =   phimin + delta_phi
        #S = (r, phi, z)
        rs          =   np.linspace(rmin, rmax, nr)
        phis        =   np.linspace(phimin, phimax, nphi)  
        
        r_0         =   rs[0] - (rs[1] - rs[0])
        r0          =   np.array([r_0])
        r_f         =   rs[-1] + (rs[-1] - rs[-2])
        rf          =   np.array([r_f])

        rs          =   np.concatenate((r0, rs, rf))

        self.ext_surf         =   np.empty((len(rs), len(phis) - 1),         dtype='object')
        
        for ir, r in enumerate(rs):
            for iphi, phi in enumerate(phis[:-1]):
                self.ext_surf[ir,iphi]           =   np.array([r, phi, 0.])
                #if ir < len(rs) - 1 and iphi < len(phis) - 1:
                #    self.d_ext_surf[ir,iphi]     =   np.array([r + delta_r(rs, ir)/2., phi + delta_phi(phis, iphi)/2., 0.])


        #phi_period          =   max(phis) - min(phis)
        self.phys_surf      =   make_periodic(np.delete(np.delete(self.ext_surf, 0, 0), -1, 0), 'surf', phi_period = delta_phi)
        self.calc_surf      =   np.delete(np.delete(self.phys_surf, -1, 0), -1, 1)
        
        '''
        for i in range(len(rs)):
            for j in range(len(phis)):
                self.surf_matrix[i,j]           =   np.array([rs[i], phis[j], 0.])
                if i < len(rs) - 1 and j < len(phis) - 1:
                    self.d_surf_matrix[i,j]     =   np.array([rs[i] + delta_r(rs, i)/2., phis[j] + delta_phi(phis, j)/2., 0.])
                if i < len(rs) - 2 and j < len(phis) - 2:
                    self.calc_surf[i,j]    =   np.array([rs[i] + delta_r(rs, i), phis[j] + delta_phi(phis, j), 0.])
        '''
                
    def get_ext_surf(self):
        return self.ext_surf 

    def get_phys_surf(self):
        return self.phys_surf 

    def get_calc_surf(self):
        return self.calc_surf 
    
    def get_all_surf(self):
        return [self.ext_surf, self.phys_surf, self.calc_surf]


def get_rs(surface, iphi):
    
    rs  = np.zeros(len(surface[iphi]))
    for ir, r_vec in enumerate(surface[:,0]):
        rs[ir]  =   r_vec[0] 
    
    return rs   
        
    
def transform_to_cart(surf):
    
    cart_surf   = np.empty(surf.shape, dtype= 'object')
    
    for ir, r_vec in np.ndenumerate(surf):
        r, phi, z           =   r_vec
        cart_surf[ir]       =   np.zeros(3)
        cart_surf[ir][0]    =   r*np.cos(phi)   
        cart_surf[ir][1]    =   r*np.sin(phi)
        cart_surf[ir][2]    =   z   
    
    return cart_surf    

def parse_surf(mat_r, mat_phi, mat_z, phi_period, sym_op, key):
    
    mat         =   np.empty(mat_r.shape, dtype = 'object')
    
    for index, val in np.ndenumerate(mat):
        mat[index]    =   np.array([mat_r[index], mat_phi[index], mat_z[index]])
    
    if key      == 'surf':
        phys_mat=   make_periodic(np.delete(np.delete(mat, 0, 0), -1, 0), \
                                        'surf', phi_period = phi_period)
    elif key    == 'umat':
        phys_mat=   make_periodic(np.delete(np.delete(mat, 0, 0), -1, 0), \
                                        'umat', sym_op = sym_op)
    #elif key    == 'normals':
    #    phys_mat     =   make_periodic(mat, 'surf', phi_period = phi_period)
    
    calc_mat    =   np.delete(np.delete(phys_mat, -1, 0), -1, 1)
        
    return mat, phys_mat, calc_mat
         
#    return surf_matrix

'''
class e_surf():
    
    def __init__(self, surf, rs, phis):
    
        form                =   surf.shape
        self.e_surf         =   np.empty((form[0] - 1, form[1] - 1) , dtype = 'object')
        
        for i in range(self.e_surf.shape[0]):
            for j in range(self.e_surf.shape[1]):
                self.e_surf[i,j]         =   np.zeros(3)
                del_vec                  =   [dr(rs, i)/2., dphi(phis, j)/2., 0]
    
                for k in range(3):
                    self.e_surf[i,j][k]  =   surf[i,j][k] + del_vec[k]
    
    def get_surf(self):            
        return self.e_surf
'''    

def deform(init_surf, u):
        
    def new_r(r): 
        #print u.u_pol(r)
        return np.array([r[0] + u.u(r, 'ur'), r[1] + u.u(r, 'uphi'), r[2] + u.u(r, 'uz')])
    
    deformed_surf = np.empty(init_surf.shape, dtype='object')
    
    for i in range(deformed_surf.shape[0]):
        for j in range(deformed_surf.shape[1]):
            deformed_surf[i,j] = new_r(init_surf[i,j])
    
    return deformed_surf 

def get_coord(surf):
    
    X   =   np.zeros(surf.shape)
    Y   =   np.zeros(surf.shape)
    Z   =   np.zeros(surf.shape)
    
    for i in range(surf.shape[0]):
        for j in range(surf.shape[1]):
            
            r       =   surf[i,j][0]
            phi     =   surf[i,j][1]
            z       =   surf[i,j][2]
            X[i,j]  =   r*np.cos(phi)
            Y[i,j]  =   r*np.sin(phi)
            Z[i,j]  =   z
    
    return X,Y,Z
