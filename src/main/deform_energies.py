'''
Created on Oct 3, 2013

@author: tohekorh
'''
import numpy as np
from plots import plot_all

class deform_energies():
    
    def __init__(self, u = None, bend_module = 1.6, strech_module = 25, sigma = 0.28):
        self.bm     =   bend_module
        self.sm     =   strech_module
        self.sigma  =   sigma
        self.u      =   None
        self.delg   =   [1e-8, 1e-8, 1e-8, 1e-8]
        if u != None:
            self.set_u(u)
 
    def calc_energies(self, consts = []):
        
        if len(consts) != 0:
            self.u.set_const(consts)
        
        E_s_surf, E_s               =   self.F_s()
        E_b_surf, E_b, normals      =   self.F_b()
            
        return E_b, E_s, E_b_surf, E_s_surf, normals 
    
    def calc_grad(self, consts, E_b, E_s, acc = 1e-8):
        
        grad                        =   np.zeros(len(consts))
        
        def add_del(index):
            new_consts  = np.zeros(len(consts))
            for ic, c in enumerate(consts):
                if ic == index:
                    new_consts[ic]   =  c + acc #self.delg[ic]
                else:
                    new_consts[ic]   =  c
            return new_consts 
        
        for i in range(len(grad)):
            self.u.set_const(add_del(i))
            grad[i]     =   (self.F_s()[1] + self.F_b()[1] - E_b - E_s) / acc #self.delg[i]
            
        return grad 
        
    def calc_energies_moldy(self, flat_umat, count = 0):
        
        self.u.parse_u(flat_umat)
            
        E_s_surf, E_s           =   self.F_s()
        E_b_surf, E_b, normals  =   self.F_b()
        
        if count % 1000 == 0 or count == 0: 
            plot_all(self.u, E_s_surf, E_b_surf, normals, show = False, \
                 save = True, path = \
                 '/space/tohekorh/Spiral/iteration/iter=%i_' %count ) 
        
        return E_b, E_s, E_b_surf, E_s_surf, normals
        
        
    def set_u(self, u):
        
        self.u  = u
        
        
    def F_s(self):
        
        calc_grid   =   self.u.calc_surf
        
        F_surf      =   np.zeros(calc_grid.shape)
        isum        =   0
        
        dudrcart, dudphicart    =   self.u.get_du()[2:]
        
        def pol_coords_der(i,j):
            
            r, phi  =   self.u.phys_surf[i, j][:2] #d_surf[i,j][:2] 
            r1, phi1=   self.u.phys_surf[i + 1, j + 1][:2] 
            
            r_c     =   (r + r1)/2.
            phi_c   =   (phi + phi1)/2.
            
            drdx    =   np.cos(phi_c)
            drdy    =   np.sin(phi_c)
            drdz    =   0
            
            dphidx  =   -np.sin(phi_c)/r_c
            dphidy  =   np.cos(phi_c)/r_c 
            dphidz  =   0
            
            dzdx    =   0
            dzdy    =   0
            dzdz    =   1
            
            return np.array([[drdx, drdy, drdz],[dphidx, dphidy, dphidz], [dzdx, dzdy, dzdz]]) 

        
        for i in range(len(F_surf)):
            for j in range(len(F_surf[0])):
                # pol_ders  =   [[drdx, drdy, drdz],[dphidx, dphidy, dphidz], [dzdx, dzdy, dzdz]]
                # u_ders    =   [[duxdr, duydr, duzdr], [duxdphi, duydphi, duzdphi]]
                pol_ders    =   pol_coords_der(i, j) 
                u_ders_x    =   np.array([(dudrcart[i, j + 1][0] +   dudrcart[i,j][0])/2. \
                                         ,(dudphicart[i + 1,j][0]+   dudphicart[i,j][0])/2.])
                u_ders_y    =   np.array([(dudrcart[i, j + 1][1] +   dudrcart[i,j][1])/2. \
                                         ,(dudphicart[i + 1,j][1]+   dudphicart[i,j][1])/2.])
                u_ders_z    =   np.array([(dudrcart[i, j + 1][2] +   dudrcart[i,j][2])/2. \
                                         ,(dudphicart[i + 1,j][2]+   dudphicart[i,j][2])/2.])
                
                
                uxx         =   np.dot(u_ders_x, pol_ders[:,0][:2])     \
                            + .5*np.dot(u_ders_x, pol_ders[:,0][:2])**2 \
                            + .5*np.dot(u_ders_y, pol_ders[:,0][:2])**2 \
                            + .5*np.dot(u_ders_z, pol_ders[:,0][:2])**2 
                
                uyy         =   np.dot(u_ders_y, pol_ders[:,1][:2])     \
                            + .5*np.dot(u_ders_y, pol_ders[:,1][:2])**2 \
                            + .5*np.dot(u_ders_x, pol_ders[:,1][:2])**2 \
                            + .5*np.dot(u_ders_z, pol_ders[:,1][:2])**2
                
                uxy         = .5*np.dot(u_ders_y, pol_ders[:,0][:2])    \
                            + .5*np.dot(u_ders_x, pol_ders[:,1][:2])    \
                            + .5*np.dot(u_ders_x, pol_ders[:,1][:2])    \
                            * np.dot(u_ders_x, pol_ders[:,0][:2])       \
                            + .5*np.dot(u_ders_y, pol_ders[:,1][:2])    \
                            * np.dot(u_ders_y, pol_ders[:,0][:2])       \
                            + .5*np.dot(u_ders_z, pol_ders[:,0][:2])    \
                            * np.dot(u_ders_z, pol_ders[:,1][:2])             
                
                F_surf[i,j] =   self.sm*(1./(2*(1-self.sigma))*(uxx + uyy)**2 - (uxx*uyy - uxy**2)) 
                
                isum       +=   F_surf[i,j]*self.u.dA_P[i,j] 
                
        return F_surf, isum
    
    def F_b(self):
        
        
        calc_grid   =   self.u.calc_surf

        dudrcart, dudphicart    =   self.u.get_du()[2:]

        
        F_b_surf    =   np.zeros(calc_grid.shape)
        normals     =   np.empty(calc_grid.shape, dtype='object')
        isum        =   0
        
        
        def S_i(i,j): 
            
            
            r, phi  =   self.u.phys_surf[i, j][:2] 
            r1, phi1=   self.u.phys_surf[i + 1, j + 1][:2] 
            
            r_c     =   (r + r1)/2.
            phi_c   =   (phi + phi1)/2.
            
            duxdr   =   (dudrcart[i, j + 1][0]   +   dudrcart[i, j][0])/2.
            duydr   =   (dudrcart[i, j + 1][1]   +   dudrcart[i, j][1])/2.
            duzdr   =   (dudrcart[i, j + 1][2]   +   dudrcart[i, j][2])/2.      
            duxdphi =   (dudphicart[i + 1, j][0] +   dudphicart[i, j][0])/2.
            duydphi =   (dudphicart[i + 1, j][1] +   dudphicart[i, j][1])/2.
            duzdphi =   (dudphicart[i + 1, j][2] +   dudphicart[i, j][2])/2.
            
            s1      =   np.array([np.cos(phi_c) + duxdr, np.sin(phi_c) + duydr, duzdr])
            s2      =   np.array([-r_c*np.sin(phi_c) + duxdphi, r*np.cos(phi_c) + duydphi, duzdphi])
            norm_s1 =   np.linalg.norm(s1)
            norm_s2 =   np.linalg.norm(s2)
            
            return s1, s2, norm_s1, norm_s2 
    
        def D_jS_i(i,j):
            
            #ddudr  =   dduxdrdr, dduydrdr, dduzdrdr, dduxdrdphi, dduydrdphi, dduzdrdphi
            #ddudphi=   dduxdphidr, dduydphidr, dduzdphidr, dduxdphidphi, dduydphidphi, dduzdphidphi
            
            u_ddr_mat, u_drdphi_mat, \
            u_dphidr_mat, u_ddphi_mat       =   self.u.get_ddu() 
            
            # use the average of second derivative in phi direction
            ddudr   =   np.array([(u_ddr_mat[i, j + 1]      +   u_ddr_mat[i, j])/2, \
                                  (u_drdphi_mat[i + 1, j]   +   u_drdphi_mat[i, j])/2])
            ddudphi =   np.array([(u_dphidr_mat[i, j + 1]   +   u_dphidr_mat[i, j])/2, \
                                  (u_ddphi_mat[i + 1, j]    +   u_ddphi_mat[i, j])/2])
            
            r, phi  =   self.u.phys_surf[i, j][:2] 
            r1, phi1=   self.u.phys_surf[i + 1, j + 1][:2] 
            
            r_c     =   (r + r1)/2.
            phi_c   =   (phi + phi1)/2.
            
            s1dr    =   np.array([ddudr[0,0], ddudr[0,1], ddudr[0,2]])
            s1dphi  =   np.array([-np.sin(phi_c) + ddudphi[0,0], np.cos(phi_c) + ddudphi[0,1], ddudphi[0,2]])
            s2dr    =   np.array([-np.sin(phi_c) + ddudr[1,0], np.cos(phi_c) + ddudr[1,1], ddudr[1,2]])
            s2dphi  =   np.array([-r_c*np.cos(phi_c) + ddudphi[1,0], -r_c*np.sin(phi_c) + ddudphi[1,1], ddudphi[1,2]])
            
            return s1dr, s1dphi, s2dr, s2dphi 
        
        for i in range(F_b_surf.shape[0]):
            for j in range(F_b_surf.shape[1]):                
                
                s1, s2, \
                ns1, ns2        =   S_i(i,j)
                s1dr, s1dphi, \
                s2dr, s2dphi    =   D_jS_i(i,j) 
                
                n               =   np.cross(s1, s2)/np.linalg.norm(np.cross(s1, s2))
                
                f1              =   s1/ns1
                f2              =   (s2 - np.dot(s2, s1) / ns1**2 * s1) \
                                  / (np.linalg.norm(s2 - np.dot(s2, s1) / ns1**2 * s1))
                
                if abs(np.dot(f1, f2)) > 1e-10:
                    print 'voi rahma!' 
                
                a               =   np.dot(s1, s2) / ns1**2
                b               =   np.linalg.norm(s2 - np.dot(s2, s1) / ns1**2 * s1)  \
                                                / (ns2**2 - np.dot(s1, s2) / ns1**2) 
                
                normals[i,j]    =   n
                
                ddkhi           =   np.zeros(3)
                ddkhi[0]        =   np.dot(n, s1dr) / ns1**2
                ddkhi[1]        =   np.dot(n, (s1dr - s1dphi - s2dr)*a*b**2 + s2dphi*b**2)
                ddkhi[2]        =   np.dot(n, -s1dr*a*b/ns1 + .5*(s1dphi + s2dr)*b / ns1)
                
                            
                F_b_surf[i,j]   =   self.bm*(1./(2*(1-self.sigma))*(ddkhi[0] + ddkhi[1])**2  \
                                                         + (ddkhi[2]**2 - ddkhi[0]*ddkhi[1]))
                
                
                isum           +=   F_b_surf[i,j]*self.u.dA_S[i,j] 
                
        return F_b_surf, isum, normals
    
    
    
    