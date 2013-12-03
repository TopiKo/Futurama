'''
Created on Oct 3, 2013

@author: tohekorh
'''
import numpy as np
from help_classes import heaviside, curve_Lr, curve_Lphi, \
         valminmax, rps, phips, make_periodic, make_periodic_der, \
         select_option, get_z_set
from surface import transform_to_cart
from tests import test_prime_maps
from plots import plot_surf
#from stress_test.tests import path_integral_r
pi = np.pi

class u():

    def __init__(self, hangle, phi_period, asurf, system = 'spiral', consts = [], moldy = False): 
        
        #self.consts                     =   consts
        self.consts                     =   None
        self.hangle                     =   hangle
        self.opt, self.sym_op           =   select_option(system, phi_period, hangle) 
        self.phi_period                 =   phi_period
        self.system                     =   system
        self.rps                        =   None
        self.set_ainit_surf(asurf)
        
        #if consts != []:
        #    print 'hoo2,' + str(consts)
        #    self.set_const(consts)
        if moldy:
            self.set_moldy(moldy)     
                
            
    def set_moldy(self, moldy):
        
        if moldy:
            self.moldy                  =   True
            self.flat_u()
    
    def set_ainit_surf(self, ainit_surf):
        
        self.ext_surf                   =   ainit_surf[0]
        self.phys_surf                  =   ainit_surf[1]
        self.calc_surf                  =   ainit_surf[2]
            
    def set_const(self, consts):
        
        if self.consts != None:
            self.consts[:len(consts)]   =   consts        
        else:
            self.consts                 =   consts
        
        
        if self.system == 'spiral_w_wave' and self.rps == None:
            self.initialize_prime_maps()
        
        if self.system == 'spiral_w_wave':
            
            self.u()
            self.du()
            
            rdev, phidev            =   [1,1] 
            
            while not phidev < 1e-10 or not rdev < 1e-10: #1e-6 < phidev or 1e-6 < rdev:
                
                self.curve_Lphi     =   curve_Lphi(self, test = False)
                
                self.phips          =   phips(self)
                
                self.u()
                self.du()
                
                rdev                =   test_prime_maps(self)[0]
                
                self.curve_Lr       =   curve_Lr(self)
                self.rps            =   rps(self) 
                
                self.u()
                self.du()
                
                phidev              =   test_prime_maps(self)[1]
                
                #print rdev, phidev
        
        self.u()  
        self.du()
        self.ddu()        
        self.dA()
        
        '''
        if self.system == 'spiral_w_wave':
            store_opt               =   self.opt
            
            self.opt                =   [1, 0, self.opt[2]] 
            n                       =   0
            rdev, phidev            =   1, 1
            while not phidev < 1e-6 or not rdev < 1e-6: #1e-6 < phidev or 1e-6 < rdev:
                
                self.u()
                self.du()
                
                if n > 0:
                    phidev          =   test_prime_maps(self)[1]
                    #path_integral_r(self)
                    
                self.curve_Lphi     =   curve_Lphi(self, test = False)
                self.phips          =   phips(self)
                #plot_surf(self)
               
                if n == 0:
                    self.opt        =   [1, 3, self.opt[2]] 
                    
                self.u()
                self.du()
                
                if n > 0:
                    rdev    =   test_prime_maps(self)[0]
                    print phidev, rdev
                
                
                
                self.curve_Lr       =   curve_Lr(self)
                self.rps            =   rps(self) 
                
                
                n                  +=   1
                self.opt            =   store_opt  
        '''
            
                   
    
        
            
            
            
    '''
    def set_u(self, umat):
        
        self.umat                   =   np.empty(3, dtype = 'object')
        #self.umat_cart              =   np.empty(3, dtype = 'object')
        
        self.ext_umat                =   umat  
        #self.umat_cart[0]           =   np.empty(umat.shape, dtype = 'object')
        
        #for index, vec in np.ndenumerate(self.ext_umat):
            #self.umat_cart[0][index]    =   np.zeros(3)
            #ur, uphi, uz                =   vec
            #r, phi                      =   self.ext_surf[index][:2]
            #self.umat_cart[0][index][0] =   (r + ur) * np.cos(phi + uphi) - r * np.cos(phi) # ux  
            #self.umat_cart[0][index][1] =   (r + ur) * np.sin(phi + uphi) - r * np.sin(phi) # uy   
            #self.umat_cart[0][index][2] =   uz                                              # uz                              # uz  
        
        
        self.phys_umat                    =   np.delete(np.delete(self.ext_umat,      0, 0), 0, 1)
        #self.umat_cart[1]               =   np.delete(np.delete(self.umat_cart[0], 0, 0), 0, 1)
        
        self.calc_umat                    =   np.delete(np.delete(self.phys_umat,      -1, 0), -1, 1)  
        #self.umat_cart[2]               =   np.delete(np.delete(self.umat_cart[1], -1, 0), -1, 1)
        
        self.du()
        self.ddu()        
        self.dA()
        self.flat_u()
    '''   
    def initialize_prime_maps(self):
        
        store_opt               =   self.opt
            
        self.opt                =   [1, 0, self.opt[2]] 
        n                       =   0
        rdev, phidev            =   1, 1
        while not phidev < 1e-10 or not rdev < 1e-10: #1e-6 < phidev or 1e-6 < rdev:
            
            self.u()
            self.du()
            
            if n > 0:
                phidev          =   test_prime_maps(self)[1]
                #path_integral_r(self)
                
            self.curve_Lphi     =   curve_Lphi(self, test = False)
            self.phips          =   phips(self)
            #plot_surf(self)
           
            if n == 0:
                self.opt        =   [1, 3, self.opt[2]] 
            else:
                rdev    =   test_prime_maps(self)[0]
                #print phidev, rdev
                
           
            self.u()
            self.du()
            
            self.curve_Lr       =   curve_Lr(self)
            self.rps            =   rps(self) 
            
            
            n                  +=   1
            self.opt            =   store_opt  
    
    
    def parse_u(self, flat_u_pol):
        
        n       =   0
        
        for index, vec in np.ndenumerate( self.ext_umat ):
            for i in range(3):
                self.ext_umat[index][i]  =   flat_u_pol[n]
                n                      +=   1
            
        self.phys_umat   =   make_periodic(np.delete(np.delete(self.ext_umat,  0, 0), \
                                                       -1, 0), 'umat', sym_op = self.sym_op)
        self.calc_umat   =   np.delete(np.delete(self.phys_umat,      -1, 0), -1, 1)  

        self.du()
        self.ddu()        
        self.dA()
        
    def u(self):
        
        self.umat       =   np.empty(3, dtype = 'object')
        
        r_mat           =   self.ext_surf
        
        if self.system == 'spiral_w_wave':
            z_set       =   get_z_set(r_mat, self.hangle, self.consts, self.phi_period)
        
        
        #print z_set
        umat_b          =   np.empty(r_mat.shape, dtype = 'object')
        umat_cart_b     =   np.empty(r_mat.shape, dtype = 'object')
        
        [rmin, rmax]    =   valminmax(r_mat)[0]
        
        for index, r_vec in np.ndenumerate( r_mat ):
            
            umat_b[index]         =   np.zeros(3)
            umat_cart_b[index]    =   np.zeros(3)
            r, phi                =   r_vec[:2] 

            
            # Map r -> r'
            if self.opt[0] == 0:
                ur     =   0
            
            elif self.opt[0] == 1:
                ur     =   self.consts[0]
            
            elif self.opt[0] == 2:
                ur     =   self.rps[index] - r
                
            else:
                print 'MIta shiivaddia' 
            

            # Map phi -> phi'
            if self.opt[1] == 0:
                uphi   =   0
                                
            elif self.opt[1] == 1:
                uphi   =   heaviside(r, self.consts[2])*((r - self.consts[2])/(rmax - self.consts[2]))**2* \
                                            self.consts[3]*np.sin(phi/self.phi_period*2*pi)
            
            elif self.opt[1] == 2:
                uphi   =   self.consts[3]
            
            elif self.opt[1] == 3:
                uphi   =   self.phips[index] - phi
                
            else:
                print 'mita shiivaddia!'
            
            
            # Map z -> z'
            if self.opt[2] == 0:
                uz     =   self.hangle*phi   

            elif self.opt[2] == 1:
                #print self.consts
                #curve_start     =   rmin + self.consts[2] # + self.consts[0]
                #uz     =   self.hangle*phi  +   heaviside(r, curve_start)*((r - curve_start)/ \
                #                        (rmax - curve_start))**2 \
                #                        *self.consts[1]*np.sin(phi/self.phi_period*2*pi*self.consts[4])
                uz      =   z_set[index]
                
            elif self.opt[2] == 2:
                uz     =   self.consts[0]/2*r**2
            
            elif self.opt[2] == 3:
                uz     =   self.hangle*phi  +   self.consts[0]*(r - rmin)  
            
            elif self.opt[2] == 4:
                uz     =   self.consts[0]   -   np.sqrt(self.consts[0]**2 - r**2) 
            
            elif self.opt[2] == 5:
                uz      =   self.hangle*phi +   \
                            (r - rmin + self.consts[0])/(rmax - rmin) \
                            *self.consts[1]*np.sin(phi/self.phi_period*2*pi)
            
            elif self.opt[2] == 6:
                phase   =   0 #-pi/2 
                uz      =   self.hangle*phi  +   heaviside(r, self.consts[2])*((r - self.consts[2])/ \
                                        (rmax - self.consts[2]))**2 \
                                        *self.consts[1]*np.cos(phi/self.phi_period*2*pi + phase)
                                        
            
            else:
                print 'mita shiivaddia!'    
            
            umat_b[index][0]    =   ur
            umat_b[index][1]    =   uphi
            umat_b[index][2]    =   uz
            
        
        #self.ext_umat           =   make_periodic(umat_b, self.sym_op) 
        self.ext_umat           =   umat_b 
        self.phys_umat          =   make_periodic(np.delete(np.delete(self.ext_umat,  0, 0), -1, 0), 'umat', sym_op = self.sym_op)
        self.calc_umat          =   np.delete(np.delete(self.phys_umat,      -1, 0), -1, 1)  

    def get_umat(self):
        
        return [self.ext_umat, self.phys_umat, self.calc_umat]

    def du(self):
        
        
#        r_mat                   =   make_periodic_surf(self.ext_surf, self.phi_period)
        r_mat                   =   make_periodic(self.ext_surf, 'surf', phi_period = self.phi_period)
 
        u_mat_pol               =   make_periodic(self.ext_umat, 'umat', sym_op = self.sym_op)
          
        dr_mat                  =   np.delete(r_mat, 0, 0)  - np.delete(r_mat, -1, 0) 
        dphi_mat                =   np.delete(r_mat, 0, 1)  - np.delete(r_mat, -1, 1)   
        
        dr_mat_r                =   np.empty(dr_mat.shape)
        dphi_mat_phi            =   np.empty(dphi_mat.shape)
        
        
        for ir,dr in np.ndenumerate(dr_mat):
            dr_mat_r[ir]        =   dr[0]
        
        for iphi, dphi in np.ndenumerate(dphi_mat):   
            dphi_mat_phi[iphi]  =   dphi[1]
        
        
        dudrpol            =   (np.delete(u_mat_pol,  0, 0)   \
                               -np.delete(u_mat_pol, -1, 0))  / dr_mat_r 
        dudphipol          =   (np.delete(u_mat_pol,  0, 1)   \
                               -np.delete(u_mat_pol, -1, 1))  / dphi_mat_phi  
        
        self.dr_mat        =   dr_mat_r
        self.dphi_mat      =   make_periodic_der(dphi_mat_phi, 1, nums = True)
        
        self.dudrpol       =   dudrpol
        
        self.dudphipol     =   make_periodic_der(dudphipol)
        
        # testi
        dudrcart           =   np.empty(self.dr_mat.shape, dtype = 'object')
        dudphicart         =   np.empty(self.dphi_mat.shape, dtype = 'object')
            
        
        for ind, vec in np.ndenumerate(self.dr_mat):
            dudrcart[ind]      =   np.zeros(3)
            
            r, phi                  =   r_mat[ind][:2] #self.ext_surf[ind][:2]
            #rp, phip                =   r + self.ext_umat[ind][0], phi + self.ext_umat[ind][1]
            rp, phip                =   r + u_mat_pol[ind][0], phi + u_mat_pol[ind][1] 
            durdr, duphidr, duzdr   =   self.dudrpol[ind] 
            
            dudrcart[ind][0]   =  -np.cos(phi) - rp*duphidr*np.sin(phip) \
                                       +(1 + durdr)*np.cos(phip) 
            dudrcart[ind][1]   =  -np.sin(phi) + rp*duphidr*np.cos(phip) \
                                       +(1 + durdr)*np.sin(phip) 
            dudrcart[ind][2]   =   duzdr

        for ind, vec in np.ndenumerate(self.dphi_mat):
            dudphicart[ind]    =   np.zeros(3)
            
            r, phi                  =   r_mat[ind][:2] #self.ext_surf[ind][:2]
            #rp, phip                =   r + self.ext_umat[ind][0], phi + self.ext_umat[ind][1] 
            rp, phip                =   r + u_mat_pol[ind][0], phi + u_mat_pol[ind][1] 
            
            durdphi, duphidphi, duzdphi \
                                    =   self.dudphipol[ind] 
            
            dudphicart[ind][0] =   r*np.sin(phi) + np.cos(phip)*durdphi \
                                       -rp*(1 + duphidphi)*np.sin(phip)   
            dudphicart[ind][1] =  -r*np.cos(phi) + np.sin(phip)*durdphi \
                                       +rp*(1 + duphidphi)*np.cos(phip)   
            dudphicart[ind][2] =   duzdphi
        
        
        self.dudrcart           =   dudrcart
        self.dudphicart         =   dudphicart
        
    def get_du(self):
        
        return  np.delete(self.dudrpol, 0, 0),                       \
                np.delete(np.delete(self.dudphipol, 0, 0), -1, 0),  \
                np.delete(self.dudrcart, 0 , 0),                    \
                np.delete(np.delete(self.dudphicart, 0, 0), - 1, 0)
    
    def ddu(self):

        dudr_mat    =   self.dudrpol
        dudphi_mat  =   self.dudphipol  

        dr_mat_r    =   self.dr_mat
        dphi_mat_phi=   np.delete(self.dphi_mat, -1, 1)
        
        ddudrdr     =   (np.delete(dudr_mat,    0, 0)   \
                        -np.delete(dudr_mat,   -1, 0))  /  np.delete(dr_mat_r, -1 ,0)
        
        ddudrdphi   =   (np.delete(dudr_mat,    0, 1)   \
                        -np.delete(dudr_mat,   -1, 1))  /  np.delete(dphi_mat_phi, -1 ,0)
        
        
        ddudphidr   =   (np.delete(dudphi_mat,  0, 0)   \
                        -np.delete(dudphi_mat, -1, 0))  /  dr_mat_r
                         
        ddudphidphi =   (np.delete(dudphi_mat,  0, 1)   \
                        -np.delete(dudphi_mat, -1, 1))  /  dphi_mat_phi
        

        self.ddudrdr        =   ddudrdr
        self.ddudrdphi      =   np.delete(make_periodic_der(ddudrdphi, -1), 0, 0)
        self.ddudphidr      =   np.delete(ddudphidr, 0, 0)
        self.ddudphidphi    =   np.delete(np.delete(make_periodic_der(ddudphidphi, -1), 0, 0), -1, 0)
            
        self.ddudrdrcart    = np.empty(self.phys_surf.shape, dtype = 'object')
        self.ddudrdphicart  = np.empty(self.phys_surf.shape, dtype = 'object')
        self.ddudphidrcart  = np.empty(self.phys_surf.shape, dtype = 'object')
        self.ddudphidphicart= np.empty(self.phys_surf.shape, dtype = 'object')
        
        dudr_mat, dudphi_mat    =   self.get_du()[:2]
        
        for ir, rvec in np.ndenumerate(self.phys_surf):
            
            r, phi              =   rvec[:2]
            ur, uphi            =   self.phys_umat[ir][:2]    
            rp, phip            =   r + ur, phi + uphi
            durdr, duphidr      =   dudr_mat[ir][0], dudr_mat[ir][1]
            durdphi, duphidphi  =   dudphi_mat[ir][0], dudphi_mat[ir][1] #[:2] # self.get_du()[1][ir]
            
            ddurdrdr, dduphidrdr, dduzdrdr          =   self.ddudrdr[ir]
            ddurdphidr, dduphidphidr, dduzdphidr    =   self.ddudphidr[ir]
            ddurdrdphi, dduphidrdphi, dduzdrdphi    =   self.ddudrdphi[ir]
            ddurdphidphi, dduphidphidphi, dduzdphidphi  =   self.ddudphidphi[ir]
            
            
            self.ddudrdrcart[ir]        = np.zeros(3)
            self.ddudphidrcart[ir]      = np.zeros(3)
            self.ddudrdphicart[ir]      = np.zeros(3)
            self.ddudphidphicart[ir]    = np.zeros(3)
            
            
            self.ddudrdrcart[ir][0]     = ddurdrdr*np.cos(phip) - 2*(1 + durdr)*duphidr*np.sin(phip) \
                                            - rp*(duphidr**2*np.cos(phip) + dduphidrdr*np.sin(phip))
            
            self.ddudrdrcart[ir][1]     = ddurdrdr*np.sin(phip) + 2*(1 + durdr)*duphidr*np.cos(phip) \
                                            - rp*(duphidr**2*np.sin(phip) - dduphidrdr*np.cos(phip))
            
            self.ddudrdrcart[ir][2]     = dduzdrdr
            
            
            self.ddudphidrcart[ir][0]   = ddurdphidr*np.cos(phip) - (1 + durdr)*(1 + duphidphi) \
                                        * np.sin(phip) - durdphi*duphidr*np.sin(phip)   \
                                        - rp*((1 + duphidphi)*duphidr*np.cos(phip)      \
                                        + dduphidphidr*np.sin(phip)) + np.sin(phi) 
            
            self.ddudphidrcart[ir][1]   = ddurdphidr*np.sin(phip) + (1 + durdr)*(1 + duphidphi) \
                                        * np.cos(phip) + durdphi*duphidr*np.cos(phip)   \
                                        - rp*(duphidr*(1 + duphidphi)*np.sin(phip)      \
                                        - dduphidphidr*np.cos(phip)) - np.cos(phi)
            
            self.ddudphidrcart[ir][2]   = dduzdphidr
            
            
            self.ddudrdphicart[ir][0]   = ddurdrdphi*np.cos(phip) - (1 + durdr)*(1 + duphidphi) \
                                        * np.sin(phip) - durdphi*duphidr*np.sin(phip)   \
                                        - rp*((1 + duphidphi)*duphidr*np.cos(phip)      \
                                        + dduphidrdphi*np.sin(phip)) + np.sin(phi) 
            self.ddudrdphicart[ir][1]   = ddurdrdphi*np.sin(phip) + (1 + durdr)*(1 + duphidphi) \
                                        * np.cos(phip) + durdphi*duphidr*np.cos(phip)   \
                                        - rp*(duphidr*(1 + duphidphi)*np.sin(phip)      \
                                        - dduphidrdphi*np.cos(phip)) - np.cos(phi)
            self.ddudrdphicart[ir][2]   = dduzdrdphi
            
            
            
            self.ddudphidphicart[ir][0] = ddurdphidphi*np.cos(phip) - durdphi*(1 + duphidphi)   \
                                        * np.sin(phip) - durdphi*(1 + duphidphi)*np.sin(phip)   \
                                        - rp*(dduphidphidphi*np.sin(phip)                       \
                                        + (1 + duphidphi)**2*np.cos(phip)) + r*np.cos(phi) 
            self.ddudphidphicart[ir][1] = ddurdphidphi*np.sin(phip) + durdphi*(1 + duphidphi)   \
                                        * np.cos(phip) + durdphi*(1 + duphidphi)*np.cos(phip)   \
                                        + rp*(dduphidphidphi*np.cos(phip)                       \
                                        - (1 + duphidphi)**2*np.sin(phip)) + r*np.sin(phi)
            self.ddudphidphicart[ir][2] = dduzdphidphi
    
    def get_ddu(self):
        
        return self.ddudrdrcart, self.ddudphidrcart, self.ddudrdphicart, self.ddudphidphicart
        
    def dA(self):
        
        self.dA_S               =   np.zeros(self.calc_surf.shape)
        self.dA_P               =   np.zeros(self.calc_surf.shape)
        
        init_surf               =   transform_to_cart(self.phys_surf)
        deformed_surf           =   transform_to_cart(self.phys_surf + self.phys_umat)
        
        for index, val in np.ndenumerate( self.dA_S ):
            i,j = index

            dr_vec_p            =   init_surf[i + 1, j]        -   init_surf[i, j]
            dphi_vec_p1         =   init_surf[i, j + 1]        -   init_surf[i, j]
            dphi_vec_p2         =   init_surf[i + 1, j + 1]    -   init_surf[i + 1, j]
            
            v1p                 =   dr_vec_p                   +   dphi_vec_p1
            v2p                 =   dphi_vec_p2                -   dr_vec_p 
            
            self.dA_P[index]    =   .5*np.linalg.norm(np.cross(v1p, v2p))
            
            dr_vec              =   deformed_surf[i + 1, j]    - deformed_surf[i, j]
            dphi_vec1           =   deformed_surf[i, j + 1]    - deformed_surf[i, j]
            dphi_vec2           =   deformed_surf[i + 1, j + 1]- deformed_surf[i + 1, j]

            v1s                 =   dr_vec        +   dphi_vec1
            v2s                 =   dphi_vec2     -   dr_vec 

            self.dA_S[index]    =   .5*np.linalg.norm(np.cross(v1s, v2s))
            
    def flat_u(self):
        
        #tmp_umat                =   self.ext_umat  
        n, m                    =   self.ext_umat.shape   #tmp_umat.shape    
        self.flat_umat          =   np.zeros(n*m*3)
        self.bounds             =   np.empty(n*m*3, dtype = 'object')
        
        n                       =   0
        dev                     =   [1.0, 0.1, 5.0]
        for index, vec in np.ndenumerate(self.ext_umat):
            j    =   index[1]
            for k, val in enumerate(vec):
                self.flat_umat[3*n + k]         =   val
                
                if j == 0:
                    if k == 1:
                        self.bounds[3*n + k]    =   (val, val)
                    elif k == 2:
                        self.bounds[3*n + k]    =   (val, val)
                    else:
                        self.bounds[3*n + k]    =   (val - dev[k], val + dev[k])   
                else:
                    #self.bounds[3*n + k]        =   (val - dev[k], val + dev[k])
                    self.bounds[3*n + k]        =   (None, None)
                
            n                  +=  1
        
    def separate_mat(self, mat):
        
        mat_r      =   np.zeros(mat.shape)
        mat_phi    =   np.zeros(mat.shape)
        mat_z      =   np.zeros(mat.shape)
        
        for ind, rvec in np.ndenumerate(mat):
            mat_r[ind]      =   rvec[0]
            mat_phi[ind]    =   rvec[1]
            mat_z[ind]      =   rvec[2]
        
        return mat_r, mat_phi, mat_z
    
def parse_u_from_file(in_file):  
    
    from bender_rw import parse_input
    from surface import surf
    
    param_set           =   parse_input(in_file)
    system              =   param_set["system"]
    nr, nphi            =   param_set["nr"],    param_set["nphi"]
    rmin, rmax          =   param_set["rmin"],  param_set["rmax"]
    phimin, phi_period  =   param_set["phimin"],param_set["phiperiod"] 
    height              =   param_set["height"]
    moldy_opm           =   param_set["moldy_opm"]
    
    hangle          =   height / 2 / pi   
    asurf           =   surf(rmin, rmax, nr, phimin, phi_period, nphi)
    ue              =   u(hangle, phi_period, asurf.get_all_surf(), \
                          system = system, moldy = moldy_opm)
    
    return ue       
   
    