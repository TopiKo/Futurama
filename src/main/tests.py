'''
Created on Oct 3, 2013

@author: tohekorh
'''
import numpy as np
from sys import exit
from surface import surf #, e_surf
from help_classes import valminmax, make_periodic #rs, phis,
from deform_energies import deform_energies
from plots import plot_all, plot_surf
from surface import transform_to_cart

        

pi = np.pi

def check(u, e_surf):
    
    screwed_up  =   False 
    
    umat_cart               =   transform_to_cart(u.phys_surf + u.phys_umat) \
                              - transform_to_cart(u.phys_surf)
    
    u_pol, u_cart           =   u.phys_umat, umat_cart #u.umat_cart[1] # u.u(e_surf[i,j], 'upol'), u.u(e_surf[i,j], 'ucart')
    dudr_pol, dudphi_pol, \
    dudr_cart, dudphi_cart =   u.get_du()  #u.du(e_surf[i,j], 'dupol'), u.du(e_surf[i,j], 'ducart')

    
    for ind, r in np.ndenumerate(e_surf):
             
        r, phi                  =   r[:2] #e_surf[i,j][:2]
        
        if abs(dudr_cart[ind][0] - ( dudr_pol[ind][0]*np.cos(phi) - u_pol[ind][1]*np.sin(phi) - r*dudr_pol[ind][1]*np.sin(phi))) > 0.01:
            print 'Alert1'
            #print r, phi, u_der_cart[0,0] - ( u_der_pol[0,0]*np.cos(phi) - u_pol[1]*np.sin(phi) - r*u_der_pol[1,0]*np.sin(phi))
            screwed_up = True
        if abs(dudr_cart[ind][1] - ( dudr_pol[ind][0]*np.sin(phi) + u_pol[ind][1]*np.cos(phi) + r*dudr_pol[ind][1]*np.cos(phi))) > 0.01:
            print 'Alert2'
            #print r, phi
            screwed_up = True               
        if abs(dudphi_cart[ind][0] - ( dudphi_pol[ind][0]*np.cos(phi) - r*dudphi_pol[ind][1]*np.sin(phi) - u_cart[ind][1])) > 0.01:
            print 'Alert3'
            #print r, phi
            screwed_up = True
        if abs(dudphi_cart[ind][1] - ( dudphi_pol[ind][0]*np.sin(phi) + r*dudphi_pol[ind][1]*np.cos(phi) + u_cart[ind][0])) > 0.01:
            print 'Alert4'
            #print r, phi
            screwed_up = True
            
            
            '''
            if abs(u_der_cart[0,0] - ( u_der_pol[0,0]*np.cos(phi) - u_pol[1]*np.sin(phi) - r*u_der_pol[0,1]*np.sin(phi))) > 0.01:
                print 'Alert1'
                #print r, phi, u_der_cart[0,0] - ( u_der_pol[0,0]*np.cos(phi) - u_pol[1]*np.sin(phi) - r*u_der_pol[1,0]*np.sin(phi))
                screwed_up = True
            if abs(u_der_cart[0,1] - ( u_der_pol[0,0]*np.sin(phi) + u_pol[1]*np.cos(phi) + r*u_der_pol[0,1]*np.cos(phi))) > 0.01:
                print 'Alert2'
                #print r, phi
                screwed_up = True               
            if abs(u_der_cart[1,0] - ( u_der_pol[1,0]*np.cos(phi) - r*u_der_pol[1,1]*np.sin(phi) - u_cart[1])) > 0.01:
                print 'Alert3'
                #print r, phi
                screwed_up = True
            if abs(u_der_cart[1,1] - ( u_der_pol[1,0]*np.sin(phi) + r*u_der_pol[1,1]*np.cos(phi) + u_cart[0])) > 0.01:
                print 'Alert4'
                #print r, phi
                screwed_up = True
            '''
    if not screwed_up:
        print 'Derivatives are ok!'
    else:
        print 'derivatives screwqed up...'
        exit(0)

class tests():
    
    
    
    def __init__(self, nr, nphi, bend_module = 1.6, strech_module = 25, sigma = 0.28):
        
        
        self.bm     =   bend_module
        self.sm     =   strech_module
        self.sigma  =   sigma
        self.nr     =   nr
        self.nphi   =   nphi
    '''
    def run(self):
        self.spiral_test()
        self.ball_test()
        self.cup_test()
        self.screw_test()
        self.compare_with_previous()
    
    def cup_test(self, param_set):
        # CUP BEGINS
        from strain import u
        
        system              =   param_set["system"]
        nr, nphi            =   param_set["nr"], param_set["nphi"]
        rmin, rmax          =   param_set["rmin"], param_set["rmax"]
        phimin, phi_period  =   param_set["phimin"], param_set["phiperiod"] 
        height              =   param_set["height"]
        consts              =   param_set["consts"]
        hangle              =   height / 2 / pi   
    
        
        
        alpha       =   consts[0]
        
        init_surf_cup                   =   surf(rmin, rmax, nr, phimin, phi_period, nphi)
        u_cup                           =   u(0, phi_period, init_surf_cup.get_all_surf(), \
                                              system = system)

        energies_cup                    =   deform_energies(u = u_cup, bend_module = self.bm,\
                                                    strech_module = self.sm, sigma = self.sigma)
        
        
        # calculate energies
        #E_s_surf, E_s                   =   energies_cup.F_s()
        #E_b_surf, E_b, normals          =   energies_cup.F_b()
        
        E_b, E_s, E_b_surf, E_s_surf, normals \
                        =   energies_cup.calc_energies(consts)
        # calc ready
        
        
        # print resutls
        print '\nThis is cup with uz = alpha/2*R**2'
        print 'code E_b                      = ' + str(E_b) 
        print 'analytical E_b (small alpha!) = ' + str(1./2.*phi_period*self.bm*alpha**2*(rmax**2 - rmin**2)*(2./(1.-self.sigma) - 1.)) # kuppi, bend
        print 'code E_s                      = ' + str(E_s) # kuppi, bend
        print 'analytical E_s (small alpha!) = ' + str(self.sm/8.*alpha**4*phi_period/(1.-self.sigma)*1./6.*(rmax**6 - rmin**6))
        # print ready
        
        plot_all(u_cup, E_s_surf, E_b_surf, normals) 
        
        # CUP ENDS
        
    def spiral_test(self, hangle = 1./2/pi, rmin = 5, rmax = 10):
        # Spiral BEGINS
        from strain import u

        phimin      =   0      
        phi_period  =   pi
        
        init_surf_spiral                    =   surf(rmin, rmax, self.nr, phimin, phi_period, self.nphi) #surf(rs_spi, phis_spi) #initialize_surf()
        u_spiral                            =   u(hangle, phi_period, init_surf_spiral.get_all_surf(), \
                                                    consts = [0], system = 'spiral')
        
        energies_spi                        =   deform_energies(u = u_spiral, bend_module = self.bm,\
                                                        strech_module = self.sm, sigma = self.sigma)
        
        E_s_surf, E_s                       =   energies_spi.F_s()
        E_b_surf, E_b, normals              =   energies_spi.F_b()
        
        # print resutls
        print '\nThis is spiral with uz = hangle*phi'
        print 'code E_b                       = ' + str(E_b) 
        print 'analytical E_b (small hangle!) = ' + str(self.bm*hangle**2*phi_period/2. \
                                                       *(1./rmin**2 - 1./rmax**2)) 
        print 'code E_s                       = ' + str(E_s)
        print 'analytical E_s (small hangle!) = ' + str(self.sm*hangle**4/(2*(1 - self.sigma))*phi_period/8. \
                                                       *(1./rmin**2 - 1./rmax**2))
        plot_all(u_spiral, E_s_surf, E_b_surf, normals) 
        # print ready
        # Spiral ENDS
    
    def screw_test(self):
        from strain import u
        # screw BEGINS
        H       =   4.0
        alpha   =   0.2
        
        
        rmin    =   5
        rmax    =   10
        phimin  =   0
        phi_period  = pi
        
        hangle  =   H/(2*pi)
      
        init_surf_screw     =   surf(rmin, rmax, self.nr, phimin, phi_period, self.nphi) 
        
        u_screw             =   u(hangle, phi_period, init_surf_screw.get_all_surf(), \
                                  system = 'screw',consts = [alpha])
       
        energies_screw      =   deform_energies(u = u_screw, bend_module = self.bm,\
                                     strech_module = self.sm, sigma = self.sigma)
        
        # calculate energies
        #energies_screw.set_u(u_screw)
        E_s_surf, E_s                       =   energies_screw.F_s()
        E_b_surf, E_b, normals              =   energies_screw.F_b()
        # calc ready

        # print resutls
        print '\nThis is screw with uz = hangle*phi + r*alpha'
        print 'code E_b                       = ' + str(E_b) 
        print 'analytical E_b (small hangle!) = ' + str(self.bm*phi_period*(hangle**2/2*(1./rmin**2 - 1./rmax**2) \
                                                            + 1./(2*(1 - self.sigma))*alpha**2*np.log(rmax/rmin)))
        print 'code E_s                       = ' + str(E_s)
        print 'analytical E_s (small hangle!) = ' + str(self.sm/8.*phi_period/(1 - self.sigma)*( \
                                                          alpha**4 / 2. * ( rmax**2  - rmin**2 )                    \
                                                        + 2 * alpha**2 * hangle**2 * np.log( rmax / rmin )    \
                                                        + hangle**4 / 2 * ( 1. / rmin**2 - 1. / rmax**2) ))
        plot_all(u_screw, E_s_surf, E_b_surf, normals) 
        # print ready
        # Screw ENDS

    def ball_test(self):
        from strain import u
        # Ball BEGINS
        R           =   10.0
        alpha       =   np.sin(pi/10)
        theta_max   =   np.arcsin(alpha)
        
        rmin        =   1
        rmax        =   alpha*R
        phimin      =   0
        phi_period  =   pi/2.
        
        init_surf_ball          =   surf(rmin, rmax, self.nr, phimin, phi_period, self.nphi) 
        u_ball                  =   u(0, phi_period, init_surf_ball.get_all_surf(), \
                                      system = 'ball', consts = [R, theta_max])
      
        energies_ball           =   deform_energies(u = u_ball, bend_module = self.bm,\
                                        strech_module = self.sm, sigma = self.sigma)
        
        
        # calculate energies
        E_s_surf, E_s           =   energies_ball.F_s()
        E_b_surf, E_b, normals  =   energies_ball.F_b()
        # calc ready
        
        alpha_max                           =   rmax/R
        alpha_min                           =   rmin/R
        theta_max                           =   np.arcsin(rmax/R)
        theta_min                           =   np.arcsin(rmin/R)
        
        # print resutls
        print '\nThis is ball with uz = R - sqrt(R**2 - r**2), from 0 to alpha*R'
        print 'code E_b_density             = ' + str(E_b/(phi_period*R**2*(np.cos(theta_min) - np.cos(theta_max)))) 
        print 'analytical E_b_density       = ' + str(self.bm / R**2 * ( 2. / (1 - self.sigma) - 1 ))
        print 'code E_s                     = ' + str(E_s)
        print 'analytical E_s               = ' + str(self.sm*R**2/(8*(1 - self.sigma))*phi_period* \
                                            (np.log(1 - alpha_max**2) + alpha_max**2/2*(alpha_max**2 - 2)/(alpha_max**2 - 1) \
                                           - np.log(1 - alpha_min**2) - alpha_min**2/2*(alpha_min**2 - 2)/(alpha_min**2 - 1)))
        
        plot_surf(u_ball)
        plot_all(u_ball, E_s_surf, E_b_surf, normals)
        # ball ENDS

    def compare_with_previous(self, rmin = 10, rmax = 20, hangle = 2./2./pi):
        from strain import u

        from calculate.energies import ebend, estrech, find_X
        
        # set sigma == 0!        
        print 'we compare with previous results: is sigma = 0? and hangle << r_min?'
        self.sigma  = 0
        
        phimin      =   0
        phi_period  =   pi/3.
        
        # make quess
        x                   =   np.sqrt(rmin**2 - hangle**2) - rmin
        quess               =   [x]
        bounds              =   [(-rmin, 3)]
        # quess ready

        
        # initialize surfaces
        init_surf_comp      =   surf(rmin, rmax, self.nr, phimin, phi_period, self.nphi) 
        u_comp              =   u(hangle, phi_period, init_surf_comp.get_all_surf(), \
                                  system = 'spiral')

        # init ready
        
        energies_comp       =   deform_energies(u = u_comp, bend_module = self.bm,\
                                   strech_module = self.sm, sigma = self.sigma)
        
        
        def optimize(u_comp, ques_consts, bounds):

            from scipy.optimize import fmin, fmin_tnc, fmin_l_bfgs_b  # bfgs, lbfgs optimointiin
            
            print 'optimizing..'
            
            def energy(consts):
                
                u_comp.set_const(consts)
                energies_comp.set_u(u_comp)
                
                E_s         =   energies_comp.F_s()[1]
                E_b         =   energies_comp.F_b()[1]
                
                print consts, E_s + E_b        
                
                return E_s + E_b
    
            
            #opm_consts  =   fmin_tnc(energy, ques_consts, bounds = bounds, approx_grad = True)[0]
            opm_consts  =   fmin_l_bfgs_b(energy, ques_consts, bounds = bounds, approx_grad = True)[0]
            
            u_comp.set_const(opm_consts)
            energies_comp.set_u(u_comp)
            #opm_consts  =   fmin(energy, ques_consts)
            
            return opm_consts
                
        opm_consts                  =   optimize(u_comp, quess, bounds)
        
        # check derivatives
        check(u_comp, u_comp.phys_surf) # init_surf_comp.get_esurf())
        # chec ready
        
        E_s_surf, E_s               =   energies_comp.F_s()
        E_b_surf, E_b, normals      =   energies_comp.F_b()
        
        x                           =   find_X(rmin, rmax - rmin, hangle)
        
        print 'x = ' + str(x)
        print 'opm _const (x) = ' + str(opm_consts[0])
        
        E2      =   ebend(hangle, rmin, rmax - rmin, x)
        E3      =   estrech(hangle, rmin, rmax - rmin, x)[0]
        
        print 'remember E_old = E_new, only for hangle << r_min and sigma ==0 !!!!'
        print 'E_old; E_b, E_s, E_b + E_s = ' + str([E2, E3, E2 + E3])
        print 'E_new; E_b, E_s, E_b + E_s = ' + str([E_b, E_s, E_b + E_s])
        plot_all(u_comp, E_s_surf, E_b_surf, normals)

    '''
    def tests(self, param_set):
        
        from strain import u
        
        system              =   param_set["system"]
        nr, nphi            =   param_set["nr"], param_set["nphi"]
        rmin, rmax          =   param_set["rmin"], param_set["rmax"]
        phimin, phi_period  =   param_set["phimin"], param_set["phiperiod"] 
        height              =   param_set["height"]
        consts              =   param_set["consts"]
        hangle              =   height / 2 / pi   
    
        init_surf           =   surf(rmin, rmax, nr, phimin, phi_period, nphi)
        
        ue                  =   u(hangle, phi_period, init_surf.get_all_surf(), \
                                              system = system)
        
        energies            =   deform_energies(u = ue, bend_module = self.bm,\
                                    strech_module = self.sm, sigma = self.sigma)
        
        
        E_b, E_s, E_b_surf, E_s_surf, normals \
                            =   energies.calc_energies(consts)
        alpha               =   consts[0]
            
        
        if system == 'cup':
            
            print '\nThis is cup with uz = alpha/2*R**2'
            print 'code E_b                      = ' + str(E_b) 
            print 'analytical E_b (small alpha!) = ' + str(1./2.*phi_period*self.bm*alpha**2*(rmax**2 - rmin**2)*(2./(1.-self.sigma) - 1.)) # kuppi, bend
            print 'code E_s                      = ' + str(E_s) 
            print 'analytical E_s (small alpha!) = ' + str(self.sm/8.*alpha**4*phi_period/(1.-self.sigma)*1./6.*(rmax**6 - rmin**6))
        
        elif system == 'spiral':
        
            print '\nThis is spiral with uz = hangle*phi'
            print 'code E_b                      = ' + str(E_b) 
            print 'analytical E_b (small hangle!)= ' + str(self.bm*hangle**2*phi_period/2. \
                                                           *(1./rmin**2 - 1./rmax**2)) 
            print 'code E_s                      = ' + str(E_s)
            print 'analytical E_s (small hangle!)= ' + str(self.sm*hangle**4/(2*(1 - self.sigma))*phi_period/8. \
                                                           *(1./rmin**2 - 1./rmax**2))
        
        elif system == 'screw':
            
            print '\nThis is screw with uz = hangle*phi + r*alpha'
            print 'code E_b                      = ' + str(E_b) 
            print 'analytical E_b (small hangle!)= ' + str(self.bm*phi_period*(hangle**2/2*(1./rmin**2 - 1./rmax**2) \
                                                                + 1./(2*(1 - self.sigma))*alpha**2*np.log(rmax/rmin)))
            print 'code E_s                      = ' + str(E_s)
            print 'analytical E_s (small hangle!)= ' + str(self.sm/8.*phi_period/(1 - self.sigma)*( \
                                                              alpha**4 / 2. * ( rmax**2  - rmin**2 )                    \
                                                            + 2 * alpha**2 * hangle**2 * np.log( rmax / rmin )    \
                                                            + hangle**4 / 2 * ( 1. / rmin**2 - 1. / rmax**2) ))
        
        elif system == 'ball':
            
            '''
            R          =   10.0
            alpha      =   np.sin(pi/10)
            theta_max  =   np.arcsin(alpha)
            rmax       =   alpha*R
            consts     =   [R, theta_max]
            '''
            
            R                   =   consts[0]
            alpha_max           =   rmax/R
            alpha_min           =   rmin/R
            theta_max           =   np.arcsin(rmax/R)
            theta_min           =   np.arcsin(rmin/R)
            
            # print resutls
            print '\nThis is ball with uz = R - sqrt(R**2 - r**2), from 0 to alpha*R'
            print 'code E_b_density              = ' + str(E_b/(phi_period*R**2*(np.cos(theta_min) - np.cos(theta_max)))) 
            print 'analytical E_b_density        = ' + str(self.bm / R**2 * ( 2. / (1 - self.sigma) - 1 ))
            print 'code E_s                      = ' + str(E_s)
            print 'analytical E_s                = ' + str(self.sm*R**2/(8*(1 - self.sigma))*phi_period* \
                                                (np.log(1 - alpha_max**2) + alpha_max**2/2*(alpha_max**2 - 2)/(alpha_max**2 - 1) \
                                               - np.log(1 - alpha_min**2) - alpha_min**2/2*(alpha_min**2 - 2)/(alpha_min**2 - 1)))
        
        
        
        plot_all(ue.phys_surf, ue.calc_surf, ue.phys_umat, ue.calc_umat, \
                 E_b_surf, E_s_surf, normals, path = '', show = True)
        #plot_all(ue, E_s_surf, E_b_surf, normals)
        
    


def test_prime_maps(u, final = False):
    
    ext_surf                        =   make_periodic(u.ext_surf, 'surf', phi_period = u.phi_period)
    ext_umat                        =   make_periodic(u.ext_umat, 'umat', sym_op = u.sym_op)
    [rmin, rmax], [phimin, phimax]  =   valminmax(ext_surf)
    
    rdev_max                        =   0
    phidev_max                      =   0
    iphidev, irdev                  =   [], []
    
    for ir, r_vec in np.ndenumerate(np.delete(ext_surf, -1, 1)):        
        
        i, j                =   ir
        C                   =   u.curve_Lphi[i]/u.phi_period #(phimax - phimin)
        
        drprimedphi         =   u.dudphipol[i,j][0]
        dphiprimedphi       =   1 + u.dudphipol[i,j][1] #(u.phips[i, j + 1] - u.phips[ir]) / u.dphi_mat[ir]  #(u.phips[i,dj] - u.phips[i,j])/dphi    
        duzdphi             =   u.dudphipol[i,j][2] # u.du(r_vec, 'duz')[1]
        
        Ct                  =   np.sqrt((r_vec[0] + ext_umat[ir][0])**2*dphiprimedphi**2 \
                                        + drprimedphi**2 + duzdphi**2)   
        
        if  phidev_max  <   abs(Ct - C):
            phidev_max  =   abs(Ct - C) 
            iphidev     =   [i,j]

    for iphi, phi_vec in np.ndenumerate(np.delete(ext_surf, -1, 0)):        
        
        i, j                =   iphi
        K                   =   u.curve_Lr[j]/(rmax - rmin)
        
        drprimedr           =   1 + u.dudrpol[i,j][0] #(u.rps[i + 1, j] - u.rps[iphi]) / u.dr_mat[iphi]   
        dphiprimedr         =   u.dudrpol[i,j][1] 
        duzdr               =   u.dudrpol[i,j][2] # u.du(r_vec, 'duz')[1]
        
        rprime              =   phi_vec[0] + ext_umat[i, j][0]
        
        Kt                  =   np.sqrt(drprimedr**2 + rprime**2*dphiprimedr**2 + duzdr**2)   
        
        if  rdev_max    <   abs(Kt - K) :
            rdev_max    =   abs(Kt - K) 
            irdev       =   [i,j]    
    
    if final and rdev_max > 1e-10:
        print 'screwed up: rdev_max = ' + str(rdev_max) 

    if final and phidev_max > 1e-10:
        print 'screwed up: phidev_max = ' + str(phidev_max) 
    
    return rdev_max, phidev_max, irdev, iphidev

def path_integral_r(u):
    
    ext_surf    =   make_periodic(u.ext_surf, 'surf', phi_period = u.phi_period)
    ext_umat    =   make_periodic(u.ext_umat, 'umat', sym_op = u.sym_op)

    inte_r      =   np.zeros(ext_umat.shape[1])
    inte_phi    =   np.zeros(ext_umat.shape[0])
    
    for ir, r_vec in np.ndenumerate(np.delete(ext_surf, -1, 0)):
        
        i,j         =   ir
        dvr         =   ext_umat[i + 1, j] + ext_surf[i + 1, j] \
                      - ext_umat[ir]       - ext_surf[ir]   
        rp          =   r_vec[0] + ext_umat[ir][0]
        
        inte_r[ir[1]] +=  np.sqrt(dvr[0]**2 + rp**2*dvr[1]**2 + dvr[2]**2) 
    
    for iphi, phi_vec in np.ndenumerate(np.delete(ext_surf, -1, 1)):

        i,j         =   iphi
        dvphi       =   ext_umat[i, j + 1]     + ext_surf[i, j + 1] \
                      - ext_umat[iphi]         - ext_surf[iphi]    
        rp          =   phi_vec[0] + ext_umat[iphi][0]
        
        inte_phi[i]+=  np.sqrt(dvphi[0]**2 + rp**2*dvphi[1]**2 + dvphi[2]**2) 
    

    for i, integ in enumerate(inte_phi):
        if abs(integ - u.curve_Lphi[i]) > 1e-10:
            print 'curve l: constant r, screwed up... dev = ' + str(abs(integ - u.curve_Lphi[i]))

    for i, integ in enumerate(inte_r):
        if abs(integ - u.curve_Lr[i]) > 1e-10:
            print 'curve l: constant phi, screwed up... dev = ' + str(abs(integ - u.curve_Lr[i]))
    
     
        
    #print inte, valminmax(ext_surf)[0][1] - valminmax(ext_surf)[0][0]
        
        
                    