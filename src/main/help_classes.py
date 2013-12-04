'''
Created on Sep 27, 2013

@author: tohekorh
'''
import numpy as np
import sys

pi = np.pi

def select_option(system, phi_period, hangle):
    
    opt     =   np.zeros(3, dtype = int)
    sym_op  =   ['', 0.]
    
    if system == 'cup':
        opt         =   [0,0,2]
        sym_op[0]   =   'rotation'
        sym_op[1]   =   phi_period 
    elif system == 'ball':
        opt         =   [0,0,4]
        sym_op[0]   =   'rotation'
        sym_op[1]   =   phi_period 

    elif system == 'spiral':
        opt         =   [1,0,0]
        sym_op[0]   =   'roto_trans'
        sym_op[1]   =   phi_period*hangle 

    elif system == 'spiral_w_wave':    
        opt     =   [2,3,1]
        sym_op[0]   =   'roto_trans'
        sym_op[1]   =   phi_period*hangle 

    elif system == 'screw':
        opt     =   [0,0,3]
        sym_op[0]   =   'roto_trans'
        sym_op[1]   =   phi_period*hangle 
    
    return opt, sym_op    
    
def heaviside(rad, rlim):
    
    if rad < rlim:
        return 0
    else:
        return 1

def get_z_set(rmat, hangle, consts, phiperiod):

    import scipy.optimize
    #print 'get_z consts = ' + str(consts) 
    w_num   =   consts[4]
    
    w_num_f =   0
    if w_num != 0:
        w_num_f = 1
    
    A       =   consts[1]*w_num_f
    z_set   =   np.zeros(rmat.shape)
    
    #print consts
    
    rs      =   strip(rmat[:,0], 'rs')
    Amps    =   np.zeros(len(rs))
    c_start =   consts[2] + rs[0] 
    
    for i, r in enumerate(rs):
        Amps[i] = heaviside(r, c_start)*((r - c_start)/ \
                                        (rs[-1] - c_start))**2*A
    
    for index, phis_vec in enumerate(rmat[:,]):
        
        phis    = strip(phis_vec, 'phis')
        r       = phis_vec[0][0]           
        phidils = phis
        beta    = np.arctan(hangle/r)

        #print r, A, beta, hangle
        
        if Amps[index] != 0:
            
            def F(phis):
                zero        =  np.sqrt(hangle**2 + r**2)*phis*np.cos(beta) \
                            - Amps[index]*np.sin(phis/phiperiod*2*np.pi*w_num)*np.sin(beta) \
                            - r*phidils
                #print zero
                #print phidils - phis
                #print 
                return zero
            
            nphi            =   scipy.optimize.broyden1(F, phis, f_tol=1e-6, \
                                                maxiter = len(phis)*1000 + 1)
            #print nphi, index, z_set
            z_set[index,:]  =   np.sqrt(hangle**2 + r**2)*nphi*np.sin(beta) \
                            + Amps[index]*np.sin(nphi/phiperiod*2*np.pi*w_num)*np.cos(beta)
        
        else:
            z_set[index,:]  = hangle*phis 
        
    return z_set
    
def strip(vec, key):
    
    new_vec = np.zeros(len(vec))
    for i, val in enumerate(vec):
        if key == 'phis':
            new_vec[i] = val[1]    
        elif key == 'rs':
            new_vec[i] = val[0]
    return new_vec    
    
def get_quess(opt, ext_surf, hangle, phiperiod, n_w):
    
    rmin,rmax           =   valminmax(ext_surf)[0]
    width               =   rmax - rmin
    dw                  =   width/30
    x                   =   np.sqrt(rmin**2 - hangle**2) - rmin
#    middle              =   rmin + width / 2 + x
    if n_w != 0:
        #print np.arctan(hangle/(rmin + x))
        #print np.arctan(3/phiperiod*2*pi)
        #print pi/2
        A_max   =   np.tan(pi/2 - pi/2/7 - np.arctan(hangle/(rmin + x)))\
                    *phiperiod/(2*pi*n_w)*(rmax + x)
        A_quess =   (hangle*2*np.pi)/50./n_w
    else:
        A_max   =   0.
        A_quess =   0.
    
    A_max       =   min(A_max, 5)
    
    if x != 0:
        middle              =   -hangle**2/(2*x) - rmin - x/2
    else:
        middle              =   0
        
    if  opt == [2,3,1]:
#        o_consts        =   [x, (hangle*2*np.pi)/340., middle - rmin, width] # 
        consts          =   [2./3.*x, A_quess, middle, width, n_w] # 
 
        bounds          =   [(x*4./3., 0.01), (0., A_max), \
                             (0.0, width), (width - dw, width + dw)]
        
    elif opt == [1,0,0]:
        consts        =   [x]
        bounds          =   [(-rmin, 0)]
    
    
    return consts, bounds

def reduce_bounds(bounds):
    
    bounds[1] = (bounds[1][0], bounds[1][1]*9/10)
    
    return bounds

def curve_Lphi(u, test = True):
    
    # Lphi = curve length along circumference with radius r
    init_surf       =   make_periodic(u.ext_surf, 'surf', phi_period = u.phi_period) 
    umat            =   make_periodic(u.ext_umat, 'umat', sym_op = u.sym_op)
    dudphipol       =   u.dudphipol
    
    Lphi            =   np.zeros(len(init_surf[:,0]))
    
    for ir, r_vec in np.ndenumerate(np.delete(init_surf, -1, 1)):
        
        dphi        =   u.dphi_mat[ir]  #init_surf[ir[0], ir[1] + 1][1] - init_surf[ir[0], ir[1]][1] 
        r           =   r_vec[0]
        
        ds          =   np.sqrt((r + umat[ir][0])**2*(1 + dudphipol[ir][1])**2 +  \
                                 + dudphipol[ir][0]**2 + dudphipol[ir][2]**2)*dphi
        #print r, ds, dudphipol[ir]
        Lphi[ir[0]]+=   ds   
    
    
    if u.consts[1] == 0 and test:
        print 'testing curve lengths...'
        scw_up          =   False
        hangle          =   u.hangle
        phimin, phimax  =   valminmax(init_surf)[1]
        for j in range(len(Lphi)):
            dev         =   Lphi[j] - (np.sqrt((init_surf[j,0][0] +  \
                            umat[j,0][0])**2 + hangle**2)*\
                            (phimax - phimin))
            if abs(dev) > 1e-12:
                print dev, init_surf[j,0][0], Lphi[j]
                scw_up  =   True
        if not scw_up:
            print 'test succesfull'
        else:
            print 'SCRWEWED UP'

        
        
    return Lphi
    
def curve_Lr(u):

    # Lr = curve length along radius with angle phi
    
    len_Lr          =   u.ext_surf.shape[1] + 1    
    Lr              =   np.zeros(len_Lr) 
    Lr[:]           =   u.consts[3]
    
    return Lr


def rps(u):
    
    # rps are the primed r coordinates
    rps                 =   np.zeros(u.ext_surf.shape)
    
    [rmin, rmax]        =   valminmax(u.ext_surf)[0]
    initial_width       =   rmax - rmin
               
    for index, r_vec in np.ndenumerate( u.ext_surf ):
        
        k, l            =   index
        K               =   u.curve_Lr[l]/initial_width
        
        if k == 0:
            rps[0, l]  =   r_vec[0]    +   u.consts[0]
        
        elif 0 < k:
            
            dr          =   u.dr_mat[k - 1, l]
            drp         =   np.sqrt(K**2 - rps[k - 1, l]**2*u.dudrpol[k - 1, l][1]**2   \
                                                          - u.dudrpol[k - 1, l][2]**2)*dr
            rps[k, l]   =   rps[k - 1, l] + drp  
                
    return rps

def phips(u):    
    
    #phips are the primed phi coordinates
    phips                   =   np.zeros(u.ext_surf.shape)
    
    r                       =   0
    for index, r_vec in np.ndenumerate( u.ext_surf ):
        k, l                =   index
        
        if l == 0:
            phips[k, 0]     =   r_vec[1] 
        
        elif 0 < l < u.ext_surf.shape[1]:
            
            C               =  u.curve_Lphi[k]/u.phi_period
            dphi            =  u.ext_surf[k, l][1] - u.ext_surf[k, l - 1][1] 
            
            if C**2 - u.dudphipol[k, l - 1][0]**2 - u.dudphipol[k, l - 1][2]**2 < 0:
                print 'pipsit kosahti'
                raise ValueError 
            dphip           =  1. / ( r + u.ext_umat[k, l - 1][0] ) \
                                *np.sqrt(C**2 - u.dudphipol[k, l - 1][0]**2 \
                                - u.dudphipol[k, l - 1][2]**2)*dphi
            
            phips[k, l]     =  phips[k, l - 1] + dphip #+ u.ext_surf[0, 0][1]
        
        r                   =   r_vec[0]
            
    return phips 


def make_periodic_der(mat, sgn = 1, nums = False):
    
    shape               =   (mat.shape[0], mat.shape[1] + 1)
    per_mat             =   np.empty(shape,dtype = 'object')
    
    if sgn == 1:
        per_mat[:,:-1]      =   mat
        for iv, v in enumerate(mat[:,0]):
            if not nums:
                per_mat[iv, -1]   =   np.zeros(3)
                per_mat[iv, -1][0]=   v[0]
                per_mat[iv, -1][1]=   v[1]
                per_mat[iv, -1][2]=   v[2]
            else:
                per_mat[iv, -1]     =   np.float
                per_mat[iv, -1]     =   v
                
    elif sgn == -1:
        per_mat[:,1:]      =   mat
        for iv, v in enumerate(mat[:,-1]):
            per_mat[iv, 0]   =   np.zeros(3)
            per_mat[iv, 0][0]=   v[0]
            per_mat[iv, 0][1]=   v[1]
            per_mat[iv, 0][2]=   v[2]
    
    return per_mat


def make_periodic(mat, key, sym_op = [], phi_period = -1):
    
    shape                           =   mat.shape[0], mat.shape[1] + 1    
    periodic_mat                   =   np.empty(shape, dtype = 'object')
    
    periodic_mat[:,:-1]            =   mat
    
    if key == 'umat':
        if sym_op[0] == 'roto_trans':
            for ind, vec in enumerate(mat[:,0]):
                periodic_mat[ind, -1]       =   np.zeros(3)
                periodic_mat[ind, -1][:2]   =   vec[:2]
                periodic_mat[ind, -1][2]    =   vec[2] + sym_op[1]
        elif sym_op[0] == 'rotation':
            for ind, vec in enumerate(mat[:,0]):
                periodic_mat[ind, -1]       =   np.zeros(3)
                periodic_mat[ind, -1][0]    =   vec[0]
                periodic_mat[ind, -1][1]    =   vec[1] #+ sym_op[1]
                periodic_mat[ind, -1][2]    =   vec[2]
        else:
            print 'EIIII'   
    elif key == 'surf':
        if phi_period == -1:
            print 'kaaaak!'
        for ind, vec in enumerate(mat[:,0]):
            periodic_mat[ind, -1]      =   np.zeros(3)
            periodic_mat[ind, -1][0]   =   vec[0]
            periodic_mat[ind, -1][1]   =   vec[1] + phi_period
            periodic_mat[ind, -1][2]   =   vec[2]

        
    return periodic_mat 
'''
def make_periodic_surf(surf, del_phi):
    
    shape                               =   surf.shape[0], surf.shape[1] + 1
    periodic_surf_mat                   =   np.empty(shape, dtype = 'object')
    
    periodic_surf_mat[:,:-1]            =   surf
    
    for ind, vec in enumerate(surf[:,0]):
        periodic_surf_mat[ind, -1]      =   np.zeros(3)
        periodic_surf_mat[ind, -1][0]   =   vec[0]
        periodic_surf_mat[ind, -1][1]   =   vec[1] + del_phi
        periodic_surf_mat[ind, -1][2]   =   vec[2]
           
    return periodic_surf_mat 
'''

def valminmax(mat):
    
    rmax    =   -1e10
    rmin    =   1e10
    phimax  =   -1e10
    phimin  =   1e10
    
    for ir, r_vec in np.ndenumerate(mat):
        r                   =   r_vec[0]
        phi                 =   r_vec[1]
        
        if r    > rmax:      rmax   =   r
        if r    < rmin:      rmin   =   r
        if phi  > phimax:    phimax = phi
        if phi  < phimin:    phimin = phi
    
    return [rmin, rmax], [phimin, phimax] 

def calc_area(u):
    
    Ap  =   0
    for ida, dA in np.ndenumerate(u.dA_P):
        Ap   += dA   
    
    As  =   0
    for ida, dA in np.ndenumerate(u.dA_S):
        As   += dA   
    
    return Ap, As

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


