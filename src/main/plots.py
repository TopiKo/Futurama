'''
Created on Oct 3, 2013

@author: tohekorh
'''
import numpy as np
from surface import get_coord 
import matplotlib.pyplot as plt
from os.path import exists
from os import makedirs
from help_classes import query_yes_no
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from bender_rw import read_bender_output

def plot_e_surfaces(in_file, show = False):
    
    asurf, umat, E_surfs, normals, path     =   read_bender_output(in_file)[1:] 
    
    plot_all(asurf[1], asurf[2], umat[1], umat[2], E_surfs[0], E_surfs[1], \
              normals, path = path, show = show)    

def plot_all(phys_surf, calc_surf, phys_umat, calc_umat, E_b, E_s, normals, path = '', show = False):
    
    phi_slice               =   0
    
    X_init, Y_init, Z_init  =   get_coord(phys_surf) #init_surf)
    X, Y, Z                 =   get_coord(phys_surf + phys_umat) #deform(init_surf, ue))
    X_ss, Y_ss, Z_ss        =   get_coord(calc_surf + calc_umat) #deform(e_surf, ue))    
    
    # PLOT BEGINS!
    
    if np.amax(E_s) != 0:
        N_s = E_s/np.amax(E_s)
    else:
        N_s = E_s
    
    if np.amax(E_b) != 0:
        N_b = E_b/np.amax(E_b)
    else:
        N_b = E_b
        
    
    x, y, nx, ny, nz, nr = [],[],[],[],[], []
    
    r = np.zeros(len(calc_surf[:,phi_slice]))
    
    for ik, k in enumerate(np.array(calc_surf[:,phi_slice] + calc_umat[:,phi_slice])): 
        r[ik]   =   k[0] 
    
    z           =   Z_ss[:,phi_slice]
    
    
    for i in X_ss:
        for j in i:
            x.append(j)
    for i in Y_ss:
        for j in i:
            y.append(j)

    n = 0
    
    for i in range(normals.shape[0]):
        for j in range(normals.shape[1]):
            normal      =   normals[i,j]
            if abs(normal[0]**2 + normal[1]**2 + normal[2]**2 - 1) > 0.00001 :
                print 'ALERT ALERT!1' 
            nx.append(normal[0])
            ny.append(normal[1])
            nnr = x[n]/np.sqrt(x[n]**2 + y[n]**2)*normal[0] + y[n]/np.sqrt(x[n]**2 + y[n]**2)*normal[1]            
            if j == phi_slice:
                nr.append(nnr)
                nz.append(normal[2]) 
            n += 1

    limits          =   [np.amin([np.amin(X),np.amin(Y),np.amin(Z)]), np.amax([np.amax(X),np.amax(Y),np.amax(Z)])]
        
    fig             =   plt.figure(figsize=plt.figaspect(0.5)*1.5)
    
    #1
    ax              =   fig.add_subplot(121, projection='3d')
    
    vv              =   max(int(len(X)/30), 1)
    
    ax.plot_surface(X_init, Y_init, Z_init, alpha = 0.2, rstride = 4*vv, cstride = 4*vv) 
    ax.plot_surface(X, Y, Z, rstride = vv, cstride = vv, \
                    alpha = 1., facecolors=cm.cool(N_b), shade=False) #OrRd
    ax.auto_scale_xyz(limits, limits, limits)
    
    
    
    #2
    ax1 = fig.add_subplot(122, projection='3d')
    ax1.plot_surface(X_init, Y_init, Z_init, alpha = 0.2, rstride = 4*vv, cstride = 4*vv) 
    ax1.plot_surface(X, Y, Z, rstride = vv, cstride = vv, \
                    alpha = 1., facecolors=cm.cool(N_s), shade=False)
    ax1.auto_scale_xyz(limits, limits, limits)
    
    if path != '':
        if not exists(path + 'pictures/'):
            makedirs(path + 'pictures/')
        
        plt.savefig(path + 'pictures/ener_surf.png')
     
    fig2 = plt.figure(figsize=plt.figaspect(0.5)*1.5)
    
    #3
    ax2 = fig2.add_subplot(131)
    ax2.quiver(x, y, nx, ny, scale=3.2) #, units='width')
    ax2.axis('equal')
    
    ax3 = fig2.add_subplot(132)
    ax3.scatter(X, Y, marker = "+") #, units='width')
    ax3.axis('equal')
    
    
    #4
    ax4 = fig2.add_subplot(133)
    ax4.plot(r,z)
    ax4.quiver(r, z, nr, nz, scale=2.8) #, units='width')
    ax4.axis('equal')

    if path != '':
        plt.savefig(path + 'pictures/normals.png')
    
    if show:
        plt.show()
    #plt.clf()
    
def plot_surf(u):
    
    from help_classes import make_periodic
    
    umat        =   make_periodic(u.ext_umat, 'umat', sym_op = u.sym_op)
    init_surf   =   make_periodic(u.ext_surf, 'surf', phi_period = u.phi_period)
    

    X_init, Y_init, Z_init  =   get_coord(init_surf)
    X, Y, Z                 =   get_coord(init_surf + umat)
    xs, ys, zs              =   get_coord(u.phys_surf + u.phys_umat)    
    xss, yss, zss           =   get_coord(u.calc_surf + u.calc_umat)    

    limits          =   [np.amin([np.amin(X),np.amin(Y),np.amin(Z)]), np.amax([np.amax(X),np.amax(Y),np.amax(Z)])]
        
    fig             =   plt.figure(figsize=plt.figaspect(0.5)*1.5)
    
    vv              =   max(int(len(X)/30), 1)
    
    
    #1
    ax              =   fig.add_subplot(131, projection='3d')
    
    ax.plot_surface(X_init, Y_init, Z_init, alpha = 0.2, rstride = vv, cstride = vv) 
    ax.plot_surface(X, Y, Z, rstride = vv, cstride = vv, alpha = 1., shade=False)
    
    ax.auto_scale_xyz(limits, limits, limits)
    
    
    
    #2
    ax1 = fig.add_subplot(132, projection='3d')
    ax1.plot_surface(xs, ys, zs, rstride = vv, cstride = vv, alpha = 1., shade=False) #, c='r', marker='o')
    ax1.auto_scale_xyz(limits, limits, limits)
    
    #3
    ax2 = fig.add_subplot(133, projection='3d')
    ax2.plot_surface(xss, yss, zss, rstride = vv, cstride = vv, alpha = 1., shade=False) #, c='r', marker='o')
    ax2.auto_scale_xyz(limits, limits, limits)
    
    
    plt.show()
    plt.clf()
    
def plot_energies(ir, rmin, rmax, heights, E_bs, E_ss, heights_hb, energies_hb, save_dir):

    
    ax1 = plt.subplots()[1]
    
    e_tot       =   E_bs + E_ss
    e_hb        =   np.zeros(len(energies_hb))
    for ie, e in enumerate(energies_hb):
        e_hb[ie] =   e - min(energies_hb) 
    

    ax1.plot(heights,       e_tot,  '-o',   label = 'code tot')
    ax1.plot(heights,       E_bs,   '--',   label = 'e_b, code')
    ax1.plot(heights,       E_ss,   '-',    label = 'e_s, code')
    ax1.plot(heights_hb,    e_hb,   '-D',   label = 'HB')
    
    ax1.set_xlabel('Height')
    ax1.set_ylabel('Total deform energy in eV')
    ax1.set_title('Deformation energy of spiral ir = ' + str(ir) \
                  + ', \n rmin, rmax = %.2f, %.2f' %(rmin, rmax))
    
    plt.legend(loc = 2)
    
    #print e_tot
    
    if len(heights) == len(heights_hb):
        ax2 = ax1.twinx()
        ax2.plot(heights_hb, 100*(e_tot/e_hb - 1), '-o', color = 'r', label = 'dev = code/hb - 1, %')
    
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        
        ax2.set_ylabel('dev between code and hb')
        plt.legend(loc = 2)
    
    
        
    plt.savefig(save_dir + 'energies_wrt_h.png')
    plt.clf()
    
def plot_relations(rad_rels, height_rels, amplitudes, n_waves, folder):
    
    ax     =   plt.subplots()[1]
    
    s   = np.zeros(len(amplitudes))
    
    for ia, a in enumerate(amplitudes):
        s = a*1000 + 10
        if n_waves[ia] == 0: color = 'b'
        elif n_waves[ia] == 1: color = 'g'
        elif n_waves[ia] == 2: color = 'r'
        elif n_waves[ia] == 3: color = 'c'
        
        ax.scatter(rad_rels[ia], height_rels[ia], c= color, s=s)   
    
    ax.set_xlabel('rmax/rmin')
    ax.set_ylabel('height/rmax')
    
    ax.set_title('Amplitude/rmax = size of the ball, \n num of waves = color, blue =0, green = 1, red = 2')
    
    plt.savefig(folder + 'amplitudes.png')
    
def plot_amplitude(ir, rmin, rmax, heights, consts, save_dir):

    
    ax1     =   plt.subplots()[1]
    
    Amps    =   np.zeros(len(consts))
    n_waves =   np.zeros(len(consts))
    
    for ic, cs in enumerate(consts):
        Amps[ic]        =   cs[1] 
        n_waves[ic]     =   cs[4]
    

    ax1.plot(heights, Amps,  '-o',   label = 'Amplitude')
    
    plt.legend()
    
    ax2 = ax1.twinx()
    ax2.plot(heights, n_waves, '-o', color = 'r', label = 'num of waves')

    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    
    ax1.set_xlabel('Height')
    ax1.set_ylabel('Amplitude angstroms')
    ax2.set_ylabel('Number of waves in u cell')
    
    ax1.set_title('Amps and n_waves ir = ' + str(ir) \
                  + ', rmin, rmax = %.2f, %.2f' %(rmin,rmax))
        
    plt.legend(loc = 2)
        
    plt.savefig(save_dir + 'Amps_waves.png')
    
    plt.clf()
    
def plot_consts_proxm(x_mesh, A_mesh, E_b_surf, E_s_surf, folder, \
                      from_file = False, acc = 100):
    
    E_surf  =   E_b_surf + E_s_surf
    
    fig     =   plt.figure(figsize=plt.figaspect(0.5)*1.5)
    ax      =   fig.add_subplot(111)
    CS      =   plt.contour(x_mesh, A_mesh, E_surf, acc,
                 colors='k')
    plt.clabel(CS, fontsize=9, inline=1)
    plt.title('Energy contours in the proximity of optimal constants')
    
    ax.set_xlabel('x')
    ax.set_ylabel('Amplitude')
    #ax.set_zlabel('Energy')

    plt.savefig(folder + 'consts_prxm.png')
    
    
    if from_file:
    
        plt.show()
        plt.clf()
    
        if query_yes_no('separate plot', "no"):
            
            if not exists(folder + 'pictures/'):
                makedirs(folder + 'pictures/')
            
            for k in range(len(A_mesh[:,0])):
                fig     =   plt.figure(figsize=plt.figaspect(0.5)*1.5)
                ax      =   fig.add_subplot(111)
                x       =   x_mesh[0, k]
                ax.plot(A_mesh[:,k], E_surf[:,k]) 
                ax.plot(A_mesh[:,k], E_b_surf[:,k], label = 'E_b') 
                ax.plot(A_mesh[:,k], E_s_surf[:,k], label = 'E_s') 
                
                ax.set_xlabel('A')
                ax.set_ylabel('Energy, x = %f.2' %x )
                plt.legend()
                plt.savefig(folder + 'pictures/x=%f.2.png' %x)
                
                plt.clf()
                
            for k in range(len(A_mesh[0,:])):
                fig     =   plt.figure(figsize=plt.figaspect(0.5)*1.5)
                ax      =   fig.add_subplot(111)
                A       =   A_mesh[k,0]
                
                ax.plot(x_mesh[k,:], E_surf[k,:]) 
                ax.plot(x_mesh[k,:], E_b_surf[k,:], label = 'E_b') 
                ax.plot(x_mesh[k,:], E_s_surf[k,:], label = 'E_s') 
                
                ax.set_xlabel('x')
                ax.set_ylabel('Energy, A = %f.2' %A )
                plt.legend()
                
                plt.savefig(folder + 'pictures/A=%f.2.png' %A)
                plt.clf()
        
    plt.clf()

def plot_curve():
    
    hangle  = 73/2/np.pi
    r       = 20 
    A       = 6.1
    acc = 100
    
    phiperiod = np.pi/3
    
    phis    = np.linspace(0, phiperiod, acc)    
    xs      = np.zeros(acc)
    cs      = np.zeros(acc)
    xsp     = np.zeros(acc)
    csp     = np.zeros(acc)
    
    beta    = np.arctan(hangle/r) 
    
    
    import scipy.optimize
    
    def F(phi):
        print phi
        return np.sqrt(hangle**2 + r**2)*phi*np.cos(beta) \
            - A*np.sin(phi/phiperiod*2*np.pi)*np.sin(beta) - r*phis
    
    
    
    for iphi, phi in enumerate(phis):
        xs[iphi] = np.sqrt(hangle**2 + r**2)*phi 
        cs[iphi] = A*np.sin(phi/phiperiod*2*np.pi)  
        xsp[iphi] = xs[iphi] * np.cos(beta) - cs[iphi]* np.sin(beta) 
        csp[iphi] = xs[iphi] * np.sin(beta) + cs[iphi]* np.cos(beta) 
        
    
    nphi    = scipy.optimize.broyden1(F, phis, f_tol=1e-14)
    ncsp    = np.sqrt(hangle**2 + r**2)*nphi*np.sin(beta) + A*np.sin(nphi/phiperiod*2*np.pi)*np.cos(beta) 
    
    print nphi
    print phis    
    print F(phis)
    
    fig     =   plt.figure(figsize=plt.figaspect(1.))
    ax      =   fig.add_subplot(111)

    ax.plot(xs, cs) 
    ax.plot(xsp, csp) 
    ax.plot(r*phis, ncsp, 'o') 
    
    ax.plot(phis*r, hangle*phis)

    #ax.plot(phis, vals) 
    plt.show()
    plt.clf()

#plot_curve()
