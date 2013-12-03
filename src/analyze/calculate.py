'''
Created on 22.11.2013

@author: tohekorh
'''


from main.bender_rw import read_bender_output, parse_path, \
        read_synopsis, read_energies, parse_input, read_total_synopsis 
from main.strain import u 
import numpy as np
from main.deform_energies import deform_energies
from numpy.lib.function_base import meshgrid
from os.path import exists
from os import makedirs
from main.plots import plot_consts_proxm, plot_energies, plot_amplitude, plot_relations

pi  =   np.pi
acc =   3

def study_height_synopsis(in_file):
    
    #in_file = '/space/tohekorh/Spiral/bender_input/calc/consts/in_8-16_H=73.0.txt'
    
    ir                        =   parse_input(in_file)[("ir")]
    data_array, syn_dir       =   read_synopsis(in_file) # heights, E_bs, E_ss, consts
    heights_hb, energies_hb   =   read_energies(ir)[:2]
    
    plot_energies(ir, data_array[("heights")], data_array[("E_bs")], data_array[("E_ss")], \
                   heights_hb, energies_hb, syn_dir)
    
    plot_amplitude(ir, data_array[("heights")], data_array[("consts")], syn_dir)

def study_total_synopsis(nr, nphi, phiperiod, system):
    
    syn_params, folder =   read_total_synopsis(nr, nphi, phiperiod, system)
    
    rad_rels        =   []
    height_rels     =   []
    amplitudes      =   []
    n_waves         =   []
        
    for params in syn_params:
        rmin        =   syn_params[(params)][("rmin")]
        rmax        =   syn_params[(params)][("rmax")]
        height      =   syn_params[(params)][("height")]
        amplitude   =   syn_params[(params)][("consts")][1]
        
        rad_rels.append(rmax/rmin)
        height_rels.append(height/rmax)
        amplitudes.append(amplitude/rmax)
        n_waves.append(syn_params[(params)][("consts")][4])
        
    print amplitudes
    print height_rels
    
    plot_relations(rad_rels, height_rels, amplitudes, n_waves, folder)
        
    
def study_consts_proximity(in_file):
    
    print 'studying consts proximity..'
    
    try:
        params, asurf =   read_bender_output(in_file)[:2] # + 'params.txt')
    except IOError as e:
        print 'there, is no file for this input... got ' + e 
                
    phi_period  =   params[("phiperiod")]
    hangle      =   params[("height")] / 2 / pi
    system      =   params[("system")]
    
    ue          =   u(hangle, phi_period, asurf, system = system, consts = params[("consts")])
    energies    =   deform_energies(ue)
    
    delx        =   abs(params[("consts")][0]) / 12.
    delA        =   abs(params[("consts")][1]) / 3.
    x, A        =   params[("consts")][:2]
    
    
    x_set       =   np.linspace(x - delx, x + delx, acc)        
    A_set       =   np.linspace(A - delA, A + delA, acc)  
    curve_start =   params[("consts")][2]
    mid         =   params[("consts")][3]
    n_wave      =   params[("consts")][4]
    
    mesh_x, mesh_A \
                =   meshgrid(x_set, A_set)  
    
    E_b_mat     =   np.zeros(mesh_x.shape)
    E_s_mat     =   np.zeros(mesh_x.shape)
    
    for index, x in np.ndenumerate(mesh_x):
        A       =   mesh_A[index]
        consts  =   [x, A, curve_start, mid, n_wave]
        E_b_mat[index], E_s_mat[index] \
                =   energies.calc_energies(consts)[:2]
        
    #print E_s_mat
    #print E_b_mat
    
    folder = parse_path(params, params[("moldy_opm")]) + 'opm_consts_proxm/'
    
    if not exists(folder):
        makedirs(folder)   
    
    np.save(folder + 'mesh_x'   , mesh_x) 
    np.save(folder + 'mesh_A'   , mesh_A) 
    np.save(folder + 'E_b_mat_consts_proxm', E_b_mat) 
    np.save(folder + 'E_s_mat_consts_proxm', E_s_mat) 
    
    plot_consts_proxm(mesh_x, mesh_A, E_b_mat, E_s_mat, folder, from_file = False, acc = acc)

def study_consts_proximity_from_file():
    
    in_file = '/space/tohekorh/Spiral/bender_input/calc/consts/in_8-16_H=73.0.txt'
    
    params  =   read_bender_output(in_file)[0] # + 'params.txt')
    folder  =   parse_path(params, params[("moldy_opm")]) + 'opm_consts_proxm/'
    
    mesh_x  =   np.load(folder + 'mesh_x.npy') 
    mesh_A  =   np.load(folder + 'mesh_A.npy') 
    E_b_mat =   np.load(folder + 'E_b_mat_consts_proxm.npy') 
    E_s_mat =   np.load(folder + 'E_s_mat_consts_proxm.npy') 
    
    #print params[("consts")][0]
    
    plot_consts_proxm(mesh_x, mesh_A, E_b_mat, E_s_mat, folder, from_file = True)

#study_height_synopsis()  
#study_consts_proximity_from_file()
#study_consts_proximity('/space/tohekorh/Spiral/bender_input/calc/consts/in_8-16_H=73.0.txt')    