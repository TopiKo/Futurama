'''
Created on Oct 3, 2013

@author: tohekorh
'''

import numpy as np

from strain import parse_u_from_file
from tests import tests
#from surface import surf 
#from deform_energies import deform_energies
#from plots import plot_all, plot_surf
#from bender_rw import write_params, read_rad, \
#                read_energies, c_logfile, \
#                 read_h_consts
from bender_rw import write_synopsis, write_bender, parse_input, \
    read_bender_output, has_been_calculated, write_total_synopsis
from help_classes import query_yes_no, get_sys_sets #, valminmax
from optimize_consts import optimize
from analyze.calculate import study_consts_proximity, \
            study_height_synopsis, study_total_synopsis
from plots import plot_e_surfaces
from os.path import exists
from os import makedirs
from datetime import datetime
import gc

pi                  =   np.pi

sigma               =   0.28    # ~ 0.2 - 0.3
bend_module         =   1.6     # = E*thickness**3/(12*(1+sigma))
strech_module       =   25 


def plot():
    
    system                          =   'spiral_w_wave_sqr'
    nr, nphi, phiperiod             =   20, 40, pi/3
    
    from os import listdir
    
    input_file_folder_consts        =   '/space/tohekorh/Spiral/bender_input/calc/consts/' \
                                           + '%s/phiper=%.2f/nr-nphi=%i-%i/'  \
                                            %(system, phiperiod, nr, nphi)
        #input_file_folder_moldy     =   '/space/tohekorh/Spiral/bender_input/calc/moldy/'

    
    for input_folder in listdir(input_file_folder_consts):
        input_file_folder   = input_file_folder_consts + input_folder + '/'
        
        if input_folder != 'store':
            for input_file in listdir(input_file_folder):
                if input_file[-3:] == 'txt':
                    if has_been_calculated(input_file_folder + input_file, 'read'):
                        in_file =  input_file_folder + input_file
                        if query_yes_no("plot this " + input_file, 'no'):
                            #study_height_synopsis(in_file)
                            #plot_e_surfaces(in_file, show = True)
                            #study_height_synopsis(in_file)
                            study_consts_proximity(in_file, acc = 30)
                        #ue              =   parse_u_from_file(input_file_folder +  input_file)
                        #params          =   read_bender_output(input_file_folder +  input_file)[0]
                        #ue.set_const(params[("consts")])
                        #energies        =   deform_energies(ue)
                        #E_b, E_s, E_b_surf, E_s_surf, normals  = energies.calc_energies()
                        
    write_total_synopsis(nr, nphi, phiperiod, system, input_file_folder_consts)
    
    study_total_synopsis(nr, nphi, phiperiod, system)
    

def run_tests():
    
    from os import listdir
    
    input_file_folder_consts    =   '/space/tohekorh/Spiral/bender_input/tests/'
        
    for input_file in listdir(input_file_folder_consts):
        if input_file[-3:] == 'txt':
            
            param_set   =   parse_input(input_file_folder_consts + input_file)
            nr, nphi    =   param_set["nr"],    param_set["nphi"]
    
            tests_a     =   tests(nr, nphi)
            tests_a.tests(param_set)


def run():
    
    from os import listdir

    system_sets  =  get_sys_sets() 
    
    for system_set in system_sets:
        
        system                          =   system_set[0] 
        [nr, nphi], phiperiod           =   system_set[1:3]
        moldy                           =   system_set[3]
        print nr,nphi, moldy
        
        gc.enable()
        
        input_file_folder_consts        =   '/space/tohekorh/Spiral/bender_input/calc/%s/' %moldy \
                                               + '%s/phiper=%.2f/nr-nphi=%i-%i/'  \
                                                %(system, phiperiod, nr, nphi)
            #input_file_folder_moldy     =   '/space/tohekorh/Spiral/bender_input/calc/moldy/'
        
        path_log    =   '/space/tohekorh/Spiral/bender_output/%s/' %system   
        
        if not exists(path_log):
            makedirs(path_log)

        
        for input_folder in listdir(input_file_folder_consts):
            input_file_folder   = input_file_folder_consts + input_folder + '/'
            
            if input_folder != 'store':
                for input_file in listdir(input_file_folder):
                    if input_file[-3:] == 'txt':
                        if not has_been_calculated(input_file_folder + input_file, 'read'):
                            logfile     =   open(path_log + 'logfile.txt', 'a')
                            in_file     =   input_file_folder + input_file
                            logfile.write(str(datetime.now()) + ' \n') 
                            logfile.write(str(system_set) + ' \n')
                            try:   
                            #if query_yes_no('system = %s, moldy = %s' %(system, moldy) + " run this " + input_file, 'no'):
                                param_set   =   parse_input(input_file_folder + input_file)
                                
                                run_bender(param_set, in_file)
                                
                                unreached   = gc.collect()
                                del gc.garbage[:]
                                print 'in collect there was %i unreached objects' %unreached
                                logfile.write(input_file + ' SUCCESS! \n\n')
                                has_been_calculated(in_file, 'write_ok')
                                
                            except Exception as e:
                                logfile.write(str(e) + ' \n')
                                logfile.write(input_file + ' ERROR! \n\n')
                                has_been_calculated(in_file, 'write_not_ok')
                                
                            logfile.close()
        
        if system != 'spiral':
            write_total_synopsis(nr, nphi, phiperiod, system, input_file_folder_consts)
            study_total_synopsis(nr, nphi, phiperiod, system)
            
        
def run_bender(param_set, in_file):

    system              =   param_set["system"]
    ir                  =   param_set["ir"]
    rmin, rmax          =   param_set["rmin"],  param_set["rmax"]
    height              =   param_set["height"]
    moldy_opm           =   param_set["moldy_opm"]
    
    print 'running system: ' + system + ', ir = ' + str(ir) 
    print 'height     = ' + str(height)
    print 'rmin, rmax = %.2f, %.2f' %(rmin, rmax)
    
    ue                  =   parse_u_from_file(in_file)
    optimize_c          =   optimize(ue, 'l_bfgs')
        
    if not moldy_opm:
            
        [E_b_surf, E_s_surf], [E_b, E_s], ue, normals   \
                        =   optimize_c.optimize_consts()
        
        write_bender(E_b_surf, E_s_surf, normals, E_b, E_s, \
                    ue, param_set, in_file)
        
        plot_e_surfaces(in_file)
        
        write_synopsis(in_file)
        
        if system != 'spiral':
            study_height_synopsis(in_file)
            try:
                study_consts_proximity(in_file, 20)
            except ValueError:
                print 'voi hoh'
        
        
        
    elif moldy_opm:
        
        #optimize_c  =   optimize(ue, 'l_bfgs')
        
        try:
            params      =   read_bender_output(in_file, read_for_moldy = True)[0]   
            ue.set_const(params[("consts")])
            ue.set_moldy(True)  
        
        except IOError:
            print 'got io error from trying to import constants for moldy opm!!'
            ue          =   optimize_c.optimize_consts()[2]
            
        
        [E_s_surf, E_b_surf], [E_b, E_s], ue    \
                    =   optimize_c.optimize_moldy()[:4]

        write_bender(E_b_surf, E_s_surf, E_b, E_s, ue, param_set)

#run_tests()
run()
#plot()
#main()  




