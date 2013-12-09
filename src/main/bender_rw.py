'''
Created on Oct 8, 2013

@author: tohekorh
'''
from os.path import exists, dirname, join
from os import makedirs, listdir
import numpy as np
from help_classes import valminmax
from numpy.core.numeric import ndarray
import pickle
#from classical import help_classes
from help_classes import select_option
from shutil import copy2
from main.help_classes import query_yes_no, get_sys_sets

    
path = '/space/tohekorh/Spiral/'


def parse_line(line, ntype = 'float'):
    
    split_line = line.split('#')[0].split(' ')
    
    data = np.zeros(len(split_line) - 1, dtype = ntype)
    
    
    for i, obj in enumerate(split_line[:-1]):
        if   ntype == 'float':
            if obj == 'pi/3':
                data[i]     =   np.pi/3
            else:
                data[i]     =   float(obj)
        elif ntype == 'int':
            data[i]         =   int(obj)
    
    return data

def parse_input(input_file):
    
    ifile = open(input_file, 'r')
    
    params  =   {}
    params["system"]                =   ifile.readline().split('#')[0].split(' ')[0]
    params["ir"]                    =   parse_line(ifile.readline(), ntype = 'int')
    params["nr"],   params["nphi"]  =   parse_line(ifile.readline(), ntype = 'int')
    params["rmin"], params["rmax"]  =   parse_line(ifile.readline(), ntype = 'float')
    params["phimin"], params["phiperiod"]   \
                                    =   parse_line(ifile.readline(), ntype = 'float')
    params["height"]                =   parse_line(ifile.readline(), ntype = 'float')[0]
    params["consts"]                =   parse_line(ifile.readline(), ntype = 'float')
    params["moldy_opm"]             =   ifile.readline().split('#')[0].split(' ')[0] in ['True', 'true']
    
    return params

class c_logfile():
    
    def __init__(self, path):
        print 'create path to the folder/create logfile to existing folder'
    
    def write(self, string):
        print' write the string to logfile with timestamp '


def parse_path(param_set, wrt_moldy):

    if wrt_moldy:
        ident   =   'opm_moldy/'
    else:
        ident   =   'opm_consts/'
    
    ir      = param_set[("ir")]
    p_path  = path + 'bender_output/' + str(param_set[("system")]) + '/' \
                     + str(ir[0]) + '-' + str(ir[1]) + \
                     '/phiper=' + str(param_set[("phiperiod")]) + \
                     '/nr-nphi=' + str(param_set[("nr")]) + '-' + str(param_set[("nphi")]) \
                     + '/H=' + str(param_set[("height")]) + '/' + ident  

    return p_path
    
def write_bender(E_b_surf, E_s_surf, normals, E_b, E_s, ue, param_set, in_file):
    
    folder = parse_path(param_set, param_set[("moldy_opm")])
    
    if not exists(folder):
        makedirs(folder)   
    
    copy2(in_file, folder + 'input.txt')

    
    # write params
    if not param_set[("moldy_opm")]:
        param_set["consts"] =   ue.consts
    param_set["E_b"]        =   E_b
    param_set["E_s"]        =   E_s

    with open(folder + 'params.txt', 'wb') as handle:
        pickle.dump(param_set, handle)
    #end write params
    
    #write surfaces
    umat_r, umat_phi, umat_z        =   ue.separate_mat(ue.ext_umat)
    init_r, init_phi, init_z        =   ue.separate_mat(ue.ext_surf)
    norm_r, norm_phi, norm_z        =   ue.separate_mat(normals)
    
    np.save(folder + 'umat_pol_r'   , umat_r) 
    np.save(folder + 'umat_pol_phi' , umat_phi) 
    np.save(folder + 'umat_pol_z'   , umat_z) 
    
    np.save(folder + 'init_pol_r'   , init_r) 
    np.save(folder + 'init_pol_phi' , init_phi) 
    np.save(folder + 'init_pol_z'   , init_z) 
    
    np.save(folder + 'norm_pol_r'   , norm_r) 
    np.save(folder + 'norm_pol_phi' , norm_phi) 
    np.save(folder + 'norm_pol_z'   , norm_z) 
    
    np.save(folder + 'E_b_surf'     , E_b_surf) 
    np.save(folder + 'E_s_surf'     , E_s_surf) 
    #end write surfaces


def write_synopsis(in_file):
    
    heights =   []
    E_bs    =   []
    E_ss    =   []
    consts  =   []
    
    params  =   parse_input(in_file)
    
    data_dir=   dirname(dirname(dirname( \
                    parse_path(params, params[("moldy_opm")]))))
    
    for sub_folder in listdir(data_dir):
        
        if sub_folder  != 'synopsis':
            
            param_dir   =   join(data_dir, sub_folder + '/opm_consts/')
            
            try:
                with open(param_dir + 'params.txt', 'rb') as handle:
                    params = pickle.loads(handle.read())
                    heights.append(params[("height")])
                    E_bs.append(params[("E_b")])
                    E_ss.append(params[("E_s")])
                    consts.append(params[("consts")])
            except IOError:
                print 'did not find params from ' + param_dir
            
    params_synopsis             =   {}
    
    params_synopsis["heights"]  =   heights
    params_synopsis["E_bs"]     =   E_bs
    params_synopsis["E_ss"]     =   E_ss
    params_synopsis["consts"]   =   consts
    
    if not exists(data_dir + '/synopsis'):
        makedirs(data_dir + '/synopsis')   
    
    with open(data_dir + '/synopsis/syn_params.txt', 'wb') as handle:
        pickle.dump(params_synopsis, handle)
        
    
    
def read_synopsis(in_file):
    
    params  =   parse_input(in_file)
    data_dir=   dirname(dirname(dirname( \
                    parse_path(params, params[("moldy_opm")])))) + '/synopsis/'
    
    with open(data_dir + 'syn_params.txt', 'rb') as handle:
        params = pickle.loads(handle.read())
    
    
    dtype   =   [('heights', float), ('E_bs', ndarray), ('E_ss', ndarray), \
                 ('consts', ndarray)]
    values  =   []
    
    for ih, h in enumerate(params[("heights")]):
        values.append((h, params[("E_bs")][ih], params[("E_ss")][ih], \
                       params[("consts")][ih]))  
    
    data    =   np.array(values, dtype=dtype)
    
    return np.sort(data, order='heights'), data_dir                        

    
            
def read_bender_output(input_file, read_for_moldy = False):
    
    from surface import parse_surf
    
    params      =   parse_input(input_file)
    
    if read_for_moldy:
        folder  =   parse_path(params, False)
    else:
        folder  =   parse_path(params, params[("moldy_opm")])

    
    #print folder
    
    with open(folder + 'params.txt', 'rb') as handle:
        params = pickle.loads(handle.read())
    
    umat_r      =   np.load(folder + 'umat_pol_r.npy') 
    umat_phi    =   np.load(folder + 'umat_pol_phi.npy') 
    umat_z      =   np.load(folder + 'umat_pol_z.npy') 
    
    i_surf_r    =   np.load(folder + 'init_pol_r.npy') 
    i_surf_phi  =   np.load(folder + 'init_pol_phi.npy') 
    i_surf_z    =   np.load(folder + 'init_pol_z.npy') 
    
    norm_r      =   np.load(folder + 'norm_pol_r.npy') 
    norm_phi    =   np.load(folder + 'norm_pol_phi.npy') 
    norm_z      =   np.load(folder + 'norm_pol_z.npy') 
    
    
    E_b_surf    =   np.load(folder + 'E_b_surf.npy') 
    E_s_surf    =   np.load(folder + 'E_s_surf.npy') 
    
    sym_op      =   select_option(params[("system")], \
                    params[("phiperiod")], params[("height")] / 2 / np.pi)[1]
    
    ext_umat, phys_umat, calc_umat \
                =   parse_surf(umat_r, umat_phi, umat_z, params[("phiperiod")], sym_op, 'umat')
    ext_surf, phys_surf, calc_surf \
                =   parse_surf(i_surf_r, i_surf_phi, i_surf_z, params[("phiperiod")], sym_op, 'surf')
    norm_surf   =   parse_surf(norm_r, norm_phi, norm_z, params[("phiperiod")], sym_op, 'surf')[0]
    
    
    return params, [ext_surf, phys_surf, calc_surf], \
                   [ext_umat, phys_umat, calc_umat], [E_b_surf, E_s_surf], norm_surf, folder  

def has_been_calculated(in_file, read_write = 'read'):

    params      =   parse_input(in_file)
    folder      =   parse_path(params, params[("moldy_opm")])
    
    
    if read_write == 'read':    
        try: 
            stat_file   =   open(folder + 'stat_file.txt', 'r')
            return stat_file.readline() in ['ready', 'error']
        except IOError:
            return False
    elif read_write == 'write_ok':
        stat_file       =   open(folder + 'stat_file.txt', 'w')
        stat_file.write('ready')
    elif read_write == 'write_not_ok':
        if not exists(folder):
            makedirs(folder) 
        stat_file       =   open(folder + 'stat_file.txt', 'w')
        stat_file.write('error')


def write_params(ue, energies, r, E_b, E_s):
    
    
    opt                             =   ue.opt
    H                               =   ue.hangle*2*np.pi
    nr, nphi                        =   ue.phys_umat.shape
    [rmin, rmax], [phimin, phimax]  =   valminmax(ue.phys_surf)
    sigma                           =   energies.sigma
    bend_module                     =   energies.bm
    strech_module                   =   energies.sm   
    
    
    path_r      =  path + 'energies/H-passivation/teor/opt=[%i,%i,%i]/' %(opt[0],opt[1], opt[2]) 
    
    # write u
    path_moldy  =   path_r + 'moldy/' 
    if not exists(path_moldy):
        makedirs(path_moldy)   
    umat_r, umat_phi, umat_z        =   ue.separate_mat(ue.ext_umat)
    init_r, init_phi, init_z        =   ue.separate_mat(ue.ext_surf)

    np.save(path_moldy + 'umat_pol_r_H=%.1f'   %H, umat_r) 
    np.save(path_moldy + 'umat_pol_phi_H=%.1f' %H, umat_phi) 
    np.save(path_moldy + 'umat_pol_z_H=%.1f'   %H, umat_z) 
    
    np.save(path_moldy + 'init_pol_r_H=%.1f'   %H, init_r) 
    np.save(path_moldy + 'init_pol_phi_H=%.1f' %H, init_phi) 
    np.save(path_moldy + 'init_pol_z_H=%.1f'   %H, init_z) 
    #
    
    def consts_string(consts):
        string_consts = ''
        for c in consts:
            string_consts += str(c) + ' '
        
        return string_consts[:-1]
    
    if not exists(path_r):
        makedirs(path_r)   
    
    if not exists(path_r + 'modules.data'):
        file_modules   =  open(path_r + 'modules.data', 'w')
        file_modules.write('# this is modules file: bend_module, strech_module, sigma \n' )
        file_modules.write(str(bend_module) + ' ' + str(strech_module) + ' ' + str(sigma))
        
    else:
        bend_module_read, strech_module_read, sigma_read = read_modules(opt)
        if bend_module_read != bend_module:
            print 'bend_modules differ we do not write!'
            return
        
        if strech_module_read != strech_module:
            print 'strech_modules differ we do not write!'
            return
        
        if sigma_read != sigma:
            print 'sigmas differ we do not write!'
            return
           
    if not exists(path_r + '%i-%i.data' %(r[0],r[1])):
        file_r                  =   open(path_r + '%i-%i.data' %(r[0],r[1]), 'a') 
        file_r.write('# Height, consts \n')
    else:
        file_r                  =   open(path_r + '%i-%i.data' %(r[0],r[1]), 'a')
    
    file_r.write(str(H) + ' [' + consts_string(ue.consts) + '] ' + 'ebs: ' +  str(E_b) + ' ' + str(E_s) \
                 + ' rs: ' + str(nr) + ' ' + str(rmin) + ' ' + str(rmax) 
                 + ' phis: ' + str(nphi) + ' ' + str(phimin) + ' ' + str(phimax) + ' \n')
    file_r.close()
    
def write_total_synopsis(nr, nphi, phiperiod, system, input_file_folder_consts):
    
    syn_dict        =   {}
    
    
    for input_folder in listdir(input_file_folder_consts):
        input_file_folder   = input_file_folder_consts + input_folder + '/'
        
        if input_folder != 'store':
            
            for input_file in listdir(input_file_folder):
                params      =   parse_input(input_file_folder + input_file)
                folder      =   parse_path(params, params[("moldy_opm")])    
                
                if input_file[-3:] == 'txt':
                    
                    try:
                        with open(folder + 'params.txt', 'rb') as handle:
                            params = pickle.loads(handle.read())
                        
                        syn_dict["%i-%i-h=%.2f" %(params[("ir")][0], \
                                           params[("ir")][1], params[("height")])]   = params
                        #print folder
                        
                    except IOError:
                        continue #print 'not ready yet, = ' + str(params[("ir")])    
                
    
    syn_path  = path + 'bender_output/' + system + '/synopsis_phiper=%.2f/' %phiperiod
    
    if not exists(syn_path):
        makedirs(syn_path)
    
    with open(syn_path + 'syn_params.txt', 'wb') as handle:
        pickle.dump(syn_dict, handle) 

def read_total_synopsis(nr, nphi, phiperiod, system):
    
    syn_path  = path + 'bender_output/' + system + '/synopsis_phiper=%.2f/' %phiperiod
    
    with open(syn_path + 'syn_params.txt', 'rb') as handle:
        syn_dict = pickle.loads(handle.read())
    
    return syn_dict, syn_path
    
def read_modules(opt):
    
    path_r      =   path + 'energies/H-passivation/teor/opt=[%i,%i,%i]/' %(opt[0],opt[1], opt[2]) 
    file_r      =   open(path_r + 'modules.data', 'r')
    
    nextline    =   file_r.readline() 
        
    while nextline != '':
        if nextline[0] != '#':
            modules     =   nextline.split(' ')
            return float(modules[0]), float(modules[1]), float(modules[2])
        nextline        = file_r.readline()
        
def read_params(r, H):
    

    path_r          =   path + 'energies/H-passivation/teor/' 
    file_r          =   open(path_r + '%i-%i.data' %(r[0],r[1]), 'r') 
        
    nextline        =   file_r.readline()
    n               =   0
    
    while nextline != '':
        data                =   nextline.split(' ')
        if nextline[0]     ==   str(H):
            consts          =   float(data[0]), float(data[1]), float(data[2])
            #rmin, rmax      =   float(data[2]), float(data[3])
            #drmin, drmax    =   float(data[4]), float(data[5])   
            return consts
        n                  +=1    
        
        nextline            =   file_r.readline()
        
    file_r.close()
        
 
def read_rad(ir):
    
    
    with open(path + 'opm_structures/H-passivation/radiuces.dict', 'rb') as handle:
        radiuces = pickle.loads(handle.read())
    
    try:
        rad = radiuces[("%i-%i" %(ir[0], ir[1]))][0]
    except KeyError as e:
        print e
    
    return rad[0], rad[1]
    
def generate_input(system_set):
    
    system              =   system_set[0]
    phimin, phiperiod   =   0, system_set[2]
    [nr, nphi]          =   system_set[1]
    if system_set[3] == 'consts':
        moldy_opm           =   'false'
    elif system_set[3] == 'moldy':
        moldy_opm           =   'True'
    
    irs                 =   [[1,6],[2,6],[4,9],[6,12],[6,17],\
                             [7,14],[8,12],[8,16],[8,18],[9,19],[11,20]] # these are set of unsuccesfull hotbit calculations
    
    for input_file in listdir(path + 'opm_structures/H-passivation/'):
        if '-' in input_file:
            ir          =   [int(input_file.split('-')[0]), int(input_file.split('-')[1]) ]
            
            if ir in irs: 
            
                if query_yes_no('make input from %i-%i, moldy %s' %(ir[0], ir[1], moldy_opm), 'yes'):
                    rmin, rmax      =   read_rad(ir)
                    try:
                        heights         =   read_energies(ir, passivated = True, k_z = 6, optimized = True)[0][1:]
                        n               =   max(1, int(len(heights)/15))
                    
                        for h in heights[::n]:
                            
                            folder = path + 'bender_input/calc/%s/%s/phiper=%.2f/nr-nphi=%i-%i/%i-%i/' \
                            %(system_set[3], system, phiperiod, nr, nphi, ir[0], ir[1])
                            
                            if not exists(folder):
                                makedirs(folder)
                            
                            input_file  =   open(folder + '/in_%i-%i_H=%.1f.txt' \
                                                  %(ir[0],ir[1],h), 'w')
                            
                            input_file.write(system + ' # system \n')
                            input_file.write(str(ir[0]) + ' ' + str(ir[1]) + ' # ir[0], ir[1] \n')
                            input_file.write(str(nr) + ' ' + str(nphi) + ' # nr, nphi \n')
                            input_file.write('%.2f %.2f # rmin, rmax \n' %(rmin, rmax))
                            input_file.write('%.2f %.2f # phimin, phiperiod \n' %(phimin, phiperiod))
                            input_file.write('%.2f # height \n' %h)
                            input_file.write('0 # consts \n')
                            input_file.write(moldy_opm + ' # moldy opm')
                            
                            input_file.close()
                    
                    except IOError:
                        print 'Has not been calculated'
                        print 
                    
 
def read_energies(r, passivated = True, k_z = 6, optimized = True):
    
    nop = ''
    if not optimized:
        nop = '_nop'
        
    if passivated:     
        epath = path + 'energies/H-passivation/kz=%i%s/' %(k_z, nop)
        H = '-H'
    else: 
        epath = path + 'energies/no-passivation/kz=%i%s/' %(k_z, nop)
        H = ''
        
    #print '%i-%i%s_kz=%i%s.data' %(r[0],r[1],H,k_z, nop)     
    e_file = open(epath + '%i-%i%s_kz=%i%s.data' %(r[0],r[1],H,k_z, nop), 'r')        
    
    nextline = e_file.readline()
    lengths_r = []
    energies_r = []
    streches_r = []
    
    while nextline != '':
        if nextline[0] != '#':
            e = nextline.split(' ')
            if len(e) == 3:
                lengths_r.append(float(e[0]))
                energies_r.append(float(e[1]))
                streches_r.append(0)
            elif len(e) == 4:
                lengths_r.append(float(e[0]))
                streches_r.append(float(e[1]))
                energies_r.append(float(e[2]))
                
        nextline = e_file.readline()
    
    e_file.close()
    
    energies    = np.zeros(len(lengths_r))
    lengths     = np.zeros(len(lengths_r))
    streches    = np.zeros(len(lengths_r))

    for i in range(len(lengths_r)):
        energies[i] = energies_r[i]
        lengths[i]  = lengths_r[i]
        streches[i] = streches_r[i]
        
    return lengths, energies, streches

def read_h_consts(ir, opt): 
    
    path_r                  =   path + 'energies/H-passivation/teor/opt=[%i,%i,%i]/' %(opt[0],opt[1], opt[2]) 
           
    file_r                  =   open(path_r + '%i-%i.data' %(ir[0],ir[1]), 'r') 
    
    nextline                =   file_r.readline()
    
    heights                 =   []
    consts                  =   []
    energies                =   []
    rdata                   =   []
    phidata                 =   []
    
    def parse_consts(line):
        cs  = []
        for c in line.split(' '):
            cs.append(float(c))
        return cs
    
    while nextline != '':
        if nextline[0] != '#':
            heights.append(nextline.split(' ')[0])
            
            consts.append(parse_consts(nextline.split(']')[0].split('[')[1]))
            
            e_str           =   nextline.split('ebs:')[1][1:].split(' ')
            
            energies.append([float(e_str[0]), float(e_str[1])])
            
            rdata_s         =   nextline.split('rs:')[1][1:].split(' ')[:3]
            phidata_s       =   nextline.split('phis:')[1][1:].split(' ')[:3]
            
            
            rdata.append([float(rdata_s[0]), float(rdata_s[1]),float(rdata_s[2])])
            phidata.append([float(phidata_s[0]), float(phidata_s[1]),float(phidata_s[2])])

        nextline    = file_r.readline(
                                      )
    file_r.close()

    dtype   =   [('height', float), ('consts', ndarray), ('energies', ndarray), \
                 ('rs', ndarray), ('phis', ndarray)]
    values  =   []
    
    
    for ih, h in enumerate(heights):
        values.append((h, consts[ih], energies[ih], rdata[ih], phidata[ih]))  
    
    data    =   np.array(values, dtype=dtype)
    
    return np.sort(data, order='height')                        


def read_umat(ir, opt, H):
    
    path_f      =   path + '/energies/H-passivation/teor/opt=[%i,%i,%i]/%i-%i/moldy/' \
                        %(opt[0],opt[1], opt[2], ir[0], ir[1])

    path_r      =   path_f + 'umat_pol_r_H=%.1f.npy'      %H 
    path_phi    =   path_f + 'umat_pol_phi_H=%.1f.npy'    %H 
    path_z      =   path_f + 'umat_pol_z_H=%.1f.npy'      %H 
    
    umat_r      =   np.load(path_r)
    umat_phi    =   np.load(path_phi)
    umat_z      =   np.load(path_z)
    
    umat        =   np.empty(umat_r.shape, dtype = 'object')
    
    for index, val in np.ndenumerate(umat):
        umat[index]     =   np.array([umat_r[index], umat_phi[index], umat_z[index]])
    
    return umat   
    
#system_sets     = get_sys_sets()
#for system_set in system_sets:
#    generate_input(system_set)     
    
    
