# Script to load, convert and save phase space data between different formats
# Can load particle data from:
#   - Impact-T phase space output files `fort.*`
#   - Plain text files in TraceWin output format `*.txt`
# Can output files as:
#   - TraceWin/PlotWin distribution file `*.dst`

import math
import numpy
import scipy.constants
import struct
import random

def get_parameters():
    """Set up the parameters needed to load, convert and write data"""
    parameters = {}
    parameters = get_file_lists(parameters)
    parameters = get_beam_parameters(parameters)
    return parameters

def get_file_lists(parameters):
    input_files = get_list('Enter filename(s) to load: ')
    if input('Do you want to merge with other input file(s)? [y/n]') == 'y':
        parameters['infiles'] = []
        parameters['infiles'].append(input_files)
        while True:
            input_files = get_list('Enter filename(s) to merge: ')
            parameters['infiles'].append(input_files)
            if input('Do you want to merge with more file(s)? [y/n]') != 'y':
                break
        parameters['merge'] = get_merge_list(parameters['infiles'])
    else:
        parameters['infiles'] = input_files
    # parameters['infiles'] = ['fort.48']
    # parameters['infiles'] = ['part_100MeV.txt']
    # parameters['infiles'] = ['fort.40', 'fort.41', 'fort.42', 'fort.43',
    #                          'fort.44', 'fort.45', 'fort.46', 'fort.47',
    #                          'fort.48']
    # parameters['infiles'] = [['part_100MeV.txt'], ['fort.48']]
    # parameters['infiles'] = [['part_100MeV.txt'],
    #                          ['fort.40', 'fort.41', 'fort.42', 'fort.43',
    #                           'fort.44', 'fort.45', 'fort.46', 'fort.47',
    #                           'fort.48']]
    # parameters['merge'] = {
    #     'x': 1, 'x`': 1, 'px': 1,
    #     'y': 1, 'y`': 1, 'py': 1,
    #     'z': 2, 'z`': 2, 'pz': 2,
    #     'phi': 2, 'gamma': 2, 'W': 2}
    parameters['outfiles'] = get_list('Enter filename(s) for output: ')
    # parameters['outfiles'] = ['particles.dst']
    return parameters

def get_list(prompt):
    print(prompt)
    print('You can enter multiple values. '
        'Press [Enter] on a blank line to finish.')
    input_list = []
    while True:
        input_value = input()
        if input_value:
            input_list.append(input_value)
        else:
            return input_list

def get_merge_list(file_lists):
    for idx, file_list in enumerate(file_lists):
        print(f'List {idx + 1}: {file_list}')
    merge_list = {}
    merge_list['x'] = get_merge_value('x', len(file_lists))
    merge_list['x`'] = merge_list['x']
    merge_list['px'] = merge_list['x']
    merge_list['y'] = get_merge_value('y', len(file_lists))
    merge_list['y`'] = merge_list['y']
    merge_list['py'] = merge_list['y']
    merge_list['z'] = get_merge_value('z', len(file_lists))
    merge_list['z`'] = merge_list['z']
    merge_list['pz'] = merge_list['z']
    merge_list['phi'] = merge_list['z']
    merge_list['gamma'] = merge_list['z']
    merge_list['W'] = merge_list['z']
    return merge_list

def get_merge_value(dim, max):
    while True:
        value = int(input(f'Which list should the {dim} values come from? '))
        if value > 0 and value <= max:
            return value
        else:
            print(f'Please enter an integer value from 1 to {max}')

def get_beam_parameters(parameters):
    if contains_text_file(parameters['infiles']):
        if input('Use beam parameters from the text file? [y/n]') == 'y':
            return parameters
    parameters['mass'] = float(input('Input the mass of the particles in MeV/c2: '))*1e6
    parameters['current'] = float(input('Input the beam current in mA: '))*1e-3
    parameters['frequency'] = float(input('Input the frequency in MHz: '))*1e6
    parameters['energy'] = float(input('Input the beam energy in MeV: '))*1e6
    # parameters['mass'] = 938.272*1e6
    # parameters['current'] = 0.0
    # parameters['frequency'] = 1.0
    # parameters['energy'] = 3.74944*1e6
    return parameters

def contains_text_file(infiles):
    """Check whether the input file includes a text file"""
    if isinstance(infiles[0], list):
        for file_list in infiles:
            if list_contains_text_file(file_list):
                return True
    else:
        if list_contains_text_file(infiles):
            return True
    return False

def list_contains_text_file(file_list):
    """If there is a txt file, may not need to enter beam parameters manually"""
    if 'txt' in [infile.split('.')[-1] for infile in file_list]:
        return True
    else:
        return False

def load_and_merge_data(parameters):
    """If more than one list of files is given, merge the data"""
    if isinstance(parameters['infiles'][0], list):
        data_list = []
        for file_list in parameters['infiles']:
            data_list.append(load_data(file_list))
        return merge_data(data_list, parameters)
    else:
        return load_data(parameters['infiles'])

def load_data(file_list):
    """Go through the list of input files and load the data"""
    data = []
    for infile in file_list:
        if infile.split('.')[0] == 'fort':
            data = data + load_impact_t_data(infile)
        elif infile.split('.')[-1] == 'txt':
            data = data + load_tracewin_text_data(infile)
        else:
            raise ValueError(f'Unrecognised input file format: {infile}')
    return data

def load_impact_t_data(infile):
    """Load data from an Impact-T output file such as `fort.40`"""
    columns = ['x', 'px', 'y', 'py', 'z', 'pz']
    print(f'Loading data from {infile}')
    with open(infile) as f:
        contents = f.readlines()
    input_data = [list(map(float, line.split())) for line in contents]
    dict_data = [dict(zip(columns, line)) for line in input_data]
    processed_data = process_impact_t_data(dict_data, parameters)
    return processed_data

def load_tracewin_text_data(infile):
    """Load data from an Impact-T output file such as `fort.40`"""
    columns = ['x', 'x`', 'y', 'y`', 'z', 'z`', 'phi', 't', 'W', 'loss']
    print(f'Loading data from {infile}')
    with open(infile) as f:
        f.readline()
        headers = list(map(float, f.readline().split()))
        npt = int(headers[0])
        parameters['mass'] = headers[1]*1e6
        parameters['energy'] = headers[2]*1e6
        parameters['frequency'] = headers[3]*1e6
        parameters['current'] = headers[4]
        parameters['charge'] = headers[5]
        f.readline()
        contents = f.readlines()
    if len(contents) != npt:
        raise ValueError(
            f'Number of rows ({len(contents)}) does not match '
            f'given number of particles ({npt})')
    input_data = [list(map(float, line.split())) for line in contents]
    dict_data = [dict(zip(columns, line)) for line in input_data]
    processed_data = process_tracewin_text_data(dict_data, parameters)
    return processed_data

def process_impact_t_data(data, parameters):
    """Calculate angles, phase and energy and add to dataset"""
    npt = len(data)
    z0 = sum([particle['z'] for particle in data])/npt
    pz0 = sum([particle['pz'] for particle in data])/npt
    for particle in data:
        particle['x`'] = particle['px'] / particle['pz']
        particle['y`'] = particle['py'] / particle['pz']
        particle['z`'] = (particle['pz'] - pz0)/particle['pz']
        particle['phi'] = - (
            (particle['z'] - z0)
            *2*scipy.constants.pi
            *parameters['frequency']
            /(pz0*scipy.constants.c))
        particle['gamma'] = 1/math.sqrt(1 - particle['pz']**2)
        particle['W'] = (particle['gamma'] - 1)*parameters['mass']
    return data

def process_tracewin_text_data(data, parameters):
    """Adjust units and calculate momenta from TraceWin text file data"""
    for particle in data:
        particle['x'] = particle['x']*1e-3
        particle['x`'] = particle['x`']*1e-3
        particle['y'] = particle['y']*1e-3
        particle['y`'] = particle['y`']*1e-3
        particle['z'] = particle['z']*1e-3
        particle['z`'] = particle['z`']*1e-3
        particle['phi'] = particle['phi']*2*math.pi/360
        particle['W'] = particle['W']*1e6
        particle['loss'] = int(particle['loss'])
        particle['gamma'] = particle['W']/parameters['mass'] + 1
        particle['pz'] = math.sqrt(1 - 1/particle['gamma']**2)
        particle['px'] = particle['x`']*particle['pz']
        particle['py'] = particle['y`']*particle['pz']
    return data

def merge_data(data_list, parameters):
    """Selectively merge the data from multiple lists into one dataset"""
    print('Merging data')
    data_list = reduce_to_smallest_npt(data_list)
    data = [{item: data_list[parameters['merge'][item]-1][i][item]
        for item in parameters['merge']}
        for i in range(0, len(data_list[0]))]
    if abs(get_energy_at_peak(data) - parameters['energy']) > 0.005:
        data = scale_energy(data, parameters['energy'])
    return data

def reduce_to_smallest_npt(data_list):
    """Reduce the size of all lists to match the smallest list"""
    npt = min([len(data) for data in data_list])
    return [random.sample(data, npt) for data in data_list]

def get_energy_at_peak(data):
    """Currently only uses the median energy, which is fast but not accurate"""
    return numpy.median([particle['W'] for particle in data])

def scale_energy(data, new_peak):
    """Scale the energy distribution so that peak is at given energy value"""
    old_peak = get_energy_at_peak(data)
    for particle in data:
        particle['W'] = particle['W'] * new_peak / old_peak
    return data

def write_data(data, parameters):
    """Go through the list of output files and save the data"""
    for outfile in parameters['outfiles']:
        if outfile.split('.')[-1] == 'dst':
            write_dst_file(data, outfile, parameters)
        else:
            raise ValueError(f'Unrecognised output file format: {outfile}')

def write_dst_file(data, outfile, parameters):
    """Write data as a TraceWin `.dst` file"""
    print(f'Writing data to {outfile}')
    npt = len(data)
    with open(outfile, 'w+b') as f:
        f.write(b'  ')
        f.write(struct.pack('i', npt))
        f.write(struct.pack('d', parameters['current']*1e3))
        f.write(struct.pack('d', parameters['frequency']*1e-6))
        f.write(b' ')
        for particle in data:
            f.write(struct.pack('6d',
                                particle['x']*100,
                                particle['x`'],
                                particle['y']*100,
                                particle['y`'],
                                particle['phi'],
                                particle['W']*1e-6))
        f.write(struct.pack('d', parameters['mass']*1e-6))

# What to do when run as a script
if __name__ == '__main__':
    parameters = get_parameters()
    data = load_and_merge_data(parameters)
    write_data(data, parameters)
