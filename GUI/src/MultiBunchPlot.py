import matplotlib.pyplot
import numpy
import sys
import pathlib

def plot_all(bunch_count):
    """Run and save all plots consecutively."""
    print('Loading experimental data...')
    try:
        experimental_results = load_experimental_results('experimental_data.txt')
    except FileNotFoundError:
        print('File not found: experimental_data.txt. '
              'Continuing without experimental data.')
        experimental_results = None
    print('Loading statistical data...')
    xdata, ydata = load_statistics_data(bunch_count)
    combined_xdata = combine_bunch_values(xdata)
    combined_ydata = combine_bunch_values(ydata)
    # Beam size
    print('Plotting beam size...')
    figure, axes = matplotlib.pyplot.subplots(dpi=300)
    plot_beam_size(axes, xdata, [], combined_xdata)
    figure.savefig('beam-size')
    if experimental_results:
        figure, axes = matplotlib.pyplot.subplots(dpi=300)
        plot_beam_size(axes, xdata, experimental_results, combined_xdata)
        figure.savefig('beam-size-vs-experiment')
    # Emittance
    print('Plotting emittance...')
    figure, axes = matplotlib.pyplot.subplots(dpi=300)
    plot_emittance(axes, xdata, ydata, combined_xdata, combined_ydata)
    figure.savefig('emittance')
    figure, axes = matplotlib.pyplot.subplots(dpi=300)
    plot_emittance_growth(axes, xdata, ydata, combined_xdata, combined_ydata)
    figure.savefig('emittance-growth')
    # Particle plots: initial
    print('Loading initial phase space data...')
    phase_space_data = load_phase_space_data(40, bunch_count)
    print('Plotting initial phase space data...')
    figure, axes = matplotlib.pyplot.subplots(nrows=2, ncols=2, dpi=300)
    plot_phase_spaces(axes, phase_space_data, title='Initial phase space',
                      bunch_count=bunch_count, grid_size=300)
    figure.savefig('phase-space-initial')
    # Particle plots: final
    print('Loading final phase space data...')
    phase_space_data = load_phase_space_data(50, bunch_count)
    print('Plotting final phase space data...')
    figure, axes = matplotlib.pyplot.subplots(nrows=2, ncols=2, dpi=300)
    plot_phase_spaces(axes, phase_space_data, title='Final phase space',
                      bunch_count=bunch_count, grid_size=300)
    figure.savefig('phase-space-final')
    # Particle plots: BPMs
    lattice = get_lattice()
    bpm_list = get_bpms(lattice)
    for location, filenumber in bpm_list:
        print(f'Loading BPM {filenumber} phase space data...')
        phase_space_data = load_phase_space_data(filenumber, bunch_count)
        print(f'Plotting BPM {filenumber} phase space data...')
        figure, axes = matplotlib.pyplot.subplots(nrows=2, ncols=2, dpi=300)
        plot_phase_spaces(axes, phase_space_data,
                          title=f'Phase space at z = {location}',
                          bunch_count=bunch_count, grid_size=300)
        figure.savefig(f'phase-space-{filenumber}')

def get_input_filename(bunch):
    """Return the filename of the input file for a particular bunch."""
    if bunch == 1:
        filename = 'ImpactT.in'
    else:
        filename = f'ImpactT{bunch}.in'
    if pathlib.Path(filename+'.rendered').is_file():
        filename += '.rendered'
    return filename

def get_bunch_count():
    """Get the number of bunches from the first input file."""
    input = read_input_file(get_input_filename(1))
    return int(input[1].split()[2])

def get_lattice():
    """Read the lattice from the first input file."""
    input = read_input_file(get_input_filename(1))
    return [line.split() for line in input[9:]]

def get_bpms(lattice):
    """Return the location and file number of all BPMs in the lattice."""
    return [(str(float(elem[4])*1000) + ' mm', int(elem[2]))
            for elem in lattice if elem[3]=='-2']

def get_bunch_counts(bunch_count):
    """Return a list of the particle counts for each bunch."""
    input_filenames = [get_input_filename(i+1) for i in range(bunch_count)]
    return [get_particle_count(filename) for filename in input_filenames]

def get_particle_count(filename):
    """Get the macroparticle count from a given input file."""
    input = read_input_file(filename)
    return int(input[2].split()[1])

def get_mass(filename):
    """Get the particle mass from a given input file."""
    input = read_input_file(filename)
    return float(input[8].split()[2])

def is_mass_matched(bunch_count):
    """Check whether the particle mass values are the same for given bunches."""
    mass = [get_mass(get_input_filename(i+1)) for i in range(bunch_count)]
    return len(set(mass)) == 1

def read_input_file(filename):
    """Read input file and return list of input lines."""
    with open(filename, 'r') as f:
        return [line.strip() for line in f.readlines() if line[0] != '!']

def load_experimental_results(filename):
    """Return z, rms beam size and error values from experimental data (m)."""
    with open(filename, 'r') as f:
        data = [line.split() for line in f.readlines()[1:]]
    return [[float(datum) for datum in row] for row in data]

def load_statistics_data(bunch_count):
    """Load simulation data per bunch from statistics files."""
    xdata = []
    ydata = []
    for i in range(bunch_count):
        with open(f'fort.{i+1}024', 'r') as f:
            xdata.append([line.split() for line in f.readlines()])
        with open(f'fort.{i+1}025', 'r') as f:
            ydata.append([line.split() for line in f.readlines()])
    return xdata, ydata

def load_phase_space_data(filenumber, bunch_count):
    """Load phase space data per bunch from particle output files."""
    data = numpy.array(numpy.loadtxt(f'fort.{filenumber}'))
    for i in range(1, bunch_count):
        data = numpy.concatenate((data, numpy.loadtxt(f'fort.{filenumber+i}')))
    return data

def combine_bunch_values(data_in):
    """Combine values of separate bunches into a single summary dataset."""
    if len(data_in) == 1:
        return data_in[0]
    data = numpy.array(data_in, dtype=float)
    npt = numpy.array(get_bunch_counts(len(data)))
    # Decompose data
    t_data = data[:,:,0]
    z0_data = data[:,:,1]
    x0_data = data[:,:,2]
    xrms_data = data[:,:,3]
    px0_data = data[:,:,4]
    pxrms_data = data[:,:,5]
    xpx_data = -data[:,:,6]
    epx_data = data[:,:,7]
    # All t data should be the same
    t = t_data[0]
    if not all([row == t.tolist() for row in t_data.tolist()]):
        raise ValueError('Time step values not in sync across bunches.')
    # Centroids can be combined by a weighted mean
    z0 = numpy.sum(z0_data.T.dot(numpy.diag(npt)).T, 0)/npt.sum()
    x0 = numpy.sum(x0_data.T.dot(numpy.diag(npt)).T, 0)/npt.sum()
    px0 = numpy.sum(px0_data.T.dot(numpy.diag(npt)).T, 0)/npt.sum()
    # To combine rms values we need the weighted sum of squares
    sumsqx = numpy.sum(numpy.square(xrms_data).T.dot(numpy.diag(npt)).T
                       + numpy.square(x0_data).T.dot(numpy.diag(npt)).T, 0)
    sumsqpx = numpy.sum(numpy.square(pxrms_data).T.dot(numpy.diag(npt)).T
                        + numpy.square(px0_data).T.dot(numpy.diag(npt)).T, 0)
    sumxpx = numpy.sum((xpx_data + x0_data*px0_data).T.dot(numpy.diag(npt)).T, 0)
    xrms = numpy.sqrt(sumsqx/npt.sum() - numpy.square(x0))
    pxrms = numpy.sqrt(sumsqpx/npt.sum() - numpy.square(px0))
    xpx = sumxpx/npt.sum() - (x0 * px0)
    epx = numpy.sqrt(xrms*xrms * pxrms*pxrms - xpx*xpx)
    # Return combined data as a standard list
    return numpy.array([t, z0, x0, xrms, px0, pxrms, -xpx, epx]).T.tolist()

def plot_beam_size(axes, data, experiment_data=[], combined_data=[]):
    """Create a plot of beam size per bunch."""
    axes.set_xlabel('z-location (mm)')
    axes.set_ylabel('Beam size (mm)')
    axes.set_title('Beam size')
    if len(experiment_data) > 0:
        plot_beam_size_experimental(axes, experiment_data)
    for i in range(len(data)):
        plot_beam_size_single(axes, data[i], '--', label=f'Bunch {i+1} rms')
    if len(combined_data) == 0:
        combined_data = combine_bunch_values(data)
    plot_beam_size_single(axes, combined_data, 'r-', label=f'Combined rms')
    axes.legend(fontsize='x-small')

def plot_beam_size_experimental(axes, data):
    """Plot experimental data points with error bars"""
    z = [float(row[0])*1.0e3 for row in data]
    rms = [float(row[1])*1.0e3 for row in data]
    error = [float(row[2])*1.0e3 for row in data]
    axes.errorbar(z, rms, yerr=error,
                  fmt='ko', markersize=2.0, elinewidth=0.5, capsize=1.0,
                  label='Experimental rms')

def plot_beam_size_single(axes, data, fmt, label):
    """Plot the rms beam size for the given data onto the given axes."""
    z = [float(row[1])*1.0e3 for row in data]
    rms = [float(row[3])*1.0e3 for row in data]
    axes.plot(z, rms, fmt, linewidth=1, label=label)

def plot_emittance(axes, xdata, ydata, combined_xdata=[], combined_ydata=[]):
    """Create a plot of average x and y emittance per bunch."""
    axes.set_xlabel('Time (ns)')
    axes.set_ylabel('Normalised rms emittance (π mm mrad)')
    axes.set_title('Emittance')
    for i in range(len(xdata)):
        plot_emittance_single(axes, xdata[i], ydata[i],
                              '--', label=f'Bunch {i+1}')
    if len(combined_xdata) == 0:
        combined_xdata = combine_bunch_values(xdata)
    if len(combined_ydata) == 0:
        combined_ydata = combine_bunch_values(ydata)
    plot_emittance_single(axes, combined_xdata, combined_ydata,
                          'r-', label=f'Combined')
    axes.legend(fontsize='x-small')

def plot_emittance_single(axes, xdata, ydata, fmt, label):
    """Plot the average x and y emittance onto the given axes."""
    t = [float(row[0])*1.0e9 for row in xdata]
    e = [(float(xdata[k][7]) + float(ydata[k][7])) / 2 * 1.0e6
         for k in range(len(xdata))]
    axes.plot(t, e, fmt, linewidth=1, label=label)

def plot_emittance_growth(axes, xdata, ydata,
                          combined_xdata=[], combined_ydata=[]):
    """Create a plot of average x and y emittance growth per bunch."""
    axes.set_xlabel('Time (ns)')
    axes.set_ylabel('Average emittance growth in x and y (relative)')
    axes.set_title('Emittance growth')
    for i in range(len(xdata)):
        plot_emittance_growth_single(axes, xdata[i], ydata[i],
                              '--', label=f'Bunch {i+1}')
    if len(combined_xdata) == 0:
        combined_xdata = combine_bunch_values(xdata)
    if len(combined_ydata) == 0:
        combined_ydata = combine_bunch_values(ydata)
    plot_emittance_growth_single(axes, combined_xdata, combined_ydata,
                                 'r-', label=f'Combined')
    axes.legend(fontsize='x-small')

def plot_emittance_growth_single(axes, xdata, ydata, fmt, label):
    """Plot the average x and y emittance growth onto the given axes."""
    t = [float(row[0])*1.0e9 for row in xdata]
    e = [(float(xdata[k][7]) + float(ydata[k][7])) / 2 * 1.0e6
         for k in range(len(xdata))]
    initial_emittance = max(e[0], 1.0e-10)
    growth = [emittance/initial_emittance - 1 for emittance in e]
    axes.plot(t, growth, fmt, linewidth=1, label=label)

def plot_phase_space(axes, x, y, xlabel, ylabel, grid_size=100):
    """Plot a single phase space onto the given axes."""
    if grid_size < 10:
        grid_size = 10
    colour_map = matplotlib.cm.get_cmap('jet')
    colour_map.set_under('white', 0.)
    # Plot 2D histogram
    hist, *_ = axes.hist2d(x, y, bins=grid_size, cmap=colour_map, cmin=1)
    axes.set_xlabel(xlabel, fontsize='x-small')
    axes.set_ylabel(ylabel, fontsize='x-small')
    axes.tick_params(labelsize='xx-small')
    # Project to 1D histograms
    xmin, xmax = axes.get_xbound()
    ymin, ymax = axes.get_ybound()
    xscale = numpy.array(range(grid_size)) / grid_size * (xmax - xmin) + xmin
    yscale = numpy.array(range(grid_size)) / grid_size * (ymax - ymin) + ymin
    xhist = numpy.nansum(hist, 1)
    yhist = numpy.nansum(hist, 0)
    xhist_scaled = xhist / xhist.max() * (ymax - ymin) * 0.2 + ymin
    yhist_scaled = yhist / yhist.max() * (xmax - xmin) * 0.2 + xmin
    # Plot 1D histograms
    axes.plot(xscale, xhist_scaled, color='green', linewidth=0.75)
    axes.plot(yhist_scaled, yscale, color='green', linewidth=0.75)

def plot_phase_spaces(axes, data, title=None, bunch_count=None, grid_size=100):
    """Plot four phase spaces onto the given array of axes."""
    if not title:
        title = 'Phase space'
    figure = axes[0,0].figure
    figure.suptitle(title)
    x = data.T[0]
    px = data.T[1]
    y = data.T[2]
    py = data.T[3]
    z = data.T[4]
    pz = data.T[5]
    xp = px/pz
    yp = py/pz
    plot_phase_space(axes[0,0], x*1e3, xp*1e3, 'x (mm)', 'x` (mrad)', grid_size)
    plot_phase_space(axes[0,1], y*1e3, yp*1e3, 'y (mm)', 'y` (mrad)', grid_size)
    plot_phase_space(axes[1,1], x*1e3, y*1e3, 'x (mm)', 'y (mm)', grid_size)
    if not bunch_count:
        bunch_count = get_bunch_count()
    if not is_mass_matched(bunch_count):
        print(f'Masses for {bunch_count} bunches do not match, '
              + 'so cannot plot combined energy distribution. '
              + 'Falling back to dimensionless momentum plot.')
        plot_phase_space(
            axes[1,0], z*1e3, pz, 'z (mm)', 'pz (dimensionless βγ)', grid_size)
    else:
        mass = get_mass(get_input_filename(1))
        gamma = 1/numpy.sqrt(1 - numpy.square(pz))
        W = (gamma - 1)*mass
        plot_phase_space(
            axes[1,0], z*1e3, W/1e6, 'z (mm)', 'Energy (MeV)', grid_size)

if __name__ == '__main__':
    matplotlib.use('agg') # Use the AGG renderer to produce PNG output
    if len(sys.argv) > 1:
        bunch_count = int(sys.argv[1])
    else:
        bunch_count = get_bunch_count()
    plot_all(bunch_count)
