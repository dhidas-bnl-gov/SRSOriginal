from matplotlib.colors import LogNorm
import numpy as np
import matplotlib.pyplot as plt


def write_power_density_csv (P, fileName) :

    x = []
    y = []
    z = []

    x = [item[0][0] for item in P]
    y = [item[0][1] for item in P]
    z = [item[1]    for item in P]

    x2 = list(set(x))
    y2 = list(set(y))

    nx = len(x2)
    ny = len(y2)

    x2.sort()
    y2.sort()

    with open(fileName, 'w') as f:

        f.write(' ,' + ', '.join(map(str, x2)) + '\n')
        for i in xrange(ny):
            f.write(str(y2[i]) + ', ')
            for j in xrange(nx):
                if j < nx - 1:
                    f.write(str(z[i * nx + j]) + ', ')
                else:
                    f.write(str(z[i * nx + j]) + '\n')




def add_power_densities(A, B):
    """Add two power density lists assuming same mesh order"""
    
    new_list = []
    
    for i in range(len(A)):
        new_list.append([A[i][0], A[i][1]+B[i][1]])
        
    return new_list



def plot_trajectory_position(trajectory, show=True, ofile=''):
    """Plot the trajectory"""

    # Get coordinate lists
    X  = [item[0][0] for item in trajectory]
    Y  = [item[0][1] for item in trajectory]
    Z  = [item[0][2] for item in trajectory]

    # Plot X and Y vs. Z
    plt.figure(1, figsize=(18, 4.5))
    plt.subplot(131)
    plt.plot(Z, X)
    plt.xlabel('Z [m]')
    plt.ylabel('X [m]')
    plt.title('Particle Trajectory')

    plt.subplot(132)
    plt.plot(Z, Y)
    plt.xlabel('Z [m]')
    plt.ylabel('Y [m]')
    plt.title('Particle Trajectory')

    plt.subplot(133)
    plt.plot(X, Y)
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.title('Particle Trajectory')


    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt


def plot_trajectory_velocity(trajectory, show=True, ofile=''):
    """Plot the trajectory"""

    # Get coordinate lists
    VX = [item[1][0] for item in trajectory]
    VY = [item[1][1] for item in trajectory]
    VZ = [item[1][2] for item in trajectory]
    T = range(len(VX))

    # Plot VX, VY, VZ vs. T
    plt.figure(1, figsize=(18, 4.5))
    plt.subplot(131)
    plt.plot(T, VX)
    plt.xlabel('T [step]')
    plt.ylabel('BX []')
    plt.title('Particle Beta X')

    plt.subplot(132)
    plt.plot(T, VY)
    plt.xlabel('T [step]')
    plt.ylabel('BY []')
    plt.title('Particle Beta Y')

    plt.subplot(133)
    plt.plot(T, VZ)
    plt.xlabel('T [step]')
    plt.ylabel('BZ []')
    plt.title('Particle Beta Z')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt
    
    
def plot_power_density(V, title='Power Density [$W / mm^2$]', xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', show=True, ofile=''):
    """Plot a 2D histogram with equal spacing"""
        
    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    P = [item[1]    for item in V]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    plt.figure(1, figsize=(18, 9))
    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    plt.colorbar()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt


def plot_flux(V, title='Flux [$\gamma / mm^2 / 0.1\%bw / s]$', xlabel='X1 Axis [$m$]', ylabel='X2 Axis [$m$]', show=True, ofile=''):
    """Plot a 2D histogram with equal spacing"""
        
    X = [item[0][0] for item in V]
    Y = [item[0][1] for item in V]
    P = [item[1]    for item in V]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))

    plt.figure(1, figsize=(18, 9))
    #plt.figure(1)
    plt.hist2d(X, Y, bins=[NX, NY], weights=P)
    plt.colorbar()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt


def plot_spectrum(S, log=False, show=True, ofile=''):
    """Plot the spectrum"""
    
    X = [item[0] for item in S]
    Y = [item[1] for item in S]
    plt.plot(X, Y)
    if log:
        plt.yscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('$\gamma / mm^2 / 0.1\%bw / s$')
    plt.title('Spectrum')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt



def plot_spectra(S, L, show=True, ofile='', title='', loc='upper left', log=False):

    plt.figure(1, figsize=(18, 9))
    for i in range(len(S)):
        s = S[i]
        
        X = [item[0] for item in s]
        Y = [item[1] for item in s]
        plt.plot(X, Y, label=L[i])

    if log:
        plt.yscale('log')
        
    plt.legend(loc=loc)
    plt.xlabel('Energy [eV]')
    plt.ylabel('[$\gamma / mm^2 / 0.1\%bw / s$]')
    if title == '':
        title='Spectrum'
    plt.title(title)
    
    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt





def plot_magnetic_field(srs, zmin, zmax, show=True, ofile=''):
    """Plot the magnetic field as a function of Z"""

    Z = np.linspace(zmin, zmax, 100000)
    Bx = [srs.get_bfield([0, 0, z])[0] for z in Z]
    By = [srs.get_bfield([0, 0, z])[1] for z in Z]
    Bz = [srs.get_bfield([0, 0, z])[2] for z in Z]

    plt.figure(1, figsize=(18, 4.5))
    plt.subplot(131)
    plt.plot(Z, Bx)
    plt.xlabel('Z [m]')
    plt.ylabel('Bx [T]')
    plt.title('Horizontal Magnetic Field')

    plt.subplot(132)
    plt.plot(Z, By)
    plt.xlabel('Z [m]')
    plt.ylabel('By [T]')
    plt.title('Vertical Magnetic Field')

    plt.subplot(133)
    plt.plot(Z, Bz)
    plt.xlabel('Z [m]')
    plt.ylabel('Bz [T]')
    plt.title('Longitudinal Magnetic Field')

    if ofile != '':
        plt.savefig(ofile, bbox_inches='tight')

    if show == True:
        plt.show()

    return plt
    
    
def total_power(pd):
    """Calculate the total power in a uniform grid.
    
    This will not work for a non-uniform grid.  Different NX and NY are ok."""
    
    X = [item[0][0] for item in pd]
    Y = [item[0][1] for item in pd]
    P = [item[1]    for item in pd]

    NX = len(np.unique(X))
    NY = len(np.unique(Y))
    
    dx = (max(X) - min(X)) / (NX - 1)
    dy = (max(Y) - min(Y)) / (NY - 1)
    
    return dx * dy * sum(P) * 1e6




def plot_electric_field_vs_time(efield, show=True, ofile=''):
    """Plot the electric field as a function of time"""

    T  = [item[0]    for item in efield]
    Ex = [item[1][0] for item in efield]
    Ey = [item[1][1] for item in efield]
    Ez = [item[1][2] for item in efield]

    plt.figure(1, figsize=(18, 9))
    plt.subplot(131)
    plt.plot(T, Ex)
    plt.xlabel('T [s]')
    plt.ylabel('Ex')
    plt.title('Electric Field (Ex)')

    plt.figure(1, figsize=(18, 9))
    plt.subplot(132)
    plt.plot(T, Ey)
    plt.xlabel('T [s]')
    plt.ylabel('Ey')
    plt.title('Electric Field (Ey)')

    plt.figure(1, figsize=(18, 9))
    plt.subplot(133)
    plt.plot(T, Ez)
    plt.xlabel('T [s]')
    plt.ylabel('Ez')
    plt.title('Electric Field (Ez)')

    return plt
    
    
    
