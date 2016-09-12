import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def power_density_3d(srs, surface,
                    normal=1, rotations=[0, 0, 0], translation=[0, 0, 0], nparticles=0, gpu=0):
    """calculate power density for and plot a parametric surface in 3d"""

    points = []


    for u in np.linspace(surface.ustart, surface.ustop, surface.nu):
        for v in np.linspace(surface.vstart, surface.vstop, surface.nv):
            points.append([surface.position(u, v), surface.normal(u, v)])


    power_density = srs.calculate_power_density(points=points, normal=normal, rotations=rotations, translation=translation, nparticles=nparticles, gpu=gpu)
    P = [item[1] for item in power_density]

    X2 = []
    Y2 = []
    Z2 = []
    for i in range(surface.nu):
        tX = []
        tY = []
        tZ = []
        for j in range(surface.nv):
            tX.append(power_density[i * surface.nv + j][0][0])
            tY.append(power_density[i * surface.nv + j][0][1])
            tZ.append(power_density[i * surface.nv + j][0][2])
        X2.append(tX)
        Y2.append(tY)
        Z2.append(tZ)

    colors =[]
    MAXP = max(P)
    PP = []
    for i in range(surface.nu):
        tmpP = []
        tmpPP = []
        for j in range(surface.nv):
            tmpP.append(P[i * surface.nv + j] / MAXP)
            tmpPP.append(P[i * surface.nv + j])

        colors.append(tmpP)
        PP.append(tmpPP)
  


    fig = plt.figure(figsize=(20, 15))
    ax = fig.gca(projection = '3d')
    #ax.view_init(30, -40)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Z [m]')
    ax.set_zlabel('Y [m]')
    #ax.set_ylim(translation_[2] - 0.1, translation_[2] + 0.1)
    ax.plot_surface(X2, Z2, Y2, facecolors=cm.jet(colors), cmap=cm.jet, rstride=1, cstride=1, alpha=1)
    ax.invert_xaxis()

    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(PP)
    plt.colorbar(m)

    plt.show()
    
    return plt
