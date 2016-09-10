from math import sin, cos, sqrt, pi

class PSTorus:
    """A Parametric surface - torus"""

    # This shape specific parameters
    R = 1
    r = 1
    

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = 0
    ustop  = 2 * pi
    vstart = 0
    vstop  = 2 * pi
    nu = 10
    nv = 10
 
    def __init__ (self, R=1, r=1, ustart=0, ustop=2*pi, nu=10, vstart=0, vstop=2*pi, nv=10):
        self.R = R
        self.r = r
        self.ustart = ustart
        self.ustop = ustop
        self.nu = nu
        self.vstart = vstart
        self.vstop = vstop
        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = (self.R + self.r * cos(u)) * cos(v)
        y = (self.R + self.r * cos(u)) * sin(v)
        z = self.r * sin(u)
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        xn = -self.r * cos(u) * (self.R + self.r * cos(u)) * cos(v)
        yn = -self.r * cos(u) * (self.R + self.r * cos(u)) * sin(v)
        zn = -self.r * (self.R + self.r * cos(u)) * sin(u)
        
        mag = sqrt(xn*xn + yn*yn + zn*zn)
        
        return [xn / mag, yn / mag, zn / mag]




class PSCylinder:
    """A Parametric surface - cylinder with no top or bottom"""

    # This shape specific parameters
    R = 1
    L = 1

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = -L/2.
    ustop  = +L/2.
    vstart = 0
    vstop  = 2 * pi
    nu = 10
    nv = 10
 
    def __init__ (self, R=1, L=1, ustart=None, ustop=None, nu=10, vstart=None, vstop=None, nv=10):
        self.R = R
        self.L = L
        if ustart is not None: self.ustart = ustart
        if ustop  is not None: self.ustop = ustop
        self.nu = nu
        if vstart is not None: self.vstart = vstart
        if vstop  is not None: self.vstop = vstop
        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = self.R * cos(v)
        y = self.R * sin(v)
        z = u
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        xn = cos(v)
        yn = sin(v)
        zn = 0
        
        mag = sqrt(xn*xn + yn*yn + zn*zn)
        
        return [xn / mag, yn / mag, zn / mag]


