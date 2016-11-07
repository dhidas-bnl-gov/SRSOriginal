from distutils.core import setup, Extension
import os

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"
os.environ["LD"] = "g++"
os.environ["LN"] = "g++"

moduleOSCARS = Extension('OSCARS',
                      include_dirs = ['include'],
                      sources = ['src/OSCARS.cc',
                                 'src/OSCARS_Python.cc',
                                 'src/T3DScalarContainer.cc',
                                 'src/TField3D_Grid.cc',
                                 'src/TField3D_Gaussian.cc',
                                 'src/TFieldContainer.cc',
                                 'src/TField3D_1D.cc',
                                 'src/TField3D_1DRegularized.cc',
                                 'src/TField3D_IdealUndulator.cc',
                                 'src/TField3D_UniformBox.cc',
                                 'src/TBFieldIdeal1D.cc',
                                 'src/TFieldPythonFunction.cc',
                                 'src/TBFieldSquareWave.cc',
                                 'src/TParticleA.cc',
                                 'src/TParticleBeam.cc',
                                 'src/TParticleBeamContainer.cc',
                                 'src/TParticleTrajectoryPoints.cc',
                                 'src/TRandomA.cc',
                                 'src/TSpectrumContainer.cc',
                                 'src/TSurfaceElement_Rectangle.cc',
                                 'src/TSurfaceOfPoints.cc',
                                 'src/TSurfacePoint.cc',
                                 'src/TSurfacePoints_3D.cc',
                                 'src/TSurfacePoints_Rectangle.cc',
                                 'src/TVector2D.cc',
                                 'src/TVector3D.cc',
                                 'src/TVector3DC.cc',
                                 'src/TVector4D.cc'],
                      extra_compile_args=['-std=c++11', '-Wno-write-strings', '-Wall', '-O3', '-pedantic', '-fPIC', '-pthread']
                     )




setup(
  version="1.0",
  description = 'This is an example of how to create a new type in a python extension',
  author = 'Dean Andrew Hidas',
  author_email = 'dhidas@bnl.gov',
  url = 'http://oscars.bnl.gov/',
  long_description = '''This creates a new python type using python-C extensions.''',
  ext_modules = [moduleOSCARS]
)
