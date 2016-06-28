from distutils.core import setup, Extension


moduleSRS = Extension('SRS',
                      include_dirs = ['include'],
                      sources = ['src/SRS.cc',
                                 'src/SRS_Python.cc',
                                 'src/TBField.cc',
                                 'src/TBField1DZ.cc',
                                 'src/TBField1DZRegularized.cc',
                                 'src/TBField3D.cc',
                                 'src/TBField3DZ.cc',
                                 'src/TBField3DZRegularized.cc',
                                 'src/TBFieldContainer.cc',
                                 'src/TBFieldIdeal1D.cc',
                                 'src/TBFieldSquareWave.cc',
                                 'src/TBFieldUniformB.cc',
                                 'src/TParticleA.cc',
                                 'src/TParticleBeam.cc',
                                 'src/TParticleBeamContainer.cc',
                                 'src/TParticleTrajectoryPoints.cc',
                                 'src/TSpectrumContainer.cc',
                                 'src/TSurfaceElement_Rectangle.cc',
                                 'src/TSurfaceOfPoints.cc',
                                 'src/TSurfacePoint.cc',
                                 'src/TSurfacePoints_BBoxSimple.cc',
                                 'src/TSurfacePoints_Rectangle.cc',
                                 'src/TSurfacePoints_RectangleSimple.cc',
                                 'src/TTwiss.cc',
                                 'src/TVector3D.cc',
                                 'src/TVector3DC.cc'
                                 ],
                      extra_compile_args=['-std=c++17', '-ansi']
                     )




setup(
  version="1.0",
  description = 'This is an example of how to create a new type in a python extension',
  author = 'Dean Andrew Hidas',
  author_email = 'dhidas@bnl.gov',
  url = 'https://github.com/dhidas',
  long_description = '''This creates a new pythong type using the python-C extensions.''',
  ext_modules = [moduleSRS]
)
