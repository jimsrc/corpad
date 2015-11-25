# setup.py

from distutils.core import setup, Extension

setup(name="read_Bmodel",
      py_modules=['read_Bmodel'], 
      ext_modules=[Extension("_read_Bmodel",
                     ["read_Bmodel.i","read_Bmodel.cc"],
                     swig_opts=['-c++', '-dlatex'],
                  )]  
    
)
