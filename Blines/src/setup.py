# setup.py

from distutils.core import setup, Extension

setup(name="Bline",
      py_modules=['Bline'], 
      ext_modules=[Extension("_Bline",
                     ["Bline.i","Bline.cc"],
                     swig_opts=['-c++'],
                  )]  
    
)
