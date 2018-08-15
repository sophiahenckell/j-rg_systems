from distutils.core import setup
from Cython.Build import cythonize

#from distutils.extension import Extension
#from Cython.Distutils import build_ext
from Cython.Distutils import build_ext, Extension 
import numpy as np

ext_modules = [
    Extension(
        name="bedsim.cython2.cconstraints",
        sources=["bedsim/cython2/cconstraints.pyx", "bedsim/cython2/constraints.c"],
            # extra_objects=["fc.o"],  # if you compile fc.cpp separately
        include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
        language="c",
        extra_compile_args = "-O3 -Wall".split(),
            # libraries=
            # extra_compile_args = "...".split(),
            # extra_link_args = "...".split()
    ),
    Extension(
        name="bedsim.cython2.ccsht",
        sources=["bedsim/cython2/ccsht.pyx", "bedsim/cython2/csht.c"],
        include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
        language="c",
        extra_compile_args = "-O3 -Wall".split(),
    ),
    Extension(
        name="bedsim.cython.constraints",
        sources = ["bedsim/cython/constraints.pyx"],
        include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
        language="c",
    ),
    Extension(
        name="bedsim.cython2.coverlap",
        sources = ["bedsim/cython2/coverlap.pyx", "bedsim/cython2/overlap.c", "bedsim/cython2/cfAB.c", "bedsim/cython2/c_aAB.c", "bedsim/cython2/cAAB.c", "bedsim/cython2/crC.c", "bedsim/cython2/cvC.c", "bedsim/cython2/cnorm.c"],
        include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
        language="c",
        #extra_compile_args = "-lm -include stdlib.h -O3 -Wall".split(), # -lm flag links against math.h
        extra_compile_args = "-lm -O3 -Wall".split(), # -lm flag links against math.h
        libraries = ["m"],
        compiler_directives={'boundscheck': False, 'wraparound': False, 'initializedcheck':False, 'language_level':3},
    ),
    #Extension(
    #    name="bedsim.boundary",
    #    sources = ["bedsim/boundary.py"],
    #    include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
    #),
    #Extension(
    #    name="bedsim.particle",
    #    sources = ["bedsim/particle.py"],
    #    include_dirs = [np.get_include()],  # .../site-packages/numpy/core/include
    #),
    ]

setup(
    name = 'bedsim',
    version = '1.0.0',
    description = 'Brownian event driven simulation framework',
    author = 'Markus Heinrich',
    author_email = 'Markus.Heinrich@uni-konstanz.de',
    license = 'EUPL',
    packages = ['bedsim'],
    #install_requires = ['numpy', 'scipy', 'simpy', 'matplotlib', 'h5py', 'moviepy'], 
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
