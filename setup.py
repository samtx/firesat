from setuptools import setup, find_packages
# from setuptools.extension import Extension
# from distutils.core import setup
# from distutils.extension import Extension
from Cython.Build import cythonize

description = (
    "A multiphysics fire detection satellite model originally "
    "described in Wertz 'Space Mission Analysis and Design', 1999"
)

# ext = cythonize(['firesat/_sgp4.pyx'], language_level = 3)
ext = []

# ext = [Extension('sgp4lib',
#                 sources=['firesat/sgp4lib.pyx', 'firesat/SGP4.cpp'],
#                 include_dirs=['./firesat/'],
#                 compiler_directives={'language_level' : 3},
#                 extra_compile_args=['-O3'],
#                 language='c++',
#                 annotate=True),]
# ext = cythonize(ext)

setup(
    name="firesat",
    author='Sam Friedman',
    author_email="samfriedman@tamu.edu",
    packages=find_packages(),
    description=description,
    python_requires='>=3',
    install_requires=[
        'numpy>=1.14',
        'sgp4>=2.8',
    ],
    test_suite='nose.collector',
    tests_require=[
        'nose'
    ],
    docs_require=[
        'sphinx'
    ],
    zip_safe=False,
    ext_modules=ext,
)
