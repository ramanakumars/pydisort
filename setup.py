from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

ext_modules = [
    Extension(
        name="pydisort",
        sources=["src/pydisort.pyx", "src/locate.c", "src/cdisort.c"],  # Cython and C sources
        include_dirs=["src"],
        extra_compile_args=['-O3', '-g']
    )
]

setup(
    name="pydisort",
    version="0.1",
    ext_modules=cythonize(ext_modules),
)
