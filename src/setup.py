from setuptools import Extension, setup

module = Extension("mykmeanssp",
                  sources=[
                    'kmeansmodule.c'
                  ])
setup(name='mykmeanssp',
     version='1.0',
     description='calculating kmeans apart from the mu',
     ext_modules=[module])