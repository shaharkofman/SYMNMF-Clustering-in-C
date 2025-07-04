from setuptools import Extension, setup

module = Extension(
    'symnmf', 
    sources=[
        'symnmf.c',
        'symnmfmodule.c'
    ]
)

setup(
    name='symnmf',
    version='1.0',
    description='symmetric non negative nmf factorization of a matrix',
    ext_modules=[module]
)
