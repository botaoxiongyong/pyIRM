from setuptools import setup

NAME = "pyIRM"

setup(
    name=NAME,
    version='1.0.0',
    url='https://github.com/botaoxiongyong/pyIRM',
    author='Liu Jiabo',
    author_email='kabol.liu@gmail.com',
    licence='MIT'
    packages=['pyIRM']
    install_requires=['matplotlib',
                      'pandas',
                      'PyQt5',
                      'scipy',
                      'numpy',
                      'sklearn',
                      'lmfit']
)
