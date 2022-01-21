from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'AeroFEM Package'
LONG_DESCRIPTION = 'A Python 3 package that provides FEM-based aircraft modeling tools.'

setup(
        name="aerofem", 
        version=VERSION,
        author="Guilherme Levi",
        author_email="guilhermelevi@usp.br",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], 
        
        keywords=['python'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)