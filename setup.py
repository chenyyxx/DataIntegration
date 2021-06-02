from setuptools import setup, find_packages

VERSION = '0.0.5' 
DESCRIPTION = 'Python data processing and formatting tools for gwas summary stats'
classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Operating System :: MacOS',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="dataintegrator", 
        version=VERSION,
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type="text/markdown",
        url = 'https://github.com/chenyyxx/DataIntegration',
        author="Yuxiang Chen",
        author_email="ychen1@hsph.harvard.edu",
        license = 'MIT',
        classifiers=classifiers,
        keywords='Data Integrater',
        packages=find_packages(include=['dataintegrator']),
        install_requires=["pyBigWig","pyliftover","numpy", "pandas"] # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'

)