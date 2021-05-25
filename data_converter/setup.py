from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Python data processing and formatting tools for gwas summary stats'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="DataConverter", 
        version=VERSION,
        author="Yuxiang Chen",
        author_email="<ychen1@hsph.harvard.edu>",
        description=DESCRIPTION,
        packages=find_packages(),
        install_requires=["pyBigWig","pyliftover","numpy", "pandas"], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        zip_safe=False

)