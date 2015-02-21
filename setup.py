from setuptools import setup
import snclass.__init__ as snclass

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name = 'snclass',
      version = snclass.__version__,
      description = 'Python SN photometric classifier',
      long_description = readme(),
      url = 'https://github.com/COINtoolbox/CosmoABC',
      author = snclass.__author__,
      author_email = snclass.__email__,
      license = 'GNU Public License',
      packages = ['snclass'],
      install_requires=[
                      'numpy>=1.8.2',
                      'matplotlib>=1.3.1',   
                      'george>=0.2.1'          
      ],
      scripts=['snclass/bin/fit_lc_george.py'],
      package_dir= {'snclass': 'snclass', 'examples':'snclass/examples'},
      zip_safe=False,
      classifiers = [
        'Programming Language :: Python',
        'Natural Language :: English',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Astronomy',
        ])
