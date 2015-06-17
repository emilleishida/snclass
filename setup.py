"""Configure snclass instalation."""

from setuptools import setup
import snclass

def readme():
    """Return README.md file."""
    with open('README.md') as my_doc:
        return my_doc.read()

setup(name='snclass',
      version=snclass.__version__,
      description='Python SN photometric classifier',
      long_description=readme(),
      url='https://github.com/emilleishida/snclass',
      author=snclass.__author__,
      author_email=snclass.__email__,
      license='GPL3',
      packages=['snclass'],
      install_requires=[
                      'numpy>=1.8.2',
                      'matplotlib>=1.3.1',
                      'gptools>=0.1'
      ],
      scripts=['snclass/bin/fit_plot_lc.py', 
               'snclass/bin/build_synthetic_spec.py'],
      package_dir={'snclass': 'snclass', 'examples':'snclass/examples'},
      zip_safe=False,
      classifiers=[
        'Programming Language :: Python',
        'Natural Language :: English',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Astronomy',
        ])
