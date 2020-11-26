
# https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56

import pathlib
from setuptools import setup, find_packages

__author__ = 'Guo Jun-Lin'
__license__ = "BSD 2 Clause"
__email__ = "guojl19@mails.tsinghua.edu.cn"

# Get the long description from the README file
# with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()

# The directory containing this file.
HERE = pathlib.Path(__file__).parent
# The text of the README file
README = (HERE / 'README.md').read_text()

install_requires = ['networkx', 'numpy', 'matplotlib',
                    'tqdm', 'scipy', 'cdlib', 'future']

keywords = ['epidemic', 'vaccination', 'graph', 'networkx', 'community',
            'propagation', 'differential-equation', 'covid-19', 'cdlib']

setup(name='epidemix',
      version='1.0.2',
      license='BSD-2-Clause',
      description='Simulation of Epidemic Propagation on a Network',
      long_description=README,
      long_description_content_type='text/markdown',
      url='https://github.com/khle08/epidemix',
      author='Guo Jun-Lin',
      author_email='guojl19@mails.tsinghua.edu.cn',
      classifiers=[
          # How mature is this project? Common values are
          #   3 - Alpha
          #   4 - Beta
          #   5 - Production/Stable
          'Development Status :: 3 - Alpha',

          # Indicate who your project is intended for
          'Intended Audience :: Developers',
          'Topic :: Software Development :: Build Tools',

          # Pick your license as you wish (should match "license" above)
          'License :: OSI Approved :: BSD License',

          "Operating System :: OS Independent",

          # Specify the Python versions you support here. In particular, ensure
          # that you indicate whether you support Python 2, Python 3 or both.
          # 'Programming Language :: Python',
          # 'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3'
      ],
      keywords=keywords,
      install_requires=install_requires,
      packages=['epidemix'],
      include_package_data=True,
      entry_points={'console_scripts': ['epidemix=epidemix.__main__:main']}
      )
