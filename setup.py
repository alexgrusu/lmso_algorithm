import setuptools
from pathlib import Path

setuptools.setup(
  name = "lmso_algorithm",
  version = 1.2,
  long_description = Path("README.md").read_text(),
  packages = setuptools.find_packages(exclude = ["tests", "example"]),
  license = 'GPLv3+',
  author="Alexandru - George Rusu",
  author_email="aqrusu@gmail.com",
  description="An Optimized LMS Algorithm",
  keywords="Adaptive filters, Echo cancellation, System identification",
  url="https://github.com/alexgrusu/lmso_algorithm",
  install_requires=['numpy', 'scipy'],
  classifiers=[
        'Development Status :: 6 - Mature',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering',
        ],
)