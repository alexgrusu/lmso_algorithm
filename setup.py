import setuptools
from pathlib import Path

setuptools.setup(
  name = "lmso_algorithm",
  version = 1.1,
  long_description = Path("README.md").read_text(),
  packages = setuptools.find_packages(exclude = ["tests", "example"]),
  author="Alexandru - George Rusu",
  author_email="aqrusu@gmail.com",
  description="An Optimized LMS Algorithm",
  keywords="Adaptive filters, Echo cancellation, System identification",
  url="https://github.com/alexgrusu/lmso_algorithm",   # project home page, if any
)