from glob import glob
from setuptools import setup
from setuptools import find_packages

if __name__== "__main__":
    setup(name="dapid-project",
          version="0.0.1",
          description="DaPID analysis",
          author="Daniel Kim",
          author_email='danielskim@stanford.edu',
          url="https://github.com/vervacity/dapid-project",
          license="MIT",
          install_requires=[],
          packages=find_packages(),
          #package_data={'ggr':'data/*.json'},
          package_data={"":["data/*.json", "data/*.txt"]},
          scripts=["bin/dapid"] + glob("R/*/*.R")
    )
