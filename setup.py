"""Pypi software definition"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name='pyowc',
    version='0.1.1',
    author="pyowc Developers",
    author_email="edsmall@bodc.ac.uk",
    description="OWC salinity calibration in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/euroargodev/argodmqc_owc",
    packages=setuptools.find_packages(),
    package_dir={'pyowc': 'pyowc'},
    package_data={'pyowc': ['data']},
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha"
    ]
)
