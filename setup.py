"""If you have a problem, if no one else can help, and if you can find them, maybe you can make use of the A-Tools."""

from setuptools import setup, find_packages

setup(
    name="atools",
    version="0.1.0",
    long_description=__doc__,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
)
