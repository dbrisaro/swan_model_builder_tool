from setuptools import setup, find_packages

setup(
    name="swan_model_builder",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "geopandas>=0.9.0",
        "shapely>=1.7.0",
        "matplotlib>=3.4.0",
        "contextily>=1.2.0",
        "netCDF4>=1.5.7",
        "PyYAML>=5.4.1",
        "cdo>=1.5.4",
        "cdsapi>=0.5.1",
    ],
) 