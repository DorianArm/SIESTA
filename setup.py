from setuptools import setup, find_packages

setup(
    name="siesta",
    version="1.0",
    author="Dorian Paillon",
    author_email="dorian.paillon@unibe.ch",
    description="Simulating Interactive Echelle Spectrogram for Targeted Applications",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/DorianArm/SIESTA",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.24",
        "pandas>=2.0",
        "dash>=3.0",
        "plotly>=5.0",
        "refractiveindex>=1.0.2"
    ],
)

# Need to add 2 lines in refractiveindex package (in refractiveindex.py)
# self._coefs = None # in __init__ method of RefractiveIndexMaterial class
# self._coefs = coefficients just after coefficients variable declaration
