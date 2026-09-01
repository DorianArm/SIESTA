# SIESTA App: Simulating Interactive Echelle Spectroscopy for Targeted Applications

# Installation
```
git clone https://github.com/DorianArm/SIESTA.git
cd SIESTA/
pip install .
```
Once the required packages installed, make sure to modify the refractiveindex.py file from its package. 2 lines have to be added to be able to use prisms with SIESTA:

```
self._coefs = None # in __init__ method of RefractiveIndexMaterial class, line 283
self._coefs = coefficients just after coefficients variable declaration, line 314
```
![image_l283](change_refractiveindex_line283.png)
![image_l314](change_refractiveindex_line314.png)


# Changing parameters
Open siesta.py and change the desired parameters. Execute the file from your IDE or run
```python siesta.py```


