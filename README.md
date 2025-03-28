# pydisort
Python wrapper for DISORT

## Installation
Install the module as follows:

```
python -m pip install .
```

## Example
```python
import numpy as np
import pydisort

disort = pydisort.PyDisort()

flags = {'usrang': True, 'usrtau': False, 'ibcnd': False, 'lamber': True, 'onlyfl': False,
         'quiet': False, 'spher': False, 'old_intensity_correction': True, 'planck': False}
disort.set_flags(**flags)

bc = {'fbeam': np.pi / 0.1, 'fisot': 0, 'albedo': 0, 'umu0': 0.1, 'phi0': 0., 'temis': 0, 'btemp': 0, 'ttemp': 0}
disort.set_boundary(**bc)

nlayers = 1
nmom = 16
nstr = 16

# get the phase function moments for an isotropic atmosphere
pmom = pydisort.get_phase_function('isotropic', nmom + 1, 0)

# copy the moments for all the layers of the atmosphere
pmom = np.repeat(pmom[:, np.newaxis], nlayers, axis=1).flatten().tolist()

# set the input atmospheric structure
# tau/SSA is defined on the layer edges
tau = [0.03125, 0.03125]
ssalb = [0.2, 0.2]

# sampling incidence angles
umu = [-1, -0.5, -0.1, 0.1, 0.5, 1]
phi = [0]

# set the input values
inp = {'tau': tau, 'ssa': ssalb, 'umu': umu, 'phi': phi, 'pmom': pmom, 'nstr': nstr,
       'nphase': nstr, 'nmom': nmom}
disort.set_input(**inp)

# Run the model
output = disort.run()

# Print the output flux
print(np.asarray(output['uu']))
```

Should return:
```
Using original intensity correction, with phase moments
[[0.         0.         0.         0.11777066 0.02641704 0.01340413]
 [0.01338263 0.02633235 0.11589789 0.         0.         0.        ]]
```
