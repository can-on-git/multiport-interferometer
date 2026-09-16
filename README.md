# Multiport-Interferometer
> This library provides a set of tools for designing multiport interferometer photonic circuits to emulate unitary transformations. Includes generation of complete circuit gds layout and steady-state power and phase simulations of multiport interferometers.

[![PyPI version](https://img.shields.io/pypi/v/multiport-interferometer.svg)](https://pypi.org/project/multiport-interferometer/)
[![License](https://img.shields.io/pypi/l/multiport-interferometer.svg)](https://pypi.org/project/multiport-interferometer/)

Some of the functionalities of the library:

    1. Bell type interferometer:
        - Schematic Logic Simulation
        - Automated GDS Layout construction
        - Physical Circuit Simulation (Emulates unitary transformations)
    
    2. Clements type inteferometer:
        - Schematic Logic Simulation
        - Build Circuit Simulation

The library accepts custom files for circuit building blocks gds components and physical simulation models.  
Uses GDSfactory for gds file processing and SAX for S-parameter simulations.

This library is part of the Open Source P.I.C Design Workflow [Tutorial](https://carlabinakad.com/tutorial).  
For complete documentation check out [Tutorial](https://carlabinakad.com/tutorial).

Supporting notebooks: 
- [Simulations](https://github.com/can-on-git/multiport-interferometer/blob/main/supporting_notebooks/Simulations.ipynb)
- [MZI Calibration](https://github.com/can-on-git/multiport-interferometer/blob/main/supporting_notebooks/MZI_calibration.ipynb)
- [Building Blocks](https://github.com/can-on-git/multiport-interferometer/blob/main/supporting_notebooks/Building_blocks.ipynb)
- [Examples](https://github.com/can-on-git/multiport-interferometer/blob/main/supporting_notebooks/code_for_tutorial.ipynb)


## Installation

```sh
pip install multiport-interferometer
```

### Import:
```sh
import multiport_interferometer as mi
```
## Usage example

### Bell Interferometer Logic Simulation:

```sh
import multiport_interferometer as mi
from scipy.stats import unitary_group
import numpy as np
import sax

dimension = 58
U_target = unitary_group(dim=dimension).rvs()
results = mi.decompose_bell(U_target)
mods = mi.ideal_bell_models(results)

netlist = mi.bell_circuit_schematic(dimension)
circuit, info = sax.circuit(netlist=netlist, models = mods) 
S = circuit() #returns a sax S-dict object
S_mat = mi.sdict_to_matrix(S) #converts the sax sdict to a numpy matrix
print('Circuit Matches Target Matrix =', np.allclose(U_target, S_mat))
```
```text
Circuit Matches Target Matrix = True
```
### Bell Interferometer Physical SImulation:
```python
dim = 5
mymesh = mi.Bell(dim)
U_target = unitary_group(dim).rvs()
S = mymesh.emulate(U_target)
S_mat = mi.sdict_to_matrix(S)
corrected = mymesh.phase_corrected(S)
print('Maximum Magnitude Error % =', np.max((abs(abs(S_mat)-abs(U_target)))*100/abs(S_mat)))
print('Maximum Phase Error = % ', np.max(np.angle(corrected/U_target)*100/np.pi))
mymesh.gds().plot()
```
```text
Maximum Magnitude Error % = 9.225697420182325
Maximum Phase Error % =  8.080086688873948
```
![Mesh](https://github.com/can-on-git/multiport-interferometer/blob/main/docs/images/interferometer.png/)

### Clements Interferometer Logic Simulation
```sh
mesh = Clements(dim)
U = unitary_group.rvs(dim)
simulation_result = sdict_to_matrix(mesh.ideal_emulate(U))
print('Circuit Matches Target Matrix =', np.allclose(U, simulation_result))
```
```text
Circuit Matches Target Matrix = True
```

## Development setup

Install development dependencies:
```sh
pip install -e ".[dev]"
```

Run tests:
```sh
pytest
```

## Contributing
Contributions are welcome, please open an issue or pull request.