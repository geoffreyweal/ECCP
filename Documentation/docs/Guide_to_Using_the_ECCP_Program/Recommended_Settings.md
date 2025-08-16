# Recommended Settings

The recommended settings for this program are given below in python format:

```python title="Recommended settings for Running the ECCP program"
from ECCP import ECCP

# This is the method use to reassemble individual molecule from the crystal. 
make_molecule_method = 'component_assembly_approach'
# This dictionary include information about determining which molecules are equivalent. Required if you want to perform ATC calculations on molecules.
molecule_equivalence_method = {'method': 'invariance_method', 'type': 'combination'} 

# This is the method use to obtain dimers between molecules in the system.
dimer_method = {'method': 'nearest_atoms_method', 'max_dimer_distance': 8.0}
# This dictionary provides information for determining which dimers are equivalent
dimer_equivalence_method = {'method': 'invariance_method', 'type': 'combination'} 

# This dictionary includes info about how to treat the enivornment surrounding dimers (where applicable).  
environment_settings = {'include_environment_in_molecule_calcs': False, 'include_environment_in_dimer_calcs': False}

# This tag indicates if you want to remove solvents from the crystal. This requires the input file to have a reference to which molecules are solvents called "SolventList"
remove_solvents = False
```

!!! tip

	It is also recommended thatr you check out the crystal before hand to make sure that there are no issues with it, and to remove sidegroups from the molecules before running the ECCP method. See:

	* The [ReCrystals Program](https://geoffreyweal.github.io/ReCrystals) to help to repair molecules in your crystals, and 
	* Check out the [RSGC Program](https://geoffreyweal.github.io/RSGC/index.html#what-is-the-remove-sidegroups-from-crystals-rsgc-program) to see how to automatically remove sidegroups (like aliphatic sidechains) from the molecules in the crystal.  
