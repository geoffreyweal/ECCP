# Methods for Making and Obtaining Unique Molecules and Dimers from the Crystal

There are several options available for obtaining moelcules and dimers from the crystal, as well as for distinguishing unique molecules and dimers in the crystal. Multiple methods are given for the ``make_molecule_method``, ``molecule_equivalence_method``, ``dimer_method``, and ``dimer_equivalence_method`` dictionaries. 


## Methods available for ``make_molecule_method``

In order to obtain the individual molecules in the crystal, the molecule may be placed in a position that crosses from one side of the periodic boundary to the opposite side of the crystal lattice. While this is all fine, it does not make viewing the molecule very user-friendly. Therefore, the ECCP program will move the segments of the disconnected molecule around so that it is more user-friendly to view. There are two approaches for reconstructing the molecule available in ECCP. Ideally, both should give the same results.

### The Super Lattice Approach

This approach will create a super-lattice consisting of 26 lattices surrounding the origin unit cell. From this, a verison of the complete connected molecule will exist. All partial segments are deleted. This approach is very simple, but it take a bit of time to run.

To use this approach, set ``make_molecule_method = 'super_lattice_approach'`` in your ``Run_ECCP.py`` script. 

### The Component Assembly Approach (Recommended Approach) { #markdown data-toc-label='The Component Assembly Approach' }

This approach will construct a graph of the components of the molecule across the crystal and use this information to determine which bonds have been disconnected in the molecule. The components are then shifted into place based on the periodic boundary condition so that the two bonds are connected again. **This is the recommended approach to use**. If this approach fail, report the issue and use the ``'super_lattice_approach'``. 

To use this approach, set ``make_molecule_method = 'component_assembly_approach'`` in your ``Run_ECCP.py`` script. 


## Methods available for ``molecule_equivalence_method``

Some of the molecules that will have been obtained will be identical to each other, having the same twists and bends. It is useful to know which molecules are identical and which are not as this will reduce the number of gaussian jobs that need to be run (as gaussian jobs run on identical molecules will give the same results). The methods that are available are: 

### The Invariance Method

This method works by using the Procrustes analysis to determine if two molecules are translationally, rotationally, and reflectively invariant with each other. 

There are three versions of the invariance method that you can use:

### The Comprehensive Invariance Method

This method is the best to use for molecules that contain mainly one or a few types of elements in high abundances (not including hydrogen). To use this method, set ``dimer_equivalence_method = {'method': 'invariance_method', 'type': 'comprehensive'}``.

#### The Minimal Elemental Abundance Invariance Method

This method is the best to use for molecules that contains one or a few elements in low abundances. This method uses these low abundance elements as references to determine how dimer should be rotated onto one another. To use this method, set ``dimer_equivalence_method = {'method': 'invariance_method', 'type': 'minimal_elemental_abundance'}``.

#### The Combination Invariance Method (Recommended Method) { #markdown data-toc-label='The Combination Invariance Method' }

The combination invariance method will check the abundance of elements in your molecules (crystal) and automatically determine whether the comprehensive invariance method or the minimal elemental abundance invariance method should be used. **This is the recommended method to use**. To use this method, set ``dimer_equivalence_method = {'method': 'invariance_method', 'type': 'combination'}``.


## Methods available for ``make_dimer_method``

In order to identify dimers you need to specify the method for how you want to identify them. There are two methods available for identifying dimers in this program:

### The Centre of Mass Method

This method identifies dimers by looking at the distance between two molecules centre's of mass. If they are smaller than a predefined maximum dimer distance, then those two molecules are considered dimers. 

To use this method, set ``dimer_method = {'method': 'CMass_method', 'max_dimer_distance': 20.0}``, where ``max_dimer_distance`` is a float that is the maximum distance between two molecules centres of mass that you deem for those two molecules to be considered a dimer. 

### The Centre of Molecule Method

This method identifies dimers by looking at the distance between the centre of each monomers. If they are smaller than a predefined maximum dimer distance, then those two molecules are considered dimers. This method ignores hydrogens in the monomer when determining the centre of the monomer, as Hydrogens can not be observed in some cases in X-ray crystallography. 

To use this method, set ``dimer_method = {'method': 'CMol_method', 'max_dimer_distance': 20.0}``, where ``max_dimer_distance`` is a float that is the maximum distance between the centre of two molecules that you deem for those two molecules to be considered a dimer. 

### The Average Distance Method

This method identifies dimers by looking at the distances of the non-hydrogen atoms between two molecules. If the average distance between an atom in molecule 1 (or molecule 2) and all non-hydrogen atoms in molecule 2 (molecule 1) is smaller than a predefined maximum dimer distance, then those two molecules are considered dimers. 

To use this method, set ``dimer_method = {'method': 'average_distance_method', 'max_dimer_distance': 5.0}``, where ``max_dimer_distance`` is a float that is the maximum distance between any atom in molecule 1 (or molecule 2) and all non-hydrogen atoms in molecule 2 (molecule 1) for molecule 1 and molecule 2 to be considered a dimer. 

### The Nearest Atoms Method

This method identifies dimers by looking at the distances of the non-hydrogen atoms between two molecules. If the distance between any two non-hydrogen atoms is smaller than a predefined maximum dimer distance, then those two molecules are considered dimers. 

To use this method, set ``dimer_method = {'method': 'nearest_atoms_method', 'max_dimer_distance': 5.0}``, where ``max_dimer_distance`` is a float that is the maximum distance between any atoms in each molecules from each other for those two molecules to be considered a dimer. 


## Methods available for ``dimer_equivalence_method``

Some of the dimers that will have been obtained will be identical to each other. We would like to know which one of these dimers are equivalent and which ones are unique, as we may not necessarily need to perform Gaussian calculations on equivalent dimer systems. There are three methods available for identifying equivalent dimers in this program:

### The Atomic Distance Method

This method works by recording the distances in each dimer, from each atom in molecule 1 to every atom in molecule 2 (as well as recording the distances from each atom in molecule 2 to every atom in molecule 1). These distances are then ordered from each atom from shortest to largest distances. If the distances between atoms in dimer 1 are the same as the distance between atoms in dimer 2, the two dimers are considiered equivalent. 

To use this method, set ``dimer_equivalence_method = {'method': 'atomic_distance_method'}``.

### The Averaging Method

This method works by recording the distances in each dimer, from each atom in molecule 1 to every atom in molecule 2 (as well as recording the distances from each atom in molecule 2 to every atom in molecule 1). All these distances are then averaged to give an average distance. This average distance is then compared between dimers to 2 decimal places. If two dimers have the same average distance, they are considered equivalent. 

To use this method, set ``dimer_equivalence_method = {'method': 'averaging_method'}``.

### The Invariance Method

This method works by using the Procrustes analysis to determine if two dimers are translationally, rotationally, and reflectively invariant with each other. 

There are three versions of the invariance method that you can use:

#### The Comprehensive Invariance Method

This method is the best to use for molecules that contain mainly one or a few types of elements in high abundances (not including hydrogen). To use this method, set ``dimer_equivalence_method = {'method': 'invariance_method', 'type': 'comprehensive'}``.

#### The Minimal Elemental Abundance Invariance Method

This method is the best to use for molecules that contain a elements in low abundances. This method uses these low abundance elements as references to determine how dimer should be rotated onto one another. To use this method, set ``dimer_equivalence_method = {'method': 'invariance_method', 'type': 'minimal_elemental_abundance'}``.

#### The Combination Invariance Method (Recommended Method)

The combination invariance method will check the abundance of elements in your molecules (crystal) and automatically determine whether the comprehensive invariance method or the minimal elemental abundance invariance method should be used. **This is the recommended method to use**.  To use this method, set ``dimer_equivalence_method = {'method': 'invariance_method', 'type': 'combination'}``.





