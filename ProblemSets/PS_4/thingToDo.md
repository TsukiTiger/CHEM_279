1. Figure out the gamma implementation
    - Done
2. Change the molecule class, so that it can read the element type as a string
     - maybe add another vector as elements
     - Done
3. Figure out the difference between n_AO and n_Basis
     - They are the same
4. Figure out how to get R from AO `vec Ra = mol.nAtoms[k].mAOs[0].get_RO();`
    - Make a "Atom" class, which is the object under "Molecule"
    - Using maps to create AO (BasisFunctions) under the Atoms
    - done
5. implement Fock Matrix and Hamiltonian
    - Working on it
6. Update AO.cpp and other utils as needed

### 2023/11/26
1. Getting the 2 electron integral of s AO done
    - Should the len be 3? Since it is only 1 AO and 3D (eq 3.3)
    - The contraction coefficients were normalized
2. Working on the eval_gammamat() function in molecule class
    - Done

 - Focus on the dim of the calculation 
