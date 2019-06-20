#Anisotropic_Package_LAMMPS

S-function expansion and ellipsoid dihedral potentials

All files in src/ should be placed in the ASPHERE folder of the LAMMPS source code. 

pair style sfunction requires 12 arguments for its pair_coeffs: epsilon_0 epsilon_000 epsilon_cc2 epsilon_220 epsilon_222 epsilon_224 sigma_0 sigma_000 sigma_cc2 sigma_220 sigma_222 sigma_224

bond and angle style eldihedral (ellipsoid dihedral) allows torsional dihedral potentials to be defined between pairs of ellipsoids or an ellipsoid and two beads, respectively.
They take four coeff parameters: k1 k2 k3 k4. These parameters are similar to the dihedral_style opls

Rigid bodies are enforced using the fix rigid/nvt and fix rigid/npt commands.

A coarse-grained ethylbenzene example is provided in the example/ directory, while an input data file can be generated using generate.py in the py/ directory. 
The example shows how to use the sfunction pair potential, the eldihedral angle potential, as well as how to specifiy rigid bodies.

in.init should be run before in.production
