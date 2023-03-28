from ChemToolsWd import *

### define some chemical data

# Corresponding atomic numbers of atomic_names
allowed_atomic=[6,7,8,9,15,16,17,35]
atomic_names=["C","N","O","F","P","S","Cl","Br"]
#Have you validated the Phosphorus chemistry? I would be careful here and think about the chemistry more. Many of the molecules when i look at the chemical space graph have wild phosphorus and sulfer chemistry 
halogens=[9,17,35]
halo_names=["F","Cl","Br"]
type_i=[0,1,2,1.5]
types=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]