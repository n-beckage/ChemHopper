from ChemTools import *

### define some chemical data

# Corresponding atomic numbers of atomic_names
allowed_atomic=[6,7,8,9,15,16,17,35]
atomic_names=["C","N","O","F","P","S","Cl","Br"]
#Have you validated the Phosphorus chemistry? I would be careful here and think about the chemistry more. Many of the molecules when i look at the chemical space graph have wild phosphorus and sulfer chemistry 
halogens=[9,17,35]
halo_names=["F","Cl","Br"]
type_i=[0,1,2,1.5]
types=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]

################################################################# BEGIN SCRIPT #####################################################################################

print("Welcome to Chemical Space! Where would you like to begin?")

print("please enter a molecule in the form of a SMILE string.")

seed = input("GRAPH SEED: ")

print("How deep will the roots go?")

depth = int(input("GRAPH DEPTH: "))

yn = input("Include exhaustive connections? (y/n):")
if yn in ['y','Y']:
	cc = True
elif yn in ['n','N']:
	cc = False

print("cc =",str(cc))

buildGraph(seed,depth,cc)