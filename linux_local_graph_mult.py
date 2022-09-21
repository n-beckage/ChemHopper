# record start time
import time
t0=time.time()

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit.Chem import Draw
import subprocess as sp
import pandas as pd
from tabulate import tabulate as tab
import os
import matplotlib.pyplot as plt
import sys
import random as ran
import re
import logging

### define some chemical data

# Corresponding atomic numbers of atomic_names
allowed_atomic=[6,7,8,9,17]#,9,15,16,17,35]#implicit hydrogen
atomic_names=["C","N","O","F","Cl"]
halogens=[9,17]
type_i=[0,1,2,1.5]
types=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]
n_save=5 #save top n_save for next round

### helper functions for molecule editing

### returns the degree of the atom (or molecule) on the chemical graph.
# IOW returns the number of non-H bonds for an atom, OR one-half the sum of said numbers for all atoms in a molecule (if a mol is passed)
def degree(at):
    return sum([types.index(str(x.GetBondType()))+1 for x in at.GetBonds()])

### returns the name/element of the atom
def aname(at):
    return atomic_names[allowed_atomic.index(at.GetAtomicNum())]

### returns the type of a given bond
def btype(b):
    return str(b.GetBondType())

### changes an atom
## old_m - the mol to be mutated
## at - the atom that is being changed. Not sure if this is needed, as it is never used within the functions
## ati - the index of the atom being changed
## newat - the new atom type (what the atom is being changed to)
## newcrg - the new formal charge of the atom
def mutate(old_m,at,ati,newat,newcrg):
    # make a new molecule with the 
    edit_mol=Chem.rdchem.RWMol(old_m)
    #nm=Chem.MolFromSmiles(Chem.MolToSmiles(m))#make new mol
    new_atom=Chem.rdchem.Atom(allowed_atomic[atomic_names.index(newat)])
    #new_atom.SetAtomicNum(allowed_atomic[atomic_names.index(newat)])
    new_atom.SetFormalCharge(newcrg)
    edit_mol.ReplaceAtom(ati,new_atom,preserveProps=True)
    new_smi=Chem.MolToSmiles(edit_mol)
    #ato=edit_mol.GetAtoms()[ati]
    #ato.SetAtomicNum(allowed_atomic[atomic_names.index(newat)])
    #ato.SetFormalCharge(newcrg)
    return new_smi

### adds a single bond to a new atom and returns a new SMILE string
## m - the mol object to be changed
## at - 
## ati - the index of the atom on the old mol to which the new atom will be attached
## newat - type of new atom to be added; given as atomic symbol (i.e. 'C', 'N', etc.)
def add_single_bond(m,at,ati,newat):
    # RWMol is a subclass of the Mol object class that is basically an editable form of a Mol
    ## https://herongyang.com/Cheminformatics/RDKit-rdkit-Chem-rdchem-RWMol-Class.html#:~:text=rdkit.Chem.rdchem.RWMol%20represents%20a%20molecule%20class%20with%20additional%20read,some%20useful%20methods%20provided%20in%20the%20RWMol%20sub-class.
    edit_mol=Chem.rdchem.RWMol(m)
    # creates a new atom based on the corresponding atomic number of the the atomic symbol
    new_atom=Chem.rdchem.Atom(allowed_atomic[atomic_names.index(newat)])
    # adds the new atom to edit_mol and returns the new index of the atom (at_ind)
    at_ind=edit_mol.AddAtom(new_atom)
    # adds the new bond, returns the new number of bonds
    edit_mol.AddBond(ati,at_ind,Chem.BondType.SINGLE)
    # the new smile string generated is sometimes weird, but still encodes the molecule so I guess it's alright
    new_smi=Chem.MolToSmiles(edit_mol)
    #currently deliminates sterochemistry uhg..
    #check out https://github.com/rdkit/rdkit/issues/2627
    #exit()
    return new_smi

### adds a double bond to a new atom and returns a new SMILE string
## m - the mol object to be changed
## at - 
## ati - the index of the atom on the old mol to which the new atom will be attached
## newat - type of new atom to be added; given as atomic symbol (i.e. 'C', 'N', etc.)
def add_double_bond(m,at,ati,newat):
    # making an editable version of the mol
    edit_mol=Chem.rdchem.RWMol(m)
    # creates a new atom based on the corresponding atomic number of the the atomic symbol
    new_atom=Chem.rdchem.Atom(allowed_atomic[atomic_names.index(newat)])
    # adds the new atom to edit_mol and returns the new index of the atom (at_ind)
    at_ind=edit_mol.AddAtom(new_atom)
    # adds the new bond, returns the new number of bonds
    edit_mol.AddBond(ati,at_ind,Chem.BondType.DOUBLE)
    # the new smile string generated is sometimes weird, but still encodes the molecule so I guess it's alright
    new_smi=Chem.MolToSmiles(edit_mol)
    #currently deliminates sterochemistry uhg..
    return new_smi
# ** Note that when converting from an RWMol to SMILE back to Mol, the atom indexing changes. This is easily visible in the jupyter sandbox **

### removes an atom from the mol
## m - the mol object to be changed
## ati - the index of the atom to be removed
def remove_atom(m,at,ati):
    edit_mol=Chem.rdchem.RWMol(m)
    #mtmp=Chem.MolFromSmiles(Chem.MolToSmiles(m))#,isomericSmiles=True
    #mw = Chem.RWMol(mtmp)
    edit_mol.RemoveAtom(ati)
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi

### Checks if the molecule is connected, i.e. without ionic bonds
# ionic bonds are represented in SMILES as a '.', i.e. [Na+].[Cl-]
def check_connected(m):
    if "." in Chem.MolToSmiles(m):
        return False
    else:
        return True

### changes a bond to a single bond
## m - the mol to be changed
## b -
## bi - the index of the bond to be changed
def make_single(m,b,bi):
    edit_mol=Chem.rdchem.RWMol(m)
    edit_mol.GetBonds()[bi].SetBondType(Chem.BondType.SINGLE)
    #change hybridization of both atoms?
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi

### changes a bond to a double bond
## m - the mol to be changed
## b -
## bi - the index of the bond to be changed
def make_double(m,b,bi):
    edit_mol=Chem.rdchem.RWMol(m)
    edit_mol.GetBonds()[bi].SetBondType(Chem.BondType.DOUBLE)
    #mw.SetBondType(mw.GetBonds()[bi],Chem.BondType.DOUBLE)
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi

### changes a bond to a triple bond
## m - the mol to be changed
## b -
## bi - the index of the bond to be changed
def make_triple(m,b,bi):
    edit_mol=Chem.rdchem.RWMol(m)
    edit_mol.GetBonds()[bi].SetBondType(Chem.BondType.TRIPLE)
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi

### removes a bond; results in the molecule being split in two, as it doesn't delete any atoms
## m - the mol to be edited
## bi - index of first atom of bond (order doesn't matter)
## bj - index of second atom of bond
def remove_bond(m,bi,bj):
    edit_mol=Chem.rdchem.RWMol(m)
    edit_mol.RemoveBond(bi,bj)
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi

### Checks if a given bond b is in between double bonds AND in a ring. Returns True/False
## m - the mol in question, although not used in the function
## b - the bond in question
def is_bet_dub_and_ring(m,b):
    bond_beg=[btype(xxx) for xxx in b.GetBeginAtom().GetBonds() if xxx !=b]
    bond_end=[btype(xxx) for xxx in b.GetEndAtom().GetBonds() if xxx !=b]
    in_bet_dub=all([any([xxx == 'DOUBLE' for xxx in bond_beg]),any([xxx == 'DOUBLE' for xxx in bond_end])])
    in_ring=b.IsInRing()
    return all([in_bet_dub,in_ring])

### adds a single bond between two atoms already present with the molecule (for intramolecular bonds)
## m - the mol to be changed
## ai - index of the first atom
## aj - index of the second atom
def add_single(m,ai,aj):
    edit_mol=Chem.rdchem.RWMol(m)
    edit_mol.AddBond(ai,aj,order=Chem.rdchem.BondType.SINGLE)
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi

### adds a double bond between two atoms already present with the molecule (for intramolecular bonds)
## m - the mol to be changed
## ai - index of the first atom
## aj - index of the second atom
def add_double(m,ai,aj):
    edit_mol=Chem.rdchem.RWMol(m)
    edit_mol.AddBond(ai,aj,order=Chem.rdchem.BondType.DOUBLE)
    new_smi=Chem.MolToSmiles(edit_mol)
    return new_smi


### ensures that rings do not break the predefined ring rules; returns 1 if passes, 0 if fails
## smii - smile string
def check_cylces(smii):
    ## ring rules:
    max_n_ring=4
    min_ring_size=5
    max_ring_size=7
    test_mo=Chem.MolFromSmiles(smii)
    rinf=test_mo.GetRingInfo()
    ring_ok=1
    #check number of rings
    nring=rinf.NumRings()
    if nring > max_n_ring:
        return 0
    #check minimum ring size
    ring_sizes=[len(aring) for aring in rinf.AtomRings()]
    if min(ring_sizes)<min_ring_size:
        return 0
    if max(ring_sizes)>max_ring_size:
        return 0
    return 1

### writes the configuration file confi.txt for a vina run
## receptor - the name of the receptor file (string, .pdbqt)
## ligand - the name of the ligand file (string, .pdbqt)
## iiter - key to what molecule this configuration corresponds to (iiter because iter is a python keyword apparently)
## fname - name of the configuration file, defaults to config
## size - the size of the search grid, in Angstroms
## exhuastiveness - how exhaustive the search is
## center_x,y,z - coordinates of the grid's center
## out - default name of the output file
## cpu - number of CPUs to use; deafult is to detect the number available and use those
## num_modes  max number of binding modes to generate
## score_only - can the search space be omitted? true/false
def configure(receptor,ligand,iteration='test',fname="config",size=20,exhaustiveness=32,center_x=12.95,center_y=15.76,center_z=2.28,out='vina_outs/out',cpu=16,num_modes=9,score_only=False):
    config_i=fname+"_"+iteration+'.txt'
    out_name=out+"_"+iteration
    os.chdir('configs')
    with open(config_i,'w') as f:
        f.write('receptor = '+receptor+'\n')
        f.write('ligand = '+ligand+'\n\n')
        f.write('out = '+out_name+'.pdbqt'+'\n\n')
        f.write('center_x = '+str(center_x)+'\n')
        f.write('center_y = '+str(center_y)+'\n')
        f.write('center_z = '+str(center_z)+'\n\n')
        f.write('size_x = '+str(size)+'\n')
        f.write('size_y = '+str(size)+'\n')
        f.write('size_z = '+str(size)+'\n\n')
        f.write('exhaustiveness = '+str(exhaustiveness)+'\n\n')
        f.write('cpu = '+str(cpu)+'\n\n')
        f.write('num_modes = '+str(num_modes))
        if score_only:
            f.write('\n\nscore_only = true')
    os.chdir('../')
    return config_i,out_name

### prepares the receptor for docking with vina, returns the file name in .pdbqt format
## prot_file - the name of the pdb file of the protein/pocket/receptor
def prepare_receptor(prot_file):
    sp.call("pythonsh ~/Desktop/Summer_Research_2022/ChemHopper/docking/prepare_receptor4.py -r "+prot_file+" -v",shell=True)
    return prot_file+'qt'

### prepares the ligand for docking with vina, returns the file name in .pdbqt format
## lig_file - the name of the pdb file of the mol/ligand
def prepare_ligand(lig_file):
    os.chdir('mols')
    sp.call("pythonsh ~/Desktop/Summer_Research_2022/ChemHopper/docking/prepare_ligand4.py -l "+lig_file+" -v",shell=True)
    os.chdir('../')
    return lig_file+'qt'

# propofol in the bound position is the reference ligand
# here I will prepare the refernce like any other ligand and run vina on it to see if I get similar results with propofol created from a smile
# ref_lig="reference_ligand.pdb"
# sp.call("pythonsh /Desktop/Summer_Research_2022/ChemHopper/docking/prepare_ligand4.py -l "+ref_lig+" -v",shell=True)
# ref_lig+='qt'
# configure(prot_file,ref_lig,out='aligned_results')
# sp.call("vina --config=config.txt",shell=True)


#### MY DOCK_IT()
### attempts to dock the molecule
## To do this, there a few unique steps. Create and prepare ligand pdb, prepare the receptor pdb (may only have to do this once if working with the same pocket), run vina (remember to log), and 
## return a docking score. First, we need to create a pdb file from the smile string we are inputting. That rewuires making a mol object, addding Hs, Embedding it (which acts like a legitimacy
## test for the molecules), and then obtimizing the forcefield. Then we prepare this pdb file (convert to pdbqt).
## conf_file - the configuration file for vina
## smile - smile string of the mol to be docked
def dock_it(smile,iiter='test'):
    # try:
    #     iiter=str(gen)+"."+str(i)
    # except NameError:
    #     iiter='test'
    mol=Chem.MolFromSmiles(smile)
    mh = Chem.AddHs(mol)
    # Generates the 3D coordinates for the molecule; serves as a test for generated molecules, as RDkit will not create coordinates if the mol is not chemically viable
    global no_embed
    if AllChem.EmbedMolecule(mh,useRandomCoords=True,maxAttempts=100)!=0: # it would be a good idea to track how many molecules cannot be embedded
        no_embed.append((iiter,smile))

    # I don't think I need to align within rdkit, probably just need to run prepare_ligand4.py
    #Chem.rdMolAlign.AlignMol(mh,ref_lig)
    # Merck Molecular Force Field: https://en.wikipedia.org/wiki/Merck_molecular_force_field
    # uses MMFF to optimize a molecule's structure
        # RETURNS: 0 if the optimization converged, -1 if the forcefield could
        # not be set up, 1 if more iterations are required.
    Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(mh,maxIters=5000)
    # create the ligand file with unique name
    lig_pdb="mols/mol_"+iiter+".pdb"
    # creates a PDB file for the ligand 
    Chem.rdmolfiles.MolToPDBFile(mh,lig_pdb)
    # prepares ligand and returns lig_pdbqt
    lig_pdbqt=prepare_ligand(lig_pdb)
    configuration,out_name=configure(prot_pdbqt,lig_pdbqt,iiter,exhaustiveness=exhaust)
    # runs vina and logs results
    logfile="logs/log_"+iiter+".txt"
    with open(logfile,'w') as log:
        run=sp.run("vina_1.2.3_linux_x86_64 --config=configs/"+configuration,shell=True,stdout=log)
    # splitting output
    sp.call("vina_split --input "+out_name+'.pdbqt',shell=True)
    # deleting the original out file from vina plus all but the best modes from vina_split
    print("OUT_NAME:",out_name)
    os.remove(out_name+'.pdbqt')
    # systematically deleting all output ligands other than the best (ligand_1)
    count=2
    fname=out_name+f'_ligand_{count}.pdbqt'
    while os.path.isfile(fname):
        os.remove(fname)
        count+=1
        fname=out_name+f'_ligand_{count}.pdbqt'
    # renaming the best ligand output
    best_out=out_name+"_ligand_1.pdbqt"
    split=best_out.rsplit("_1",1)
    os.rename(best_out,''.join(split))
    # Opening the log file to read in the best affinity
    with open(logfile,"r") as results:
        lines=results.readlines()[::-1]
    best_mode=[i for i in lines if re.match('\s+1\s',i)][0]
    print("9th line from the end:\n",best_mode)
    s=0
    while len(best_mode.split(' ')[s])<=1:
        s+=1
    print("THE SPLIT:",best_mode.split(' ')[s])
    affinity=float(best_mode.split(' ')[s])
    print(smile,affinity)
    return(affinity)


### returns a list containg the grid images of the molecules and their correpsonding affinities
## mols - the list of mol objects to be printed
## affinities - the corresponding affinities of the mols - SAME INDEXING IS ASSUMED
## row_len - # of mols to be printed per grid row
## rows_per_page - # of rows of mols to be printed per page
def save_grids(mols,affinities,row_len=8,col_len=10,highlight=True):
    num_mols=len(mols)
    nrows=num_mols/row_len
    if not nrows.is_integer():
        nrows=int(nrows+1)
    npgs=nrows/col_len
    if not npgs.is_integer():
        npgs=int(npgs+1)
    mols_per_pg=row_len*col_len
    mol_chunks=[]
    affin_chunks=[]
    for i in range(0, len(mols), mols_per_pg):
        mol_chunks.append(mols[i:i+mols_per_pg])
        affin_chunks.append(affinities[i:i+mols_per_pg])
    # returns index of the top molecule in uniq_desc
    best_index=uniq_desc.index(best_path[-1][0])
    top=uniq_mols[best_index]
    grids=[]
    for i in range(int(npgs)):
        if highlight==True:
            img=Draw.MolsToGridImage(mol_chunks[i],molsPerRow=row_len,subImgSize=(350,350),legends=[str(a) for a in affin_chunks[i]],highlightAtomLists=[mol.GetSubstructMatch(top) for mol in mol_chunks[i]])
        else:
            img=Draw.MolsToGridImage(mol_chunks[i],molsPerRow=row_len,subImgSize=(350,350),legends=[str(a) for a in affin_chunks[i]])
        img.save("grids/gen_"+str(gen+1)+"_grid_"+str(i)+".png")
        grids.append(img)
    return grids

######################################################### BEGIN SCRIPT ########################################################################################

# making directories top organize our output files
os.mkdir("logs")
os.mkdir("np_objs")
os.mkdir("configs")
os.mkdir("grids")
os.mkdir("vina_outs")
os.mkdir("mols")

# setting up log format (only needs to be done once)
logging.basicConfig(filename='error_log.log', filemode='w',format='%(asctime)s -%(name)s:%(levelname)s- %(message)s')


# preparing the receptor
prot_pdb='1e7a_aligned.pdb'
prot_pdbqt=prepare_receptor(prot_pdb)
prot_pdbqt=prot_pdb+'qt'

### Propofol SMILE string - parent molecule, patient 0, generations[0][0], etc
parent_name='Propofol'
parent_0='CC(C)C1=C(C(=CC=C1)C(C)C)O'
#propofol=Chem.MolFromSmiles(propofol)
#Draw.MolToFile(propofol,"propofol.png")

# defining the exhaustiveness for Autodock Vina, now using keyboard input:
print("\n\n\nWelcome to Chemical Space!\nHow deep will you go?")

# defining the depth of our search
str_depth=input("depth = ")
depth=int(str_depth)

print("\nHow exhaustive will your search be?")

# taking in exhaustiveness hyperparameter
str_exhaust=input("exhaustiveness = ")
exhaust=int(str_exhaust)

print("\n\nGreat! Let's begin")

print('\n\n############\n DEPTH = '+str(depth)+' \n############\n')

# for keeping track of molecules that do not embed
no_embed=[]

# create generations (master) list, add the original parent results as the frist entry
generations=[[(parent_0,dock_it(parent_0,'0.0'))]]
# for testing purposes (will not call dock_it())
# generations=[[(parent_0,ran.uniform(-6,0))]]

# Keeping a set containing all UNIQUE molecules created and tested
master_set=set(parent_0)

# creates best_path list (top performers from each search), which will be of length depth+1, due to the parent molecule being the first entry
best_path=[generations[0][0]]

# Begin iterations
for gen in range(depth):
    # this is the temporarry array that will hold the molecules created throughout the editing process
    all_desc=[]
    # this defines the current parent as the best molecule from the previous iteration
    parent_smi=best_path[-1][0]
    parent_mol=Chem.MolFromSmiles(parent_smi)
    # Not sure what this line does specifically, but it is necessary so that the function degree() returns the proper bond counts.
    Chem.rdmolops.Kekulize(parent_mol,clearAromaticFlags=True)

    # Begin working through the molecular editing process - iterating through ATOMS
    for i,atom in enumerate(parent_mol.GetAtoms()):
        at_deg=degree(atom)
        atom_n=aname(atom)
        # checkpt
        #### Try to Mutate Atoms
        if at_deg == 4:#could be C or N
            if atom_n=="C":
                all_desc.append(mutate(parent_mol,atom,i,"N",1))#mutates atom in molecule to "N"
            if atom_n=="N":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))   
        if at_deg == 3:#could be a C or N
            if atom_n=="C":
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                #all_desc.append(mutate(parent_mol,atom,i,"N",1))#could be protonated?
            if atom_n=="N":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
        if at_deg <=2: #could be a C or N or O
            if atom_n=="C":
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                #all_desc.append(mutate(parent_mol,atom,i,"N",1))#could be protonated?
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
            if atom_n=="N":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
            if atom_n=="O":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                #all_desc.append(mutate(parent_mol,atom,i,"N",1))#could be protonated?
                # if the atom is a halogen, it has degree of 1
            if atom_n in ["F","Cl"]:
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
        if at_deg == 1:
            all_desc.append(mutate(parent_mol,atom,i,"F",0))
            all_desc.append(mutate(parent_mol,atom,i,"Cl",0))
        

        
        #### Try to add a single atom w/ new atom bonded to the existing atom on molecule
        if at_deg == 3:# could be a C or N, could add a single bond C,N,or O
            #assume not an "O"
            if atom_n =="C":
                #FIX TO ALLOW N(-C)4, non-heteroatom quaternary nitrogen
                all_desc.append(add_single_bond(parent_mol,atom,i,"N"))#adds single bonded "N" to atom in molecule
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"O"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"F"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"Cl"))
        if at_deg ==2:# could be a C or N or O, but dont add to O
            if atom_n == "C":
                #add single bonded atom
                all_desc.append(add_single_bond(parent_mol,atom,i,"N"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"O"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"F"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"Cl"))
                #add double bonded atom but only if no double bond on c already
                if all([btype(bond)!="DOUBLE" for bond in atom.GetBonds()]):
                    all_desc.append(add_double_bond(parent_mol,atom,i,"N"))
                    all_desc.append(add_double_bond(parent_mol,atom,i,"C"))
                    all_desc.append(add_double_bond(parent_mol,atom,i,"O"))
            elif atom_n=="N":
                all_desc.append(add_single_bond(parent_mol,atom,i,"N"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))
        if at_deg ==1:# could be a C or N or O, but dont add double bond to O
            if atom_n == "C" :
                #add single bonded atom
                all_desc.append(add_single_bond(parent_mol,atom,i,"N"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"O"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"F"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"Cl"))
                #add double bonded atom
                all_desc.append(add_double_bond(parent_mol,atom,i,"N"))
                all_desc.append(add_double_bond(parent_mol,atom,i,"C"))
                all_desc.append(add_double_bond(parent_mol,atom,i,"O"))
            elif atom_n =="N":
                all_desc.append(add_single_bond(parent_mol,atom,i,"N"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))
                all_desc.append(add_double_bond(parent_mol,atom,i,"C"))
            elif atom_n=="O":
                #add single bonded atom
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))
                all_desc.append(add_single_bond(parent_mol,atom,i,"O"))
        
        #### Try removing atom
        if at_deg==1:#no worries
            all_desc.append(remove_atom(parent_mol,atom,i))
        elif at_deg>=2:#only can remove if the resulting molecule is connected (doesn't have floating framents)
            tmp_mols=remove_atom(parent_mol,atom,i)
            tmp_mol=Chem.MolFromSmiles(tmp_mols)
            if check_connected(tmp_mol)==1:
                all_desc.append(tmp_mols)
        
    ### Begin iteration through BONDS
    idone=0
    for i,bond in enumerate(parent_mol.GetBonds()):
        #ignore bonds to halogens
        if parent_mol.GetAtoms()[bond.GetBeginAtomIdx()].GetAtomicNum() not in halogens and parent_mol.GetAtoms()[bond.GetEndAtomIdx()].GetAtomicNum() not in halogens:
            #determine elements in bond
            atom1=atomic_names[allowed_atomic.index(parent_mol.GetAtoms()[bond.GetBeginAtomIdx()].GetAtomicNum())]
            atom2=atomic_names[allowed_atomic.index(parent_mol.GetAtoms()[bond.GetEndAtomIdx()].GetAtomicNum())]
            #determine degree of atoms in bond
            degree1=degree(parent_mol.GetAtoms()[bond.GetBeginAtomIdx()])
            degree2=degree(parent_mol.GetAtoms()[bond.GetEndAtomIdx()])
            #check bond type and apply chemical rules
            if btype(bond) =="SINGLE":
                #check if not (between two double bonds and in ring)
                if not is_bet_dub_and_ring(parent_mol,bond):
                    #print(idone,bond)
                    #try making it a double bond
                    #print(i+1,is_aromatic(parent_mol,bond))
                    if (atom1 == "C") and (atom2 == "C"):#(atom1 == "C" or atom1 == "N") and (atom2 == "C" or atom2 == "N"):
                        if degree1<=3 and degree2<=3:#only two <terterary carbons can form a double bond (ignoring N+ )
                            all_desc.append(make_double(parent_mol,bond,i))
                            idone+=1
                    elif (atom1 == "O") and (atom2 == "C"):# a single bond with O-C
                        if degree1==1 and degree2<=3:#check if O is single degree and C is <terterary
                            all_desc.append(make_double(parent_mol,bond,i))# hydroxyl -> carbonyl
                    elif (atom1 == "C") and (atom2 == "O"):# a single bond with C-O
                        if degree1<=3 and degree2==1:#check if O is single degree and C is <terterary
                            all_desc.append(make_double(parent_mol,bond,i))# hydroxyl -> carbonyl      
                    elif atom1 == "C" and atom2 == "N":# a single bond with C-N
                        if degree1<=3 and degree2<=2:#check if N is 2 degree and C is <terterary
                            all_desc.append(make_double(parent_mol,bond,i))# 
                    elif atom1 == "N" and atom2 == "C":# a single bond with N-C
                        if degree1<=2 and degree2<=3:#check if N is 2 degree and C is <terterary
                            all_desc.append(make_double(parent_mol,bond,i))# 
                    #try making it a triple bond
                    #-C#N or -C#C or -C#C-
                    if (atom1 =="N" and degree1==1 and atom2=="C" and degree2<=2) or (atom1 =="C" and degree1<=2 and atom2=="N" and degree2==1) :# -C-N -> -C#N
                        all_desc.append(make_triple(parent_mol,bond,i))
                    elif atom1 =="C" and atom2=="C" and degree1<=2 and degree2<=2:
                        all_desc.append(make_triple(parent_mol,bond,i))
            elif btype(bond) =="DOUBLE":
                #try making a single bond, should be fine. but maybe diols are unlikely?
                all_desc.append(make_single(parent_mol,bond,i))
                #try making it a triple bond if not in ring
                if (atom1 == "C" or atom1 == "N") and (atom2 == "C" or atom2 == "N") and not (atom1 =="N" and atom2 == "N") and not bond.IsInRing():
                    if degree1<=2 and degree2<=2:
                        all_desc.append(make_triple(parent_mol,bond,i))
            elif btype(bond) =="TRIPLE":
                #try making a single bond
                all_desc.append(make_single(parent_mol,bond,i))
                #try making a double bond, it already was a triple, so should be fine
                all_desc.append(make_double(parent_mol,bond,i))
                #try removing the bond
            tmp_mols=remove_bond(parent_mol,bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())
            tmp_mol=Chem.MolFromSmiles(tmp_mols)
            if check_connected(tmp_mol)==1:
                all_desc.append(tmp_mols)

    ### Have yet to work on rules for intramolecular bonds - i.e. making rings
    
    # code to remove duplicate smile strings for current generation and prvious generation. Since each generation will be unique, there is no need to continuous check ALL previous generations each
    # time we go through the loop; this line of code has a recursive structure
    uniq_desc=list(set(all_desc).difference(master_set))
    master_set=master_set.union(uniq_desc)

    # time to dock each descendant in uniq_desc and record results
    results=[] # the results for current GEN ONLY

    ### CHECKPOINT

    print('\n################\n DEPTH = '+str(gen+1)+'/'+str(depth)+' \n################\n')

    print("# OF UNIQUE DESC:",str(len(uniq_desc)))
    for i,smi in enumerate(uniq_desc):
        iiter=str(gen+1)+"."+str(i)
        print("ITER =",str(iiter))
        try:
            results.append((smi,dock_it(smi,iiter)))
        except Exception as e:
            logging.exception("IITER: "+str(iiter)+"; TYPE: "+type(e).__name__)
        # # for testing purposes (will not run dock_it())
        # try:
        #     if i%7==0:
        #         ran.uniform(whupsies)
        #     else:
        #         results.append((smi,ran.uniform(-12,0)))
        # except Exception as e:
        #     logging.exception("IITER: "+str(iiter)+"; TYPE: "+type(e).__name__)

    # record this generations results to the master list:
    generations.append(results)

    # pulling just docking scores from the results list
    affinities=[x[1] for x in results]
    # getting the index of the best performer in results
    best_index=affinities.index(min(affinities))
    # adding the best perfomer to the best_path list
    best_path.append(results[best_index])

    ### produce an image of the molecules
    # first create a mol list
    uniq_mols=[]
    for smi in uniq_desc:
        mol=Chem.MolFromSmiles(smi)
        Chem.rdmolops.Kekulize(mol,clearAromaticFlags=True)
        uniq_mols.append(mol)
    # img=Draw.MolsToGridImage(uniq_mols,molsPerRow=int(len(uniq_mols)/10),subImgSize=(350,350),legends=[str(a) for a in affinities],maxMols=len(uniq_mols))

    save_grids(uniq_mols,affinities,8,10)

    print("BEST PATH:\n",best_path)

    np.save("np_objs/generations.npy",np.array(generations,dtype=object))
    np.save("np_objs/best_path.npy",np.array(best_path))
    np.save("np_objs/no_embed.npy",np.array(no_embed))

# tabulating the best path
df=pd.DataFrame(best_path)
table=tab(df,headers=["Generation","Best Molecule","Binding Affinity"])

# Getting the runtime of the program
t1=time.time()
elapsed=t1-t0
runtime=time.strftime("%H:%M:%S",time.gmtime(elapsed))

print('\nRUNTIME: '+runtime)

# writing a log file for the overall run
with open('RUN_LOG.txt','w') as log:
    log.write('STARTING MOL: '+parent_name+'\n')
    log.write('    SMILE: '+parent_0+'\n')
    log.write('RECEPTOR: 1e7a_aligned\n')
    log.write('    POCKET: S3\n')
    log.write('DEPTH: '+str(depth)+'\n')
    log.write('EXHAUSTIVENESS: '+str(exhaust)+'\n\n')
    log.write('RUNTIME: '+runtime+'\n\n')
    log.write('Number of mols that could not be embedded: '+str(len(no_embed))+'\n\n')
    log.write(table)
