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

############# ChemTools.py ##############
# Includes all the functions neccessary for the ChemHopper algorithm
# program should work if the script imports this file, assuming I didn't make a mistake
#########################################

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

### adds a triple bond to a new atom and returns a new SMILE string
## m - the mol object to be changed
## at - 
## ati - the index of the atom on the old mol to which the new atom will be attached
## newat - type of new atom to be added; given as atomic symbol (i.e. 'C', 'N', etc.)
def add_triple_bond(m,at,ati,newat):
    # making an editable version of the mol
    edit_mol=Chem.rdchem.RWMol(m)
    # creates a new atom based on the corresponding atomic number of the the atomic symbol
    new_atom=Chem.rdchem.Atom(allowed_atomic[atomic_names.index(newat)])
    # adds the new atom to edit_mol and returns the new index of the atom (at_ind)
    at_ind=edit_mol.AddAtom(new_atom)
    # adds the new bond, returns the new number of bonds
    edit_mol.AddBond(ati,at_ind,Chem.BondType.TRIPLE)
    # the new smile string generated is sometimes weird, but still encodes the molecule so I guess it's alright
    new_smi=Chem.MolToSmiles(edit_mol)
    #currently deliminates sterochemistry uhg..
    return new_smi

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
def configure(receptor,ligand,iteration='test',fname="config",size=20,exhaustiveness=8,center_x=12.95,center_y=15.76,center_z=2.28,out='vina_outs/out',cpu=16,num_modes=9,score_only=False):
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
    sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_receptor4.py -r "+prot_file+" -v",shell=True)
    return prot_file+'qt'

### prepares the ligand for docking with vina, returns the file name in .pdbqt format
## lig_file - the name of the pdb file of the mol/ligand
def prepare_ligand(lig_file):
    os.chdir('mols')
    sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_ligand4.py -l "+lig_file+" -v",shell=True)
    os.chdir('../')
    return lig_file+'qt'

    #### MY DOCK_IT()
### attempts to dock the molecule
## To do this, there a few unique steps. Create and prepare ligand pdb, prepare the receptor pdb (may only have to do this once if working with the same pocket), run vina (remember to log), and 
## return a docking score. First, we need to create a pdb file from the smile string we are inputting. That rewuires making a mol object, addding Hs, Embedding it (which acts like a legitimacy
## test for the molecules), and then obtimizing the forcefield. Then we prepare this pdb file (convert to pdbqt).
## conf_file - the configuration file for vina
## lig_smile - smile string of the mol to be docked
## prot_pdbqt - name of protein file to be docked to
def dock_it(lig_smile,prot_pdbqt,exhaustiveness=8,iiter='test'):
    # try:
    #     iiter=str(gen)+"."+str(i)
    # except NameError:
    #     iiter='test'
    mol=Chem.MolFromSmiles(lig_smile)
    mh = Chem.AddHs(mol)
    # Generates the 3D coordinates for the molecule; serves as a test for generated molecules, as RDkit will not create coordinates if the mol is not chemically viable
    global no_embed
    if AllChem.EmbedMolecule(mh,useRandomCoords=True,maxAttempts=100)!=0: # it would be a good idea to track how many molecules cannot be embedded
        no_embed.append((iiter,lig_smile))

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
    configuration,out_name=configure(prot_pdbqt,lig_pdbqt,iiter,exhaustiveness=exhaustiveness)
    # runs vina and logs results
    logfile="logs/log_"+iiter+".txt"
    with open(logfile,'w') as log:
        run=sp.run("vina --config=configs/"+configuration,shell=True,stdout=log)
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
    print(lig_smile,affinity)
    return(affinity)


### returns a list containg the grid images of the molecules and their correpsonding affinities
## mols - the list of mol objects to be printed
## top - best molecule for the generation (to be highlighted)
## affinities - the corresponding affinities of the mols - SAME INDEXING IS ASSUMED
## gen - the generation number (needed for th e)
## row_len - # of mols to be printed per grid row
## rows_per_page - # of rows of mols to be printed per page
def save_grids(mols,top,affinities,gen,row_len=8,col_len=10,highlight=True):
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
    grids=[]
    for i in range(int(npgs)):
        if highlight==True:
            img=Draw.MolsToGridImage(mol_chunks[i],molsPerRow=row_len,subImgSize=(350,350),legends=[str(a) for a in affin_chunks[i]],highlightAtomLists=[mol.GetSubstructMatch(top) for mol in mol_chunks[i]])
        else:
            img=Draw.MolsToGridImage(mol_chunks[i],molsPerRow=row_len,subImgSize=(350,350),legends=[str(a) for a in affin_chunks[i]])
        img.save("grids/gen_"+str(gen+1)+"_grid_"+str(i)+".png")
        grids.append(img)
    return grids
