from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import QED
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from faerun import Faerun
import logging
import matplotlib.pyplot as plt
from matplotlib.table import Cell
from mhfp.encoder import MHFPEncoder
import networkx as nx
import numpy as np
import networkx as nx
import os
import pandas as pd
import platform as pf
import random as ran
import re
from sklearn.manifold import TSNE
import sys
import subprocess as sp
from tabulate import tabulate as tab
from time import time


###JMR. Carefuly im disabling this but if you comment this out you can see all the warnings being generated.
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')                                                                                                                                                           

############# ChemTools.py ##############
# Includes all the functions neccessary for the ChemHopper algorithm
# program should work if the script imports this file, assuming I didn't make a mistake
#########################################

### define some chemical data

# Corresponding atomic numbers of atomic_names
allowed_atomic=[6,7,8,9,15,16,17,35]
atomic_names=["C","N","O","F","P","S","Cl","Br"]
#Have you validated the Phosphorus chemistry? I would be careful here and think about the chemistry more. Many of the molecules when i look at the chemical space graph have wild phosphorus and sulfer chemistry 
halogens=[9,17,35]
halo_names=["F","Cl","Br"]
type_i=[0,1,2,1.5]
types=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]

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

### adds a single intramolecular bond and returns a new SMILE string
## m - the mol object to be changed
## ati - the index of the first atom
## atj - the index of the second atom
def add_single_intrabond(m,ati,atj):
    # RWMol is a subclass of the Mol object class that is basically an editable form of a Mol
    ## https://herongyang.com/Cheminformatics/RDKit-rdkit-Chem-rdchem-RWMol-Class.html#:~:text=rdkit.Chem.rdchem.RWMol%20represents%20a%20molecule%20class%20with%20additional%20read,some%20useful%20methods%20provided%20in%20the%20RWMol%20sub-class.
    edit_mol=Chem.rdchem.RWMol(m)
    # adds the new bond, returns the new number of bonds
    edit_mol.AddBond(ati,atj,Chem.BondType.SINGLE)
    # the new smile string generated is sometimes weird, but still encodes the molecule so I guess it's alright
    new_smi=Chem.MolToSmiles(edit_mol)
    #currently deliminates sterochemistry uhg..
    #check out https://github.com/rdkit/rdkit/issues/2627
    #exit()
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

### def extract_intersection - gets the extracted elemenets between two lists 
##### list master - the list to be referred to for intersection function (no duplicates)
##### list target - the list that will get its intersection removed (no duplicates); should be all molecules generated from nextGen())
### returns:
##### list uniq - the list of unique elements in target (elements that are not in master)
##### list intersect -  the list of elements in target that are also in master (the intersection)
def extract_intersection(target, master):
    # this if statement will make sure master_list has no duplicates (it shouldn't at this point - this is a failsafe)
    if (len(master)!=len(set(master))):
        print("Duplicates in master_list detected!")
        master = list(dict.fromkeys(master))
    uniq = [i for i in target if i not in master]
    intersect = [i for i in target if i in master]
    return uniq, intersect

### def make_dirs(dirs)
# makes the directories needed for the Gradient Descent Program
### string list dirs - list of directories to make
def make_dirs(dirs):
    for d in dirs:
        if os.path.isdir(d):
            continue
        else: os.mkdir(d)
    return

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
def configure(receptor,ligand,iteration='test',fname="config",size=20,exhaustiveness=8,center_x=12.95,center_y=15.76,center_z=2.28,out='vina_outs/out',cpu=16,num_modes=9,seed=0,verbosity=1,score_only=False):
    config_i=fname+"_"+iteration+'.txt'
    out_name=out+"_"+iteration
    os.chdir('configs')
    f = open(config_i,'w')
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
    f.write('num_modes = '+str(num_modes)+'\n\n')
    f.write('seed = '+str(seed)+'\n\n')
    f.write('verbosity = '+str(verbosity))
    if score_only:
        f.write('\n\nscore_only = true')
    f.close()
    os.chdir('../')
    return config_i,out_name

### prepares the receptor for docking with vina, returns the file name in .pdbqt format
## prot_file - the name of the pdb file of the protein/pocket/receptor
def prepare_receptor(prot_file):
    if pf.system()=='Linux':
        sp.call("pythonsh ~/Desktop/ChemHopper/docking/prepare_receptor4.py -r "+prot_file+" -v",shell=True)
    elif pf.system()=='Windows':
        sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_receptor4.py -r "+prot_file+" -v",shell=True)
    return prot_file+'qt'

### prepares the ligand for docking with vina, returns the file name in .pdbqt format
## lig_file - the name of the pdb file of the mol/ligand
def prepare_ligand(lig_file):
    os.chdir('mols')
    if pf.system()=='Linux':
        sp.call("pythonsh ~/Desktop/ChemHopper/docking/prepare_ligand4.py -l "+lig_file+" -v",shell=True)
    elif pf.system()=='Windows':
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
    configuration,out_name=configure(prot_pdbqt,lig_pdbqt,iiter,exhaustiveness=exhaustiveness,seed=1,verbosity=2)
    # runs vina and logs results
    logfile="logs/log_"+iiter+".txt"
    log = open(logfile, 'w')
    if pf.system()=='Linux':
        run=sp.run("vina_1.2.3_linux_x86_64 --config=configs/"+configuration,shell=True,stdout=log)
    elif pf.system()=='Windows':
        run=sp.run("vina --config=configs/"+configuration,shell=True,stdout=log)
    log.close()
    # splitting output
    sp.call("vina_split --input "+out_name+'.pdbqt',shell=True)
    # deleting the original out file from vina plus all but the best modes from vina_split
    print("OUT_NAME:",out_name)
    os.remove(out_name+'.pdbqt')
    # systematically deleting all output ligands other than the best (ligand_1)
    count=2
    fname = '{}_ligand_{}.pdbqt'.format(out_name, count)
    while os.path.isfile(fname):
        os.remove(fname)
        count+=1
        fname = '{}_ligand_{}.pdbqt'.format(out_name, count)
    # renaming the best ligand output
    best_out=out_name+"_ligand_1.pdbqt"
    split=best_out.rsplit("_1",1)
    os.rename(best_out,''.join(split))
    # Opening the log file to read in the best affinity
    results = open(logfile, 'r')
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


# Ring-building chunk to add to nextGen() ***IN-PROGRESS****
def buildRings(smile):
    new_smis = []
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.rdchem.RWMol(mol)
    # kekulize so that aromatic bonds are treated as alternative double/single bonds, but keep aromatic flags on atoms
    Chem.rdmolops.Kekulize(mol,clearAromaticFlags=False) # clearAromaticFlags default param is False btw
    for i, atom in enumerate(mol.GetAtoms()):
        # use i+1 to avoid i=j; this would cause GetShortestPath() to throw an error
        # You might think it would be efficient to start the enumeration here at i+4, that way we don't waste time generating paths
        # between atoms that are less than 4 bonds away; good thinking, but doesn't work due to the way rdkit assigns atom indices.
        # These indices quickly become arbitrary as the complexity of the molecule increases, which means that the difference in
        # indices does not relate to the distance between atoms. For instance, for propofol, atoms 3 and 13 are seperated by 4 bonds,
        # but the difference in indices is 10.
        for j in range(i+1, mol.GetNumAtoms()):
            # a list of atoms to traverse to get from atom i to j (includes i and j)
            path = Chem.GetShortestPath(mol, i, j)
            # the number of bonds that seperate atom i from j
            distance = len(path)-1
            # only try to from 5 or 6 membered rings (remember distance is the # of bonds seperating 2 atoms)
            if distance==4 or distance==5: # the largest filter - goes first
                # list of path atoms, but converted to True or False depending wether or not each path atom is in a ring
                path_atoms_in_ring = [mol.GetAtomWithIdx(a).IsInRing() for a in path]
                # Rules for paths over aromatic rings: if any atoms within the path are flagged as aromatic...
                if any(mol.GetAtomWithIdx(a).GetIsAromatic()==True for a in path):
                    # getting ring info from the mol
                    rings = mol.GetRingInfo()
                    # list containing the number of rings an atom is part of for all the atoms in path
                    num_rings_atoms_are_in = [rings.NumAtomRings(a) for a in path]
                    # valid building pattern for 3 ring atoms
                    val_pat = [1,2,1]
                    if sum(path_atoms_in_ring)==3:
                        if any(num_rings_atoms_are_in[x:x+len(val_pat)] == val_pat for x in range(len(num_rings_atoms_are_in)-len(val_pat)+1)):
                            new_smis.append(add_single_intrabond(mol,i,j))
                    if sum(path_atoms_in_ring) in [1,2]:
    #                     print(path_atoms_in_ring)
                        new_smis.append(add_single_intrabond(mol,i,j))
                    if sum(path_atoms_in_ring)==4:
                        ai, aj = mol.GetAtomWithIdx(i).GetSymbol(), mol.GetAtomWithIdx(j).GetSymbol()
                        if (ai in ["N","O","S"]) ^ (aj in ["N","O","S"]):
                            new_smis.append(add_single_intrabond(mol,i,j))
                # Rules for paths over non-aromatic rings
                # ...if no atoms are aromatic, but still in a ring
                # if the number of path atoms in a non-aromatic ring are exactly 1,2,3, or 4 (these rules apply to 5 and 6-membered rings)
                elif sum(path_atoms_in_ring) in [1,2,3,4]:
                    new_smis.append(add_single_intrabond(mol,i,j))
                # Rules for paths over non-rings (atom chains)
                # ... if path atoms are neither aromatic or in a ring
                elif sum(path_atoms_in_ring)==0:
                    new_smis.append(add_single_intrabond(mol,i,j))
    # aligning smiles, removing bad smiles, removing duplicates
    return tightenSmileList(new_smis)

def tightenSmileList(new_smis):
    ## the below lines 'align' smile strings basically (see duplicate bug described in phenethylamine test)
    # converting all new smiles into mol objects
    new_mols = [Chem.MolFromSmiles(smi) for smi in new_smis]
    # converting all new mols back into smis - to fix the synonymous smiles bug
    corrected_new_smis = [Chem.MolToSmiles(mol) for mol in new_mols if mol is not None]
    # remove duplicates
    new_smis_uniq=list(dict.fromkeys(corrected_new_smis))
    return new_smis_uniq


### def nextGen(parent_smi,master_list)
##### string parent_smi - the SMILE string of the parent molecule
### returns:
##### string list all_desc - SMILE list of all immediate descendants
def nextGen(parent_smi):
    all_desc = []
    parent_mol=Chem.MolFromSmiles(parent_smi)
    # Not sure what this line does specifically, but it is necessary so that the function degree() returns the proper bond counts.  
    # JMR "kekulizing" means turning groups of aromatic bonds into alternating single and double bonds.
    if parent_mol == None:
        return None
    Chem.rdmolops.Kekulize(parent_mol,clearAromaticFlags=True)



    # Begin working through the molecular editing process - iterating through ATOMS
    for i,atom in enumerate(parent_mol.GetAtoms()):
        at_deg=degree(atom)
        atom_n=aname(atom)
        bond_types = [str(x.GetBondType()) for x in atom.GetBonds()]
        #### Try to Mutate Atoms
        if at_deg == 4:# could be C, S, or N
            if atom_n=="C":
                all_desc.append(mutate(parent_mol,atom,i,"N",1))# mutates atom in molecule to "N"
            if atom_n=="N":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
            if atom_n=="S":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",1))
            if atom_n!="S" and bond_types.count("DOUBLE")>=1: # if there is at least 1 double bond with deg=4, then can be sulfur
                all_desc.append(mutate(parent_mol,atom,i,"S",0))
        if at_deg == 3:#could be a C, N, or P
            if atom_n=="C":
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                #all_desc.append(mutate(parent_mol,atom,i,"N",1))#could be protonated?
            if atom_n=="N":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
            if atom_n=="P":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
            if atom_n!="P" and bond_types.count("TRIPLE")==0: # phosphorus doesn't form triple bonds as far as I know
                all_desc.append(mutate(parent_mol,atom,i,"P",0))
        if at_deg <=2: #could be a C or N or O
            if atom_n=="C":
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                #all_desc.append(mutate(parent_mol,atom,i,"N",1))#could be protonated?
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
                all_desc.append(mutate(parent_mol,atom,i,"S",0))
                all_desc.append(mutate(parent_mol,atom,i,"P",0))
            if atom_n=="N":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
                all_desc.append(mutate(parent_mol,atom,i,"S",0))
                all_desc.append(mutate(parent_mol,atom,i,"P",0))
            if atom_n=="O":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                all_desc.append(mutate(parent_mol,atom,i,"S",0))
                all_desc.append(mutate(parent_mol,atom,i,"P",0))
                #all_desc.append(mutate(parent_mol,atom,i,"N",1))#could be protonated?
            if atom_n=="P":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                all_desc.append(mutate(parent_mol,atom,i,"S",0))
                all_desc.append(mutate(parent_mol,atom,i,"O",0))       
            if atom_n=="S":
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
                all_desc.append(mutate(parent_mol,atom,i,"P",0))
                # if the atom is a halogen, it has degree of 1
            if atom_n in halo_names:
                all_desc.append(mutate(parent_mol,atom,i,"C",0))
                all_desc.append(mutate(parent_mol,atom,i,"N",0))
                all_desc.append(mutate(parent_mol,atom,i,"O",0))
                all_desc.append(mutate(parent_mol,atom,i,"P",0))
                all_desc.append(mutate(parent_mol,atom,i,"S",0))
        if at_deg <= 1:
            for element in atomic_names:
                if element!=atom_n:
                    all_desc.append(mutate(parent_mol,atom,i,element,0))
            if atom_n=="F":
                all_desc.append(mutate(parent_mol,atom,i,"Cl",0))
                all_desc.append(mutate(parent_mol,atom,i,"Br",0))
            if atom_n=="Cl":
                all_desc.append(mutate(parent_mol,atom,i,"F",0))
                all_desc.append(mutate(parent_mol,atom,i,"Br",0))
            if atom_n=="Br":
                all_desc.append(mutate(parent_mol,atom,i,"Cl",0))
                all_desc.append(mutate(parent_mol,atom,i,"F",0))

        
        #### Try to add a single atom w/ new atom bonded to the existing atom on molecule
        if atom_n=="S" and at_deg==4: # sulfur has deg 4, it can take a final double bond
            for element in ["C","N","O","S","P"]:
                all_desc.append(add_double_bond(parent_mol,atom,i,element))
        if at_deg == 3:# could be a C, N, or P
            if atom_n=="C": # phosphorus does not form 4 non-H bonds to my knowledge (nor does N)
                #FIX TO ALLOW N(-C)4, non-heteroatom quaternary nitrogen
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
            if atom_n=="P":
                for element in ["C","N","O","S","P"]:
                    all_desc.append(add_double_bond(parent_mol,atom,i,element))
        if at_deg ==2:# could be a C, N, O, S, or P but dont add to O
            if atom_n == "C":
                #add single bonded atom
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
                #add double bonded atom but only if the bond is not in a ring
                if all([not bond.IsInRing() for bond in atom.GetBonds()]):
                    for element in ["C","N","O","S","P"]:
                        all_desc.append(add_double_bond(parent_mol,atom,i,element))
            elif atom_n=="N":
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
            elif atom_n=="S":
                # sulfur can't have 4 single bonds, but can have 2 double or 1 double +2 single
                for element in ["C","N","O","S","P"]:
                    all_desc.append(add_double_bond(parent_mol,atom,i,element))
            elif atom_n=="P":
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
        if at_deg <=1:# could be a C, N, O, P, or S, but dont add double bond to O
            if atom_n == "C" :
                #add single bonded atom
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
                #add double bonded atom
                for element in ["C","N","O","S","P"]:
                    all_desc.append(add_double_bond(parent_mol,atom,i,element))
                #add triple bonded atom
                all_desc.append(add_triple_bond(parent_mol,atom,i,"N"))
                all_desc.append(add_triple_bond(parent_mol,atom,i,"C"))
            elif atom_n =="N":
                # add single bond to everything
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
                for element in ["C","N","O","S","P"]:
                    all_desc.append(add_double_bond(parent_mol,atom,i,element))
                for element in ["C","N"]:
                    all_desc.append(add_triple_bond(parent_mol,atom,i,element))
            elif atom_n=="O":
                #add single bonded atom
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
                # add double bond to carbon
                all_desc.append(add_double_bond(parent_mol,atom,i,"C"))
            elif atom_n=="P":
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
                for element in ["C","N","O","S","P"]:
                    all_desc.append(add_double_bond(parent_mol,atom,i,element))
            elif atom_n=="S":
                for element in atomic_names:
                    all_desc.append(add_single_bond(parent_mol,atom,i,element))
                for element in ["C","N","O","S","P"]:
                    all_desc.append(add_double_bond(parent_mol,atom,i,element))
            elif atom_n in halo_names:
                all_desc.append(add_single_bond(parent_mol,atom,i,"C"))      
        
        #### Try removing atom
        if at_deg==1:#no worries
            all_desc.append(remove_atom(parent_mol,atom,i))
        elif at_deg>=2:#only can remove if the resulting molecule is connected (doesn't have floating framents)
            tmp_mols=remove_atom(parent_mol,atom,i)
            tmp_mol=Chem.MolFromSmiles(tmp_mols)
            if check_connected(tmp_mol)==1:
                all_desc.append(tmp_mols)
        
    ### Begin iteration through BONDS
    for i,bond in enumerate(parent_mol.GetBonds()):
        # get the indices of the beginning and end atoms within the bond
        aidx1 = bond.GetBeginAtomIdx()
        aidx2 = bond.GetEndAtomIdx()
        #ignore bonds to halogens
        if parent_mol.GetAtoms()[aidx1].GetAtomicNum() not in halogens and parent_mol.GetAtoms()[aidx2].GetAtomicNum() not in halogens:
            #determine the names of elements in bond
            atom1=parent_mol.GetAtoms()[aidx1]
            atom2=parent_mol.GetAtoms()[aidx2]
            aname1=atomic_names[allowed_atomic.index(atom1.GetAtomicNum())]
            aname2=atomic_names[allowed_atomic.index(atom2.GetAtomicNum())]
            #determine degree of atoms in bond
            degree1=degree(parent_mol.GetAtoms()[aidx1])
            degree2=degree(parent_mol.GetAtoms()[aidx2])
            #check bond type and apply chemical rules
            if btype(bond) =="SINGLE":
                #check if not (between two double bonds and in ring)
                if not is_bet_dub_and_ring(parent_mol,bond):
                    if(aname1=="C" and aname2=="C" and degree1<=3 and degree2<=3):
                        all_desc.append(make_double(parent_mol,bond,i))
                    elif(aname1=="C" and degree1<=3) or (aname2=="C" and degree2<=3):
                        if(aname2 in ["N","P"] and degree2<=2) or (aname1 in ["N","P"] and degree1<=2):
                            all_desc.append(make_double(parent_mol,bond,i))
                        elif(aname2 in ["O","S"] and degree2==1) or (aname1 in ["O","S"] and degree1==1):
                            all_desc.append(make_double(parent_mol,bond,i))
                    elif(aname1 in ["N","P"] and degree1<=2) or (aname2 in ["N","P"] and degree2<=2):
                        if(aname1 in ["N","P"] and degree1<=2) and (aname2 in ["N","P"] and degree2<=2):
                            all_desc.append(make_double(parent_mol,bond,i))
                        elif(aname1 in ["O","S"] and degree1==1) or (aname2 in ["O","S"] and degree2==1):
                            all_desc.append(make_double(parent_mol,bond,i))
                    elif(degree1==1 and degree2==1):
                        all_desc.append(make_double(parent_mol,bond,i))
                    #try making it a triple bond    
                    #-C#N or -C#C or -C#C-
                    if (aname1 =="N" and degree1==1 and aname2=="C" and degree2<=2) or (aname1 =="C" and degree1<=2 and aname2=="N" and degree2==1) :# -C-N -> -C#N
                        all_desc.append(make_triple(parent_mol,bond,i))
                    elif aname1 =="C" and aname2=="C" and degree1<=2 and degree2<=2:
                        all_desc.append(make_triple(parent_mol,bond,i))
            elif btype(bond) =="DOUBLE":
                #try making a single bond, should be fine. but maybe diols are unlikely?
                all_desc.append(make_single(parent_mol,bond,i))
                #try making it a triple bond if bond is not in ring
                if not bond.IsInRing():
                    if (aname1=="C" and aname2=="C" and degree1<=3 and degree2<=3):
                        all_desc.append(make_triple(parent_mol,bond,i))
                    if (aname1=="C" or aname2=="C") and (aname1=="N" or aname2=="N"): #if we have a C=N or N=C
                        if (aname1=="C"):
                            if(degree1<=3 and degree2<3):
                                all_desc.append(make_triple(parent_mol,bond,i))
                        elif (degree2<=3 and degree1<3):
                            all_desc.append(make_triple(parent_mol,bond,i))
                    if (aname1=="N" and aname2=="N" and degree1<3 and degree2<3):
                        all_desc.append(make_triple(parent_mol,bond,i))
                    if (aname1=="P" and aname2=="P" and degree1<3 and degree2<3):
                        all_desc.append(make_triple(parent_mol,bond,i))
                        
            elif btype(bond) =="TRIPLE":
                #try making a single bond
                all_desc.append(make_single(parent_mol,bond,i))
                #try making a double bond, it already was a triple, so should be fine
                all_desc.append(make_double(parent_mol,bond,i))
                #try removing the bond
            tmp_mols=remove_bond(parent_mol,aidx1,aidx2)
            tmp_mol=Chem.MolFromSmiles(tmp_mols)
            if check_connected(tmp_mol)==1:
                all_desc.append(tmp_mols)
    # remove bad smiles, duplicates, align smiles
    return tightenSmileList(all_desc)

### JMR returns some atom properties based on your GPCRLigNet.frame_data.py
def get_mol_props(smiles):
    mol=Chem.MolFromSmiles(smiles)
    NHD=Lipinski.NumHDonors(mol)
    NHA=Lipinski.NumHAcceptors(mol)
    MWT=Descriptors.ExactMolWt(mol)
    MLP=Descriptors.MolLogP(mol)
    MMR=Crippen.MolMR(mol)
    NAT=mol.GetNumAtoms()#notice the mol object generated from the smiles has this method
    PSA=QED.properties(mol)[4]#a bit redundant sincewe already calculated some of these
    #QEDproperties(MW=180.15899999999996, ALOGP=1.3101, HBA=4, HBD=1, PSA=63.60000000000001, ROTB=2, AROM=1, ALERTS=2)
    return NHD,NHA,MWT,MLP,MMR,NAT,PSA

### JMR checks if the molecule is real, returns 1 if it is not
def check_if_not_real(smiles):
    #try to add hydrogens and give it 3d coordinates
    try:
        mol=Chem.MolFromSmiles(smiles)
        m2=Chem.AddHs(mol)
        check_val=AllChem.EmbedMolecule(m2)
    except:
        check_val=1
        
    return check_val


### def buildGraph:
##### string seed - the SMILE string of the starting molecule to seed the graph
##### int depth - the number of generations to explore
##### boolean complete_connections - flag to determine wether or not to add the remaing connections to outermost nodes post-loop
##### boolean write_to_log - flag todetermine if log file for grpah build should be written. Default True.
### returns nx.Graph chemical_space_graph - the completed chemical space graph
def buildGraph(seed, depth, complete_connections = False, write_to_log = True):

    log_file_name = "log_"+seed+"_d"+str(depth)+"_ec"+str(complete_connections)+".txt"
    if write_to_log:
        logfile = open(log_file_name, "w")
        logfile.write(tab([[seed, depth, complete_connections]],headers=['Seed','Depth','Complete Connections'])+'\n')

    start_graph_time = time()
    print('GRAPH PARAMETERS')
    print(tab([[seed, depth, complete_connections]],headers=['Seed','Depth','Complete Connections']),'\n')
    
    leafs=[seed]#current 'leafs' (molecules whose neighbors we have not calculated yet)
    chemical_space_graph = nx.Graph()# regular graph does not allow duplicate edges. Digraph would add directionality to edges. Regular graph is fine for now 

    chemical_space_graph.add_node(seed)#I think adding the smiles string as the name of the node may simplify things

    #now build the graph
    for gen in range(1,depth+1):
        #list of leafs for the next generation
        new_leafs=[]
        #list of times it took to process a leaf
        leaf_times=[]
        branch_start = time()
        #loop over all leafs - leafs is the list of molecules whose neighbors have not been generated yet
        for leaf in leafs:
            #get the wall time
            ti=time()
            #comptue neighbors for the current leaf in the loop
            all_neigh = nextGen(leaf)
            ##Possibly speed this up by using 
            for i,mol_smiles in enumerate(all_neigh): # iterates through the generated neighbors of the current leaf
                #check if we have computed found this node yet
                if mol_smiles not in chemical_space_graph.nodes:
                    #add it to the graph
                    chemical_space_graph.add_node(mol_smiles)
                    #since its a new smiles add it to the new leaf list
                    new_leafs.append(mol_smiles)
                #add the edge to the graph
                chemical_space_graph.add_edge(leaf,mol_smiles)
            #get the difference in wall time
            iter_time = time()-ti
            leaf_times.append(iter_time)
        #set the leaf list to the new one
        leafs=new_leafs
        print("actual time for last iter:",reportTime(time() - branch_start))
        print("number of leafs for next iter",len(leafs))
        if write_to_log:
            logfile.write("\nnumber of leafs for next iter: "+str(len(leafs)))
        # estimating and reporting time for next generation of nodes to be added
        expected_time = len(leafs)*np.average(leaf_times)
        print("expected time for next iter:",reportTime(expected_time))
        if write_to_log:
            logfile.write("\nexpected time for next iter: "+reportTime(expected_time))
    
    if complete_connections:
        print("all nodes created; now adding remaing edges (expected time above)")
        if write_to_log:
            logfile.write("\nall nodes created; now adding remaing edges (expected time above)")
        # it will be most efficient to just add the last remaining connections post-loop, once all nodes are created
        ti = time()
        for leaf in leafs:
            all_neigh=nextGen(leaf)
            for neigh in all_neigh:
                if neigh in chemical_space_graph.nodes:
                    chemical_space_graph.add_edge(leaf,neigh)
        tf = time() - ti
        print("adding final edges actually took:",reportTime(tf))
        if write_to_log:
            logfile.write("\nadding final edges actually took: "+reportTime(tf))

    build_time = time() - start_graph_time
    print("total build time for this graph:",reportTime(build_time))

    print("total number of nodes:",chemical_space_graph.number_of_nodes())
    print("total number of edges",chemical_space_graph.number_of_edges())

    if write_to_log:
        logfile.write("\ntotal build time for this graph: "+reportTime(build_time))
        logfile.write("\ntotal number of nodes: "+str(chemical_space_graph.number_of_nodes()))
        logfile.write("\ntotal number of edges: "+str(chemical_space_graph.number_of_edges()))
    logfile.close()

    return chemical_space_graph

### def reportTime
##### float t - total time in seconds
### returns a string reporting the time in hours:minutes:seconds
def reportTime(t):
    hours, remainder = divmod(t, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:.0f}:{:02.0f}:{:02}".format(hours,minutes,round(seconds))

### def get_node_labels_fps()
##### nx.Graph csg - the chemical space graph
### returns list - the SMILE strings of each node [list] and their corresponding fingerprints [list]
def get_node_labels_fps(csg):
    enc = MHFPEncoder(1024)
    #the min hash fingerprint from the mhfp package
    fps=[]
    print('encoding molecules using MHFP')
    node_labels = []#Faerun has something internal that turns the smiles into a molecule image (sweeeet)
    for i, s in enumerate(csg.nodes):
        node_labels.append(s)
        ###instead of the minimal spanning tree of the chemical space graph we will use a more traditional method for determinine node positions
        fps.append(np.array(enc.encode_mol(Chem.MolFromSmiles(s))))
    fps=np.array(fps)
    return node_labels, fps

### def getPlotList()
##### list node_labels - all of the smiles from the graph
### returns prop_list - 8 molecular properties for each of the smiles; a single list for each property, 8 lists total-> [8x[]]
### dependencies: get_mol_props()
def getPropList(node_labels):
    #properties to include
    NHD=[]
    NHA=[]
    MWT=[]
    MLP=[]
    MMR=[]
    NAT=[]
    PSA=[]
    # adding a continuous druglikeness metric
    qed=[]

    # prop_list = [NHD,NHA,MWT,MLP,MMR,NAT,PSA,qed]
    prop_list = [[],[],[],[],[],[],[],[]]
    for node_smiles in node_labels:
        # mol_props is a list of all the molecular properties for a given smile
        mol_props=get_mol_props(node_smiles)
        # for i, prop in enumerate(mol_props):
        #     prop_list[i].append(mol_props[i])
        #     # from rdkit.Chem import QED
        # prop_list[7] = QED.default(Chem.MolFromSmiles(node_smiles))
        NHD.append(mol_props[0])
        NHA.append(mol_props[1])
        MWT.append(mol_props[2])
        MLP.append(mol_props[3])
        MMR.append(mol_props[4])
        NAT.append(mol_props[5])
        PSA.append(mol_props[6])
        qed.append(QED.default(Chem.MolFromSmiles(node_smiles)))
    return NHD, NHA, MWT, MLP, MMR, NAT, PSA, qed


# mol_props - list: a list of properties for a given smile from get_mol_props
# node_labels - list of all the smiles in the graph
# fps - corresponding list of all the fingerprints from the graph
# prop_list - list of all the molecular properties of all the molecules in the graph (corresponds with above). 
def faerunPlot(chemical_space_graph, scatter_name,node_labels,fps):
    # enc = MHFPEncoder(1024)
    # #the min hash fingerprint from the mhfp package
    # fps=[]
    # print('encoding molecules using MHFP')
    # node_labels = []#Faerun has something internal that turns the smiles into a molecule image (sweeeet)
    # for i, s in enumerate(chemical_space_graph.nodes):
    #     node_labels.append(s)
    #     ###instead of the minimal spanning tree of the chemical space graph we will use a more traditional method for determinine node positions
    #     fps.append(np.array(enc.encode_mol(Chem.MolFromSmiles(s))))

    #lf = tm.LSHForest(1024, 64)
    #x_posit, y_posit, start_tmap, termini_tmap, _ = tm.layout_from_lsh_forest(lf)# this is also not supported on my macOS
    # fps=np.array(fps)
    print('fingerprint array shape',fps.shape)
    #use TSNE to visualize the chemical space
    time_init=time()
    print('computing the t-distributed stochastic neighbor embedding')
    X_embedded = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=3).fit_transform(fps)
    x_posit=X_embedded[:,0]
    y_posit=X_embedded[:,1]
    print('computing the t-distributed stochastic neighbor embedding took',(time()-time_init),'seconds')
    ### I tried to directly add the chemical space graph, but my macOS didn't support layout_from_edge_list (smh)
    ### we can still add the links from the chemical space graph
    print('extracting edges from the CSG')
    start_nodes=[]
    end_nodes=[]
    for edges in chemical_space_graph.edges:
        start_nodes.append(node_labels.index(edges[0]))
        end_nodes.append(node_labels.index(edges[1]))
    # ## Turns out we dont need tmap for a minimal spanning tree. tmap's thing is speed, so maybe we want to move to it in the future
    # ## I ran into prioblems specific to macOS that i dont want to debug https://github.com/reymond-group/tmap/issues/12
    # chemical_space_tree=nx.minimum_spanning_tree(chemical_space_graph)
    # ## now compute an xy representation of the tree (https://stackoverflow.com/questions/57512155/how-to-draw-a-tree-more-beautifully-in-networkx)
    # ##this was pretty slow for me
    # pos_dict = nx.nx_agraph.graphviz_layout(chemical_space_tree, prog="twopi")
    # x_posit=[]
    # y_posit=[]
    # for node in node_labels:
    #    x_posit.append(pos_dict[node][0])
    #    y_posit.append(pos_dict[node][1])
    # x_posit=(np.array(x_posit)-np.min(x_posit))/(np.max(x_posit)-np.min(x_posit))
    # y_posit=(np.array(y_posit)-np.min(y_posit))/(np.max(y_posit)-np.min(y_posit))
    print('computing some molecule properties')
    #use a physchem prop for coloring
    color_values=[]

    NHD, NHA, MWT, MLP, MMR, NAT, PSA, qed = getPropList(node_labels)

    # plt.scatter(x_posit,y_posit)
    # plt.show()
    print('generating faerun .html')
    faerun = Faerun(view="front", coords=False)
#     scatter_name="Aspirin_Chemical_Space"
    tree_name="Aspirin_Chemical_Space_Tree"
    faerun.add_scatter(
        scatter_name,
        {"x": x_posit, "y": y_posit, "c":[NHD,
                    NHA,
                    MWT,
                    MLP,
                    MMR,
                    NAT,
                    PSA,
                    qed], "labels": node_labels},
        shader="smoothCircle",
        point_scale=2.0,
        max_point_size=20,

        categorical=[True, True, False, False, False, False,False],
        colormap=["tab10", "tab10", "Blues", "Blues", "Blues", "Blues","Blues","plasma"],
        series_title=['NHD','NHA','MWT','MLP','MMR','NAT','PSA','QED'],
        has_legend=True,
    )

    faerun.add_tree(tree_name, {"from": start_nodes, "to": end_nodes}, point_helper=scatter_name, color="#222222")


    faerun.plot(scatter_name,template="smiles")
