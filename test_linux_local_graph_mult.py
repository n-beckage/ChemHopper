# record start time
import time
t0=time.time()
from ChemTools import *

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

######################################################### BEGIN SCRIPT ########################################################################################

# making directories to organize our output files
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
# parent_0='CC(C)C1=C(C(=CC=C1)C(C)C)O'
# using acetone as an example
parent_name='Ethane'
parent_0='CC'
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
generations=[[(parent_0,dock_it(lig_smile=parent_0,prot_pdbqt=prot_pdbqt,exhaustiveness=exhaust,iiter='0.0'))]]
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
            results.append((smi,dock_it(smi,prot_pdbqt,exhaust,iiter)))
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
