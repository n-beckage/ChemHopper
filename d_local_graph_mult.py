from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit.Chem import Draw
import subprocess as sp
import os
import matplotlib.pyplot as plt
import sys
import shutil
import glob
from time import time

### Propofol SMILE string - used as a sample molecule for functionality testing
propofol='CC(C)C1=C(C(=CC=C1)C(C)C)O'
propofol=Chem.MolFromSmiles(propofol)
#Draw.MolToFile(propofol,"propofol.png")

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
## ring rules:
max_n_ring=4
min_ring_size=5
max_ring_size=7
## smii - smile string
def check_cylces(smii):
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
def configure(receptor,ligand,iteration='test',fname="config",size=20,exhaustiveness=8,center_x=12.95,center_y=15.76,center_z=2.28,out='out',cpu=1,num_modes=9,score_only=False):
    config_i=fname+"_"+iteration+'.txt'
    out_name=out+"_"+iteration
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
    return config_i,out_name

### prepares the receptor for docking with vina
## prot_file - the name of the pdb file of the protein/pocket/receptor
prot_pdb='1e7a_aligned.pdb'
def prepare_receptor(prot_file):
    sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_receptor4.py -r "+prot_file+" -v",shell=True)
prot_pdbqt=prot_pdb+'qt'

# propofol in the bound position is the reference ligand
# here I will prepare the refernce like any other ligand and run vina on it to see if I get similar results with propofol created from a smile
# ref_lig="reference_ligand.pdb"
# sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_ligand4.py -l "+ref_lig+" -v",shell=True)
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
def dock_it(smile):
    try:
        iiter=str(i)+"."+str(j)
    except NameError:
        iiter='test'

    
    mol=Chem.MolFromSmiles(smile)
    mh = Chem.AddHs(mol)
    # Generates the 3D coordinates for the molecule; serves as a test for generated molecules, as RDkit will not create coordinates if the mol is not chemically viable
    AllChem.EmbedMolecule(mh,useRandomCoords=True,maxAttempts=100) # it would be a good idea to track how many molecules cannot be embedded
    # I don't think I need to align within rdkit, probably just need to run prepare_ligand4.py
    #Chem.rdMolAlign.AlignMol(mh,ref_lig)
    # Merck Molecular Force Field: https://en.wikipedia.org/wiki/Merck_molecular_force_field
    # uses MMFF to optimize a molecule's structure
        # RETURNS: 0 if the optimization converged, -1 if the forcefield could
        # not be set up, 1 if more iterations are required.
    Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(mh,maxIters=5000)
    # create the ligand file with unique name
    lig_pdb="mol_"+iiter+".pdb"
    # creates a PDB file for the ligand 
    Chem.rdmolfiles.MolToPDBFile(mh,lig_pdb)
    # calls prepare_ligand4.py to prepare lig_file (convert to .pdbqt)
    sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_ligand4.py -l "+lig_pdb+" -v",shell=True)
    # add qt to the end of lig_pdb so that we have ligand.pdbqt
    lig_pdbqt=lig_pdb+"qt"
    configuration,out_name=configure(prot_pdbqt,lig_pdbqt,iiter)
    # runs vina and logs results
    logfile="log_"+iiter+".txt"
    with open(logfile,'w') as log:
        run=sp.run("vina --config="+configuration,shell=True,stdout=log)
    # splitting output
    sp.call("vina_split --input "+out_name+'.pdbqt',shell=True)
    # deleting the original out file from vina plus all but the best modes from vina_split
    os.remove(out_name+'.pdbqt')
    for i in range(8): # range should be num_modes-1
        os.remove(out_name+"_ligand_"+str(i+2)+".pdbqt")
    # renaming the output
    out=out_name+"_ligand_1.pdbqt"
    split=out.rsplit("_1",1)
    os.rename(out,''.join(split))
    # Opening the log file to read in the best affinity
    with open(logfile,"r") as results:
        lines=results.readlines()
        best_mode=lines[-9]
        affinity=float(best_mode.split(' ')[10])
    print(smile,affinity)
    return(affinity)



#### JACOB'S DOCK_IT()
### attempts to dock the molecule
## To do this, there a few unique steps. Create and prepare ligand pdb, prepare the receptor pdb (may only have to do this once if working with the same pocket), run vina (remember to log), and 
## return a docking score. First, we need to create a pdb file from the smile string we are inputting. That rewuires making a mol object, addding Hs, Embedding it (which acts like a legitimacy
## test for the molecules), and then obtimizing the forcefield. Then we prepare this pdb file (convert to pdbqt).
## sidx -
## sijx -
## smile - smile string of the mol to be docked
## lig_file - name of the PDB file to be created; needs to have .pdb ending
## prot_file - name of the receptor pdb file
## exhaustiveness - how exhaustive of a search you want
lig_file="reference_ligand.pdb"
def dock_it(sidx,sijx,smile,lig_file):
    try:
        mol=Chem.MolFromSmiles(smile)
        mh = Chem.AddHs(mol)
        # Generates the 3D coordinates for the molecule; serves as a test for generated molecules, as RDkit will not create coordinates if the mol is not chemically viable
        AllChem.EmbedMolecule(mh,useRandomCoords=True,maxAttempts=100) # it would be a good idea to track how many molecules cannot be embedded
        # Merck Molecular Force Field: https://en.wikipedia.org/wiki/Merck_molecular_force_field
        # uses MMFF to optimize a molecule's structure
            # RETURNS: 0 if the optimization converged, -1 if the forcefield could
            # not be set up, 1 if more iterations are required.
        Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(mh,maxIters=5000)
        # creates a PDB file for the ligand 
        Chem.rdmolfiles.MolToPDBFile(mh,lig_file)
        # calls prepare_ligand4.py to prepare lig_file (convert to .pdbqt)
        sp.call("python2 C:/Users/nbeck/Desktop/Summer_Research_2022/docking/prepare_ligand4.py -l "+lig_file+" -v",shell=True)
        # add qt to the end of lig_file so that we have ligand.pdbqt
        lig_file+="qt"
        # renaming path of lif file for some reason?
        os.rename("/home/jacob/More_Data/Ligand_NN/albumin_2022/"+lig_file.split("/")[-1],lig_file)
        # string is the command that will be called to run vina
        stringg="/home/jacob/software/autodock_vina_1_1_2_linux_x86/bin/vina --receptor "+prot_pdb+" --ligand "+lig_file
        stringg+=" --center_x 70.885 --center_y 14.481 --center_z 80.641 --size_x 20.0 --size_y 20.0  --size_z 20.0 --out "+lig_file.split(".")[0]+"_dockedt.pdbqt"
        stringg+=" --log mylog.log --cpu 16 --num_modes=1 --exhaustiveness "+str(exhaustiveness)
        
        # calling vina
        sp.call(stringg,shell=True,stdout=sp.DEVNULL)

        # calling vina_split
        sp.run("/home/jacob/software/autodock_vina_1_1_2_linux_x86/bin/vina_split --input "+lig_file.split(".")[0]+"_dockedt.pdbqt",shell=True,stdout=sp.DEVNULL)

        # not sure what lig_2 is or why it is here...
        # but it seems that vina is ran again on lig_2 and logged
        lig_2=lig_file.split(".")[0]+"_dockedt_ligand_1.pdbqt"#lig_0_0_dockedt_ligand_1.pdbqt
        stringg="/home/jacob/software/autodock_vina_1_1_2_linux_x86/bin/vina --receptor "+prot_pdb+" --ligand "+lig_2
        stringg+=" --center_x 70.885 --center_y 14.481 --center_z 80.641 --size_x 20.0 --size_y 20.0  --size_z 20.0 --score_only --out "+lig_file.split(".")[0]+"_docked.pdbqt"#--weight_gauss2 -0.01042 # --weight_rot 0.02423
        stringg+=" --log mylog2.log --cpu 1 --num_modes=1 --exhaustiveness 1"
        sp.run(stringg,shell=True,stdout=sp.DEVNULL)

        # looks like he's pulling the affinity out of the log file
        fo=open("mylog2.log",'r')
        tmp=0
        for line in fo:
            stufz=line.strip().split()
            if "Affinity:" in stufz:
                energies=float(stufz[1])
                break
        fo.close()
        os.remove("mylog.log")
        os.remove("mylog2.log")
        os.remove(lig_file)
        os.remove(lig_file[:-2])
        #os.remove(lig_file.split(".")[0]+"_dockedt.pdbqt")
        os.remove(lig_2)
    except:
        energies=0.
    print(sidx,sijx,smile,energies)
    return energies



##### BEGIN JACOB'S SCRIPTING #####

rname="attempt_3_exhaust_1_semigreed_nohol_noring"

rdir="/home/jacob/More_Data/Ligand_NN/albumin_2022/"
#sdf_file=rdir+"1e7a_lig.sdf"

# the initial ligand file
int_lig_file=rdir+"1e7a_lig.pdb"
# the receptor file
prot_pdb=rdir+"1e7a_prot.pdb"
#molecules = Chem.rdmolfiles.SDMolSupplier(sdf_file)
#for mol in molecules:
#  min_smiles=Chem.MolFromSmiles(Chem.MolToSmiles(mol))
#Chem.rdmolops.Kekulize(min_smiles)

# preparing receptor
sp.call("/home/jacob/software/mgltools_x86_64Linux2_1.5.6/bin/python /home/jacob/software/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py  -r "+prot_pdb,shell=True)
prot_pdb+="qt"

# will iterate 100 times
iter_do=range(100)#10
exhaustiveness=1
#dock_it(0,min_smiles,"/home/jacob/More_Data/Ligand_NN/albumin_2022/1e7a_vina/lig_"+str(0)+"_"+str(0)+".pdb")
#exit()

#define some chemical data
allowed_atomic=[6,7,8,9,17]#,9,15,16,17,35]#implicit hydrogen
atomic_names=["C","N","O","F","Cl"]
halogens=[9,17]
type_i=[0,1,2,1.5]
types=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]
n_save=5 #save top n_save for next round

#disallow C=C=C
#apply pains filter
#apply drug likeness filter from noah
#apply smarts pattern removals


cwd=os.getcwd()
for x in iter_do:
    #list of smiles strings for new molecules
    new_mol_full=[]
    print("ON ITERATION "+str(x))
    #loop over nodes from previous iteration
    for smiles_opt in min_smiles:
        new_mol=[smiles_opt]
        
        test_mol=Chem.MolFromSmiles(smiles_opt)

        # Not sure what this line does specifically, but it is necessary so that the function degree() returns the proper bond counts.
        Chem.rdmolops.Kekulize(test_mol,clearAromaticFlags=True)
        
        ### loop over atoms in the molecule to apply atom editing rules
        for i,atom in enumerate(test_mol.GetAtoms()):
            at_deg=degree(atom)
            atom_n=aname(atom)
    
            #### Try to Mutate Atoms
            if at_deg == 4:#could be C or N
                if atom_n=="C":
                    new_mol.append(mutate(test_mol,atom,i,"N",1))#mutates atom in molecule to "N"
                if atom_n=="N":
                    new_mol.append(mutate(test_mol,atom,i,"C",0))   
            if at_deg == 3:#could be a C or N
                if atom_n=="C":
                    new_mol.append(mutate(test_mol,atom,i,"N",0))
                    #new_mol.append(mutate(test_mol,atom,i,"N",1))#could be protonated?
                if atom_n=="N":
                    new_mol.append(mutate(test_mol,atom,i,"C",0))
            if at_deg <=2: #could be a C or N or O
                if atom_n=="C":
                    new_mol.append(mutate(test_mol,atom,i,"N",0))
                    #new_mol.append(mutate(test_mol,atom,i,"N",1))#could be protonated?
                    new_mol.append(mutate(test_mol,atom,i,"O",0))
                if atom_n=="N":
                    new_mol.append(mutate(test_mol,atom,i,"C",0))
                    new_mol.append(mutate(test_mol,atom,i,"O",0))
                if atom_n=="O":
                    new_mol.append(mutate(test_mol,atom,i,"C",0))
                    new_mol.append(mutate(test_mol,atom,i,"N",0))
                    #new_mol.append(mutate(test_mol,atom,i,"N",1))#could be protonated?
                # if the atom is a halogen, it has degree of 1
                if atom_n in ["F","Cl"]:
                    new_mol.append(mutate(test_mol,atom,i,"C",0))
                    new_mol.append(mutate(test_mol,atom,i,"N",0))
                    new_mol.append(mutate(test_mol,atom,i,"O",0))
            if at_deg == 1
                new_mol.append(mutate(test_mol,atom,i,"F",0))
                new_mol.append(mutate(test_mol,atom,i,"Cl",0))
    
            #### Try to add a single atom w/ new atom bonded to the existing atom on molecule
            if at_deg == 3:# could be a C or N, could add a single bond C,N,or O
               #assume not an "O"
               if atom_n =="C":
                   #FIX TO ALLOW N(-C)4, non-heteroatom quaternary nitrogen
                   new_mol.append(add_single_bond(test_mol,atom,i,"N"))#adds single bonded "N" to atom in molecule
                   new_mol.append(add_single_bond(test_mol,atom,i,"C"))
                   new_mol.append(add_single_bond(test_mol,atom,i,"O"))
                   #new_mol.append(add_single_bond(test_mol,atom,i,"F"))
                   #new_mol.append(add_single_bond(test_mol,atom,i,"Cl"))
            if at_deg ==2:# could be a C or N or O, but dont add to O
               if atom_n == "C":
                    #add single bonded atom
                    new_mol.append(add_single_bond(test_mol,atom,i,"N"))
                    new_mol.append(add_single_bond(test_mol,atom,i,"C"))
                    new_mol.append(add_single_bond(test_mol,atom,i,"O"))
                    #new_mol.append(add_single_bond(test_mol,atom,i,"F"))
                    #new_mol.append(add_single_bond(test_mol,atom,i,"Cl"))
                    #add double bonded atom but only if no double bond on c already
                    if all([btype(bond)!="DOUBLE" for bond in atom.GetBonds()]):
                        new_mol.append(add_double_bond(test_mol,atom,i,"N"))
                        new_mol.append(add_double_bond(test_mol,atom,i,"C"))
                        new_mol.append(add_double_bond(test_mol,atom,i,"O"))
               elif atom_n=="N":
                    new_mol.append(add_single_bond(test_mol,atom,i,"N"))
                    new_mol.append(add_single_bond(test_mol,atom,i,"C"))
            if at_deg ==1:# could be a C or N or O, but dont add double bond to O
               if atom_n == "C" :
                   #add single bonded atom
                   new_mol.append(add_single_bond(test_mol,atom,i,"N"))
                   new_mol.append(add_single_bond(test_mol,atom,i,"C"))
                   new_mol.append(add_single_bond(test_mol,atom,i,"O"))
                   #new_mol.append(add_single_bond(test_mol,atom,i,"F"))
                   #new_mol.append(add_single_bond(test_mol,atom,i,"Cl"))
                   #add double bonded atom
                   new_mol.append(add_double_bond(test_mol,atom,i,"N"))
                   new_mol.append(add_double_bond(test_mol,atom,i,"C"))
                   new_mol.append(add_double_bond(test_mol,atom,i,"O"))
               elif atom_n =="N":
                   new_mol.append(add_single_bond(test_mol,atom,i,"N"))
                   new_mol.append(add_single_bond(test_mol,atom,i,"C"))
                   new_mol.append(add_double_bond(test_mol,atom,i,"C"))
               elif atom_n=="O":
                   #add single bonded atom
                   new_mol.append(add_single_bond(test_mol,atom,i,"C"))
                   new_mol.append(add_single_bond(test_mol,atom,i,"O"))
    
            #### Try removing atom
            if at_deg==1:#no worries
               new_mol.append(remove_atom(test_mol,atom,i))
            elif at_deg>=2:#only can remove if the resulting molecule is connected (doesn't have floating framents)
               tmp_mols=remove_atom(test_mol,atom,i)
               tmp_mol=Chem.MolFromSmiles(tmp_mols)
               if check_connected(tmp_mol)==1:
                   new_mol.append(tmp_mols)
    
        #### TRY bond manipulations
        idone=0
        for i,bond in enumerate(test_mol.GetBonds()):
            #ignore bonds to halogens
            if test_mol.GetAtoms()[bond.GetBeginAtomIdx()].GetAtomicNum() not in halogens and test_mol.GetAtoms()[bond.GetEndAtomIdx()].GetAtomicNum() not in halogens:
                #determine elements in bond
                atom1=atomic_names[allowed_atomic.index(test_mol.GetAtoms()[bond.GetBeginAtomIdx()].GetAtomicNum())]
                atom2=atomic_names[allowed_atomic.index(test_mol.GetAtoms()[bond.GetEndAtomIdx()].GetAtomicNum())]
                #determine degree of atoms in bond
                degree1=degree(test_mol.GetAtoms()[bond.GetBeginAtomIdx()])
                degree2=degree(test_mol.GetAtoms()[bond.GetEndAtomIdx()])
                #check bond type and apply chemical rules
                if btype(bond) =="SINGLE":
                    #check if not (between two double bonds and in ring)
                    if not is_bet_dub_and_ring(test_mol,bond):
                        #print(idone,bond)
                        #try making it a double bond
                        #print(i+1,is_aromatic(test_mol,bond))
                        if (atom1 == "C") and (atom2 == "C"):#(atom1 == "C" or atom1 == "N") and (atom2 == "C" or atom2 == "N"):
                            if degree1<=3 and degree2<=3:#only two <terterary carbons can form a double bond (ignoring N+ )
                                new_mol.append(make_double(test_mol,bond,i))
                                idone+=1
                        elif (atom1 == "O") and (atom2 == "C"):# a single bond with O-C
                            if degree1==1 and degree2<=3:#check if O is single degree and C is <terterary
                                new_mol.append(make_double(test_mol,bond,i))# hydroxyl -> carbonyl
                        elif (atom1 == "C") and (atom2 == "O"):# a single bond with C-O
                            if degree1<=3 and degree2==1:#check if O is single degree and C is <terterary
                                new_mol.append(make_double(test_mol,bond,i))# hydroxyl -> carbonyl      
                        elif atom1 == "C" and atom2 == "N":# a single bond with C-N
                            if degree1<=3 and degree2<=2:#check if N is 2 degree and C is <terterary
                                new_mol.append(make_double(test_mol,bond,i))# 
                        elif atom1 == "N" and atom2 == "C":# a single bond with N-C
                            if degree1<=2 and degree2<=3:#check if N is 2 degree and C is <terterary
                                new_mol.append(make_double(test_mol,bond,i))# 
                        #try making it a triple bond
                        #-C#N or -C#C or -C#C-
                        if (atom1 =="N" and degree1==1 and atom2=="C" and degree2<=2) or (atom1 =="C" and degree1<=2 and atom2=="N" and degree2==1) :# -C-N -> -C#N
                            new_mol.append(make_triple(test_mol,bond,i))
                        elif atom1 =="C" and atom2=="C" and degree1<=2 and degree2<=2:
                            new_mol.append(make_triple(test_mol,bond,i))
                elif btype(bond) =="DOUBLE":
                    #try making a single bond, should be fine. but maybe diols are unlikely?
                    new_mol.append(make_single(test_mol,bond,i))
                    #try making it a triple bond if not in ring
                    if (atom1 == "C" or atom1 == "N") and (atom2 == "C" or atom2 == "N") and not (atom1 =="N" and atom2 == "N") and not bond.IsInRing():
                      if degree1<=2 and degree2<=2:
                          new_mol.append(make_triple(test_mol,bond,i))
                elif btype(bond) =="TRIPLE":
                    #try making a single bond
                    new_mol.append(make_single(test_mol,bond,i))
                    #try making a double bond, it already was a triple, so should be fine
                    new_mol.append(make_double(test_mol,bond,i))
                #try removing the bond
                tmp_mols=remove_bond(test_mol,bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())
                tmp_mol=Chem.MolFromSmiles(tmp_mols)
                if check_connected(tmp_mol)==1:
                   new_mol.append(tmp_mols)
    
        #### Try adding new bonds (rings) in the molecule
        # creating cycle rules, find pairs of atoms that a bond could be added to make a cycle. Cycle mainpulation? (6<->5)
        bond_inds=[set([bnd.GetBeginAtomIdx(),bnd.GetEndAtomIdx()]) for bnd in test_mol.GetBonds()]
        for aii,atomi in enumerate(test_mol.GetAtoms()):
           for ajj,atomj in enumerate(test_mol.GetAtoms()):
               if ajj > aii:
                   if set([aii,ajj]) not in bond_inds:
                       #determine elements in new bond
                       atom1=atomic_names[allowed_atomic.index(atomi.GetAtomicNum())]
                       atom2=atomic_names[allowed_atomic.index(atomj.GetAtomicNum())]
                       #determine degree of atoms in new bond
                       degree1=degree(atomi)
                       degree2=degree(atomj)
                       #try creating a cycle and checking cycle info is ok
                       if (atom1 == "C") and (atom2 == "C"):
                           if degree1==3 and degree2==3:
                               #can add a single
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                           elif degree1<=2 and degree2<=2:
                               #can add single or double
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                               tst_s=add_double(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                       elif (atom1 == "C") and (atom2 == "N"):
                           if degree1==3 and degree2==1:
                               #can add a single
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                           if degree1<=2 and degree2==1:
                               #can add a single or double
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                               tst_s=add_double(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                       elif (atom1 == "N") and (atom2 == "C"):
                           if degree1==1 and degree2==3:
                               #can add a single
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                           if degree1==1 and degree2<=2:
                               #can add a single or double
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                               tst_s=add_double(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                       elif (atom1 == "C") and (atom2 == "O"):
                           if degree1<=3 and degree2==1:
                               #can add a single
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
                       elif (atom1 == "O") and (atom2 == "C"):
                           if degree1==1 and degree2<=3:
                               #can add a single
                               tst_s=add_single(test_mol,aii,ajj)
                               if check_cylces(tst_s):
                                   new_mol.append(tst_s)
        ### try expanding/contracting cycles by removing/added atoms
        #save list of list of smiles strings
        new_mol_full.append(new_mol)
    #remove repetitive smiles which can arise from symmetry
    #print(len(new_mol))
    trim_mol=[]
    trim_mol_flat=[]
    new_mol_full_flat=[smix for smi_list in new_mol_full for smix in smi_list]
    for smi_list in new_mol_full:
        trim_list=[]
        for smi in smi_list:
            if trim_mol_flat.count(smi)==0:
                trim_list.append(smi)
                trim_mol_flat.append(smi)
        trim_mol.append(trim_list)
    print(len(trim_mol_flat),[len(xxx) for xxx in trim_mol])
    ## for visualizing a test
    #mols=[]
    #for x in trim_mol:
    #    print(x)
    #    x_mol=Chem.MolFromSmiles(x)
    #    Chem.rdmolops.Kekulize(x_mol,clearAromaticFlags=True)
    #    mols.append(x_mol)
    ##print(len(mols))
    #img=Draw.MolsToGridImage(mols,molsPerRow=int(len(mols)/10),subImgSize=(350,350))#,legends=[x.GetProp("_Name") for x in new_mol])   
    #img.save('mol_test.png')
    flat_energies=[]
    flat_smilss=[]
    energies=[]
    smilss=[]
    for ij,ij_list in enumerate(trim_mol):
        energy_list=[]
        smilss_list=[]
        for it,moel in enumerate(ij_list):
            #ti=time()
            try:
                smilss_list.append(moel)
                energy_list.append(dock_it(it,ij,moel,"/home/jacob/More_Data/Ligand_NN/albumin_2022/"+rname+"/lig_"+str(x)+"_"+str(ij)+"_"+str(it)+".pdb"))
            except:
                smilss_list.append("bad")
                energy_list.append(0.0)
        #exit()
        smilss.append(smilss_list)
        energies.append(energy_list)
        flat_energies.extend(energy_list)
        flat_smilss.extend(smilss_list)
    np.save(rname+"/smilss_"+str(x),smilss)
    np.save(rname+"/energies_"+str(x),energies)
    
    for ij,ij_list in enumerate(trim_mol):
        mols=[]
        for it,moel in enumerate(ij_list):
            x_mol=Chem.MolFromSmiles(moel)
            Chem.rdmolops.Kekulize(x_mol,clearAromaticFlags=True)
            mols.append(x_mol)
        img=Draw.MolsToGridImage(mols,molsPerRow=int(len(ij_list)/10),subImgSize=(350,350),legends=[str(it)+" "+str(xx) for it,xx in enumerate(energies[ij])])#,legends=[x.GetProp("_Name") for x in new_mol])   
        img.save(rname+'/atom_mute_'+str(x)+'_'+str(ij)+'.png')
    #get indexes of top n_save scores
    top_flat_inds=np.argsort(flat_energies)[:n_save]
    top_smiles=np.array(flat_smilss)[top_flat_inds]


    #if x==0:
    #    #use the first as old
    #    old_energy=[energies[0][0]]
    #else:
    #    old_energy=[min(min_energy) for xxxx in min_energy]
    ##min_smiles is a list of all the smiles from the previous iteration that beat the previous best energy
    min_smiles=[]
    min_energy=[]
    graph_data=[]#graph_data[i]=[it,ij] of the molecule that made it to the next generation
    #loop over old energies
    #print("check",len(old_energy),len(trim_mol))
    #exit()

    for ixx in range(len(trim_mol)):
        for eni in range(len(energies[ixx])):
            #if energies[ixx][eni]<old_energy[ixx]:
            if smilss[ixx][eni] in top_smiles:
                min_smiles.append(smilss[ixx][eni])
                min_energy.append(energies[ixx][eni])
                graph_data.append([ixx,eni])
                print(str(ixx)+" "+smilss[ixx][eni]+" "+str(energies[ixx][eni])+" ginload lig_"+str(x)+"_"+str(ixx)+"_"+str(eni)+"_dockedt.pdbqt")
                
    print("restarting with "+str(len(min_smiles))+" new molecules ",str(min(min_energy)))
    np.save(rname+"/min_smiles_"+str(x),min_smiles)
    np.save(rname+"/min_energy_"+str(x),min_energy)
    np.save(rname+"/graph_data_"+str(x),graph_data)
    #fo=open("best_smile_"+str(x)+".txt",'w')
    #fo.write(min_smiles+"\n")
    #fo.write(str(min_energy))
    #fo.close()
    #print(energies.index(min_energy),min_energy)

## exhaustiveness of 10
#i indsmiles  docking_score
#-1 0 CC(C)C1=C(O)C(=CC=C1)C(C)C -7.04906
#0 78 CC1(C)C2=CC=CC(=C2O)C1(C)C -7.68937
#1 39 CC1=C2C=CC=C1C(C)(C)C2(C)C -8.12077
#2 26 CC1=C2C=C(F)C=C1C(C)(C)C2(C)C -8.59493
#3 29 CC1=C2C=C(F)C(C)=C1C(C)(C)C2(C)C -8.78752
#4 19 CC1=C2C(C)=C(C(F)=C1F)C(C)(C)C2(C)C -9.09168
#5 69 CC1=C2C(F)C(F)C(C)=C1C(C)(C)C2(C)C -9.10516
#6 30 CC1=C2C(F)C(F)(Cl)C(C)=C1C(C)(C)C2(C)C -9.29423
#7 64 CCC1(C)C2=C(C)C(=C(C)C(F)(Cl)C2F)C1(C)C -9.80458
#8 20 C=C(C)C1(C)C2=C(C)C(=C(C)C(F)(Cl)C2F)C1(C)C -9.73556
#9 6  CC(=CF)C1(C)C2=C(C)C(=C(C)C(F)(Cl)C2F)C1(C)C -9.95064

#mkdir attempt_0_exhaust_10_greedy
#mv atom_mute_* attempt_0_exhaust_10_greedy/.
#mv energies_* attempt_0_exhaust_10_greedy/.
#mv smilss_* attempt_0_exhaust_10_greedy/.
#mkdir attempt_0_exhaust_10_greedy/1e7a_vina
#mv 1e7a_vina/* attempt_0_exhaust_10_greedy/1e7a_vina/.
#mv best_smile_* attempt_0_exhaust_10_greedy/.


#"attempt_1_exhaust_1_semigreed"
#0 0 CC(C)C1=C(O)C(=CC=C1)C(C)C -7.04683

#0 CC1=C(C(C)C)C(O)=C(C(C)C)C=C1 -7.37789 lig_0_0_28.pdb 
#0 CC(C)C1=C(O)C(C(C)C)=C(F)C=C1 -7.45173 lig_0_0_30.pdb 
#0 CC1=CC=C(C(C)C)C(O)=C1C(C)C -7.39873 lig_0_0_42.pdb 
#0 CC(C)C1=CC=CC2=C1OCC2C -7.48011 lig_0_0_75.pdb 
#0 CC(C)C1=C2OCC(C)C2=CC=C1 -7.67585 lig_0_0_80.pdb 

#1 CC(C)C1=C2OCC(C)C2=C(F)C=C1 -7.90311 lig_1_1_72.pdb 
#2 CC1=CCC(C(C)C)C(O)=C1C(C)C -8.17568 lig_1_2_69.pdb 
#3 CC(C)C1=CC=C(F)C2=C1OCC2C -8.03137 lig_1_3_34.pdb 
#3 C=C1OC2=C(C=CC=C2C(C)C)C1C -7.86888 lig_1_3_52.pdb 
#4 CC(C)C1=C2OC(C)C(C)C2=CC=C1 -7.95822 lig_1_4_26.pdb 

#1 CC1=CCC2C(O)=C1C(C)CCC2C -8.61477 lig_2_1_96.pdb 
#3 CC(C)C1=CC=CC2=C1OC(=CF)C2C -8.26832 lig_2_3_6.pdb 
#3 C=C=C1OC2=C(C=CC=C2C(C)C)C1C -8.45065 lig_2_3_9.pdb 
#3 C=C1OC2=C(C=CC(F)=C2C(C)C)C1C -8.28066 lig_2_3_36.pdb 
#4 CC(C)C1=C2OC(C)(C)C(C)C2=CC=C1 -8.35973 lig_2_4_26.pdb

#0 CC1=CCC2C(C)=C1C(C)CCC2C -8.84246 lig_3_0_39.pdb 
#2 C=C=C1OC2=C(C=CC(C)=C2C(C)C)C1C -8.99927 lig_3_2_36.pdb 
#2 C=C=C1OC2=C(C=CC=C2C(C)(C)C)C1C -9.11159 lig_3_2_44.pdb 
#2 C=C=C1OC2=C(C=CC=C2C(C)(C)F)C1C -8.82773 lig_3_2_46.pdb 
#2 C=C=C1OC2=C(C=CC=C2C(C)(C)Cl)C1C -8.90615 lig_3_2_47.pdb

#0 CC1=CCC2C(C)=C1C(C)CCC2(C)C -9.44863 lig_4_0_93.pdb 
#1 C=C=C1OC2=C(C=CC(C)=C2C(C)(C)C)C1C -9.62981 lig_4_1_50.pdb 
#1 C=C=C1OC2=C(C=CC(C)=C2C(C)(C)Cl)C1C -9.46982 lig_4_1_53.pdb 
#3 C=C=C1OC2=C(C=CC(F)=C2C(C)(C)F)C1C -9.48256 lig_4_3_38.pdb 
#4 C=C=C1OC2=C(C=CC(F)=C2C(C)(C)Cl)C1C -9.34176 lig_4_4_38.pdb 

#0 CC1=C(F)CC2C(C)=C1C(C)CCC2(C)C -9.57683 lig_5_0_17.pdb 
#0 CC1=CCC2C(C)=C1C(C)(F)CCC2(C)C -9.77622 lig_5_0_56.pdb 
#1 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C)C1C -9.93398 lig_5_1_32.pdb 
#2 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)Cl)C1C -9.71741 lig_5_2_32.pdb 
#3 C=C=C1OC2=C(C=CC(F)=C2C(C)(C)C)C1C -9.5647 lig_5_3_53.pdb 

#0 CC1=C(F)C(=O)C2C(C)=C1C(C)CCC2(C)C -9.96768 lig_6_0_27.pdb 
#0 CC1=C(F)CC2C(C)=C1C(C)(F)CCC2(C)C -10.03102 lig_6_0_54.pdb 
#2 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C)C1C -9.92719 lig_6_2_0.pdb 
#2 C=C=C1NC2=C(C=C(F)C(C)=C2C(C)(C)C)C1C -9.96595 lig_6_2_15.pdb 
#2 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)CC)C1C -9.92132 lig_6_2_50.pdb 

#0 C=C1CC(C)C2=C(C)C(C(=O)C(F)=C2C)C1(C)C -10.25806 lig_7_0_79.pdb 
#3 C=C=C1NC2=C(C(F)=C(F)C(C)=C2C(C)(C)C)C1C -10.00501 lig_7_3_23.pdb 
#3 C=C=C1NC2=C(C=C(F)C(C)=C2C(C)(C)C)C1=C -10.22954 lig_7_3_86.pdb 
#4 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)F)C1C -10.10238 lig_7_4_61.pdb 
#4 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)=N)C1C -10.06236 lig_7_4_63.pdb 

#0 C=C1CC(C)C2=C(C)C(C(=O)C(F)=C2C)C1(C)C -10.24776 lig_8_0_0.pdb 
#2 C=C=C1NC2=C(C=C(F)C(C)=C2C(C)(C)C)C1=CF -10.21696 lig_8_2_65.pdb 
#3 C=C=C1CC2=C(C=C(F)C(C)=C2C(C)(C)C(C)F)C1C -10.25346 lig_8_3_14.pdb 
#3 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)(F)F)C1C -10.66535 lig_8_3_62.pdb 
#3 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)(F)Cl)C1C -10.40524 lig_8_3_63.pdb

#0 C=C1CC(C)C2=C(C)C(C(=O)C(F)=C2C)C1(C)C -10.24776 lig_8_0_0.pdb 
#2 C=C=C1NC2=C(C=C(F)C(C)=C2C(C)(C)C)C1=CF -10.21696 lig_8_2_65.pdb 
#3 C=C=C1CC2=C(C=C(F)C(C)=C2C(C)(C)C(C)F)C1C -10.25346 lig_8_3_14.pdb 
#3 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)(F)F)C1C -10.66535 lig_8_3_62.pdb 
#3 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)(F)Cl)C1C -10.40524 lig_8_3_63.pdb 

#1 C=C=C1NC2=C(C1=CF)C(F)=C(F)C(C)=C2C(C)(C)C -10.51968 lig_9_1_27.pdb 
#3 C=C=C1O[N+]2=C(C=C(F)C(C)=C2C(C)(C)C(C)(F)F)C1C -10.41162 lig_9_3_16.pdb 
#3 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(F)(F)CCl)C1C -10.53416 lig_9_3_64.pdb 
#3 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(F)(F)C=C)C1C -10.59971 lig_9_3_66.pdb 
#4 C=C=C1OC2=C(C=C(F)C(C)=C2C(C)(C)C(C)(F)Cl)C1C -10.51531 lig_9_4_0.pdb 