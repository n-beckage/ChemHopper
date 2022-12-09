# record start time
import time
t0=time.time()
from ChemToolsLn import *

### define some chemical data

# Corresponding atomic numbers of atomic_names
allowed_atomic=[6,7,8,9,15,16,17,35]
atomic_names=["C","N","O","F","P","S","Cl","Br"]
halogens=[9,17,35]
type_i=[0,1,2,1.5]
types=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]
######################################################### BEGIN SCRIPT ########################################################################################

# making directories to organize our output files
make_dirs(["logs","np_objs","configs","grids","vina_outs","mols"])

# setting up log format (only needs to be done once)
logging.basicConfig(filename='error_log.log', filemode='w',level=logging.INFO,format='%(asctime)s -%(name)s:%(levelname)s- %(message)s')


# preparing the receptor
prot_pdb='1e7a_aligned.pdb'
# prot_pdbqt=prepare_receptor(prot_pdb)
prot_pdbqt=prot_pdb+'qt'

### Propofol SMILE string - parent molecule, patient 0, generations[0][0], etc
# parent_0='CC(C)C1=C(C(=CC=C1)C(C)C)O'
parent_name='Ethane'
parent_0='CC'
#propofol=Chem.MolFromSmiles(propofol)
#Draw.MolToFile(propofol,"propofol.png")

# defining the exhaustiveness for Autodock Vina, now using keyboard input:
print("\n\n\nWelcome to Chemical Space!")

print("\nPlease enter the name and SMILE string of your starting ligand:")
parent_name=input("Name = ")
parent_0=input("SMILE = ")

print("\nHow deep will you go?")

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

# DOCKING - FOR REAL RUNS
generations=[[(parent_0,dock_it(lig_smile=parent_0,prot_pdbqt=prot_pdbqt,exhaustiveness=exhaust,iiter='0.0'))]]

# NO DOCKING - USE FOR TESTING
# generations=[[(parent_0,ran.uniform(-6,0))]]

# Keeping a list containing all UNIQUE molecules created and tested
master_list = [parent_0]

# creates best_path list (top performers from each search), which will be of length depth+1, due to the parent molecule being the first entry
best_path=[generations[0][0]]

# Begin iterations
parent = parent_0
for gen in range(depth):
    # call the chemical exploration function
    all_desc = nextGen(parent)

    # extract uniq descendants and previously discovered descendants (intersect)
    uniq_desc, intersect_desc = extract_intersection(all_desc, master_list)
    
    # add novel uniq descendants to master list
    master_list = master_list + uniq_desc

    # time to dock each descendant in uniq_desc and record results
    results=[] # the results for current GEN ONLY

    ### CHECKPOINT

    print('\n################\n DEPTH = '+str(gen+1)+'/'+str(depth)+' \n################')

    print("# OF UNIQUE DESC:",str(len(uniq_desc)),"\n")
    for i,smi in enumerate(uniq_desc):
        iiter=str(gen+1)+"."+str(i)
        print("ITER =",str(iiter))
        try:
            # DOCKING - FOR REAL RUNS
            results.append((smi,dock_it(smi,prot_pdbqt,exhaust,iiter)))

            # NO DOCKING - FOR TEST RUNS
            # results.append((smi,ran.uniform(-6,0)))
        except Exception as e:
            logging.exception("IITER: "+str(iiter)+"; TYPE: "+type(e).__name__)
        
    # record this generations results to the master list:
    generations.append(results)

    # pulling just docking scores from the results list
    affinities=[x[1] for x in results]
    # getting the index of the best performer in results
    best_index=affinities.index(min(affinities))
    # adding the best perfomer to the best_path list
    best_path.append(results[best_index])
    # Assigning the top molecule to be the next generation's parent
    parent = results[best_index][0]

    ### produce an image of the molecules
    # first create a mol list
    uniq_mols=[]
    bad_smis=[]
    for smi in uniq_desc:
        mol=Chem.MolFromSmiles(smi)
        if mol is not None:
            Chem.rdmolops.Kekulize(mol,clearAromaticFlags=True)
            uniq_mols.append(mol)
        else:
            bad_smis.append(smi)
            logging.info("BAD SMILE: "+smi)
            print("BAD SMILE: ",smi)

    # returns index of the top molecule in uniq_desc
    best_index=uniq_desc.index(best_path[-1][0])
    top=uniq_mols[best_index]
    save_grids(uniq_mols,top,affinities,gen,8,10)

    # print("Generations list is:")
    # print(generations)

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