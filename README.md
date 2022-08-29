# ChemHopper Project README

### 8/29/22

Oftentimes the best molecule chosen from generation 1 is not consistent across runs with the following exceptions:

- run 7 and 11 have the same best gen 1

- 13 and 5 as well



#### Questions for Jacob

- Is stereoisomers accounted for when generating new molecules? If so are they considered "unique" mols?
- Where do I find other ligands for serum albumin that target the same pocket as propofol?
- Congruency bug - does it make sense to go backward in our search of chemical space? 
  - Allowing our search to move backwards seems like a good way to get trapped in local minima, doesn't it?
- What Vina parameters did you/can you use in order to get more reliable docking scores?
  - I know Vina has a monte carlo algorithm, so there is a degree of stochasticity (IOW the algorithm is not really consistent; the same docking parameters will give different scores over and over again)
- Is there a reason why in your script you commented out halogen transformations? If so why?
  - In my runs adding fluorines or chlorines seems to be a pretty favored choice for the algorithm and I'm not sure how realistic that is

### 8/15/22

#### Congruency Bug - bring this up with Jacob

**Description**: the molecular editing gauntlet allows for adding an atom, mutating an atom/bond, and removing an atom/bond. These later two transformations represent potential backwards moves in chemical space, potentially leading to congruent molecules across multiple generations (not just immediately neighboring generations either). This is a highly inefficient bug as it wastes resources docking molecules that have already been tested. So the question remains, how do we ensure that no duplicate molecules inadvertently recreated through the gauntlet get docked more than once?

**Solution**: The answer to this question lies within sets. I created a `master_set` that initially contains just the smile string of the parent molecule, `parent_0`. Each generation, we go through the entire gauntlet as per usual. However, we define the list `uniq_desc` as the *difference* between the master set and the current generation of descendants. In other words, `uniq_desc` will only contain molecules that are not in the master set, or *have not been created before*. The master set is then joined with the uniq list of novel descendants `uniq_desc` via the `.union()` method, yielding an updated master set containing all uniq descendants to date.

**Note: I haven only tested this patch in jupyter notebooks, not yet with an entire program run**

### 8/14/22

- added error logging to script when running dock_it iteratively (ln 586); produces `error_log.txt`
- Using regex, fixed bug within `dock_it()` that had prevented the best docking score from the autodock vina output (log files) from being read in correctly (ln 329)
- Ran entire program with `depth=5` and propofol as `parent_0`. This run was called `run_5_depth_5`

#### Bugs to fix

- in the grid images produced each loop. sometimes more than one molecule is highlighted. This is because if the best molecule is the result of a mutation, not an addition, then at least one other molecule created via addition will share the substructure of the best mutated molecule (for example, a molecule might gain a bond to chlorine (addition), whereas another may mutate a carbon to a chlorine. These two molecules would share a common substructure)

- congruent molecules are made across generations via mutation. Say the best mol from the previous generation is butane, and the best mol of the current generation is 1-chloropropane (mutating a terminal carbon to a chlorine). Given that the first parent was propane (butane being a descendant via adding a carbon to the end of propane), 1-chloropropane was already created and docked in the first generation. Now in the second gen, it is created and tested again via mutation of butane. This is a big waste of computational resources and should be avoided.

  - **This has been fixed by implementing sets** 
  - Also, it speaks to the large influence stochasticity is having within autodock vina right now. What ends up being the best  1st-gen descendant of propofol is inconsistent so far. Look at the gen 1 grids for different runs and compare highlighted performers to see for your self 

-  **Done** - organizing folders to contain log files, config files, pdbs and pdbqts, etc

  

### 8/12/22

**Tasks**

- work on intramolecular code chunks
- write basic parameters like `depth`, and important output like`generations` and `results` to a text file at the end of the script - **done**, called `RUN_LOG.txt`

### 8/11/22

Ran the script completely and fully in the folder `firs_full_test`, ensured everything worked on my machine. Took a little over 30 minutes to run starting with propofol with a depth of 1. Saved all results as expected. **Next step should be to get it running on linux workstation, attempt validation**

### 8/10/22

#### Meeting with Jianing

- Validation for the algorithm: starting with one ligand for an Albumin pocket that is similar to another ligand that fits the same pocket; if the algorithm is working right, it should traverse through chemical space and converge on the similar ligand. 
- Report - write a report for what you've done this summer, what results I have so far, and where to go next
- **Next steps:** A big problem with the algorithm as-is is that it builds molecules in 2D space and then converts them to 3D (with rdkit). Autodock Vina utilizes a Monte Carlo algorithm, thus there is a component of randomness to it. That means that the binding score is not the absolute best - it is not an exhaustive docking algorithm. The problem with this is that we know the crystal structure of some bound ligands - say propofol for example. Ideally we would like to convert just one atom/bond at a time and perserve the conformation of the molecule in the pocket, and see how that one change affects delta G (free energy). However, currently rdkit doesn't have the ability to mutate molecules in 3D space, and thus we have to redock each time we create a new molecule. Apparently **Schrodinger has a Python interface that we can use to make script 3D changes**. This should be the next step focus of the project
- For now, adjust Autodock Vina parameters to make improve the certainty of the dokcing positions. This may mean increasing the number of modes (attempts) the algorithm tries, increase exhaustiveness, reduce the grid size, etc. **Might be a good idea to consult Jacob about this.**

### 8/7/22

Manually tested the molecular editing code chunks in jupyter notebook - works insofar as it is tested (with acetone as starting molecule).

**To Continue**: work on defining `iter` within the script and saving the grid drawings in an orderly manner

### 8/6/22

Another note on the **naming system** - each file, no matter how many generations far down, should only have two numbers to identify it. The first number denotes what search the molecules comes from (2) would mean it is part of the second search, aka it descends from the parent of the second search, which naturally is the best-performing molecule from the first search. The second number is more straightforward. It is simply an index of the descendant in the order it was created.

Created a new script that will be my organized, official version of Jacobs. Created a new jupyter notebook to visually test the molecular editing algorithm. Continue working on this next. 

### 8/5/22

Worked through molecular editing code with notebook. Developed tentative naming scheme for files: `fname+str(i)+"."+str(j)+fextension` where `i` is the generation, and `j` is the index of the descendant. `i` also corresponds to the best descendant of the previous generation. For example `1.1` means the second descendant of the parent molecule (1 denotes the generation; indeed the descendants of the parent (0) would be the 1st gen). `5.5` would indicate the 6th descendant of the 5th generation.

Next steps are to read through the pseudocode written for the algorithm and start writing it out in code, continuing to psudeocode and problem-solve as necessary. 

#### Attempt to run one iteration of docking with the script

This job begins at line 221, where the receptor is prepared ONCE. After prepare_receptor is called, it does not need to be called again. Use `prot_pdbqt` from this point on as the object containing the receptor filename.

At this point, it's becoming clear that I need to develop a naming system for each created ligand and their respective output, log, and confi files.

### 7/28/22

Continued work on the `dock_it()` function and pipeline.

#### To be done once in the script:

- prepare the receptor with `prepare_receptor4.py`, which can be done with the function `prepare_receptor()`
  - This is already done; the receptor file that we will be using for vina from here on it is called `1e7a_aligned.pdbqt`

#### To be done iteratively in `dock_it()`:

- create a mol from the given smile string, convert that mol into a systematically named `.pdb` file and prepare it for docking with `prepare_ligand4.py`
- create a systematically named configuration file using `configure()`; the only parameters that should be changing for each iteration is the name of the ligand and name if the output (and maybe the name of the configuration file itself if you want to track configurations for each run of vina).
  - might be worth it to create a special function that only changes the names of the ligand  and out files in `config.txt` instead of completing rewriting the entire configuration file each time. This would be more efficient, especially if no other parameters in `config.txt	` change.
- create a systematically named log file to track the written output from vina
- run `vina_split` on the vina results and discard all output except the best mode (conformation) from each vina run
- extract the binding affinity from the best mode (in the log file) and have the function print the smile string and binding affinity; return binding affinity
- there will probably be other things worth tracking for each run -- stay tuned.

##### Questions for jacob

- why do many of the mol-editing functions have parameters for both the atom object and the atom index? only the atom index is ever used; it seems redundant
- Does it make sense to ignore protonation? i think so, as amine protonation is dependent on pH, so not really much of a solid molecular change
