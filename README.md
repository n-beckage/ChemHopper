# ChemHopper Project Shop Log

### 10/11/22

Successfully separated functions into separate script, `ChemTools.py`, from `test_linux`.

- Noticed that chemical transformations are not as robust as I thought - see the grid produced by  the test run with ethane. When carbon is mutated to oxygen or nitrogen, bond mutations should be possible but are not explored. Neither is the adding a triple bont to nitrogen, which should be possible.

### 10/10/22 - Back into the garage

**Chemical Neighborhood Exploration**

1. Enumerate all *k*-nearest chemical neighbors, where *k* = `depth`
   1. No duplicates
   2. Graph structure (IOW edges) should be preserved
2. Complete molecular editing gauntlet
   1. Add rules for intermolecular bonds (rings)

For **Ligand Optimization**, will need to limit the number of halogenations

- i.e if current parent has 1 or 2 halogenations, stop adding more
- *Foresight Problem*: what if the best ligand is not found by constantly going with the best child of a given generation? imagine that the nth descendant of a so-so child is better than the nth descendant of the "best child" lineage? IOW, what if the path to the best ligand is not smooth like a well? what if it is more like a volcano?

### *BUG -10/9/22*

**Problem**: the program is not timing itself accurately - the displayed runtime is significantly less than what it actually was.

This came up with the runs 21/22, where each program took at least 48 hours, however the the runtime in the `RUN_LOG` file indicated a total runtime less than 24 hours. This is clearly not an accurate runtime. 

**1st Hypothesis**: I suspect the timer in the python script is only recording the time that the python code itself takes to run; IOW it is not keeping time when Vina is running. Since Vina is the lengthiest component of the program, especially at high exhaustiveness (runs 21/22 were at e512), this would explain why the recorded runtime is so much less than the actual runtime.

**1st Solution**: If the above hypothesis is correct, then I need to find a way to time the program outside of python, maybe through some external timer that can be run independently on the Linux  terminal.

**2nd Hypothesis**: With the latest runs (21/22), exhaustiveness was too high and maxed out the cpus on the local linux machine. Each was held at constant 100% usage. The computer didn't crash or freeze, however it really couldn't do much else but run the programs; it was very laggy when I tried to browse the web or do anything. I suppose it's possible that the python timer itself was lagging, especially when Vina was running (which would pretty much, between the two programs, be all the time). This would explain why the recorded runtime was less than the actual runtime.

**2nd Solution**: Try a run with high exhaustiveness run one at time and see if it is accurate. 1st solution could also confirm/reject this hypothesis.

**Immediate Goals:**

- **Explore the local chemical space around a ligand**
- **Start with a ligand that has a chemical series associated with it**
- **let's go n steps away from a ligand and characteize that chemical space with the tools we have**
- **Preserve the graph structure** - make graph with tmap
- Abstract for this!! submit for the ACS meeting Spring 2023

### 9/27/22

Meeting With Jacob and Jianing

- Start with a crystal structure and build in 3D, ideally within the same pocket.  That way the binding conformation would be conserved
- Maybe shrodinger has reaction data types to store each transformation as a reaction
- Matched molecular pair analysis 
  - MMPDB - open sourced 
  - MPA is often used for QSAR - seeing how an incremental change affects
- Keep track of changes through a dictionary of possible changes.
- Find cases of proteins where and pockets with multiple ligands for the same pocket
- Out method is unique in it's sampling 
- Our algorithm is only as good as whatever scoring function you use
- At the end of the day we need experimental data to validate any novel scaffolds
- Steps to take to grow confidence in the method:
  - Design a function to assess how close to a known conformation a given pose
- How many molecules are 


### 8/29/22

Oftentimes the best molecule chosen from generation 1 is not consistent across runs with the following exceptions:

- run 7 and 11 have the same best gen 1

- 13 and 5 as well



#### Questions for Jacob

- Is stereoisomers accounted for when generating new molecules? If so are they considered "unique" mols?

  - RDkit has functions to enumerate through stereocenters  - conisder enumerating through the

  - RA score - how synthesizable is a molcule

- Where do I find other ligands for serum albumin that target the same pocket as propofol?

- Congruency bug - does it make sense to go backward in our search of chemical space? (removing atoms)

- What Vina parameters did you/can you use in order to get more reliable docking scores?
  - I know Vina has a monte carlo algorithm, so there is a degree of stochasticity (IOW the algorithm is not really consistent; the same docking parameters will give different scores over and over again)
  - Exhaustiveness INCREASE

- Is there a reason why in your script you commented out halogen transformations? If so why?
  - In my runs adding fluorines or chlorines seems to be a pretty favored choice for the algorithm and I'm not sure how realistic that is
  - Limit to maybe 2 or 3 halogenations

  **Where to go from here**

  1. First thing first - *need to make vina calculate more consistent dockings scores - or else use a different method (free energy?) entirely*
     1. Shouldn't runs with identical parameters give at least near-identical results? Right now it's rare for vina to consistently decide which molecule one gen removed from the original parent is the best, let alone molecules at a greater depth
        1. **Solution: increase exhaustiveness *even more***
           1. Exhaustiveness has a roughly linear relationship with run time - doubling exhaustiveness will double run time
  2. *Changing the structure of the algorithm*
     1. right now the program sort of 'marches forward' until it reaches the required depth. It does pick the molecule with the best docking score from each generation, but said molecule isn't necessarily better than it's parent
        1. This has the effect of skirting potential local minima; it may give a better picture of the overall energy landscape (assuming vina results can be trusted), but it will not directly seek out local minima
     2. My intuition is to restructure the algorithm such that it only moves on from a given generation if a child molecule has a better docking score than it's parent. Given the margin of error with Vina, it may make sense to re-dock a given generation (or maybe a portion - top 20%?) a number of times to ensure that none of the molecules are actually better than the parent. In a perfect world where vina is 100% accurate, there would be no need to re-test but alas
        1. Re-docking of course is computationally expensive and will likely greatly lengthen the run time of the program. Will probably move away from `depth` as the guiding parameter and instead just have a `do while` loop to run the algorithm until a local minimum is found
  3. Spoke with Jianing about *validating the algorithm* with another S3 ligand - traversing from propofol to the second ligand with the algorithm would validate the method
     1. Where to find a ligand like this
        1. Try aligning pdb files to 1e7a
  4. *Future Direction* - moving away from vina and rdkit to using a schrodinger's python API to directly manipulate ligands whilst still in the pocket and calculating the effect a given change has on the free energy
     1. This method will be more accurate than constantly redocking with vina, which has to try docking molecules starting at a random conformation
  5. *AutoGrow4* - sounds a lot like our algorithm. You mentioned it in your own research proposal claiming that it had a vague continuous definition of chemical space, whereas our algorithm defines CHSP discretely. Can you speak more to that difference?
     1. It uses a GA to optimize leads, combine fragments, etc; no stepwise chemical transformations
     2. Seems like the major difference with our algorithm vs others is an exhaustive, stepwise-exploration of chemical space

### 8/15/22

#### Congruency Bug - FIXED

**Description**: the molecular editing gauntlet allows for adding an atom, mutating an atom/bond, and removing an atom/bond. These later two transformations represent potential backwards moves in chemical space, potentially leading to congruent molecules across multiple generations (not just immediately neighboring generations either). This is a highly inefficient bug as it wastes resources docking molecules that have already been tested. So the question remains, how do we ensure that no duplicate molecules inadvertently recreated through the gauntlet get docked more than once?

**Solution**: The answer to this question lies within sets. I created a `master_set` that initially contains just the smile string of the parent molecule, `parent_0`. Each generation, we go through the entire gauntlet as per usual. However, we define the list `uniq_desc` as the *difference* between the master set and the current generation of descendants. In other words, `uniq_desc` will only contain molecules that are not in the master set, or *have not been created before*. The master set is then joined with the uniq list of novel descendants `uniq_desc` via the `.union()` method, yielding an updated master set containing all uniq descendants to date.

**Note: I haven only tested this patch in jupyter notebooks, not yet with an entire program run**

### 8/14/22

- added error logging to script when running dock_it iteratively (ln 586); produces `error_log.txt`
- Using regex, fixed bug within `dock_it()` that had prevented the best docking score from the autodock vina output (log files) from being read in correctly (ln 329)
- Ran entire program with `depth=5` and propofol as `parent_0`. This run was called `run_5_depth_5`

#### Bugs to fix

- **ACTIVE ** the grid images produced each loop. sometimes more than one molecule is highlighted. This is because if the best molecule is the result of a mutation, not an addition, then at least one other molecule created via addition will share the substructure of the best mutated molecule (for example, a molecule might gain a bond to chlorine (addition), whereas another may mutate a carbon to a chlorine. These two molecules would share a common substructure)

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
