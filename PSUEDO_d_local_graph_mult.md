### PSUEDOCODE for the Script

before anything, the receptor can be prepared

`prepare receptor`

obviously we start with a given mol, which happens to be propofol in this case. Let's refer to this molecule as our parent molecule, as all mutants will be 'descendants' of this molecule in a sense.

`parent_mol='CC(C)C1=C(C(=CC=C1)C(C)C)O'`

Next we should define the generational depth of our graph exploration. Generational depth can be thought of as the extent of our exploration of the chemical space surrounding our parent molecule. For example, a generational depth of 5 means that we will explore descendants that are within 5 steps of our parent on the chemical space graph. This parameter defines the exhaustiveness of our search. This exhaustiveness can be **"wide"** or **"narrow"**. **Wide** would mean that we explore all the descendants (within our depth) of every descendant for each generation; in this way, our pool of descendants to explore grows exponentially. This knd of search is more a more thorough (dare I say complete) exploration of the chemical space surrounding our parent, however it is much more computationally expensive. A **Narrow** exploration, by contrast, means we explore the offspring of *only the best descendant* for a given generation. This is a computationally efficient path to a local minimum, but the tradeoff is that it is a limited search that ultimately leads us to one local minimum.

`depth=5`

Now we can put the parent molecule through the molecular editing gauntlet. Each descendent molecule will be saved in an array, call it something like `all_descendants` (this is equivalent to the list `new_mol` in Jacob's script). it's likely that `all_desc` contains identical molecules through symmetry. We need to remove those identical smile strings. That will leave us with a clean array containing SMILEs of only the all the unique descendants; we can call this `uniq_desc` or simply `descendants`


`for i in range(depth):
	# this is the arra
	all_desc=[]
	# begin molecular edits #`

Now that we have an array of unique descendants, it's time to dock each one and save the results to a list (`results`). It's a good idea to have one array that contains the SMILE for each molecule assessed and it's corresponding affinity, and also one that contains this information for only the best-performing molecule from each generation (let's call it `best_path`). That way, when the algorithm is over, we will end up with an array of length equal to depth+1, containing our parent molecule and its affinity are at index 0, and the best 5th-generation molecule and its afinnity at index 5.

*this stores all the results from a given generation as a list of tuples, where each tuple contains the result for a single descendant *
`results=[]
for j in uniq_desc:
	results.append((j,dock_it(j))`

this adds all results from a given generation into one list, which is then added to the master list `generations` The parent molecules is added to the beginning of `generation`, and will have index `[0][0]`. The numeric name for each molecule and its respective files will correspond to its index in `generations`. Thus, `generations` is basically the master list, containing all descendants explored throughout the algorithm and their docking scores

`# this initializeds generations such that genations[0][0] will give the string and affinity of the parent molecule`
`generations=[[(parent_mol,dock_it(parent_mol))]]`
`# this records all the results from each iteration
generations.append(results)`

Now let's write a line to get the best molecule and it's results from each generation
`affinities=[x[1] for x in results]`
`best_index=affinities.index(min(affinities))`
`best_path.append(all_results[best_index])`

Now make the best molecule from the previous generation the new parent, and begin the loop again:
`parent_mol=best_path[-1][0]

** At this point, we have saved a total list of all the results, and a list with length=depth+1 (because it will include the parent mol, gen 0) of the best results from each generation **
Later on it might be a code idea to print images of each molecule and it's corresponding affinity like jacob did.

#### END OF ITERATION
repeat now for depth of 2
