import argparse as ap
import cProfile
import pstats
from ChemTools import *
import csv
import time

# bgRun=cProfile.Profile()
# bgRun.run('buildGraph("C",3,True)')
# stats = pstats.Stats(bgRun)

# stats.strip_dirs()
# stats.sort_stats('cumtime')
# stats.print_stats()
# print("Total time: ", stats.total_tt)


# new_smis = ['C','CC','CO','CJK']
# new_mols = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in new_smis if Chem.MolFromSmiles(smi) is not None]

# print(new_mols)

runData = []
seed = 'C'
for depth in range(1, 4):
	for exhaustive in [True, False]:
		bgRun=cProfile.Profile()
		run_params = "buildGraph("+"'"+seed+"'"+","+str(depth)+","+str(exhaustive)+")"
		start_time = time.time()
		graph = buildGraph(seed,depth,exhaustive)
		print(run_params)
		fae_params = "faerunPlot(graph,'test_graph')"
		bgRun.run(fae_params)
		end_time = time.time()
		stats = pstats.Stats(bgRun)
		stats.strip_dirs()
		stats.sort_stats('cumtime')
		log_file_name = f"faerunPlot_profile_{seed}_d{depth}_ec{exhaustive}.txt"
		with open(log_file_name, "w") as f:
			f.write(f"Seed: {seed}\n")
			f.write(f"Depth: {depth}\n")
			f.write(f"Exhaustive: {exhaustive}\n\n")
			stats.stream = f
			stats.print_stats()
			f.write(f"Total faerunPlot time: {stats.total_tt}\n")
			f.write("Total buildGraph+faerunPlot time: "+reportTime(end_time-start_time))
		runData.append([seed,depth,exhaustive,stats.total_tt])

csv_file_name = "Linux_C_plot_Data.csv"
with open(csv_file_name, "w", newline='') as csvfile:
	csvfile = open(csv_file_name, "w", newline='')
	writer = csv.writer(csvfile)
	writer.writerow(["Seed", "Depth", "Exhaustive", "Total Time"])
	for run in runData:
		writer.writerow(run)