import argparse as ap
import cProfile
from ChemTools import *

bgRun=cProfile.Profile()
bgRun.run('buildGraph("C",2,True)')
stats = pstats.Stats(bgRun)

stats.strip_dirs()
stats.sort_stats('cumtime')
stats.print_stats()
print("Total time: ", stats.total_tt)

