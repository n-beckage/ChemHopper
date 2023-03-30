#!/usr/bin/env python

import argparse as ap
from ChemTools import *

################################################################# BEGIN SCRIPT #####################################################################################

# start by creating command-line interface
csBuilder = ap.ArgumentParser(prog="CSGraphExplorer 1.0",description="a program that builds a local chemical space graph starting at a given seed molecule. Returns an HTML file containing a faerun visualization of the graph")
# SMILE seed is the only positional argument
csBuilder.add_argument('SMILE',type=str,help="the SMILE string of the molecule to seed the graph with")
# depth and cc have defaults
csBuilder.add_argument('--depth','-d',type=int,metavar=None,default=1,help="integer specifying depth of the graph. default is 1.")
csBuilder.add_argument('--exhaustive_connections','-ec',action="store_true",help="option to add exhaustive connections to the graph. default is false.")
# parses args from terminal
args = csBuilder.parse_args()

# print("\nWelcome to Chemical Space! Where would you like to begin?\n")

# print("please enter a molecule in the form of a SMILE string.\n")

# bad_smi = True
# while bad_smi:
# 	seed = input("GRAPH SEED: ")
# 	if Chem.MolFromSmiles(seed) is None:
# 		print(seed,"is not a valid SMILE string. Please try again.")
# 	elif Chem.MolFromSmiles(seed) is not None:
# 		bad_smi = False

# print("\nHow deep will the roots go?\n")

# bad_depth = True
# while bad_depth:
# 	depth = input("GRAPH DEPTH: ")
# 	try:
# 		depth = int(depth)
# 		bad_depth = False
# 	except:
# 		print("Invalid input. Please enter an integer.")

# yn=''
# while yn not in ['y','Y','n','N']:
# 	yn = input("\nInclude exhaustive connections? (y/n):")
# 	if yn in ['y','Y']:
# 		cc = True
# 	elif yn in ['n','N']:
# 		cc = False
# 	else:
# 		print(yn,"is invalid input. Please enter y/n")

# print()

# running buildGraph() with args
csg = buildGraph(args.SMILE,args.depth,args.exhaustive_connections)

gname = input("Please enter a filename for the graph: ")
csg_name = gname+"_d"+str(depth)+"_cc"+str(cc)
faerunPlot(csg, csg_name)