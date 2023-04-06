#!/usr/bin/env python3

import argparse
from ChemTools import *
import datetime
import time

################################################################# BEGIN SCRIPT #####################################################################################

# start by creating command-line interface
csBuilder = argparse.ArgumentParser(prog="CSGraphExplorer 1.0",description="a program that builds a local chemical space graph starting at a given seed molecule. Returns an HTML file containing a faerun visualization of the graph")
# SMILE seed is the only positional argument
csBuilder.add_argument('SMILE',type=str,help="the SMILE string of the molecule to seed the graph with")
# depth and cc have defaults
csBuilder.add_argument('--depth','-d',type=int,metavar=None,default=1,help="integer specifying depth of the graph. default is 1.")
csBuilder.add_argument('--exhaustive_connections','-ec',action="store_true",help="option to add exhaustive connections to the graph. default is false.")
csBuilder.add_argument('--name','-n',type=str,metavar=None, default = None,help="option to provide alternamtive name of the seed molecule (common name, abbreviation, etc). Default name of the seed is simply the SMILE string.")
csBuilder.add_argument('--plot','-pl',action="store_false",help="option to plot graph. Default is true")

# parses args from terminal
args = csBuilder.parse_args()

seed = args.SMILE
depth = args.depth
cc = args.exhaustive_connections
name = args.name
plot = args.plot


####################### chatGPT ################


###### Build the graph ######
# Record start time
start_time = time.time()

# Run buildGraph function
csg = buildGraph(seed,depth,cc)

# Record time taken for buildGraph function
buildGraph_time = time.time() - start_time

### name graph if provided
if name is None:
	csg_name = seed+"_d"+str(depth)+"_cc"+str(cc)
else: csg_name = name+"_d"+str(depth)+"_cc"+str(cc)

##### Plot the Graph #####
plot_stamp = time.time()
if plot:
	# get list of smile strings from graph, get corresponding fingerprints list
	node_labels, fps = get_node_labels_fps(csg)

	# Run faerunPlot function:
	faerunPlot(csg, csg_name,node_labels,fps)
# Record time taken for faerunPlot function
faerunPlot_time = time.time() - plot_stamp

####################### END chatGPT ################



##### Save molecular properties data #####
data_stamp = time.time()

# if the plot argument was given, then plot is false and node_labels must be created
if plot is False:
	node_labels = []
	for smile in csg.nodes():
		node_labels.append(smile)

# get the molecular properties nested lists
NHD, NHA, MWT, MLP, MMR, NAT, PSA, qed = getPropList(node_labels)

data_time = time.time() - data_stamp

prop_cols = {"SMILE":node_labels, "nhd":NHD, "nha":NHA, "mwt:":MWT, "mlp":MLP, "mmr":MMR, "nat":NAT, "psa":PSA, "qed":qed}
df = pd.DataFrame(prop_cols)
fname = seed + "_d" + str(depth) + "_ec" + str(cc) + ".csv"
df.to_csv(fname, index=False)

##### Logging & Results ######
# Calculate total time taken for entire program
total_time = buildGraph_time + faerunPlot_time

# Create filename for text file
log_file_name = "log_"+seed+"_d"+str(depth)+"_ec"+str(cc)+".txt"

# Write data to text file
# with open(filename, "w") as f:
f = open(log_file_name, "a")
f.write("\n\nBuildGraph time: " + reportTime(buildGraph_time) + "\n")
if plot:
	f.write("faerunPlot time: " + reportTime(faerunPlot_time) + "\n")
f.write("save data time: "+reportTime(data_time)+"\n")
f.write("total time: " + reportTime(total_time) + "\n")
f.close()
