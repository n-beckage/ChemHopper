import argparse as ap
import cProfile
import pstats
from ChemTools import *
import csv
import pandas as pd
import time


node_labels = ["C","CC","CCC"]

NHD, NHA, MWT, MLP, MMR, NAT, PSA, qed = getPropList(node_labels)

print(qed)

prop_cols = {"SMILE":node_labels, "nhd":NHD, "nha":NHA, "mwt:":MWT, "mlp":MLP, "mmr":MMR, "nat":NAT, "psa":PSA, "qed":qed}
df = pd.DataFrame(prop_cols)
fname = "testData.csv"
df.to_csv(fname, index=False)