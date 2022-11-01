import argparse

import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('var')
parser.add_argument('file1')
parser.add_argument('file2')

args = parser.parse_args()

df1 = pd.read_csv(args.file1)
df1_index = df1.columns[0]
df1.set_index(df1_index)

df2 = pd.read_csv(args.file2)
df2_index = df2.columns[0]
df2.set_index(df2_index)

fig, ax = plt.subplots(figsize=(30, 10))

df1.plot(ax=ax, x=df1_index, y=args.var)
df2.plot(ax=ax, x=df2_index, y=args.var)

plt.show()
