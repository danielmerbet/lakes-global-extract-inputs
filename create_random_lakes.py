import argparse
import random
import csv

parser = argparse.ArgumentParser()

parser.add_argument('count', type=int)

args = parser.parse_args()

id_format = 'lake_{{:0{}d}}'.format(len(str(args.count)))

with open('lakes.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['id', 'lat', 'lon'])
    for i in range(args.count):
        writer.writerow([id_format.format(i+1), random.uniform(-90, 90), random.uniform(-180, 180)])
