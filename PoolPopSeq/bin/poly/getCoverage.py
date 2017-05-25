from __future__ import division
import os, math, decimal, sys, argparse

#coverage = str(sys.argv[1])
#coverage_merged = str(sys.argv[2])
#evidence = str(sys.argv[3])

def condenseCoverage(IN, OUT_path):
    OUT = open(OUT_path, 'w')
    with open(IN) as f:
        coverage_max = 0
        position_start = 0
        count = 0
        for x ,line in enumerate(f):

            line_list = line.split()
            coverage = line_list[3]
            contig = line_list[0]
            position = line_list[1]
            if x == 0:
                coverage_max = coverage
                position_start = position
                count == 0
            else:
                if coverage == coverage_max:
                    count += 1
                else:
                    print>> OUT, contig, position_start, int(position_start) + int(count), coverage_max

                    coverage_max = coverage
                    position_start = position
                    count = 0

        print>> OUT, contig, position_start, int(position_start) + int(count), coverage_max

    OUT.close()

def selectSites(OUT1_path, OUT2_path):
    poly_dict = {}
    #OUT_subset =
    OUT_groups = OUT2_path.split('/')
    OUT_subset = open('/'.join(OUT_groups[:-1]) + '/coverage_merged_subset.txt', 'w')
    with open(OUT2_path) as f:
        for x ,line in enumerate(f):
            line_split = line.split()
            if len(line_split) < 3:
                continue
            seq_id = line_split[3]
            position = line_split[4]
            if seq_id in poly_dict:
                poly_dict[seq_id].append(int(position))
            else:
                poly_dict[seq_id] = [int(position)]
    with open(OUT1_path) as g:
        # g = gd file
        contig = ''
        L = 0
        for x ,line in enumerate(g):
            if line_split[0] != contig:
                contig = line_split[0]
            line_split = line.split()

            contig = line_split[0]
            sites = poly_dict[contig]
            l2 = [i for i in sites if i >= int(line_split[1]) and i <= int(line_split[2])]
            if len(l2) > 0:
                for l in l2:
                    sites.remove(l)
                    print int(line_split[1]), l, int(line_split[2])
                    print>> OUT_subset, L, contig, line_split[1], line_split[2], line_split[3]
                    print len(sites)
                L += 1


    OUT_subset.close()

#selectSites(coverage_merged, evidence)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "condense coverage")
    parser.add_argument('-c', action='store_true', default=False)
    parser.add_argument('-v', action='store_true', default=False)

    parser.add_argument('-i', type = str, default = "", help = "in file")
    parser.add_argument('-o', type = str, default = "", help = "out file")

    params = parser.parse_args()

    if params.v == False and params.c == True:
        condenseCoverage(params.i, params.o)
    else:
        print "No argument provided\nExiting program"
