#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
import random
from random import randint
import argparse
import re
import os
import sys
from operator import itemgetter
random.seed(9001)
import statistics
import networkx as nx
import matplotlib

__author__ = "Moreau Baptiste"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Moreau Baptiste"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Moreau Baptiste"
__email__ = "moreaubapt@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()



def read_fastq(fastq_file):
    #Read the sequences from a file and assert they are well formed otherwise raise an exception
    with open(fastq_file, 'r') as f:
        while True:
            line = f.readline()
            if len(line) ==0 :
                break
            line2 = f.readline()
            f.readline()
            f.readline()
            yield line2[:-1]



def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    dic = {}
    for i in read_fastq(fastq_file):
        for j in cut_kmer(i,kmer_size):
            try:
                dic[j] +=1
            except KeyError:
                dic[j] = 1
    return dic

def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for key in kmer_dict.keys():
        graph.add_edge(key[:-1],key[1:],weight=kmer_dict[key])
    return graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    pass

def get_sink_nodes(graph):
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    pass

def save_contigs(contigs_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    for i in read_fastq(args.fastq_file):
        print(i)
    #dic = build_kmer_dict(args.fastq_file,args.kmer_size)
    #print(dic)
if __name__ == '__main__':
    main()
