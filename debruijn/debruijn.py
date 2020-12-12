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
random.seed(9001)
import statistics
import networkx as nx
import matplotlib
import argparse
import os
import sys
from operator import itemgetter


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
    starting_nodes = []
    for n in graph.nodes():
        a=0
        for i,j in enumerate(graph.predecessors(n)):
            a+=1
        if a == 0:
            starting_nodes.append(n)
    return starting_nodes

def get_sink_nodes(graph):
    sink_nodes = []
    for n in graph.nodes():
        a = 0
        for i,j in enumerate(graph.successors(n)):
            a+=1
        if a == 0:
            sink_nodes.append(n)
    return sink_nodes

def get_contig(graph, starting_node, ending):
    c = starting_node
    if starting_node == ending:
        return starting_node
    a = 0
    for n in graph.successors(starting_node):
        a+=1
        #print(starting_node,ending,"get contig with starting",n)
        r = get_contig(graph,n,ending)
        if r != 'FALSE' :
            c+=r[1:]
            #print(starting_node,ending,c,"get contig with ending")
            return c
    if a == 0:
        #print(starting_node,ending,"pas de successors")
        return "FALSE"
    else:
        return 'FALSE'

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for sN in starting_nodes:
        for eN in ending_nodes:
            #print(sN,eN,"starting")
            c = get_contig(graph,sN,eN)
            #print(c)
            if c != 'FALSE':
                contigs.append((c,len(c)))
    return contigs

def save_contigs(contigs_list, output_file):
    def fill(text, width=80):
        """Split text with a line return to respect fasta format"""
        return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

    with open(output_file,"w") as file:
        for i, contig in enumerate(contigs_list):
            seq, t = contig
            file.write(">contig_{0} len={1}\n".format(i,t))
            file.write(fill(seq+"\n"))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    pass
    # Get arguments
    #args = get_arguments()
    #graph = nx.DiGraph()
    #graph.add_edges_from([("TC", "CA"), ("AC", "CA"), ("CA", "AG"), ("AG", "GC"), ("GC", "CG"), ("CG", "GA"), ("GA", "AT"), ("GA", "AA")])
    #contig_list = get_contigs(graph, ["TC", "AC"], ["AT" , "AA"])
    #save_contigs(contig_list,"oui")
    #print(contig_list)
    #dic = build_kmer_dict(args.fastq_file,args.kmer_size)
    #print(dic)
if __name__ == '__main__':
    main()
