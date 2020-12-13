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
    """Generator reading sequences contained in the file.
      :Parameters:
         fastq_file : Path to the file
    """
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
    """Generator cuting kmer contained in the sequence.
      :Parameters:
         read : sequence
         kmer_size : size of the kmer
    """
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    """Create a kmer dictionnary based on the sequences of a file with the specified size.
      :Parameters:
         fastq_file : Path of the file
         kmer_size : size of the kmer
    """
    dic = {}
    for i in read_fastq(fastq_file):
        for j in cut_kmer(i,kmer_size):
            try:
                dic[j] +=1
            except KeyError:
                dic[j] = 1
    return dic

def build_graph(kmer_dict):
    """Return the corresponding oriented graph of a kmer dictionnary.
      :Parameters:
         kmer_dict : kmer dictionnary
    """
    graph = nx.DiGraph()
    for key in kmer_dict.keys():
        graph.add_edge(key[:-1],key[1:],weight=kmer_dict[key])
    return graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """qui prend un graphe et une liste de chemin,
    la variable booléenne delete_entry_node pour indiquer si les noeuds d’entrée
    seront supprimés et
    la variable booléenne delete_sink_node pour indiquer si les noeuds de sortie
    seront supprimés et retourne un graphe nettoyé des chemins indésirables."""
    for path in path_list:
        for i in range(len(path)-1):

            try:
                graph.remove_edge(path[i],path[i+1])
                print("removing edge({},{})".format(path[i],path[i+1]))
            except nx.exception.NetworkXError:
                pass
        if delete_entry_node:
            graph.remove_node(path[0])
        if delete_sink_node:
            graph.remove_node(path[-1])
        for k in path:
            if graph.has_node(k):
                a=0
                for i,j in enumerate(graph.predecessors(k)):
                    a+=1
                for i,j in enumerate(graph.successors(k)):
                    a+=1
                if a == 0:
                    graph.remove_node(k)
    return graph


def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """prend un graphe, une liste de chemin, une liste donnant la longueur
    de chaque chemin, une liste donnant le poids moyen de chaque chemin,
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés
    et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés
     et retourne un graphe nettoyé des chemins indésirables.
    Par défaut, delete_entry_node et delete_sink_node seront ici à False"""
    best_weights=[0]
    for i,weight in enumerate(weight_avg_list):
        if i != best_weights[0] :
            if weight > weight_avg_list[best_weights[0]]:
                best_weights = best_weights[0:1]
                best_weights[0] = i
            if weight == weight_avg_list[best_weights[0]] and i != best_weights[0]:
                best_weights.append(i)
    best_length=[best_weights[0]]
    for i in best_weights:
        if i != best_length[0] :
            if path_length[i] > path_length[best_length[0]]:
                best_length = best_length[0:1]
                best_length[0] = i
            if path_length[i] == path_length[best_length[0]] and i != best_length[0]:
                best_length.append(i)
    best_path_index = best_length[random.randint(0,len(best_length)-1)]
    path_list.pop(best_path_index)
    return remove_paths(graph, path_list, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    """Return the avergae weight of a path.
      :Parameters:
         graph : graph
         path : path within the graph
    """
    data = []
    for i in range(len(path)-1):
        if graph.has_edge(path[i],path[i+1]) :
            data.append(graph.get_edge_data(path[i],path[i+1])["weight"])
    return statistics.mean(data)

def get_ancestors(graph,node):
    a=0
    if graph.has_node(node):
        for n in graph.predecessors(node):
            yield(a,n)
            a+=1

def get_bubble_paths(graph,starting_node, ending):
    paths=[]
    valids=[]
    for (a,n) in get_ancestors(graph,ending):
        l = [n,ending]
        paths.append(l)
    def constructPaths(graph,starting_node,paths,valids):
        prt_paths=[]
        for path in paths:
            j=0
            for (a,n) in get_ancestors(graph,path[0]):
                if j==0:
                    path.insert(0,n)
                else:
                    prt_paths.append(path[1:].insert(0,n))
                if n == starting_node:
                    valids.append(path)
                    paths.remove(path)
                j+=1
            if j == 0:
                paths.remove(path)
        paths = paths+prt_paths
        return (paths,valids)

    while len(paths) > 0:
        (paths,valids) = constructPaths(graph,starting_node,paths,valids)
    return valids



def solve_bubble(graph, ancestor_node, descendant_node):
    """qui prend un graphe, un nœud ancêtre, un nœud descendant et
    retourne un graph nettoyé de la bulle se trouvant
     entre ces deux nœuds en utilisant les fonctions précédemment développées"""
    paths = get_bubble_paths(graph, ancestor_node, descendant_node)
    paths_length = []
    paths_weight = []
    for path in paths:
        paths_length.append(len(path))
        paths_weight.append(path_average_weight(graph,path))
        print(path,len(path),path_average_weight(graph,path))

    return select_best_path(graph,paths,paths_length,paths_weight)

def simplify_bubbles(graph):
    nodes = list(graph.nodes)
    def solver(graph,ancestor_from,ending):
        for (a,starting_node) in get_ancestors(graph,ancestor_from):
            paths = get_bubble_paths(graph,starting_node, ending)
            if len(paths) > 1 :
                return solve_bubble(graph,starting_node,ending)
        return False
    for i in range(len(nodes)):
        for j in range(i,len(nodes)):
            x = solver(graph,nodes[i],nodes[j])
            y = x
            while x != False:
                y = x
                x = solver(x,nodes[i],nodes[j])
            if y != False:
                graph = y
    return graph


def get_successor(graph,node):
    a=0
    if graph.has_node(node):
        for n in graph.successors(node):
            yield(a,n)
            a+=1


def solve_entry_tips(graph, starting_nodes):
    ancestor = next(get_successor(graph,starting_nodes[0]))
    def asCommonSuccesor(graph,starting_nodes,ancestor):
        future_succ = []
        for node in starting_nodes:
            successors = next(get_successor(graph,node))
            future_succ.append(successors)
            if successors == ancestor:
                return True
        else:
            return future_succ
    x = asCommonSuccesor(graph,starting_nodes,ancestor)
    while x != True:
        while x != []:
             x = asCommonSuccesor(graph,x,ancestor)
             if x == True:
                 break
        if x != True:
            ancestor = next(get_successor(graph,starting_nodes[0]))
            x = asCommonSuccesor(graph,starting_nodes,ancestor)
    print(ancestor[1])
    for node in starting_nodes:
        print(get_bubble_paths(graph,node,ancestor[1]))





def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    """Return nodes of the graph with no ancestor.
      :Parameters:
         graph : the graph
    """
    starting_nodes = []
    for n in graph.nodes():
        a=0
        for i,j in enumerate(graph.predecessors(n)):
            a+=1
        if a == 0:
            starting_nodes.append(n)
    return starting_nodes

def get_sink_nodes(graph):
    """Return nodes of the graph with no successors.
      :Parameters:
         graph : the graph
    """
    sink_nodes = []
    for n in graph.nodes():
        a = 0
        for i,j in enumerate(graph.successors(n)):
            a+=1
        if a == 0:
            sink_nodes.append(n)
    return sink_nodes

def get_contig(graph, starting_node, ending):
    """Return the contig of the graph for the specified starting node and ending node
       or 'FALSE' if theree is no such contig.
      :Parameters:
         graph : the graph
         starting_node : the node that will be used as start position in the graph to
         find contig
         ending : the node that will be used as end position in the graph to find contig
    """
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
    """Return the contigs of the graph as tuple formed as
    (contig,len of the contig) for the specified starting nodes and ending nodes.
      :Parameters:
         graph : the graph
         starting_node : the nodes that will be used as start position in
         the graph to find contigs
         ending_nodes : the nodes that will be used as end position in
         the graph to find contigs
    """
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
    """
    :Parameters:
         contigs_list : the contig list
         output_file : Path of the file
    """
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
    #graph.add_edges_from([("TC", "CA"), ("AC", "CA"), ("CA", "AG"), ("AG", "GC"),
    # ("GC", "CG"), ("CG", "GA"), ("GA", "AT"), ("GA", "AA")])
    #contig_list = get_contigs(graph, ["TC", "AC"], ["AT" , "AA"])
    #save_contigs(contig_list,"oui")
    #print(contig_list)
    #dic = build_kmer_dict(args.fastq_file,args.kmer_size)
    #print(dic)
if __name__ == '__main__':
    main()
