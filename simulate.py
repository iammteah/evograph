from asymmetree.datastructures  import  PhyloTree
from asymmetree.datastructures.Tree import LCA

import asymmetree.treeevolve  as te
import asymmetree.hgt as hgt
from asymmetree.cograph  import  Cotree

import asymmetree.tools.GraphTools  as gt
import itertools
import matplotlib.pyplot as plt
import random

import networkx as nx
from networkx.algorithms import approximation
from networkx.algorithms import operators

# GROUP PARAMETERS:
HGT_REPLACE = 1.0
DEATH_RATE_INCR = [0.0, 0.25, 0.5, 0.75, 0.9]
BIRTH_RATE_CONST = 1.0

# OTHER PARAMETERS:
NON_BINARY = random.uniform(0.0,0.5)
NUM_SPC_GENE_TREE_PAIRS = 100
NUM_LEAFS = 10 #is set
RNG_LEAFS = random.uniform(10,50) #set finally
DUPL_LOSS_HGT = [
    (0.25,0.25,0.25),
    (0.5, 0.5, 0.5),
    (0.5, 0.5, 1.0),
    (0.5, 0.5, 1.5),
    (1.0, 1.0, 0.5),
    (1.0, 1.0, 1.0),
    (1.5, 1.5, 1.5)
]

# metrics for analysis 
def add_metrics(ct):
    try:
        ct["recall"] = ct["tp"]/(ct["tp"]+ct["fn"])
        ct["precision"] = ct["tp"]/(ct["tp"]+ct["fp"])
        ct["accuracy"] = (ct["tp"] + ct["tn"])/(ct["tp"]+ct["tn"]+ct["fp"]+ct["fn"])
    except:
        pass
    return ct

# graph tools
def complete_graph_from_list(L, create_using=None):
    nodes = L
    edges = itertools.combinations(nodes, 2)
    g = nx.Graph()
    g.add_nodes_from(nodes)
    g.add_edges_from(edges)
    return g

def get_color_of_node(_id, _tree):
    nodes = _tree.sorted_nodes()
    for n in nodes:
        if n.label == id:
            return n.color

#Cluser Deletion
def get_clique_list(cograph):
    clique_list = []
    
    c = approximation.max_clique(cograph)
    while c != set():
        clique_list.append(complete_graph_from_list(c))
        for n in c:
            cograph.remove_node(n)
        c = approximation.max_clique(cograph)
        
    return clique_list

skipped = 0
tree_pair_counter = 0
for _rate in DEATH_RATE_INCR:
    for _dupl, _loss, _hgt in DUPL_LOSS_HGT:
        
        print ("Results for dupl:{} loss:{} hgt:{}".format(_dupl, _loss, _hgt))
        fraction = 0
        for i in range(NUM_SPC_GENE_TREE_PAIRS):
            
            #CREATE SIMULATED GRAPHS AND TREES
            species_tree_S = te.simulate_species_tree (NUM_LEAFS, model= "BDP", non_binary_prob = NON_BINARY, planted = 1, birth_rate = BIRTH_RATE_CONST, death_rate = _rate, rescale_to_height = 1.0) # Are these the correct parameters?
            observable_species_tree_S = te.observable_tree(species_tree_S)
            gene_tree_G = te.simulate_dated_gene_tree(species_tree_S, dupl_rate = _dupl, loss_rate = _loss, hgt_rate = _hgt, prohibit_extinction='per_familiy', replaceprob = HGT_REPLACE)
            observable_gene_tree_OG = te.observable_tree(gene_tree_G)

            ldt_graph_LDT = hgt.ldt_graph(observable_gene_tree_OG , species_tree_S)
            transfer_edges_TE = hgt.true_transfer_edges(observable_gene_tree_OG)
            fitch_FG = hgt.undirected_fitch(observable_gene_tree_OG , transfer_edges_TE)
            
            cotree_CT = Cotree.cotree(ldt_graph_LDT)

            # STORE ORIGINAL TREES
            # BUILD IDs
            tree_pair_counter =+ 1
            S_ID=tree_pair_counter + 'S' + str(_rate +_dupl + _loss +_hgt)
            T_ID=tree_pair_counter + 'T' + str(_rate +_dupl + _loss +_hgt)
            #species_tree_S.serialize('./data/species_tree_' + S_ID + '.pickle')
            #gene_tree_G.serialize('./data/gene_tree_' + T_ID + '.pickle')

            # ---
            # SPECIFIC QUESTION 1
            # ---
            
            matches = 0
            for t_edge in fitch_FG.edges:
                if ldt_graph_LDT.has_edge(*t_edge):
                    matches += 1
            try:
                fraction += matches / len(fitch_FG.edges)
            except:
                skipped += 1
                
            # ---
            # SPECIFIC QUESTION 2
            # ---
            
            # - CLUSTER DELETION METHOD FOR CONSTRUCTING FITCH GRAPH
            
            complement_ldt = operators.complement(ldt_graph_LDT)
            clique_graph_list = get_clique_list(complement_ldt)
            complete_graph = nx.Graph()
            for _g in clique_graph_list:
                complete_graph = nx.compose(complete_graph,_g)
            cluster_deletion_fitch_graph = operators.complement(complete_graph)

            # - METHOD FOR CONSTRUCTING RS FITCH GRAPH
            
            constructor = hgt.RsScenarioConstructor(ldt_graph_LDT)
            
            if constructor.run():
                S2 = constructor.S
                T2 = constructor.T
                transfer_edges2 = hgt.rs_transfer_edges(T2, S2)
                rs_fitch = hgt.undirected_fitch(T2 , transfer_edges2)
            
            
            ct_cluster_deletion = gt.contingency_table(fitch_FG , cluster_deletion_fitch_graph)
            ct_cluster_deletion = add_metrics(ct_cluster_deletion)
            ct_rs = gt.contingency_table(fitch_FG , rs_fitch)
            ct_rs = add_metrics(ct_rs)
            
            print("Cluster Deletion: ")
            print(ct_cluster_deletion)
            print("\nRS: ")
            print(ct_rs)
            
            # ---
            # SPECIFIC QUESTION 3
            # ---
            
            species_triples = observable_species_tree_S.get_triples(id_only=True)
            gene_triples = observable_gene_tree_OG.get_triples(id_only=True)
            
            T_G = []
            S_G_o = []
            
            ldt_triples = itertools.combinations(ldt_graph_LDT, 3)
            ldt_triples_cleaned = []
            
            for triple in ldt_triples: # add combinations not generated by itertools.combinations
                ldt_triples_cleaned.append(triple)
                ldt_triples_cleaned.append((triple[2],triple[1],triple[0]))
                ldt_triples_cleaned.append((triple[0],triple[2],triple[1]))
            
            ldt_colors = nx.get_node_attributes(ldt_graph_LDT, 'color')
            for triple in ldt_triples_cleaned:
                
                if triple[0] != triple[1] and triple[0] != triple[2] and triple[1] != triple[2]:
                    if ldt_graph_LDT.has_edge(triple[0], triple[1]) and not ldt_graph_LDT.has_edge(triple[0], triple[2]) and not ldt_graph_LDT.has_edge(triple[1], triple[2]):
                        T_G.append(triple)
                        
                if ldt_colors[triple[0]] != ldt_colors[triple[1]] and ldt_colors[triple[1]] != ldt_colors[triple[2]] and ldt_colors[triple[0]] != ldt_colors[triple[2]]:
                    if ldt_graph_LDT.has_edge(triple[0], triple[2]) and ldt_graph_LDT.has_edge(triple[1], triple[2]) and not ldt_graph_LDT.has_edge(triple[0], triple[1]):
                        S_G_o.append(triple)
                        
            print(len(T_G)/len(gene_triples))
            print(len(S_G_o)/len(species_triples))
            
            break
            
        fraction = fraction / NUM_SPC_GENE_TREE_PAIRS
        print("Xenolog Fraction visible in LDT Graphs: " + str(fraction))
        break
    break
    print("Skipped " + str(skipped) + " fitch graphs due to zero edges.")
