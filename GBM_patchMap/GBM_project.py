"""
                                    ---------------------
                                     P A T C H  -  M A P
                                    ---------------------
                        GBM PROJECT - STRUCTURAL MAPPING OF THE NETWORKS
                                          Analysis

                                        Cansu Dincer
                                    Nurcan Tuncbag, PhD
"""
########################################################################################################################
########################################################################################################################
# IMPORTING #

from GBM_curation import Proteome, Int3D_Protein, Interfaces_II, Interfaces_PRISM, \
    Protein_Interfaces, Protein_Mutations, Mutations, Domain
from GBM_curation import read_patient_networks, human_proteome, interactome_3d_protein, \
    interactome_insider_interactions, PRISM_interactions
from GBM_curation import domain, protein_interface, protein_mutations, mutations
from GBM_curation import all_nodes, all_edges, gene_unip, all_node_uniprot_dict, all_edge_uniprot_dict

import pandas, pickle, numpy, os, networkx, re
import seaborn as sns
import matplotlib.pyplot as plt

########################################################################################################################
########################################################################################################################
# Numbers of nodes/edges & structural nodes/edges with sources

# Read patient network dictionary
x = print(domain_list1)
patient_networks = read_patient_networks()

patient_network_info_df = pandas.DataFrame(0, index=list(patient_networks.keys()),
                                           columns=["#Nodes", "#Edges", "#Str_Nodes", "#Str_Edges",
                                                    "#PDB_Nodes", "#ModBase_Nodes",
                                                    "#PDB_Edges", "#I3D_Edges", "#PRISM_Edges",
                                                    "#ECLAIR_Edges", "#NonStr_Nodes", "#NonStr_Edges"])

i3d = interactome_3d_protein()
ii = interactome_insider_interactions()
prism = PRISM_interactions()
interface = protein_interface()
proteome = human_proteome()


def check_edge_source(edge):
    sources = []
    if edge[0] in interface.keys() and interface[edge[0]].interface_source != {}:
        if edge[1] in interface[edge[0]].interface_source.keys():
            for d in interface[edge[0]].interface_source[edge[1]]:
                for k in d.keys():
                    if k not in sources: sources.append(k)
    if edge[1] in interface.keys() and interface[edge[1]].interface_source != {}:
        if edge[0] in interface[edge[1]].interface_source.keys():
            for d in interface[edge[1]].interface_source[edge[0]]:
                for k in d.keys():
                    if k not in sources: sources.append(k)
    if "PDB" in sources:
        return "PDB"
    elif "I3D" in sources:
        return "I3D"
    elif "PRISM" in sources:
        return "PRISM"
    elif "ECLAIR" in sources:
        return "ECLAIR"
    else:
        return None


def check_node_source(node):
    if node in i3d.keys() and i3d[node].coverageFile != dict():
        for l in i3d[node].coverageFile[1]:
            if l[2] == 1:
                if l[3] == "PDB":
                    return "PDB"
                elif l[3] == "MODBASE":
                    return "MODBASE"

    else:
        return None


for patient, ntw in patient_networks.items():
    patient_network_info_df.at[patient, "#Nodes"] = len(ntw.nodes)
    patient_network_info_df.at[patient, "#Edges"] = len(ntw.edges)
    patient_network_info_df.at[patient, "#Str_Nodes"] = len(
        [node for node in ntw.nodes if node in i3d.keys() and i3d[node].coveredResidues != []])
    patient_network_info_df.at[patient, "#Str_Edges"] = len(
        [edge for edge in ntw.edges if check_edge_source(edge=edge) != None])
    # Only Rank 1-1 will be considered
    patient_network_info_df.at[patient, "#PDB_Nodes"] = len(
        [node for node in ntw.nodes if check_node_source(node=node) == "PDB"])
    patient_network_info_df.at[patient, "#ModBase_Nodes"] = len(
        [node for node in ntw.nodes if check_node_source(node=node) == "MODBASE"])
    patient_network_info_df.at[patient, "#PDB_Edges"] = len(
        [edge for edge in ntw.edges if check_edge_source(edge=edge) == "PDB"])
    patient_network_info_df.at[patient, "#I3D_Edges"] = len(
        [edge for edge in ntw.edges if check_edge_source(edge=edge) == "I3D"])
    patient_network_info_df.at[patient, "#PRISM_Edges"] = len(
        [edge for edge in ntw.edges if check_edge_source(edge=edge) == "PRISM"])
    patient_network_info_df.at[patient, "#ECLAIR_Edges"] = len(
        [edge for edge in ntw.edges if check_edge_source(edge=edge) == "ECLAIR"])
    patient_network_info_df.at[patient, "#NonStr_Nodes"] = len(
        [node for node in ntw.nodes if node in i3d.keys() and i3d[node].coveredResidues == []])
    patient_network_info_df.at[patient, "#NonStr_Edges"] = len(
        [edge for edge in ntw.edges if check_edge_source(edge=edge) == None])

patient_network_info_df.to_csv(result_path + "Table_1.csv")


def take_percentages(all, x):
    p = (x * 100.0) / float(all)
    return p


percentages_patient_network_info_df = pandas.DataFrame(0, index=patient_networks.keys(),
                                                       columns=["%PDB_Nodes", "%MODBASE_Nodes",
                                                                "%NonStr_Nodes", "%PDB_Edges",
                                                                "%I3D_Edges", "%PRISM_Edges",
                                                                "%ECLAIR_Edges", "%NonStr_Edges"])

for patient, ntw in patient_networks.items():
    percentages_patient_network_info_df.at[patient, "%PDB_Nodes"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Nodes"], x=patient_network_info_df.loc[patient, "#PDB_Nodes"])
    percentages_patient_network_info_df.at[patient, "%MODBASE_Nodes"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Nodes"], x=patient_network_info_df.loc[patient, "#ModBase_Nodes"])
    percentages_patient_network_info_df.at[patient, "%NonStr_Nodes"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Nodes"], x=patient_network_info_df.loc[patient, "#NonStr_Nodes"])
    percentages_patient_network_info_df.at[patient, "%PDB_Edges"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Edges"], x=patient_network_info_df.loc[patient, "#PDB_Edges"])
    percentages_patient_network_info_df.at[patient, "%I3D_Edges"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Edges"], x=patient_network_info_df.loc[patient, "#I3D_Edges"])
    percentages_patient_network_info_df.at[patient, "%PRISM_Edges"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Edges"], x=patient_network_info_df.loc[patient, "#PRISM_Edges"])
    percentages_patient_network_info_df.at[patient, "%ECLAIR_Edges"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Edges"], x=patient_network_info_df.loc[patient, "#ECLAIR_Edges"])
    percentages_patient_network_info_df.at[patient, "%NonStr_Edges"] = take_percentages(
        all=patient_network_info_df.loc[patient, "#Edges"], x=patient_network_info_df.loc[patient, "#NonStr_Edges"])

node_percentages = percentages_patient_network_info_df.copy()
node_percentages = node_percentages.reset_index()
node_percentages = node_percentages.melt(id_vars=["index"],
                                         value_vars=["%PDB_Nodes", "%MODBASE_Nodes", "%NonStr_Nodes"])
node_percentages.columns = ["Patient", "Type", "Percentages"]
node_percentages["Method"] = node_percentages.apply(lambda x: "EXP" if x["Type"][1:].split("_")[0] == "PDB" else (
    "NO METHOD" if x.Type[1:].split("_")[0] == "NonStr" else "MODEL"), axis=1)
node_percentages["Column"] = "NODE"

edge_percentages = percentages_patient_network_info_df.copy()
edge_percentages = edge_percentages.reset_index()
edge_percentages = edge_percentages.melt(id_vars=["index"],
                                         value_vars=["%PDB_Edges", "%I3D_Edges", "%PRISM_Edges",
                                                     "%ECLAIR_Edges", "%NonStr_Edges"])
edge_percentages.columns = ["Patient", "Type", "Percentages"]
edge_percentages["Method"] = edge_percentages.apply(lambda x: "EXP" if x["Type"][1:].split("_")[0] == "PDB" else (
    "NO METHOD" if x.Type[1:].split("_")[0] == "NonStr" else "MODEL"), axis=1)
edge_percentages["Column"] = "EDGE"

percentages = pandas.concat([node_percentages, edge_percentages], axis=0)

flierprops = {"color": "black", "marker": "o", "markersize": 1.5}
colors = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
sns.set_palette(colors)
sns.boxplot(data=node_percentages, x="Type", y="Percentages", hue="Method",
            flierprops=flierprops)
plt.tight_layout()
plt.savefig(result_path + "Node_structure_percentages.pdf")
plt.savefig(result_path + "Node_structure_percentages.png")
plt.close()

sns.boxplot(data=edge_percentages, x="Type", y="Percentages", hue="Method",
            flierprops=flierprops)
plt.tight_layout()
plt.savefig(result_path + "Edge_structure_percentages.pdf")
plt.savefig(result_path + "Edge_structure_percentages.png")
plt.close()

########################################################################################################################
########################################################################################################################
# Huge Network

# Gene - Gene Pair
patient_networks = read_patient_networks()

# Gene -Domain Pair

domains = domain()
all_domains = [i + "-" + dom for i,l in domains.items() if l.domain_list != [] for dom in l.domain_list]
domain_gene_network = networkx.DiGraph()

for gene, obj in domains.items():
    for dom in obj.domain_list:
        dom_name = gene + "-" + dom
        domain_gene_network.add_edge(gene, dom_name)

# Domain - Interface Residue Pair

domain_interface_residue_network = networkx.DiGraph()

for gene, obj in domains.items():
    for dom, res_list in obj.interface_domains.items():
        dom_name = gene + "-" + dom
        for res in res_list:
            domain_interface_residue_network.add_edge(dom_name, str(res))

# Domain -Domain Pair
protein_muts = protein_mutations()
"""
gene_domain_res_network = networkx.Graph()
for edge in all_edges:
    edge1, edge2 = edge.split("-")[0], edge.split("-")[1]

    # Take the domain lists

    if edge1 in list(domain_gene_network.nodes):
        domain_list1 = [i for i in networkx.all_neighbors(domain_gene_network, edge1)]
    if edge2 in list(domain_gene_network.nodes):
        domain_list2 = [i for i in networkx.all_neighbors(domain_gene_network, edge2)]

    # Check the interface residues inside and if they match with the other partners' residues

    if domain_list1 != []:
        for dom1 in domain_list1:
            res_check = False
            mutated_res1 = []
            if dom1 in list(domain_interface_residue_network.nodes):
                for res1 in networkx.all_neighbors(domain_interface_residue_network, dom1):
                    if edge1 in interface.keys():
                        if edge2 in interface[edge1].interface_dict.keys():
                            if int(res1) in interface[edge1].interface_dict[edge2]:
                                res_check = True
                                if edge1 in protein_muts.keys():
                                    if res1 in protein_muts[edge1].mutatedResidues: mutated_res1.append(res1)
            if res_check:
                gene_domain_res_network.add_edge(edge1, dom1, node1 = "gene", node2 = "domain", edge_type = "GD")
                for mr in mutated_res1: gene_domain_res_network.add_edge(dom1, edge + "__" + mr,
                                                                         node1 = "domain", node2 = "Residue",
                                                                         edge_type = "DR")


    if domain_list2 != []:
        for dom2 in domain_list2:
            res_check = False
            mutated_res2 = []
            if dom2 in list(domain_interface_residue_network.nodes):
                for res2 in networkx.all_neighbors(domain_interface_residue_network, dom2):
                    if edge2 in interface.keys():
                        if edge1 in interface[edge2].interface_dict.keys():
                            if int(res2) in interface[edge2].interface_dict[edge1]:
                                res_check = True
                                if edge2 in protein_muts.keys():
                                    if res2 in protein_muts[edge2].mutatedResidues: mutated_res2.append(res2)
            if res_check:
                gene_domain_res_network.add_edge(edge2, dom2, node1 = "gene", node2 = "domain", edge_type = "GD")
                for mr in mutated_res2: gene_domain_res_network.add_edge(dom2, edge + "__" +mr,
                                                                         node1 = "domain", node2 = "Residue",
                                                                         edge_type="DR")

    if domain_list1 != [] and domain_list2 != []:
        one, two = False, False
        for dom1 in domain_list1:
            for dom2 in domain_list2:
                if dom1 in list(domain_interface_residue_network.nodes) and dom2 in list(domain_interface_residue_network.nodes):
                    for res1 in networkx.all_neighbors(domain_interface_residue_network, dom1):
                        if edge1 in interface.keys() and edge2 in interface[edge1].interface_dict.keys():
                            if int(res1) in interface[edge1].interface_dict[edge2]:
                                one = True
                                for res2 in networkx.all_neighbors(domain_interface_residue_network, dom2):
                                    if edge2 in interface.keys() and edge1 in interface[edge2].interface_dict.keys():
                                        if int(res2) in interface[edge2].interface_dict[edge1]:
                                            two = True

                if one and two:
                    gene_domain_res_network.add_edge(dom1, dom2, node1= "domain", node2 ="domain",
                                                     edge_type = "DD")
                    gene_domain_res_network.add_edge(edge1, dom1, node1= "gene", node2 ="domain",
                                                     edge_type = "GD")
                    gene_domain_res_network.add_edge(edge2, dom2, node1= "gene", node2 ="domain",
                                                     edge_type = "GD")

                elif one is False and two:gene_domain_res_network.add_edge(dom1, edge2, node1= "gene", node2 ="domain",
                                                                           edge_type = "GD")
                elif one and two is False: gene_domain_res_network.add_edge(edge1, dom2, node1= "domain", node2 ="gene",
                                                                            edge_type = "GD")
                else: gene_domain_res_network.add_edge(edge1, edge2, node1 = "gene", node2 = "gene", edge_type = "GG")
"""

def networkx_to_node_attributes(network_object, label, group):

    global all_nodes, all_domains

    node_att_df = pandas.DataFrame("none", index=list(network_object.nodes), columns=["Type"])

    for i in network_object.nodes:
        if i in all_nodes: node_att_df.at[i, "Type"] = "Gene"
        elif i in all_domains: node_att_df.at[i, "Type"] = "Domain"
        else: node_att_df.at[i, "Type"] = "Residue"

    name = label + "_" + group if group is not None else label
    node_att_df.to_csv(result_path + name + '_mega_network_node_attributes.csv')
    return 1


def mega_network(label, edges, group):
    global domain_gene_network

    network = networkx.Graph()

    for edge in edges:
        edge1, edge2 = edge.split("-")[0], edge.split("-")[1]
        domain_list1, domain_list2 = [], []
        # Take the domain lists

        if edge1 in list(domain_gene_network.nodes):
            domain_list1 = [i for i in networkx.all_neighbors(domain_gene_network, edge1)]
        if edge2 in list(domain_gene_network.nodes):
            domain_list2 = [i for i in networkx.all_neighbors(domain_gene_network, edge2)]

        # Check the interface residues inside and if they match with the other partners' residues

        if domain_list1 != []:
            for dom1 in domain_list1:
                res_check = False
                mutated_res1 = []
                if dom1 in list(domain_interface_residue_network.nodes):
                    for res1 in networkx.all_neighbors(domain_interface_residue_network, dom1):
                        if edge1 in interface.keys():
                            if edge2 in interface[edge1].interface_dict.keys():
                                if int(res1) in interface[edge1].interface_dict[edge2]:
                                    res_check = True
                                    if edge1 in protein_muts.keys():
                                        if res1 in protein_muts[edge1].mutatedResidues: mutated_res1.append(res1)
                if res_check:
                    network.add_edge(edge1, dom1, node1="gene", node2="domain", edge_type="GD")
                    for mr in mutated_res1: network.add_edge(dom1, edge + "__" + mr, node1="domain",
                                                             node2="Residue", edge_type="DR")

        if domain_list2 != []:
            for dom2 in domain_list2:
                res_check = False
                mutated_res2 = []
                if dom2 in list(domain_interface_residue_network.nodes):
                    for res2 in networkx.all_neighbors(domain_interface_residue_network, dom2):
                        if edge2 in interface.keys():
                            if edge1 in interface[edge2].interface_dict.keys():
                                if int(res2) in interface[edge2].interface_dict[edge1]:
                                    res_check = True
                                    if edge2 in protein_muts.keys():
                                        if res2 in protein_muts[edge2].mutatedResidues: mutated_res2.append(res2)
                if res_check:
                    network.add_edge(edge2, dom2, node1="gene", node2="domain", edge_type="GD")
                    for mr in mutated_res2: network.add_edge(dom2, edge + "__" + mr, node1="domain",
                                                             node2="Residue", edge_type="DR")

        if domain_list1 != [] and domain_list2 != []:
            one, two = False, False
            for dom1 in domain_list1:
                for dom2 in domain_list2:
                    if dom1 in list(domain_interface_residue_network.nodes) and dom2 in list(
                            domain_interface_residue_network.nodes):
                        for res1 in networkx.all_neighbors(domain_interface_residue_network, dom1):
                            if edge1 in interface.keys() and edge2 in interface[edge1].interface_dict.keys():
                                if int(res1) in interface[edge1].interface_dict[edge2]:
                                    one = True
                                    for res2 in networkx.all_neighbors(domain_interface_residue_network, dom2):
                                        if edge2 in interface.keys() and edge1 in interface[edge2].interface_dict.keys():
                                            if int(res2) in interface[edge2].interface_dict[edge1]:
                                                two = True

                    if one and two:
                        network.add_edge(dom1, dom2, node1="domain", node2="domain", edge_type="DD")
                        network.add_edge(edge1, dom1, node1="gene", node2="domain", edge_type="GD")
                        network.add_edge(edge2, dom2, node1="gene", node2="domain", edge_type="GD")

                    elif one is False and two:
                        network.add_edge(dom1, edge2, node1="gene", node2="domain", edge_type="GD")
                        network.add_edge(edge1, dom1, node1="gene", node2="domain", edge_type="GD")

                    elif one and two is False:
                        network.add_edge(edge1, dom2, node1="domain", node2="gene", edge_type="GD")
                        network.add_edge(edge2, dom2, node1="gene", node2="domain", edge_type="GD")

                    else:
                        network.add_edge(edge1, edge2, node1="gene", node2="gene", edge_type="GG")

    networkx.write_edgelist(network, result_path + label + '_mega_network.csv', data=['node1', 'node2', 'edge_type'])

    networkx_to_node_attributes(network_object = network, label= label, group= group)

    return network


mega_net = mega_network(label="whole", edges=all_edges, group = None)


# Patient Networks

patient_mega_nets = {}
for patient, network_object in patient_networks.items():
    pat_mega_net = mega_network(label=patient, edges=[edge[0] + "-" + edge[1] for edge in list(network_object.edges)],
                                group = None)
    patient_mega_nets[patient] = pat_mega_net


def wout_GG_mega_net(mega_net, label, group = "GGless"):

    new_GGless_mega = networkx.Graph()

    for i in mega_net.edges(data=True):
        if i[2]["edge_type"] != "GG":
            new_GGless_mega.add_edge(i[0], i[1])
            attrs = {(i[0], i[1]): i[2]}
            networkx.set_edge_attributes(new_GGless_mega, attrs)

    group_name = "_" + group if group is not None else None

    print(len(new_GGless_mega.nodes))
    networkx.write_edgelist(new_GGless_mega, result_path + label + group_name + '_mega_network.csv',
                            data=['node1', 'node2', 'edge_type'])

    networkx_to_node_attributes(network_object=new_GGless_mega, label = label, group = group)
    return new_GGless_mega

wout_GG_mega_net(mega_net= mega_net, label = "whole")

patient_mega_GGless_nets = {}
for patient, network_object in patient_mega_nets.items():
    print(len(network_object.nodes))
    wout_GG_mega_net(mega_net=network_object, label=patient)