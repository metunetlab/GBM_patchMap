"""
                                    ---------------------
                                     P A T C H  -  M A P
                                    ---------------------
                        GBM PROJECT - STRUCTURAL MAPPING OF THE NETWORKS
                                   Data Curation & Objects

                                        Cansu Dincer
                                    Nurcan Tuncbag, PhD
"""
########################################################################################################################
########################################################################################################################
# IMPORTING #

import pandas, pickle, numpy, re, networkx
from paths import data_path, raw_path, curated_path, object_path, network_path, prism_path, result_path


########################################################################################################################
########################################################################################################################
# NETWORKS #

# Take patient network data and create network objects

"""
network_dict = {}
removed_nodes = ["icrogid:4115704", "icrogid:2178303", "icrogid:5343788", "DUMMY"]
for file in os.listdir(network_path):
    if file.split("_")[-1] == "multipledouble.txt":
        ntw = pandas.read_csv(network_path + file, sep = "\t")
        ntw.columns = ["Node1", "Node2", "Score"]
        G = networkx.Graph()
        for ind, row in ntw.iterrows():
            if row.Node1 not in removed_nodes and row.Node2 not in removed_nodes:
                G.add_edge(row.Node1, row.Node2, weight= row.Score)
        network_dict[file.split("_")[0]] = G

pickle.dump(network_dict, open(result_path + "ALL_PATIENT_NETWORKS.p", "wb"))
"""

def read_patient_networks():
    network_dict = pickle.load(open(result_path + "ALL_PATIENT_NETWORKS.p", "rb"))
    return network_dict


# Initial network information - # edge and nodes

"""
network_info_df = pandas.DataFrame(0, index = network_dict.keys(),
                                   columns=["Node_Number", "Edge_Number"])
network_info_df["Node_Number"] = network_info_df.apply(lambda x: len(network_dict[x.name].nodes), axis=1)
network_info_df["Edge_Number"] = network_info_df.apply(lambda x: len(network_dict[x.name].edges), axis=1)

network_info_df.to_csv(result_path + "first_network_info_df.csv")
"""


########################################################################################################################
########################################################################################################################
# PROTEOME #

class Proteome(object):
    def __init__(self, uniprot, status, gene, synonyms, sequence):
        self.uniprot = uniprot
        self.status = status
        self.gene = gene.split("; ")
        self.sequence, self.seqLength = sequence, len(sequence)
        self.synonyms = synonyms

"""
proteome_dict = {}
proteome = pandas.read_csv(raw_path + "proteome.csv", index_col=0)
for ind, row in proteome.iterrows():
    if type(row.synonyms) != float:
        synonym_list = row.synonyms.split(" ")
    else: synonym_list = None
    if type(row.primary_gene) != float:
        obj = Proteome(uniprot=ind, status=row.status, gene=row.primary_gene,
                       synonyms=synonym_list,
                       sequence=row.sequence)
        proteome_dict[ind] = obj


pickle.dump(proteome_dict, open(object_path + "proteome.p", "wb"))
"""

def human_proteome():
    proteome_dict = pickle.load(open(object_path + "proteome.p", "rb"))
    return proteome_dict


########################################################################################################################
########################################################################################################################
# NODE & EDGES #

"""
# All nodes
all_nodes = []
for patient, ntw in network_dict.items():
    for node in ntw.nodes:
        if node not in all_nodes: all_nodes.append(node)


# Dictionary for Hugo Symbols of the genes to Uniprot IDs

gene_unip = {}
for node in all_nodes:
    for gene, obj in proteome_dict.items():
        if obj.status == "reviewed":
            if node in obj.gene: gene_unip[node] = obj.uniprot
            else:
                if obj.synonyms != None and node in obj.synonyms:
                    gene_unip[node] = obj.uniprot

# All uniprot nodes
all_node_uniprot_dict = dict()
for patient, ntw in network_dict.items():
    for node in ntw.nodes:
        if node not in all_node_uniprot_dict.keys(): all_node_uniprot_dict[node] = gene_unip[node]


# All edges & All uniprot edges
all_edges = []
all_edge_uniprot_dict = dict()
for patient, ntw in network_dict.items():
    for edge in ntw.edges:
        if edge[0]+"-"+edge[1] not in all_edges and edge[1]+"-"+edge[0] not in all_edges:
            all_edges.append(edge[0]+"-"+edge[1])
        unip1, unip2 = gene_unip[edge[0]], gene_unip[edge[1]]
        if edge[0] + "-" + edge[1] not in all_edge_uniprot_dict.keys() and edge[1] + "-" + edge[0] not in all_edge_uniprot_dict.keys():
            all_edge_uniprot_dict[edge[0] + "-" + edge[1]] = unip1+"-"+unip2

pickle.dump(all_nodes, open(curated_path + "All_Nodes.p", "wb"))
pickle.dump(all_edges, open(curated_path + "All_Edges.p", "wb"))
pickle.dump(gene_unip, open(curated_path + "Gene_to_UniprotID.p", "wb"))
pickle.dump(all_node_uniprot_dict, open(curated_path + "All_Node_UniprotID_dict.p", "wb"))
pickle.dump(all_edge_uniprot_dict, open(curated_path + "All_Edge_UniprotID_dict.p", "wb"))
"""

all_nodes = pickle.load(open(curated_path + "All_Nodes.p", "rb"))
all_edges = pickle.load(open(curated_path + "All_Edges.p", "rb"))
gene_unip = pickle.load(open(curated_path + "Gene_to_UniprotID.p", "rb"))
all_node_uniprot_dict = pickle.load(open(curated_path + "All_Node_UniprotID_dict.p", "rb"))
all_edge_uniprot_dict =pickle.load(open(curated_path + "All_Edge_UniprotID_dict.p", "rb"))

########################################################################################################################
########################################################################################################################
# STRUCTURES #

##################################################################################
# Interactome3D #

# Interactome3D has both PDB for experimental structures and  ModBase
# for computationally predictive models

class Int3D_Protein(object):
    def __init__(self, uniprot, gene):
        self.gene = gene
        self.uniprot = uniprot
        self.threeD = dict()
        self.pdb = list()
        self.coverage = dict()
        self.coverageFile = dict()
        self.coveredResidues = list()
    def addStructure(self, type, pdb, chain, mpqs, zdope, ga431, rank1, rank2):
        self.pdb.append(pdb+ "_" + chain)
        self.threeD[pdb+ "_" + chain] = {"type": type, "rank": (rank1,rank2)}
        self.threeD[pdb + "_" + chain]["reliable"] = True if mpqs >= 1.1 and ga431 >= 0.7 and zdope < 0 else False if type == "model" else True
    def addCoverage(self, rank1, rank2, seqbegin, seqend, file, source, reliable = True):
        for pdb, property in self.threeD.items():
            if (rank1, rank2) == property["rank"] and property["reliable"] == True:
                if rank1 not in self.coverage:
                    self.coverage[rank1] = [(seqbegin, seqend, rank2)]
                else:
                    self.coverage[rank1].append((seqbegin, seqend, rank2))
                if rank1 not in self.coverageFile.keys():
                    self.coverageFile[rank1] = [(seqbegin, seqend, rank2, source, file)]
                else:
                    self.coverageFile[rank1].append((seqbegin, seqend, rank2, source, file))
                for res in list(range(seqbegin, seqend+1)):
                    if res not in self.coveredResidues: self.coveredResidues.append(res)

"""
int3D_proteins_df = pandas.read_csv(raw_path + "proteins.dat", sep = "\t", index_col = 0)
int3d_proteins = {}
for gene, uniprot_ind in all_node_uniprot_dict.items():
    df = int3D_proteins_df[int3D_proteins_df.index == uniprot_ind]
    obj = Int3D_Protein(uniprot_ind, gene)
    for ind, row in df.iterrows():
        obj.addStructure(type=row.TYPE, pdb=row.PDB_ID, chain=row.CHAIN,
                         mpqs=row.MPQS, zdope=row.ZDOPE, ga431=row.GA431,
                         rank1=row.RANK_MAJOR, rank2=row.RANK_MINOR)
        obj.addCoverage(seqbegin=row.SEQ_BEGIN, seqend=row.SEQ_END,
                        file=row.FILENAME, rank1=row.RANK_MAJOR, rank2=row.RANK_MINOR,
                        source = "PDB" if "EXP" in row.FILENAME.split("-") else "MODBASE")
    int3d_proteins[gene] = obj

pickle.dump(int3d_proteins, open(object_path + "int3d_proteins.p", "wb"))
"""

def interactome_3d_protein():
    int3d_proteins = pickle.load(open(object_path + "int3d_proteins.p", "rb"))
    return int3d_proteins

##################################################################################
# Interactome Insider #

# Interactome Insider gas both PDB for experimental interface residues,
# and I3D & ECLAIR interface residue indices (not structures itself)

class Interfaces_II(object):
    def __init__(self, uniprot1, uniprot2, interaction):
        self.partners = interaction
        self.partner1, self.partner2 = uniprot1, uniprot2
        self.bs1, self.bs2 = list(), list()
        self.interfaces = dict()

    def addInterfaces(self, res1, res2, source):
        exp = r"\d+-\d"
        if res1 != "[]":
            res1 = res1.replace("[","")
            res1 = res1.replace("]","")
            res1_list = res1.split(",")
            for res in res1_list:
                if re.match(exp, res):
                    start, end = res.split("-")
                    self.bs1 += range(int(start), int(end)+1)
                else:
                    self.bs1.append(int(res))
        if res2 != "[]":
            res2 = res2.replace("[","")
            res2 = res2.replace("]","")
            res2_list = res2.split(",")
            for res in res2_list:
                if re.match(exp, res):
                    start, end = res.split("-")
                    self.bs2 += range(int(start), int(end)+1)
                else:
                    self.bs2.append(int(res))
        if source not in self.interfaces.keys(): self.interfaces[source] = [(self.bs1, self.bs2)]
        else: self.interfaces[source].append((self.bs1, self.bs2))

"""
interactome_insider = pandas.read_csv(raw_path + "Interactome_Insider_Interfaces.csv", index_col=0)
ii_interactions = {}
for gene_partner, uniprot_partner in all_edge_uniprot_dict.items():
    dfs = []
    if uniprot_partner in interactome_insider.index:
        df1 = interactome_insider[interactome_insider.index == uniprot_partner]
        dfs.append(df1)
    if uniprot_partner.split("-")[1] + "-" + uniprot_partner.split("-")[0] in interactome_insider.index:
        df2 = interactome_insider[interactome_insider.index == uniprot_partner.split("-")[1] + "-" + uniprot_partner.split("-")[0]]
        dfs.append(df2)

    if dfs != list():
        df = pandas.concat(dfs)

        obj = Interfaces_II(uniprot1=uniprot_partner.split("-")[0],
                            uniprot2=uniprot_partner.split("-")[1], interaction=uniprot_partner)

        for interaction, row in df.iterrows():
            if interaction == uniprot_partner:
                obj.addInterfaces(res1=row.P1_IRES, res2=row.P2_IRES, source=row.Source)
            elif interaction.split("-")[1]+"-"+interaction.split("-")[0] == uniprot_partner:
                obj.addInterfaces(res1=row.P2_IRES, res2=row.P1_IRES, source=row.Source)

            if gene_partner not in ii_interactions.keys() and \
                    gene_partner.split("-")[1]+"-"+gene_partner.split("-")[0] not in ii_interactions.keys():
                ii_interactions[gene_partner] = obj

            else: print(gene_partner)

pickle.dump(ii_interactions, open(object_path + "ii_interactions.p", "wb"))
"""

def interactome_insider_interactions():
    ii_interactions = pickle.load(open(object_path + "ii_interactions.p", "rb"))
    return ii_interactions

##################################################################################
# PRISM #

class Interfaces_PRISM(object):
    def __init__(self, uniprot1, uniprot2):
        self.partners = uniprot1 + "-" + uniprot2
        self.partner1, self.partner2 = uniprot1, uniprot2
        self.bs1, self.bs2 = list(), list()
        self.interfaces = dict()

    def addInterfaces(self, int_df):
        if set(int_df[1]) != set(): self.bs1 = list(set(int_df[1]))
        if set(int_df[3]) != set(): self.bs2 = list(set(int_df[3]))
        source = "PRISM"
        if source not in self.interfaces.keys(): self.interfaces[source] = [(self.bs1, self.bs2)]
        else: self.interfaces[source].append((self.bs1, self.bs2))

"""
prism_interaction_files = [file.split(".")[0] for file in os.listdir(prism_path)]
prism_interactions = {}

for gene_partner, uniprot_partner in all_edge_uniprot_dict.items():

    if uniprot_partner in prism_interaction_files:
        df = pandas.read_csv(prism_path+ uniprot_partner + ".txt",
                             sep = "\t", header = None)
        print(df)
        u1, u2 = list(set(df[0]))[0], list(set(df[2]))[0]
        obj = Interfaces_PRISM(uniprot1=u1, uniprot2=u2)
        obj.addInterfaces(int_df = df)
        if u1+"-"+u2 == uniprot_partner: prism_interactions[gene_partner] = obj

    elif uniprot_partner.split("-")[1] + "-" + uniprot_partner.split("-")[0] in prism_interaction_files:
        df = pandas.read_csv(prism_path + uniprot_partner.split("-")[1] + "-" +
                             uniprot_partner.split("-")[0] + ".txt", sep="\t", header=None)

        u1, u2 = list(set(df[0]))[0], list(set(df[2]))[0]
        obj = Interfaces_PRISM(uniprot1=u1, uniprot2=u2)
        obj.addInterfaces(int_df=df)
        if u2+"-"+u1 == uniprot_partner:
            prism_interactions[gene_partner.split("-")[1]+"-"+gene_partner.split("-")[0]] = obj

pickle.dump(prism_interactions, open(object_path + "prism_interactions.p", "wb"))
"""

def PRISM_interactions():
    prism_interactions = pickle.load(open(object_path + "prism_interactions.p", "rb"))
    return prism_interactions


########################################################################################################################
########################################################################################################################
# INTERFACES #

def checkPartner(uniprot, partner, obj):

    if uniprot == partner.split("-")[0]:
        return {"Partner": obj.partner2, "bs": obj.bs1, "Source": obj.interfaces}
    elif uniprot == partner.split("-")[1]:
        return {"Partner": obj.partner1, "bs": obj.bs2, "Source": obj.interfaces}
    else:
        return False

class Protein_Interfaces(object):
    def __init__(self, uniprot, gene):
        self.uniprot, self.gene = uniprot, gene
        self.interface_dict = dict()
        self.interface_source = dict()
        self.all_interface_pos = list()

    def interfaceInfo(self, ii, prism, proteome_dict):

        global checkPartner

        for partner, obj in ii.items():
            if checkPartner(self.uniprot, obj.partners, obj):
                d = checkPartner(self.uniprot, obj.partners, obj)
                genes = [g for g in proteome_dict[d["Partner"]].gene if g in all_nodes]
                if genes != []:
                    gene_partner=genes[0]

                    if d["bs"] != []:
                        for res in d["bs"]:
                            if self.all_interface_pos == list():
                                self.all_interface_pos.append(res)
                            else:
                                if res not in self.all_interface_pos: self.all_interface_pos.append(res)
                    if gene_partner not in self.interface_dict.keys():
                        self.interface_dict[gene_partner] = d["bs"]
                    else:
                        for res in d["bs"]:
                            if res not in self.interface_dict[gene_partner]:
                                self.interface_dict[gene_partner].append(res)
                    if gene_partner not in self.interface_source.keys():
                        self.interface_source[gene_partner] = [d["Source"]]
                    else:
                        if d["Source"] not in self.interface_source[gene_partner]:
                            self.interface_source[gene_partner].append(d["Source"])


        for partner, obj in prism.items():
            if checkPartner(self.uniprot, obj.partners, obj):
                d = checkPartner(self.uniprot, obj.partners, obj)
                genes = [g for g in proteome_dict[d["Partner"]].gene if g in all_nodes]
                if genes != []:
                    gene_partner=genes[0]
                    if d["bs"] != []:
                        for res in d["bs"]:
                            if self.all_interface_pos == list():
                                self.all_interface_pos.append(res)
                            else:
                                if res not in self.all_interface_pos: self.all_interface_pos.append(res)
                    if gene_partner not in self.interface_dict.keys():
                        self.interface_dict[gene_partner] = d["bs"]
                    else:
                        for res in d["bs"]:
                            if res not in self.interface_dict[gene_partner]:
                                self.interface_dict[gene_partner].append(res)
                    if gene_partner not in self.interface_source.keys():
                        self.interface_source[gene_partner] = [d["Source"]]
                    else:
                        if d["Source"] not in self.interface_source[gene_partner]:
                            self.interface_source[gene_partner].append(d["Source"])

"""
protein_interfaces = {}
for gene, uniprotID in all_node_uniprot_dict.items():
    obj = Protein_Interfaces(uniprot=uniprotID, gene= gene)
    obj.interfaceInfo(ii = interactome_insider_interactions(), prism = PRISM_interactions(),
                      proteome_dict=human_proteome())
    protein_interfaces[gene] = obj

pickle.dump(protein_interfaces, open(object_path + "protein_interfaces.p", "wb"))

"""

def protein_interface():
    protein_interface = pickle.load(open(object_path + "protein_interfaces.p", "rb"))
    return protein_interface


########################################################################################################################
########################################################################################################################
# DOMAIN #

class Domain(object):
    def __init__(self, gene, uniprot):
        self.uniprot, self.gene = uniprot, gene
        self.domain_list = list()
        self.all_domain_res = list()
        self.domain_seq = dict()
        self.structured_domains = dict()
        self.interface_domains = dict()

    def domainInfo(self, seqstrat, seqend, domain_name):

        if domain_name not in self.domain_list: self.domain_list.append(domain_name)
        for i in range(seqstrat, seqend+1):
            if i not in self.all_domain_res: self.all_domain_res.append(i)
        if domain_name not in self.domain_seq.keys():
            self.domain_seq[domain_name] = list(range(seqstrat, seqend+1))
        else:
            for i in range(seqstrat, seqend + 1):
                if i not in self.domain_seq[domain_name]:
                    self.domain_seq[domain_name].append(i)

    def structureInfo(self, structure_object, interface_object, seqstrat, seqend, domain_name):

        if self.gene in structure_object.keys():
            for i in range(seqstrat, seqend + 1):
                if i in structure_object[self.gene].coveredResidues:
                    if domain_name not in self.structured_domains.keys():
                        self.structured_domains[domain_name] = [i]
                    else: self.structured_domains[domain_name].append(i)

        if self.gene in interface_object.keys():
            for i in range(seqstrat, seqend + 1):
                if i in interface_object[self.gene].all_interface_pos:
                    if domain_name not in self.interface_domains.keys():
                        self.interface_domains[domain_name] = [i]
                    else: self.interface_domains[domain_name].append(i)

"""
# PFAM data

pfam_df = pandas.read_csv(raw_path + "PFAM_9606_v.33_1.tsv", sep = "\t", index_col=0)
pfam_df = pfam_df[pfam_df.type == "Domain"]

domain_dict = {}
for gene, uniprotID in all_node_uniprot_dict.items():
    print(uniprotID)
    if uniprotID in pfam_df.index:
        df = pfam_df[pfam_df.index == uniprotID]
        obj = Domain(uniprot=uniprotID, gene= gene)
        for ind, row in df.iterrows():
            obj.domainInfo(seqstrat = row["alignment start"], seqend = row["alignment end"],
                           domain_name = row["hmm name"])
            obj.structureInfo(seqstrat = row["alignment start"], seqend = row["alignment end"],
                              domain_name = row["hmm name"],structure_object = interactome_3d_protein(),
                              interface_object= protein_interface())
        domain_dict[gene] = obj

pickle.dump(domain_dict, open(object_path + "pfam_domains.p", "wb"))
"""

def domain():
    domain_dict = pickle.load(open(object_path + "pfam_domains.p", "rb"))
    return domain_dict


########################################################################################################################
########################################################################################################################
# MUTATIONS #

class Protein_Mutations(object):
    def __init__(self, gene, uniprot):
        self.gene, self.uniprot = gene, uniprot
        self.mutations = list()
        self.mutatedResidues = list()
        self.interface_mutations = list()
        self.structured_mutations = list()
        self.domain_mutations = list()

    def add_mutationInfo(self, mutation, position):

        if mutation not in self.mutations: self.mutations.append(mutation)
        if position not in self.mutatedResidues: self.mutatedResidues.append(position)

    def add_structureInfo(self, mutation, position, structure_object, interface_object):

        if mutation.split(" p.")[0] in structure_object.keys():
            if int(position) in structure_object[mutation.split(" p.")[0]].coveredResidues:
                if mutation not in self.structured_mutations:
                    self.structured_mutations.append(mutation)

        if mutation.split(" p.")[0] in interface_object.keys():
            if int(position) in interface_object[mutation.split(" p.")[0]].all_interface_pos:
                if mutation not in self.interface_mutations:
                    self.interface_mutations.append(mutation)

    def add_domainInfo(self, mutation, position, domain_object):

        if mutation.split(" p.")[0] in domain_object.keys():
            if int(position) in domain_object[mutation.split(" p.")[0]].all_domain_res:
                if mutation not in self.domain_mutations:
                    self.domain_mutations.append(mutation)

class Mutations(object):
    def __init__(self, mutation, isDel, isTCGA, isCOSMIC):
        self.mutation = mutation
        self.gene = None
        self.isDeleterious, self.TCGAhotspot, self.COSMIChotspot = isDel, isTCGA, isCOSMIC
        self.DepMaps = list()
        self.onInterface = False
        self.onDomain = False
        self.haveStructure = False
        self.onPatch = None

    def add_geneInfo(self, proteome):
        if self.mutation.split(" p.")[0] in proteome.keys():
            self.gene = self.mutation.split(" p.")[0]

    def add_cellLines(self, depMap):
        if depMap not in self.DepMaps: self.DepMaps.append(depMap)

    def add_structureInfo(self, position, structure_object, interface_object):

        if self.gene in structure_object.keys():
            if int(position) in structure_object[self.gene].coveredResidues:
                self.haveStructure = True

        if self.gene in interface_object.keys():
            if int(position) in interface_object[self.gene].all_interface_pos:
                self.onInterface = True

    def add_domainInfo(self, position, domain_object):

        if self.gene in domain_object.keys():
            if int(position) in interface_object[self.gene].all_domain_res:
                self.onDomain = True

    def add_patchInfo(self, position, patchResidues):
        if position in patchResidues: self.onPatch = True

"""
cl_info = pandas.read_csv(raw_path + "sample_info_20Q2.csv", index_col=0)
gbm_cl_info = cl_info[cl_info.Subtype == "Glioblastoma"]
gbm_cls = set(gbm_cl_info.index)

mutation_raw = pandas.read_csv(raw_path + "CCLE_mutations_20Q2.csv", sep = "\t", index_col=0)
gbm_mutation = mutation_raw[mutation_raw.DepMap_ID.isin(list(gbm_cls))]
gbm_mutation = gbm_mutation[~pandas.isna(gbm_mutation.Protein_Change)]
gbm_mutation = gbm_mutation[gbm_mutation.Variant_Classification.isin(['Nonsense_Mutation', 'Missense_Mutation'])]
gbm_mutation["Check_mutation"] = gbm_mutation.apply(
    lambda x: True if len(x.Protein_Change.split("p.")[1].split(">")) == 1 or
                       len(x.Protein_Change.split("p.")[1].split("_")) == 1 else False, axis =1)
gbm_mutation = gbm_mutation[gbm_mutation.Check_mutation]
gbm_mutation = gbm_mutation[["Protein_Change", "isDeleterious", "isTCGAhotspot", "isCOSMIChotspot", "DepMap_ID"]]
gbm_mutation["Mutation"] = gbm_mutation.apply(lambda x: str(x.name) + " " + str(x.Protein_Change), axis=1)
gbm_mutation["Position"] = gbm_mutation.apply(lambda x: x.Protein_Change.split(".")[1][1:-1], axis=1)
gbm_mutation = gbm_mutation[gbm_mutation.index.isin(all_nodes)]
gbm_mutation = gbm_mutation.drop_duplicates()
gbm_mutation.to_csv(curated_path + "Mutations_20Q2_CCLE.csv")

protein_mutations_dict = {}
for gene, uniprotID in all_node_uniprot_dict.items():
    if gene in gbm_mutation.index:
        df = gbm_mutation[gbm_mutation.index == gene]
        obj = Protein_Mutations(uniprot=uniprotID, gene= gene)
        for ind, row in df.iterrows():
            obj.add_mutationInfo(mutation=row.Mutation, position = row.Position)
            obj.add_structureInfo(mutation=row.Mutation, position=row.Position,
                                  structure_object = interactome_3d_protein(),
                                  interface_object= protein_interface())
            obj.add_domainInfo(mutation =row.Mutation, position = row.Position,
                               domain_object = domain())
        protein_mutations_dict[gene] = obj

pickle.dump(protein_mutations_dict, open(object_path + "protein_mutations.p", "wb"))

mutations_dict = {}
for mutation in set(gbm_mutation.Mutation):
    df = gbm_mutation[gbm_mutation.Mutation == mutation]
    obj = Mutations(mutation= mutation, isDel= df.isDeleterious.values[0],
                    isTCGA = df.isTCGAhotspot.values[0], isCOSMIC= df.isCOSMIChotspot.values[0])
    obj.add_geneInfo(proteome = human_proteome())
    for depmap in df.DepMap_ID.values: obj.add_cellLines(depMap=depmap)
    obj.add_structureInfo(position = df.Position.values[0],
                          structure_object = interactome_3d_protein(),
                          interface_object = protein_interface())
    obj.add_domainInfo(position = df.Position.values[0], domain_object = domain())

    mutations_dict[gene] = obj

pickle.dump(mutations_dict, open(object_path + "gbm_mutations.p", "wb"))
"""

def protein_mutations():
    protein_mutations_dict = pickle.load(open(object_path + "protein_mutations.p", "rb"))
    return protein_mutations_dict


def mutations():
    mutations_dict = pickle.load(open(object_path + "gbm_mutations.p", "rb"))
    return mutations_dict


########################################################################################################################
########################################################################################################################
