import xml.etree.ElementTree as ET
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# Analyzing KGML 
def parse_kgml(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    entries = {}
    relations = []
    
    # Analyzing entries
    for entry in root.findall('entry'):
        entry_id = entry.get('id')
        entry_type = entry.get('type')
        if entry_type == 'gene':
            entries[entry_id] = entry.get('name')
    
    # Analyzing relations
    for relation in root.findall('relation'):
        entry1 = relation.get('entry1')
        entry2 = relation.get('entry2')
        if entry1 in entries and entry2 in entries:
            relations.append((entry1, entry2))
    
    print('relations:', relations)
    print('entry_id:', entries)
    return entries, relations

def create_network(entries, relations):
    G = nx.DiGraph()

    for entry_id, name in entries.items():
        G.add_node(entry_id, label=name, type='gene')
    
    G.add_edges_from(relations)
    
    return G

def remove_redundant_edges(G):
    for node in G.nodes:
        lengths = nx.single_source_shortest_path_length(G, node)
        for target, length in lengths.items():
            if length > 1 and G.has_edge(node, target):
                G.remove_edge(node, target)

def draw_network(G):
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, seed=42) 
    labels = nx.get_node_attributes(G, 'label')
    nx.draw(G, pos, with_labels=True,  node_size=20, node_color='skyblue', edge_color='gray', font_size=10, font_weight='bold')
    plt.title("Gene Interaction Network")
    plt.show()

def main():
    kgml_file = "D:\\ADLAF\\ADLAF\\ADLAF_main\\KEGG.xml"  
    entries, relations = parse_kgml(kgml_file)
    G = create_network(entries, relations)
    remove_redundant_edges(G)
    draw_network(G)

    #with open('entry.txt', 'w') as f:
    #    for entry_id, name in entries.items():
    #        f.write(f"{entry_id}: {name}\n")
    #print("Entry IDs in entry.txt")
    
    #relations_df = pd.DataFrame(relations, columns=['Entry 1', 'Entry 2'])
    #relations_df.to_excel('relations.xlsx', index=False)
    #print("Relations in relations.xlsx")

if __name__ == "__main__":
    main()


    



