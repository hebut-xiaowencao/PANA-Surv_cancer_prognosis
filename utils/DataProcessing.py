import os
import numpy as np
from typing import Union, List, Tuple, Optional, Callable
from scipy.stats import boxcox
import torch
from torch_geometric.data import InMemoryDataset, Data
from config import RELATIONSHIPS_PATH

def get_omics_data(features):
    duration = np.array(features.iloc[1, :].values.tolist()[2:])
    event = np.array(features.iloc[0, :].values.tolist()[2:])

    gene_exp = features.loc[features['Platform'] == 'geneExp']
    gene_exp = gene_exp.drop(columns=['Platform'], axis=1)

    meth = features.loc[features['Platform'] == 'methylation']
    meth = meth.drop(columns=['Platform'], axis=1)

    omic_data = pd.merge(gene_exp, meth, on=['GeneSymbol'], how='inner')
    return omic_data, duration, event

def get_edge(omic_data, gene_relationship):

    gene_name = omic_data['GeneSymbol']
    gene_relationship.columns = ['gene_x', 'gene_y']

    gene_idx_dic = {'gene_name': gene_name, 'node_idx': list(np.arange(len(gene_name)))}
    gene_idx = pd.DataFrame(gene_idx_dic)

    # Merge gene relationship data with gene indexes to generate connections between nodes
    tmp1 = gene_idx.rename(columns={'gene_name': 'gene_x'})
    tmp_nodex = pd.merge(gene_relationship, tmp1, on='gene_x').drop_duplicates().reset_index(drop=True)

    tmp2 = gene_idx.rename(columns={'gene_name': 'gene_y'})
    adj_df = pd.merge(tmp_nodex, tmp2, on='gene_y').drop_duplicates().reset_index(drop=True)

    adj_df = adj_df[['node_idx_x', 'node_idx_y']]

    return adj_df

def read_data(raw_dir, raw_file_names):
    import pandas as pd
  
    #edge
    feature_path = os.path.join(raw_dir, raw_file_names[0])
    
    edge_path = RELATIONSHIPS_PATH
    if not edge_path.exists():
        raise FileNotFoundError(f"Gene relationship file not found: {edge_path}")
  
    dataset_original = pd.read_excel(feature_path)
    gene_relationship = pd.read_excel(edge_path)
    omic_data, time, status = get_omics_data(dataset_original)
    
    features = np.transpose(omic_data.drop(columns=['GeneSymbol'], axis=1).values)

    for i in range(len(time)):
        if status[i] == 0:
            time[i] = -time[i]
    
    data_time = time.reshape(-1, 1)

    samples_num = features.shape[0] // 2

    if data_time.shape[0] == status.shape[0] == samples_num:
        print("There are %d samples" % (samples_num))
    
    #Splitting histology data into gene expression methylation data
    gene_exp,  meth = features[0:samples_num], features[samples_num:]
    edges = get_edge(omic_data, gene_relationship)

    return gene_exp, meth, data_time, status, edges

class CancerDataset(InMemoryDataset):

    def __init__(self, root: Optional[str] = None, transform: Optional[Callable] = None,
                 pre_transform: Optional[Callable] = None, pre_filter: Optional[Callable] = None):
        super().__init__(root, transform, pre_transform, pre_filter)
        self.data, self.slices = torch.load(self.processed_paths[0],weights_only=False)

    @property
    def raw_dir(self) -> str:
        return os.path.join(self.root, 'raw')
    
    @property
    def processed_dir(self) -> str:
        return os.path.join(self.root, 'processed')

    @property
    def raw_file_names(self) -> Union[str, List[str], Tuple]:
        names = os.listdir(self.raw_dir)
        return names

    @property
    def processed_file_names(self) -> Union[str, List[str], Tuple]:
        names = ['{}.pt'.format(os.path.basename(self.root))]
        return names

    def download(self):
        pass

    def save_clinic(self, duration, event):
        torch.save(obj={'duration': duration, 'event': event},
                   f=os.path.join(self.processed_dir, '{}_info.pkl'.format(os.path.basename(self.root))))

    def process(self):
        data_list = []
      
        gene_exp, meth, duration, event, edges = read_data(raw_dir=self.raw_dir,
                                                                     raw_file_names=self.raw_file_names)
        edge_index = torch.stack([torch.tensor(edges['node_idx_x'].to_list(), dtype=torch.long),
                                  torch.tensor(edges['node_idx_y'].to_list(), dtype=torch.long)])

        print('Converting samples to PyG graphs')

        for index, gene in enumerate(gene_exp):

            feature = np.vstack((gene, meth[index]))
            feature = torch.tensor(data=feature, dtype=torch.float32).t()
            data_list.append(Data(x=feature, edge_index=edge_index))


        data, slices = self.collate(data_list=data_list)
        torch.save((data, slices), self.processed_paths[0])

        self.save_clinic(duration=duration, event=event)
