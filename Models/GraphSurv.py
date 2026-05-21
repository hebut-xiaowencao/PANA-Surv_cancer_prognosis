import torch
from Models.appGNN import GNN
from Models.CoxNN import Coxnn


class GraphSurv(torch.nn.Module):
    def __init__(self, input_dim, n_layer: int = 2):
        super(GraphSurv, self).__init__()
        torch.manual_seed(1234)
        torch.set_printoptions(profile='full')

        self.alpha_raw = torch.nn.Parameter(torch.tensor(0.0))
        self.GNN = GNN(gnn='gcn', n_layer=int(n_layer), feature_len=2, dim=input_dim)
        self.CoxNN = Coxnn(input_dim=input_dim)

    def get_alpha(self):
        return torch.sigmoid(self.alpha_raw)

    def forward(self, x, edge_index):
        GNN_out = self.GNN(x, edge_index)
        risk = self.CoxNN(GNN_out)
        return risk, GNN_out

    def forward_fusion(self, x1, x2, edge_index):
        h1 = self.GNN(x1, edge_index)
        h2 = self.GNN(x2, edge_index)
        alpha = self.get_alpha()
        h_fuse = alpha * h1 + (1 - alpha) * h2
        risk = self.CoxNN(h_fuse)
        return risk, h_fuse
