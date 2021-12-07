# =============================================================================
#Capsnet consists of two parts: encoder and decoder. The first three layers are encoders and the last three layers are decoders:
#Layer 1: convolution layer
#Second layer: primarycaps layer
#Layer 3: digitacaps layer
#Layer 4: the first fully connected layer
#Layer 5: the second fully connected layer
#Layer 6: the third fully connected layer
# =============================================================================

from utils.myexptorch import *
#Construction of capsule neural network
#1. Capsule activation function

def squash(x, dim=-1, scale_factor=0.5, eplison=1e-7):
    x_norm = torch.norm(x, p=2, dim=dim, keepdim=True) 
    scale = x_norm**2 / (scale_factor + x_norm**2) 
    v = scale * x / (x_norm + eplison) 
    return v

# 2.sampleToMatrix
class sampleToMatrix(nn.Module):
    # =============================================================================
    # Due to the natural defects of scExpr, or there are not many genes expressed by scExpr themselves, which will lead to the great sparsity of expression matrix
    # For example, the original expression matrix is n * m, in which the value of n is about 30000 ~ 50000 without any filtering.
    # After preliminary threshold screening, the value of n is about 10000. However, when accurate to each cell, the expressed genes may be less than 500 ~ 1000
    # Therefore, if this sparse vector (a single sample is a vector) is directly used for neural network training, it will lead to difficulty or even failure in training
    #
    # Therefore, I developed a feature extraction model, which uses different shallow neural networks to extract the corresponding initialization features,
    # Multi dimensional gene expression is mapped into the expression of low latitude gene clusters
    # =============================================================================
    # =============================================================================
    # A total of 28 feature extraction layers are defined
    # Each layer consists of only one full connection for feature extraction, and the output of each full connection is 28
    # In this way, the long vector sample can be transformed into a 28 * 28 matrix
    # Because the initialization parameters of each layer are different, and there is only one network layer
    # In this way, the feature can be extracted according to the weight of each layer (similar to the principal component)
    #
    # n_ Genes is the total number of genes
    # out_ Featur defaults to 28. In principle, the following reshape {i} layer should be added synchronously if adjustment is required
    # =============================================================================
    
    def __init__(self, n_genes,out_feature):
        self.n_genes = n_genes
        self.out_feature = out_feature
        super(sampleToMatrix,self).__init__()
        self.reshapeLayer =nn.Sequential(
            #nn.BatchNorm1d(n_genes),
            nn.Linear(self.n_genes,  self.out_feature * self.out_feature),nn.ReLU(inplace=True)
            )
    
    def forward(self, sample):
        x = self.reshapeLayer(sample)
        x = x.view(-1, 1, self.out_feature, self.out_feature) # (batch, 通道, 长, 宽)
        return x                                                                                                                                                                             

# 3. PrimaryCaps
class PrimaryCaps(nn.Module):

    
    def __init__(self, input_channels, output_caps, output_dim, kernel_size, stride):
        super(PrimaryCaps, self).__init__()
        self.conv = nn.Conv2d(input_channels, output_caps * output_dim, kernel_size=kernel_size, stride=stride)
        self.input_channels = input_channels
        self.output_caps = output_caps
        self.output_dim = output_dim

    def forward(self, x):
        out = self.conv(x)
        N, C, H, W = out.size()
        out = out.view(N, self.output_caps, self.output_dim, H, W)

        out = out.permute(0, 1, 3, 4, 2).contiguous()
        out = out.view(out.size(0), -1, out.size(4))
        out = squash(out)
        return out

# 4. DigitCaps
class DigitCaps(nn.Module):
    def __init__(self, in_dim, in_caps, num_caps, dim_caps, num_routing):
        super(DigitCaps, self).__init__()
        self.in_dim = in_dim
        self.in_caps = in_caps 
        self.num_caps = num_caps
        self.dim_caps = dim_caps
        self.num_routing = num_routing
        self.device = device
        self.W = nn.Parameter(0.01 * torch.randn(1, num_caps, in_caps, dim_caps, in_dim),
                              requires_grad=True)

    def forward(self, x):
        batch_size = x.size(0)
        x = x.unsqueeze(1).unsqueeze(4)
        u_hat = torch.matmul(self.W, x)
        u_hat = u_hat.squeeze(-1)
        temp_u_hat = u_hat.detach()
        b = torch.zeros(batch_size, self.num_caps, self.in_caps, 1).to(self.device)
        for route_iter in range(self.num_routing - 1):
            c = b.softmax(dim=1)
            s = (c * temp_u_hat).sum(dim=2)
            v = squash(s)
            uv = torch.matmul(temp_u_hat, v.unsqueeze(-1))
            b += uv
            
        c = b.softmax(dim=1)
        s = (c * u_hat).sum(dim=2)
        v = squash(s)
        return v

# 5.ReconstructionNet
class ReconstructionNet(nn.Module):
    def __init__(self, n_dim, n_classes, n_genes):
        super(ReconstructionNet, self).__init__()
        self.fc1 = nn.Linear(n_dim * n_classes, 512)
        self.fc2 = nn.Linear(512, 1024)
        self.fc3 = nn.Linear(1024, n_genes)
        self.n_dim = n_dim
        self.n_genes = n_genes
        self.n_classes = n_classes
    def forward(self, x):
        x = x.view(-1,  self.n_dim * self.n_classes)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x)) 
        return x


class CapsLoss(nn.Module):
    def __init__(self, m_pos=0.9, m_neg=0.1, lamb=0.5):
        super(CapsLoss, self).__init__()
        self.m_pos = m_pos
        self.m_neg = m_neg
        self.lamb = lamb
    def forward(self, vector, targets, size_average=False):
        losses = targets.float() * F.relu(self.m_pos - vector).pow(2) + \
                 self.lamb * (1. - targets.float()) * F.relu(vector - self.m_neg).pow(2)
        return losses.mean() if size_average else losses.sum()