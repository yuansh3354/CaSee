# -*- coding: utf-8 -*-
import os
import math
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib import cm

import torch
import torchvision
import torch.utils.data as Data
import torch.nn as nn
import torch.nn.functional as F 
import pytorch_lightning as pl
from sklearn.preprocessing import StandardScaler,MinMaxScaler
from sklearn.model_selection import train_test_split

import warnings
from utils.myexptorch import *

warnings.filterwarnings("ignore")

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') 
def deviceChange(x):
    if (torch.cuda.is_available() & x == 'cuda'):
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    return device


def expToTorch(x):
    return torch.from_numpy(np.array(x)).float()


def labelTotorch(y):
    return torch.LongTensor(y)


def makeDataiter(x,y,batch_size,shuffle=True):
    return Data.DataLoader(Data.TensorDataset(x, y), batch_size, shuffle=shuffle)

def DataFrame_normalize(my_data,std_method):
    if std_method == 'minmax':
        method = MinMaxScaler()
    if std_method == 'std':
        method = StandardScaler()
    my_data = pd.DataFrame(method.fit(my_data).transform(
        my_data),columns=my_data.columns,index=my_data.index)
    return my_data



def toOneHot(ylabels,n_class):
    onehot = torch.zeros(ylabels.shape[0],n_class)
    index = torch.LongTensor(ylabels).view(-1,1)
    onehot.scatter_(dim=1, index=index, value=1)
    return onehot

def toLabel(ylabels):
    return torch.topk(ylabels, 1)[1].squeeze(1)


def show_tensor_pictures(x, y):
    plt.figure(figsize=(12,12))
    fig, axes = plt.subplots(1, len(x))
    for ax, image, label in zip(axes, x, y):
        ax.imshow(image.view((28, 28)).numpy())
        ax.set_title(label)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
    plt.show()