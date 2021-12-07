# -*- coding: utf-8 -*-

from utils.capsnetutils import *
import torchmetrics

                 
class CapsuleNet(pl.LightningModule):
    train_epoch_loss = []
    train_epoch_acc = []
    train_epoch_aucroc = []
    
    val_epoch_loss = []
    val_epoch_acc = []
    val_epoch_aucroc = []
    
    test_predict = []
    test_sample_label = []
    test_decoder = []
    
    test_digitcaps = []
    test_matrix = []
    test_conv = []
    test_primary = []
    test_decoder = []
    def __init__(self, n_class=10, n_genes=2000, loss_func=None, 
                 out_feature = 28,need_to_matrix=True,
                 cov_in=1, cov_out=256, cov_kernel_size=9,
                 pc_input_channels=256,pc_output_caps=32,pc_output_dim=8, pc_kernel_size=9,pc_stride=2,
                 rc_in_dim=8,rc_in_caps=32 * 6 * 6,rc_dim_caps=16,rc_num_routing=3,lr=1e-4):
        super(CapsuleNet, self).__init__()
        self.n_class = n_class
        self.n_genes = n_genes
        self.reconstruction_alpha = 1e-4
        self.myloss = loss_func
        self.lr = lr
        self.pl_accuracy = torchmetrics.Accuracy()
        if need_to_matrix:
            self.to_matrix = sampleToMatrix(n_genes = n_genes,
            out_feature = out_feature)
        
        self.conv = nn.Conv2d(cov_in, cov_out, kernel_size = cov_kernel_size)
        self.relu = nn.ReLU(inplace=True)
        self.primary_caps = PrimaryCaps(input_channels=pc_input_channels, 
                                        output_caps=pc_output_caps,
                                        output_dim=pc_output_dim, 
                                        kernel_size=pc_kernel_size,
                                        stride=pc_stride)
        
        self.digit_caps = DigitCaps(in_dim=rc_in_dim,
                                    in_caps=rc_in_caps,
                                    num_caps=n_class,
                                    dim_caps=rc_dim_caps,
                                    num_routing=rc_num_routing)
        self.decoder = ReconstructionNet(n_dim=rc_dim_caps,
                                         n_classes=n_class,
                                        n_genes=n_genes)
        
    def forward(self,vector_sample):
        
        x = self.to_matrix(vector_sample)
        x = self.relu(self.conv(x))
        x = self.primary_caps(x)
        x = self.digit_caps(x)
        probs = x.pow(2).sum(dim=2).sqrt()    
        reconstruction = self.decoder(x)
        return reconstruction, probs

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        return optimizer
    
    def training_step(self, train_batch, batch_ix):
        vector_sample, sample_label = train_batch
        decoder, probs= self.forward(vector_sample)
        reconstruction_loss = F.mse_loss(decoder, vector_sample)
        margin_loss  = self.myloss(probs, sample_label)
        loss = self.reconstruction_alpha * reconstruction_loss + margin_loss
        acc = self.pl_accuracy(toLabel(probs), toLabel(sample_label))
        
        mylogdict = {'loss': loss, 'log': {'train_loss': loss, 'train_acc': acc}}
        return mylogdict
    
    def validation_step(self, validation_batch, batch_ix):
        vector_sample, sample_label = validation_batch
        decoder, probs= self.forward(vector_sample)
        reconstruction_loss = F.mse_loss(decoder, vector_sample)
        margin_loss  = self.myloss(probs, sample_label)
        loss = self.reconstruction_alpha * reconstruction_loss + margin_loss
        val_acc = self.pl_accuracy(toLabel(probs), toLabel(sample_label))
        self.log_dict({'val_loss': loss, 'val_acc': val_acc})
        mylogdict = {'log': {'val_loss': loss, 'val_acc': val_acc}}
        return mylogdict
    
    def test_step(self, test_batch, batch_ix):
        vector_sample, sample_label = test_batch
        x = self.to_matrix(vector_sample) 
        test_matrix = x.cpu()
        
        x = self.relu(self.conv(x))
        test_conv = x.cpu()
        
        x = self.primary_caps(x)
        test_primary = x.cpu()
        
        x = self.digit_caps(x)
        test_digitcaps = x.cpu()
        
        probs = x.pow(2).sum(dim=2).sqrt()    
        reconstruction = self.decoder(x)
        
        self.test_predict.append(probs.cpu())
        self.test_sample_label.append(sample_label.cpu())
        self.test_decoder.append(reconstruction.cpu())
        self.test_matrix.append(test_matrix)
        self.test_conv.append(test_conv)
        self.test_primary.append(test_primary)
        self.test_digitcaps.append(test_digitcaps)
        return {'test': 'test epoch finish ....'}
    
    def training_epoch_end(self, output):
        train_loss = sum([out['log']['train_loss'].item() for out in output]) / len(output)
        self.train_epoch_loss.append(train_loss)
        
        train_acc = sum([out['log']['train_acc'].item() for out in output]) / len(output)
        self.train_epoch_acc.append(train_acc)
    
    def validation_epoch_end(self, output):
        val_loss = sum([out['log']['val_loss'].item() for out in output]) / len(output)
        self.val_epoch_loss.append(val_loss)
        
        val_acc = sum([out['log']['val_acc'].item() for out in output]) / len(output)
        self.val_epoch_acc.append(val_acc)
        print('mean_val_loss: ', val_loss, '\t', 'mean_val_acc: ', val_acc)
        