import numpy as np
import pandas as pd 

#specify file names and modification types
def GetFileName():
    file_names = ['GlyGlyAla (K)Sites_op.csv','GlyGlyCamCys (K)Sites_op.csv','GlyGlyAsp (K)Sites_op.csv','GlyGlyGlu (K)Sites_op.csv','GlyGlyPhe (K)Sites_op.csv'
                    ,'GlyGlyGly (K)Sites_op.csv','GlyGlyHis (K)Sites_op.csv','GlyGlyLys (K)Sites_op.csv','GlyGlyMet (K)Sites_op.csv','GlyGlyAsn (K)Sites_op.csv'
                    ,'GlyGlyPro (K)Sites_op.csv','GlyGlyGln (K)Sites_op.csv','GlyGlyArg (K)Sites_op.csv','GlyGlySer (K)Sites_op.csv','GlyGlyThr (K)Sites_op.csv'
                    ,'GlyGlyVal (K)Sites_op.csv','GlyGlyTrp (K)Sites_op.csv','GlyGlyXle (K)Sites_op.csv','GlyGlyTyr (K)Sites_op.csv','GlyGlyoxMet (K)Sites_op.csv'
                    ]
    xgg_types = ['AGG','CamCGG','DGG','EGG','FGG','GGG','HGG','KGG','MGG','NGG','PGG','QGG','RGG','SGG','TGG','VGG','WGG','XleGG','YGG','oxMGG']
    return [file_names,xgg_types]

#specify file path and file names
file_names = GetFileName()[0]
xgg_types = GetFileName()[1]
file_path = '..\\GGX\\Sites\\'
file_output = 'Character\\All_combined.csv'

#loop for modification datasets combining
for file_name, type_xgg in zip(file_names,xgg_types):
    file_input = file_path + file_name
    print(file_name,type_xgg)

    if type_xgg == 'AGG':
        dataset = pd.read_table(file_input,sep=',')
        print('AGG')
    else:
        dataset1 = pd.read_table(file_input,sep=',')
        dataset = dataset.append(dataset1)
        print(type_xgg)

dataset.to_csv(file_output,index=False)

