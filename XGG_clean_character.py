import numpy as np
import pandas as pd 
from pymsfilereader import MSFileReader

#specify files and modification types
def GetFileName():
    file_names = ['GlyGlyAla (K)Sites.txt','GlyGlyCamCys (K)Sites.txt','GlyGlyAsp (K)Sites.txt','GlyGlyGlu (K)Sites.txt','GlyGlyPhe (K)Sites.txt'
                    ,'GlyGlyGly (K)Sites.txt','GlyGlyHis (K)Sites.txt','GlyGlyLys (K)Sites.txt','GlyGlyMet (K)Sites.txt','GlyGlyAsn (K)Sites.txt'
                    ,'GlyGlyPro (K)Sites.txt','GlyGlyGln (K)Sites.txt','GlyGlyArg (K)Sites.txt','GlyGlySer (K)Sites.txt','GlyGlyThr (K)Sites.txt'
                    ,'GlyGlyVal (K)Sites.txt','GlyGlyTrp (K)Sites.txt','GlyGlyXle (K)Sites.txt','GlyGlyTyr (K)Sites.txt','GlyGlyoxMet (K)Sites.txt'
                    ,'GlyGly (K)Sites.txt']
    xgg_types = ['AGG','CamCGG','DGG','EGG','FGG','GGG','HGG','KGG','MGG','NGG','PGG','QGG','RGG','SGG','TGG','VGG','WGG','XleGG','YGG','oxMGG','GG']
    return [file_names,xgg_types]

#specifiy modification characteristic fragments m/z
def Character_mz(type_xgg):
    xgg_mz_dict = {
        'MGG': 246.0834,
        'oxMGG': 262.0783,
        'AGG': 186.0800,
        'CamCGG':275.0736,
        'DGG': 230.0699,
        'EGG': 244.0855,
        'FGG': 262.1113,
        'GGG': 172.0644,
        'HGG': 252.1018,
        'XleGG': 228.1270,
        'KGG': 243.1379,
        'NGG': 229.0859,
        'PGG': 212.0957,
        'QGG': 243.1015,
        'RGG': 271.1440,
        'SGG': 202.0750,
        'TGG': 216.0906,
        'VGG': 214.1113,
        'WGG': 301.1222,
        'YGG': 278.1063,
        'GG': 115.0429
    }
    if type_xgg in xgg_mz_dict:
        return xgg_mz_dict[type_xgg]
    else:
        print('No such modification!')
        return 0

#characteristic fragment detection with 0.1 tolerance in Best score raw file and scan number
def CharacterIdentify(datalist,raw_file_path,character_type,delta_character=0.1):
    character_peak = Character_mz(character_type)
    for index,row in datalist.iterrows():
        rawfile = MSFileReader(raw_file_path+'\\'+row['Best score raw file']+'.raw')
        mz = list(rawfile.GetMassListFromScanNum(row['Best score scan number'])[0][0])
        rawfile.Close()
        mznp = np.array(mz)
        mznp = abs(mznp - character_peak)
        mznp.sort()  
        if mznp[0] <= delta_character:
            datalist.loc[index,'Characteristic Fragment'] = 1
        print(character_type, index)
    return datalist

#data clean by potential contaminant and reverse decoy, and protein site annotation
def DataClean(dataset):
    dataset.replace('#DIV/0!', np.nan, inplace=True)
    dataset_clean = dataset.loc[((dataset['Potential contaminant'].apply(str) != '+') & (dataset['Reverse'].apply(str) != '+'))]
    dataset_clean.loc[:,'ACCID'] = dataset_clean.loc[:,'Protein'].str.extract('sp\|(.*?)\|')
    dataset_clean.loc[:,'PROTEINS'] = dataset_clean.loc[:,'Protein'].str.extract('sp\|.*?\|(.*?)_')
    dataset_clean.loc[:,'ORGANISM'] = dataset_clean.loc[:,'Fasta headers'].str.extract('_(.*?)\s')
    dataset_clean.loc[:,'FullNames'] = dataset_clean.loc[:,'Fasta headers'].str.extract('\s(.*)\sOS=')
    return dataset_clean

#specify file path and file names
file_path = '..\\GGX\\Sites\\'
output_path = 'Character\\'
ubisite_raw_path = '..\\ubisite\\'
file_names = GetFileName()[0]
xgg_types = GetFileName()[1]

#loop for characteristic fragment detection and data clean
for file_name, type_xgg in zip(file_names,xgg_types):
    file_input = file_path + file_name
    file_output = output_path + file_name[:-4] + '_op.csv'
    print(file_name,type_xgg)

    dataset = pd.read_table(file_input)
    dataset_clean = DataClean(dataset)
    dataset_char = CharacterIdentify(dataset_clean, ubisite_raw_path, type_xgg,)
    dataset_op = dataset_char
    dataset_op.loc[:,'Mod'] = type_xgg
    dataset_op.to_csv(file_output,index=False)