

import pandas as pd
from sklearn import preprocessing
import chem
import numpy as np

pd.options.mode.chained_assignment = None

class PreProcessor(object):
    def __init__(self,df,cols=None):
        dfcols = list(df.columns)
        self.signal_cols = cols
        self.df = df
        self.scaled = False
        self.filtered =False
        if self.signal_cols ==None:
            print("Choose numerical columns to use for analysis (signal cols)")
            self.signal_cols = self._column_picker(self.df)
        elif isinstance(self.signal_cols,list):
            self.signal_cols = [dfcols[int(i)] for i in cols]

    def _scale_values(self,df2scale):
        '''scale the values in cols and append to the end of the 
        data frame with <orig_colname>_scaled naming convention'''
        scaler = preprocessing.MinMaxScaler()
        scaled_vals = scaler.fit_transform(df2scale[self.signal_cols].values)
        new_colnames = ['{}_scaled'.format(c) for c in self.signal_cols]
        for i,ncn in enumerate(new_colnames):
            df2scale[ncn] = scaled_vals[:,i]
        return df2scale

    def _column_picker(self,df2pick):
        columns = list(df2pick.columns)
        for i,c in enumerate(columns):
            print('[{}] : {}'.format(i,c))
        selCs = input("\nwhich columns would you like to choose? > ")
        print('\n')
        selCs = [int(c) for c in selCs.split(',')] 
        return [columns[i] for i in selCs]

    def scale(self):
        '''groupby "SET" and then scale the columns of import'''
        dfs = []
        cols2build = {}
        for col in ["LIN_MEAN","CYC_MEAN"]:
            print("Pick Columns for {}".format(col))
            cols2build[col] = self._column_picker(self.df)
    
        for s,g in self.df.groupby('SET'):
            scaled_df =self._scale_values(g)
            for col,cols in cols2build.items():
                scaled_df[col] = np.mean([scaled_df['{}_scaled'.format(cols[0])].values,
                    scaled_df['{}_scaled'.format(cols[1])].values],axis=0)
            
            scaled_df['CYC_RATIO'] = scaled_df['{}_scaled'.format(cols2build["CYC_MEAN"][0])] / scaled_df['{}_scaled'.format(cols2build["CYC_MEAN"][1])]
            
            scaled_df["CYC2LIN_RATIO"] = scaled_df['CYC_MEAN'] / scaled_df['LIN_MEAN']
            dfs.append(scaled_df)   
        self.df = pd.concat(dfs,ignore_index=True)
        self.scaled = True

    def filt(self,intensity=99.9,ratio=99):
        if not self.scaled:
            self.scale
        
        print("*"*40)
        fdf = self.df[np.log2(self.df.CYC_RATIO) < 2 ]
        fdf = fdf[np.log2(fdf.CYC_RATIO) > -2]
        print("Sequences after Variance filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))
        print("*"*40)

        intp = np.percentile(fdf.LIN_MEAN,intensity)
        fdf = fdf[fdf.LIN_MEAN > intp]
        print("Sequences after Intensity filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))
        print("*"*40)

        ratp = np.percentile(fdf.CYC2LIN_RATIO, ratio)
        fdf = fdf[fdf.CYC2LIN_RATIO > ratp]
        print("Sequences after Ratio filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))
        print("*"*40)

        self.filtered = True
        self.fdf = fdf

    def rfilt(self,intensity=99.9,ratio=99):
        if not self.scaled:
            self.scale
        
        print("*"*40)
        fdf = self.df[np.log2(self.df.CYC_RATIO) < 2 ]
        fdf = fdf[np.log2(fdf.CYC_RATIO) > -2]
        print("Sequences after Variance filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))
        print("*"*40)

        ratp = np.percentile(fdf.CYC2LIN_RATIO, ratio)
        fdf = fdf[fdf.CYC2LIN_RATIO > ratp]
        print("Sequences after Ratio filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))
        print("*"*40)

        intp = np.percentile(fdf.LIN_MEAN,intensity)
        fdf = fdf[fdf.LIN_MEAN > intp]
        print("Sequences after Intensity filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))
        print("*"*40)

        
        self.filtered = True
        self.fdf = fdf

    def gen_sequences(self):
        if not self.filtered:
            print("need to filter first use self.filt or self.rfilt")
            return False
        self.hit_seqs = []
        new_sequence_names = []
        for probe in self.fdf.PROBE_SEQUENCE.values:
            mons = [mond["_"].residue]
            mons.extend([mond[x].residue for x in probe])
            new_seq = chem.seq.from_dict(mons)
            new_seq.cyclize()
            new_seq.deprotect()
            self.hit_seqs.append(new_seq)
            new_sequence_names.append(str(new_seq))
        self.fdf['new_seq_name'] = new_sequence_names



crude_monomers = pd.read_excel('C:\\Users\\cpye\\repos\\HitExtractor\\crude_monomers.xlsx')
sing_letter_codes = pd.read_excel('C:\\Users\\cpye\\repos\\HitExtractor\\single_letter_codes.xlsx')
print(list(sing_letter_codes.columns))
crude_monomers.head()

sing2abbrv = {str(sing).replace(u'\xa0',''):str(abbrv).replace(u'\xa0','') for sing,abbrv in sing_letter_codes.values}
abbrv2sing = {str(abbrv).replace(u'\xa0',''):str(sing).replace(u'\xa0','') for sing,abbrv in sing_letter_codes.values}
abbrv2smi = {str(abbrv).replace(u'\xa0',''):smi for abbrv,smi in crude_monomers[['abbrv','smiles']].values}

mons = []
mond = {}
for used_mon in sing_letter_codes.abbrv.values:
    used_mon = used_mon.replace(u'\xa0','')
    
    try:
        mon = chem.monomer(abbrv=used_mon,smiles=abbrv2smi[used_mon])
        mons.append(mon)
        mond[abbrv2sing[used_mon]] = mon
    except KeyError:
        print("bad mon: {}".format(used_mon))
