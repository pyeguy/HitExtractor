

import pandas as pd
from sklearn import preprocessing
import chem
import numpy as np
from collections import OrderedDict

class PreProcessor(object):
    def __init__(self,df,cols=None):
        self.cols = cols
        self.df = df
        self.scaled = False
        self.filtered =False
        if self.cols ==None:
            print("Choose numerical columns to use for analysis (signal cols)")
            self.cols = self._column_picker(self.df)
    def _scale_values(self,df2scale):
        '''scale the values in cols and append to the end of the 
        data frame with <orig_colname>_scaled naming convention'''
        scaler = preprocessing.MinMaxScaler()
        scaled_vals = scaler.fit_transform(df2scale[self.cols].values)
        new_colnames = ['{}_scaled'.format(c) for c in self.cols]
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
        
        fdf = self.df[np.log2(self.df.CYC_RATIO) < 2 ]
        fdf = fdf[np.log2(self.df.CYC_RATIO) > -2]
        print("Sequences after Variance filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))

        intp = np.percentile(fdf.LIN_MEAN,intensity)
        fdf = fdf[fdf.LIN_MEAN > intp]
        print("Sequences after Intensity filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))

        ratp = np.percentile(fdf.CYC2LIN_RATIO, ratio)
        fdf = fdf[fdf.CYC2LIN_RATIO > ratp]
        print("Sequences after Ratio filter: {}".format(len(set(fdf.PROBE_SEQUENCE))))

        self.filtered = True
        self.fdf = fdf

 



