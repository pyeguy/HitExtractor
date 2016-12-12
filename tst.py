from HitExtractor import PreProcessor
import pandas as pd
from IPython import embed

df = pd.read_csv('./Data/GST_Biotin.csv',index_col=0)
pp = PreProcessor(df)
pp.scale()
pp.filt()
embed()