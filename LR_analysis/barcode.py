import pandas as pd
import os,sys
import os.path
import glob

abspath = os. getcwd()
print(abspath)
inFilePath = sys.argv[1]
inFilePath2 = sys.argv[2]
outFilePath = sys.argv[3]

df1 = pd.read_csv(inFilePath2,sep='\t',encoding = "ISO-8859-1",header = (0))

try:
   try:
      df = pd.read_csv(inFilePath,sep='\t',encoding = "ISO-8859-1",header = (0))
   except pd.errors.EmptyDataError as e:
      print("file error, I'll try again")
except IndexError as e1:
   print("No data")

df2=pd.merge(df,df1,how='inner',on=["barcode"],right_index=True)

df2.to_csv(outFilePath, index = False,sep='\t')
