import pandas as pd
import numpy

def plothistogram(csvfile):
   df=pd.read_csv(csvfile,comment='#')
   print(df)
   zbin=numpy.zeros(int(len(df)))
   snbin=numpy.zeros(len(df),dtype=numpy.int32)
   magbin=numpy.zeros(len(df),dtype=numpy.int32)
   snmag=numpy.zeros(5383)
   snz=numpy.zeros(5383)
   for i in range(len(df)):
      snbin[i]=df['SNIaTotal'].iloc[i]
      magbin[i]=df['i-mag'].iloc[i]
      zbin[i]=df['z'].iloc[i]

   count=0 ; dz=0.05
   print(zbin)
   print(snbin)
   for i in range(len(df)-1):
      dmag=(magbin[i+1]-magbin[i])/float(snbin[i])
      dz=(zbin[i+1]-zbin[i])/float(snbin[i])
      for j in range(snbin[i]):
         print(i,j,count)
         snmag[count]=magbin[i]+dmag*float(j)
         snz[count]=zbin[i]+dz*float(j)
         count+=1

   for k in range(5383):
        print(k,snmag[k],snz[k])

   #   for j in range(.0)):
   #      count+=1
   #   print(i,num,z,imag)

csvfile='../csvfiles/RomanSN_WhitePaper.csv'
plothistogram(csvfile)

