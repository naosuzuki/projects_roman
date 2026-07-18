import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy
import sys

# From SSP Overview Paper
# 34.56500  -4.85000
# 35.81388  -5.30852
# 35.59346  -4.06236
# 36.83275  -4.46317

sxdsx=[34.565,35.5934582,36.83275,35.813875]
sxdsy=[-4.85,-4.0623556,-4.4631722,-5.3085194]

def plot_transient(master_sn,master_agn,hstlist):
   dfsn=pd.read_csv(master_sn,usecols=[0,2,3],skiprows=1,names=['snname','ra','dec'])
   print(dfsn)
   dfagn=pd.read_csv(master_agn,usecols=[1,4,5],skiprows=1,names=['snname','ra','dec'])
   print(dfagn)
   dfhst=pd.read_csv(hstlist,usecols=[0,1,2],delim_whitespace=True,skiprows=1,names=['snname','ra','dec'],comment='#')
   print(dfhst)

   ra=dfsn.loc[:,'ra']    ; dec=dfsn.loc[:,'dec']
   hstra=dfhst.loc[:,'ra']; hstdec=dfhst.loc[:,'dec']
   agnra=dfagn.loc[:,'ra']; agndec=dfagn.loc[:,'dec']

   print(dfsn)	
   fig,ax=plt.subplots()
   # Font to be Times
   plt.rc('font',family='serif')
   #ymin=-1.0*ymax*0.05  ; ymax=ymax*1.5; dy=ymax-ymin
   #plt.title(sn.snname)

   xmin=37.9  ; xmax=33.5  ; dx=xmax-xmin
   ymin=-6.2  ; ymax=-3.1  ; dx=xmax-xmin

   plt.xlim([xmin,xmax])
   plt.ylim([ymin,ymax])
   plt.xlabel('RA')
   plt.ylabel('DEC')

   #Ia Good
   plt.scatter(hstra,hstdec,c='r',s=30.0,label='HST SNIa')
   plt.scatter(ra,dec,facecolors='none',edgecolors='r',s=10.0,label='Good SNIa')
   plt.scatter(agnra,agndec,c='b',s=2.0,label='Transients')
   radius=numpy.sqrt(2.00/numpy.pi)

   for i in range(4):
      r=patches.Circle((sxdsx[i],sxdsy[i]),radius,color='k',fill=False)
      ax.add_artist(r)

   plt.text(34.9,-5.9,r'SXDS Ultra Deep')
   outputfilename='HSCtransients_sxds.png'
   ax.legend()
   plt.savefig(outputfilename)
   plt.clf()
   plt.close()

   for i in range(4):
      r=patches.Circle((sxdsx[i],sxdsy[i]),radius,color='k',fill=False)
      ax.add_artist(r)

#plot_transient(master_sn,master_agn,hstlist)
plot_desi(master_sn,master_agn,hstlist,desifile)
