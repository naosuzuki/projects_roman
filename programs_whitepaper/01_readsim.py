import pandas as pd
import numpy
import matplotlib.pyplot as plt
import sys

def read_sim(csvfile):
   df=pd.read_csv(csvfile,comment='#',delim_whitespace=True)
   print(df)
#   snIa_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRMAX']>=5.0) & (df['SNRSUM']>5.)]
#   snIa_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRMAX']>=5.0) & (df['SNRSUM']>5.)]
#   snIa_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRMAX']>=5.0)]
#   snIa_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRMAX']>=5.0)]
   snIa_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRSUM']>=40.0)]
   snIa_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRSUM']>=40.0)]
#   ccsn_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']>0) & (df['SNRSUM']>=1.0)]
#   ccsn_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']>0) & (df['SNRSUM']>=1.0)]
   ccsn_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']!=0)]
   ccsn_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']!=0)]
   #print(ccsn_wide)
   #for i in range(len(ccsn_deep)):
   #    print(i,ccsn_deep['SIM_TEMPLATE_INDEX'].iloc[i])
   #sys.exit(1)

   zcmb_wide=snIa_wide['SIM_ZCMB'].to_numpy()
   zcmb_deep=snIa_deep['SIM_ZCMB'].to_numpy()
   zccsn_wide=ccsn_wide['SIM_ZCMB'].to_numpy()
   zccsn_deep=ccsn_deep['SIM_ZCMB'].to_numpy()

   #plot_z_vs_sn(zcmb_wide,zcmb_deep,zccsn_wide,zccsn_deep)

   host_zmag_wide=snIa_wide['HOST_MAG_Z'].to_numpy()
   host_zmag_deep=snIa_deep['HOST_MAG_Z'].to_numpy()
   hostccsn_zmag_wide=ccsn_wide['HOST_MAG_Z'].to_numpy()
   hostccsn_zmag_deep=ccsn_deep['HOST_MAG_Z'].to_numpy()
   #plot_mag_vs_hostgalaxy(host_zmag_wide,host_zmag_deep,hostccsn_zmag_wide,hostccsn_zmag_deep)
   plot_mag_and_z_vs_hostgalaxy(host_zmag_wide,host_zmag_deep,hostccsn_zmag_wide,hostccsn_zmag_deep)


def plot_z_vs_sn_north():
   zbin=numpy.arange(30)*0.1+0.05
   snbin=numpy.zeros(30)
# North Factor
   nmax_wide=int(len(host_zmag_wide)*38./(38.+27.))
   N_host_zmag_wide=host_zmag_wide[0:nmax_wide]
   nmax_deep=int(len(host_zmag_deep)*7./(7.+16.))
   N_host_zmag_deep=host_zmag_deep[0:nmax_deep]
   print('SNIa Wide=',nmax_wide,'Deep=',nmax_deep)

#   for i in range(30): 
   
#   for i in range(len(N_host_zmag_wide)):
       
        


def plot_z_vs_sn(zcmb_wide,zcmb_deep,zccsn_wide,zccsn_deep):
#  zbin=numpy.arange(30)*0.1
   zbin=numpy.arange(30)*0.1+0.05
   plt.hist(zcmb_wide,bins=zbin)
   plt.hist(zcmb_deep,bins=zbin,alpha=0.5)
   plt.hist(zccsn_wide,bins=zbin,alpha=0.3)
   plt.hist(zccsn_deep,bins=zbin,alpha=0.3)
   plt.savefig('z_vs_SNIa_v12.png')
   plt.clf()
   plt.close()
   print('SNIa WIDE total=',len(zcmb_wide),'North=',len(zcmb_wide)/65.*38.)
   print('SNIa DEEP total=',len(zcmb_deep),'North=',len(zcmb_deep)/23.*7.)
   print('CCSN WIDE total=',len(zccsn_wide),'North=',len(zccsn_wide)/65.*38.)
   print('CCSN DEEP total=',len(zccsn_deep),'North=',len(zccsn_deep)/23.*7.)

def plot_mag_vs_hostgalaxy(host_zmag_wide,host_zmag_deep,hostccsn_zmag_wide,hostccsn_zmag_deep):
   magbin=numpy.arange(24)*0.5+18.0
   plt.hist(host_zmag_wide,bins=magbin)
   plt.hist(host_zmag_deep,bins=magbin,alpha=0.5)
   plt.hist(hostccsn_zmag_wide,bins=magbin,alpha=0.3)
   plt.hist(hostccsn_zmag_deep,bins=magbin,alpha=0.3)
   plt.savefig('z_vs_hostzmag_v12.png')
   plt.clf()
   plt.close()

def plot_mag_and_z_vs_hostgalaxy(host_zmag_wide,host_zmag_deep,hostccsn_zmag_wide,hostccsn_zmag_deep):
# North Factor
   nmax_wide=int(len(host_zmag_wide)*38./(38.+27.))
   N_host_zmag_wide=host_zmag_wide[0:nmax_wide]
   nmax_deep=int(len(host_zmag_deep)*7./(7.+16.))
   N_host_zmag_deep=host_zmag_deep[0:nmax_deep]
   print('SNIa Wide=',nmax_wide,'Deep=',nmax_deep)

# For CCSN
   nmax_ccsn_wide=int(len(hostccsn_zmag_wide)*38./(38.+27.))
   N_hostccsn_zmag_wide=hostccsn_zmag_wide[0:nmax_ccsn_wide]
   nmax_ccsn_deep=int(len(hostccsn_zmag_deep)*7./(7.+16.))
   N_hostccsn_zmag_deep=hostccsn_zmag_deep[0:nmax_ccsn_deep]
   print('CCSN Wide=',nmax_ccsn_wide,'Deep=',nmax_ccsn_deep)

   zmaghost=numpy.concatenate((N_host_zmag_wide,N_host_zmag_deep),axis=0)
   print('SNIa Total=',len(zmaghost))

   zmaghost_ccsn=numpy.concatenate((N_hostccsn_zmag_wide,N_hostccsn_zmag_deep),axis=0)
   print('CCSN Total=',len(zmaghost_ccsn))

   magbin=numpy.arange(24)*0.5+18.0
   speczbin=numpy.zeros(24)
   speczbin_ccsn=numpy.zeros(24)
   zmaghist,zmagbin=numpy.histogram(zmaghost,bins=magbin)
   zmaghist_ccsn,zmagbin=numpy.histogram(zmaghost_ccsn,bins=magbin)
   for i in range(23):
      if(magbin[i]<20.0):
          factor=1.0
          speczbin[i]=zmaghist[i] 
          speczbin_ccsn[i]=zmaghist_ccsn[i] 
      elif(magbin[i]<26.0):
          x=magbin[i]
          factor=1.0+0.87884*x-0.140627*x**2+0.00751911*x**3-0.000133365*x**4
# PFS Configuration Factor
      factor*=0.7
      speczbin[i]=zmaghist[i]*factor
      speczbin_ccsn[i]=zmaghist_ccsn[i]*factor
   plt.hist(zmaghost,bins=magbin,color='b',align='left')
   plt.bar(magbin,speczbin,color='r')
   print('SNIa Total targets=',len(zmaghost))
   print('SNIa Total specz=',numpy.sum(speczbin))
   plt.savefig('z_vs_hostzmag_SNIa_v18.png')
   plt.clf()
   plt.close()

## CCSN
   plt.hist(zmaghost_ccsn,bins=magbin,color='b',align='left')
   plt.bar(magbin,speczbin_ccsn,color='r')
   print('CCSN Total targets=',len(zmaghost_ccsn))
   print('CCSN Total specz=',numpy.sum(speczbin_ccsn))
   plt.savefig('z_vs_hostzmag_CCSN_v18.png')
   plt.clf()
   plt.close()

   for i in range(23):
     print(magbin[i],zmaghist[i],"%5.2f"%(speczbin[i]),"%5.2f"%(speczbin[i]/zmaghist[i]),\
                                 "%6.2f"%(speczbin_ccsn[i]),"%5.2f"%(speczbin_ccsn[i]/zmaghist_ccsn[i]))

csvfile='../data/OUT_ROMAN_WIDE+DEEP.TEXT'
csvfile='../data/OUT_ROMAN_PASSCUTS.SNANA.TEXT'
csvfile='../data/OUT_ROMAN_TRIGGER.SNANA.TEXT'
read_sim(csvfile)
