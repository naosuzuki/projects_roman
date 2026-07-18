import pandas as pd
import numpy
import matplotlib.pyplot as plt
import sys

def read_sim(csvfile):
   df=pd.read_csv(csvfile,comment='#',delim_whitespace=True)
   print(df)
   snIa_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRSUM']>=40.0)]
   snIa_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']==0) & (df['SNRSUM']>=40.0)]
   ccsn_wide=df[(df['FIELD']=='WIDE') & (df['SIM_TEMPLATE_INDEX']!=0)]
   ccsn_deep=df[(df['FIELD']=='DEEP') & (df['SIM_TEMPLATE_INDEX']!=0)]

   mjd_wide=snIa_wide['SIM_PEAKMJD'].to_numpy()
   mjd_deep=snIa_deep['SIM_PEAKMJD'].to_numpy()
   ccsnmjd_wide=ccsn_wide['SIM_PEAKMJD'].to_numpy()
   ccsnmjd_deep=ccsn_deep['SIM_PEAKMJD'].to_numpy()

   zcmb_wide=snIa_wide['SIM_ZCMB'].to_numpy()
   zcmb_deep=snIa_deep['SIM_ZCMB'].to_numpy()
   zccsn_wide=ccsn_wide['SIM_ZCMB'].to_numpy()
   zccsn_deep=ccsn_deep['SIM_ZCMB'].to_numpy()

   host_zmag_wide=snIa_wide['HOST_MAG_Z'].to_numpy()
   host_zmag_deep=snIa_deep['HOST_MAG_Z'].to_numpy()
   hostccsn_zmag_wide=ccsn_wide['HOST_MAG_Z'].to_numpy()
   hostccsn_zmag_deep=ccsn_deep['HOST_MAG_Z'].to_numpy()

   mb_wide=snIa_wide['SIM_DLMAG'].to_numpy()
   mb_deep=snIa_deep['SIM_DLMAG'].to_numpy()

   #mb_wide=snIa_wide['SIM_mB'].to_numpy()
   #mb_deep=snIa_deep['SIM_mB'].to_numpy()
   mb_wide=snIa_wide['PEAKMAG_Z'].to_numpy()
   mb_deep=snIa_deep['PEAKMAG_Z'].to_numpy()
   find_activeSNIa_num(mb_wide,mb_deep)

   mbccsn_wide=ccsn_wide['PEAKMAG_Z'].to_numpy()
   mbccsn_deep=ccsn_deep['PEAKMAG_Z'].to_numpy()
   find_activeCCSN_num(mbccsn_wide,mbccsn_deep)
   sys.exit(1)
   #plot_mjd_vs_mag(mjd_wide,mjd_deep,mb_wide,mb_deep)
   plot_mjd_vs_mag(ccsnmjd_wide,ccsnmjd_deep,mbccsn_wide,mbccsn_deep)

   ra_wide=snIa_wide['HOST_RA'].to_numpy()
   dec_wide=snIa_wide['HOST_DEC'].to_numpy()
   ra_deep=snIa_deep['HOST_RA'].to_numpy()
   dec_deep=snIa_deep['HOST_DEC'].to_numpy()
   plot_radec(ra_wide,dec_wide,ra_deep,dec_deep)

def find_activeCCSN_num(mb_wide,mb_deep):
# North Factor
   nmax_wide=int(len(mb_wide)*38./(38.+27.))
   N_mb_wide=mb_wide[0:nmax_wide]
   nmax_deep=int(len(mb_deep)*7./(7.+16.))
   N_mb_deep=mb_deep[0:nmax_deep]
   print('CCSN Wide=',nmax_wide,'Deep=',nmax_deep,'Total=',nmax_wide+nmax_deep)
   flag_wide=numpy.where(N_mb_wide<24.0,1,0)
   flag_deep=numpy.where(N_mb_deep<24.0,1,0)
   print('Total CCSN < 24','wide=',numpy.sum(flag_wide),'deep=',numpy.sum(flag_deep),'total=',numpy.sum(flag_wide)+numpy.sum(flag_deep))
   # Visibility Factor
   print('Subaru Visible Total CCSN < 24','wide=',numpy.sum(flag_wide)*10./24.,'deep=',numpy.sum(flag_deep)*10./24.,'total=',(numpy.sum(flag_wide)+numpy.sum(flag_deep))*10./24.)
   print('Subaru specz CCSN < 24','wide=',numpy.sum(flag_wide)*10./24.*0.7,'deep=',numpy.sum(flag_deep)*10./24.*0.7,'total=',(numpy.sum(flag_wide)+numpy.sum(flag_deep))*10./24.*0.7)

def find_activeSNIa_num(mb_wide,mb_deep):
# North Factor
   nmax_wide=int(len(mb_wide)*38./(38.+27.))
   N_mb_wide=mb_wide[0:nmax_wide]
   nmax_deep=int(len(mb_deep)*7./(7.+16.))
   N_mb_deep=mb_deep[0:nmax_deep]
   print('SNIa Wide=',nmax_wide,'Deep=',nmax_deep,'Total=',nmax_wide+nmax_deep)
   flag_wide=numpy.where(N_mb_wide<24.0,1,0)
   flag_deep=numpy.where(N_mb_deep<24.0,1,0)
   print('Total SNIa < 24','wide=',numpy.sum(flag_wide),'deep=',numpy.sum(flag_deep),'total=',numpy.sum(flag_wide)+numpy.sum(flag_deep))
   # Visibility Factor
   print('Subaru Visible Total SNIa < 24','wide=',numpy.sum(flag_wide)*10./24.,'deep=',numpy.sum(flag_deep)*10./24.,'total=',(numpy.sum(flag_wide)+numpy.sum(flag_deep))*10./24.)
   print('Subaru Specz SNIa < 24 : 70% success','wide=',numpy.sum(flag_wide)*10./24.*0.7,'deep=',numpy.sum(flag_deep)*10./24.*0.7,'total=',(numpy.sum(flag_wide)+numpy.sum(flag_deep))*10./24.*0.7)


def plot_mjd_vs_mag(mjd_wide,mjd_deep,mb_wide,mb_deep):
#  zbin=numpy.arange(30)*0.1
   plt.scatter(mjd_wide,mb_wide,c='r',s=0.1)
   plt.scatter(mjd_deep,mb_deep,c='b',s=0.1)
   plt.savefig('time_vs_mag_v04_ccsn.png')
   plt.clf()
   plt.close()

def plot_radec(ra_wide,dec_wide,ra_deep,dec_deep):
   plt.scatter(ra_wide,dec_wide,c='b',s=0.1)
   plt.scatter(ra_deep,dec_deep,c='r',s=0.1)
   plt.savefig('snIamap_v01.png')
   plt.clf()
   plt.close()


csvfile='../data/OUT_ROMAN_TRIGGER.SNANA.TEXT'
csvfile='../data/OUT_ROMAN_TRIGGER+PEAKMAGS.SNANA.TEXT'
read_sim(csvfile)
