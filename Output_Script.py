###############################################################################
# Author = Nathan van der Wielen
# Date = 2022
# Purpose = This script is used to analyse the output files of the SDA analysis
# This script works in tandem with SDA.py, Reflector.py, ABAQUS, and ESATAN.
###############################################################################

#------------------------------------------------------------------------------
# LIBRARIES
#------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time

t = time.process_time()

#------------------------------------------------------------------------------
# GLOBAL VARIABLES
#------------------------------------------------------------------------------
C = 299792458
SaveFig = True

#------------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Calculates the Root mean square
# -----------------------------------------------------------------------------
def rms(data, shape):  # Function returns the RMS values
    dsquare = np.power(data,2)
    dsum = dsquare.sum(axis=1)
    rms = np.sqrt(dsum/shape[1])*1000
    return rms

# -----------------------------------------------------------------------------
# Checks last time CSV file was accessed. Used to limit number of times CSV files are read
# -----------------------------------------------------------------------------
def last_opened(filepath):
    data = pd.read_csv(open('D:/1-Master Thesis/ABAQUS/Data Analysis/Database_Time.csv'),
                       sep=',',
                       header=0,
                       index_col=(0))
    new_time = os.path.getmtime(filepath)
    if filepath in data.index:
        old_time = data.at[filepath,'Time']
    else:
        old_time = 0.0
    if int(new_time) != int(old_time):
        data.at[filepath,'Time'] = new_time
        data.to_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/Database_Time.csv')
        check = True  # Means the data is new and need to be updated
    else:
        check = False  # Means the data in not new and doesn't need to be updated  
    # check = True # Overide to update anyways
    return check
    
#------------------------------------------------------------------------------
# CLASSES
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Abaqus object to access and analyse Abaqus results
#------------------------------------------------------------------------------
class ABAQUS:
    def __init__(self, filepath):
        print(str(filepath)+' loading...')
        # ---------------------------------------------------------------------
        # Find the Case name of Object from Abaqus file
        # ---------------------------------------------------------------------
        self.address = filepath.split(sep='/') 
        self.case = self.address[-1]
        
        # ---------------------------------------------------------------------
        # To not repeatdly access Displacement CSV file, output to data analysis file only if CSV data is new
        # ---------------------------------------------------------------------
        if os.path.exists(filepath+'_deflection.csv'):
            self.check_u = last_opened(filepath=filepath+'_deflection.csv')
            
            # -----------------------------------------------------------------
            # Output new displacement data analysis to data analysis file
            # -----------------------------------------------------------------
            if self.check_u == True:
                print('Updating deflection...')
                self.file_u = pd.read_csv(open(filepath+'_deflection.csv','r'),
                                          sep=',',
                                          header=0) # Read the Deflection file create DataFrame
                
                self.file_u.columns = self.file_u.columns.str.replace(' ', '') # Remove the blank spaces in the columns
                self.data_u = self.file_u.pivot_table(index='Frame',
                                                      columns='NodeLabel',
                                                      values='U-Magnitude')
                
                self.time = pd.DataFrame(self.data_u.index)
                self.time[['Increment','Time']] = self.time['Frame'].str.split('=', 1, expand=True) 
                self.time['Time'] = self.time['Time'].astype(float)
                # self.time = self.time.iloc[:-1]
                self.time.set_index('Frame',
                                    inplace=True)
                self.time.drop(['Increment'],
                               axis=1,
                               inplace=True)
                self.time = self.time.squeeze(axis=1)
                
                self.shape_u = self.data_u.shape
                self.data_u = self.data_u.iloc[ 0:self.shape_u[0] , 0:self.shape_u[1]]
                
                self.mean_u = self.data_u.mean(axis=1)
                self.max_u = self.data_u.max(axis=1)
                self.min_u = self.data_u.min(axis=1)
                self.rms = rms(data=self.data_u, shape=self.shape_u)
                
                if os.path.exists('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv'):
                    os.remove('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv')
                
                self.analysis = pd.concat([self.mean_u, self.max_u, self.min_u, self.rms, self.time],
                                          axis=1,
                                          ignore_index=False)
                self.analysis.columns = ['mean_u','max_u','min_u','RMS','time']
                self.analysis.to_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv')
            
            # -----------------------------------------------------------------
            # Get the displacement data from data analysis file if data hasn't changed
            # -----------------------------------------------------------------
            else:
                print('Fetching deflection data...')
                self.time = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv',
                                        usecols=(['time']),
                                        dtype=float)
                self.mean_u = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv',
                                          usecols=(['mean_u']),
                                          dtype=float)
                self.max_u = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv',
                                         usecols=(['max_u']),
                                         dtype=float)
                self.min_u = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv',
                                         usecols=(['min_u']),
                                         dtype=float)
                self.rms = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_u.csv',
                                       usecols=(['RMS']),
                                       dtype=float)
        else:
            print('File does not exist')    
        
        # ---------------------------------------------------------------------
        # To not repeatdly access temperature CSV file, output to data analysis file only if CSV data is new
        # ---------------------------------------------------------------------
        if os.path.exists(filepath+'_thermal.csv'):
            self.check_t = last_opened(filepath=filepath+'_thermal.csv')
            
            # -----------------------------------------------------------------
            # Output new temperature data analysis to data analysis file
            # -----------------------------------------------------------------
            if self.check_t == True:
                print('Updating thermal...')
                self.file_t = pd.read_csv(open(filepath+'_thermal.csv', 'r'),
                                          sep=',',
                                          header=0) # Read the Deflection file create DataFrame

                self.file_t.columns = self.file_t.columns.str.replace(' ', '') # Remove the blank spaces in the columns

                self.data_t = self.file_t.pivot_table(index='Frame',
                                                      columns='NodeLabel',
                                                      values='TEMP')
                self.shape_t = self.data_t.shape
                self.data_t = self.data_t.iloc[0:self.shape_t[0], 0:self.shape_t[1]]
        
                self.mean_t = self.data_t.mean(axis=1)
                self.max_t = self.data_t.max(axis=1)
                self.min_t = self.data_t.min(axis=1)
                self.delta_t = self.max_t - self.min_t
                
                if os.path.exists('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv'):
                    os.remove('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv')
                
                self.analysis = pd.concat([self.mean_t, self.max_t, self.min_t,self.delta_t],
                                          axis=1)
                self.analysis.columns = ['mean_t','max_t','min_t','delta_t']
                self.analysis.to_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv')
            
            # -----------------------------------------------------------------
            # Get the temperature data from data analysis file if data hasn't changed
            # -----------------------------------------------------------------
            else:
                print('Fetching thermal data...')
                self.mean_t = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv',
                                          usecols=(['mean_t']),
                                          dtype=float)
                self.max_t = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv',
                                         usecols=(['max_t']),
                                         dtype=float)
                self.min_t = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv',
                                         usecols=(['min_t']),
                                         dtype=float)
                self.delta_t = pd.read_csv('D:/1-Master Thesis/ABAQUS/Data Analysis/'+self.case+'_analysis_t.csv',
                                           usecols=(['delta_t']),
                                           dtype=float)
        else:
            print('File does not exist')
        print(str(filepath)+' done loading')
    
    # -------------------------------------------------------------------------
    # Plot the temperature of the reflector for object
    # -------------------------------------------------------------------------  
    def plot_t(self):
        self.fig, self.ax = plt.subplots()
        self.ax.plot(self.time.astype(float),
                     self.mean_t, marker='.',
                     label='Average Temperature')
        self.ax.plot(self.time.astype(float),
                     self.max_t, marker='.',
                     label='Max Temperature')
        self.ax.plot(self.time.astype(float),
                     self.min_t,
                     marker='.',
                     label='Min Temperature')
        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_title(label=self.case+' Temperature Study')
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Temperature (K)')
        self.fig.legend(bbox_to_anchor=(1, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
        if SaveFig:
            self.fig.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/'+self.case+'_Temperature_Study.png',
                             bbox_inches='tight' )
    
    # -------------------------------------------------------------------------
    # Plot the displacement of the reflector for object
    # -------------------------------------------------------------------------  
    def plot_u(self):
        self.fig, self.ax = plt.subplots()
        self.ax.plot(self.time.astype(float),
                     self.mean_u,
                     marker='.',
                     label='Average Displacement')
        self.ax.plot(self.time.astype(float),
                     self.max_u,
                     marker='.',
                     label='Max Displacement')
        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_title(label=self.case+' Deflection Study')
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Displacement (m)')
        self.fig.legend(bbox_to_anchor=(1, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
        if SaveFig:
            self.fig.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/'+self.case+'_Displacement_Study.png',
                             bbox_inches='tight' )
    
    # -------------------------------------------------------------------------
    # Plot the RMS of surface displacement for object
    # -------------------------------------------------------------------------  
    def plot_RMS(self):        
        self.fig, self.ax = plt.subplots()
        self.ax.plot(self.time.astype(float),
                     self.rms,
                     marker='.',
                     label='Average RMS')
        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_title(label=self.case+' RMS Study')
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('RMS')
        self.fig.legend(bbox_to_anchor=(1, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
    
    # -------------------------------------------------------------------------
    # Get max, min, max delta, and eclipse times for object
    # -------------------------------------------------------------------------  
    def get_info(self):
        print('Max temperature:', self.max_t.max(),'(K)      @ Time:', self.max_t.idxmax(),'(s)')
        print('Min temperature:', self.min_t.min(),'(K)      @ Time:', self.min_t.idxmin(),'(s)')
        print('Max temperature delta:', self.delta_t.max(),'(K)      @ Time:', self.delta_t.idxmax(),'(s)')
        print('Max deflection:', self.max_u.max(),'(mm)      @ Time:', self.max_u.idxmax(),'(s)')


#------------------------------------------------------------------------------
# Esatan object to access and analyse Esatan results
#------------------------------------------------------------------------------
class ESATAN:
    def __init__(self, filepath):
        # ---------------------------------------------------------------------
        # Find the Case name of Object from ESATAN file
        # ---------------------------------------------------------------------
        self.address = filepath.split(sep='/')
        self.case = self.address[3]+'_'+self.address[-1]
        
        # ---------------------------------------------------------------------
        # Import data resulting from ESATAN analysis
        # ---------------------------------------------------------------------
        self.r_top_temp = pd.read_csv(filepath + '/Top_Temperature.csv',
                                      sep=',',
                                      header=1,
                                      index_col=('TIME'))
        self.r_top_qa = pd.read_csv(filepath +'/Top_QA.csv',
                                    sep=',',
                                    header=1,
                                    index_col=('TIME'))
        self.r_Top_QE = pd.read_csv(filepath +'/Top_QE.csv',
                                    sep=',',
                                    header=1,
                                    index_col=('TIME'))
        self.r_Top_QS = pd.read_csv(filepath +'/Top_QS.csv',
                                    sep=',',
                                    header=1,
                                    index_col=('TIME'))
        self.r_Bot_QA = pd.read_csv(filepath +'/Bottom_QA.csv',
                                    sep=',',
                                    header=1,
                                    index_col=('TIME'))
        self.r_Bot_QE = pd.read_csv(filepath +'/Bottom_QE.csv',
                                    sep=',',
                                    header=1,
                                    index_col=('TIME'))
        self.r_Bot_QS = pd.read_csv(filepath +'/Bottom_QS.csv',
                                    sep=',',
                                    header=1,
                                    index_col=('TIME'))
        
        # ---------------------------------------------------------------------
        # Find the shape of the ESATAN csv files
        # ---------------------------------------------------------------------
        self.shape = self.r_top_temp.shape
        
        # ---------------------------------------------------------------------
        # Get the time from results
        # ---------------------------------------------------------------------
        self.time = pd.DataFrame(data=self.r_top_temp.index,
                                 columns=['TIME'],
                                 dtype= float)
        self.time = self.time.iloc[:-1]
        self.time = self.time.squeeze(axis=1)
        
        # ---------------------------------------------------------------------
        # Get the temperature and energy values from results 
        # ---------------------------------------------------------------------
        self.Top_Temp = self.r_top_temp.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        self.Top_QA = self.r_top_qa.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        self.Top_QE = self.r_Top_QE.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        self.Top_QS = self.r_Top_QS.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        self.Bottom_QA = self.r_Bot_QA.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        self.Bottom_QE = self.r_Bot_QE.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        self.Bottom_QS = self.r_Bot_QS.iloc[ 0:self.shape[0]-1 , 1:self.shape[1]-1 ]
        
        # ---------------------------------------------------------------------
        # Translate time to angle of orbit
        # ---------------------------------------------------------------------
        self.angle = (self.time-self.time.min())/(self.time.max()-self.time.min())*360
    
        # ---------------------------------------------------------------------
        # Find the mean, max, min, and delta of the tempertures
        # ---------------------------------------------------------------------
        self.mean = self.Top_Temp.mean(axis=1)
        self.max = self.Top_Temp.max(axis=1)
        self.min = self.Top_Temp.min(axis=1)
        self.delta = self.max - self.min
    
    
        # ---------------------------------------------------------------------
        # If Report Exists, fetch eclipse data by parsing ESATAN report
        # ---------------------------------------------------------------------
        if os.path.exists(filepath + '/Report.rpt'):
            self.RTP = True
            self.report = pd.read_fwf(open(filepath + '/Report.rpt','r'),
                                      header = None)
            self.report = self.report[self.report[0].str.contains('Eclipse')]
            self.report[['Parameter', 'Value']] = self.report[0].str.split('    ', 1, expand=True)
            self.report.drop(self.report.columns[[0]], 
                             axis=1,
                             inplace=True)
            self.report = self.report.astype({'Value': float})
            self.eclipse_entry = self.report.iloc[0,1]
            self.eclipse_exit = self.report.iloc[1,1]
            self.eclipse_percent = self.report.iloc[2,1]
        else:
            self.RTP = False
            
        print(filepath+' done loading')
    
    # -------------------------------------------------------------------------
    # Plot the max, min, and mean temperature of object
    # -------------------------------------------------------------------------
    def plot(self):
        self.fig, self.ax = plt.subplots()
               
        self.ax.plot(self.mean,
                     marker='.',
                     label='Average Temperature')
        self.ax.plot(self.max,
                     marker='.',
                     label='Max Temperature')
        self.ax.plot(self.min,
                     marker='.',
                     label='Min Temperature')
        self.ax.fill_between(x=self.time,
                             y1=self.min,
                             y2=self.max,
                             alpha=0.5)
        
        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Temperature (K)')
         
        self.secax = self.ax.twiny()
        self.secax.set_xticks(range(0,396,36))
        if self.RTP == True:
            self.secax.axvline(x= self.eclipse_entry,
                               color='blue',
                               linestyle='dotted',
                               label='Eclipse Entry')
            self.secax.axvline(x= self.eclipse_exit,
                               color='red',
                               linestyle='dotted',
                               label='Eclipse Exit')
        self.secax.set_xlabel('Orbital angle')
        
        self.fig.suptitle(t=self.case+' Temperature Study',
                          y = 1.05)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
        if SaveFig:
            self.fig.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/'+self.case+'_Temperature_Study.png',
                             bbox_inches='tight' )
    
    # -------------------------------------------------------------------------
    # Plot the different energies absorbed by the reflector on both sides
    # -------------------------------------------------------------------------
    def plot_temp_energy(self):
        self.fig, self.ax = plt.subplots()
               
        self.ax.plot(self.mean,
                     marker='.',
                     label='Average Temperature')
        
        self.ax2 = self.ax.twinx()
        self.ax2.plot(self.Top_QA.mean(axis=1),
                      marker='.',
                      label='Top QA')
        self.ax2.plot(self.Top_QE.mean(axis=1),
                      marker='.',
                      label='Top QE')
        self.ax2.plot(self.Top_QS.mean(axis=1),
                      marker='.',
                      label='Top QS')
        
        self.ax2.plot(self.Bottom_QA.mean(axis=1),
                      marker='.',
                      label='Bottom QA')
        self.ax2.plot(self.Bottom_QE.mean(axis=1),
                      marker='.',
                      label='Bottom QE')
        self.ax2.plot(self.Bottom_QS.mean(axis=1),
                      marker='.',
                      label='Bottom QS')
        
        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Temperature (K)')
        self.ax2.set_ylabel('Energy (W)')
        
        self.secax = self.ax.twiny()
        self.secax.set_xticks(range(0,396,36))
        self.secax.axvline(x= self.eclipse_entry,
                           color='blue',
                           linestyle='dotted',
                           label='Eclipse Entry')
        self.secax.axvline(x= self.eclipse_exit,
                           color='red',
                           linestyle='dotted',
                           label='Eclipse Exit')
        self.secax.set_xlabel('Orbital angle')
        
        self.fig.suptitle(t=self.case+' Energy Study',
                          y = 1.05)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
        if SaveFig:
            self.fig.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/'+self.case+'_Energy_Study.png',
                             bbox_inches='tight')
    
    # -------------------------------------------------------------------------
    # plot the change in temperature for object
    # -------------------------------------------------------------------------
    def plot_diff(self):
        self.fig, self.ax = plt.subplots()
               
        self.ax.plot(self.mean.diff(),
                     marker='.',
                     label='Change in Average Temperature')
        self.ax.plot(self.max.diff(),
                     marker='.',
                     label='Change in Max Temperature')
        self.ax.plot(self.min.diff(),
                     marker='.',
                     label='Change in Min Temperature')
        
        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Temperature (K)')
         
        self.secax = self.ax.twiny()
        self.secax.set_xticks(range(0,396,36))
        if self.RTP == True:
            self.secax.axvline(x= self.eclipse_entry,
                               color='blue',
                               linestyle='dotted',
                               label='Eclipse Entry')
            self.secax.axvline(x= self.eclipse_exit,
                               color='red',
                               linestyle='dotted',
                               label='Eclipse Exit')
        self.secax.set_xlabel('Orbital angle')
        
        self.fig.suptitle(t=self.case+' Change in Temperature Study',
                          y = 1.05)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
        if SaveFig:
            self.fig.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/'+self.case+'_Temperature_Study_Derivative.png',
                             bbox_inches='tight' )
    
    # -------------------------------------------------------------------------
    # plot the change in temperature for object
    # -------------------------------------------------------------------------
    def plot_delta(self):
        self.fig, self.ax = plt.subplots()
               
        self.ax.plot(self.max - self.min,
                     marker='.',
                     label='Delta in Temperature across Reflector')

        self.ax.autoscale(axis='x',
                          tight=True)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Temperature (K)')
         
        self.secax = self.ax.twiny()
        self.secax.set_xticks(range(0,396,36))
        if self.RTP == True:
            self.secax.axvline(x= self.eclipse_entry,
                               color='blue',
                               linestyle='dotted',
                               label='Eclipse Entry')
            self.secax.axvline(x= self.eclipse_exit,
                               color='red',
                               linestyle='dotted',
                               label='Eclipse Exit')
        self.secax.set_xlabel('Orbital angle')
        
        self.fig.suptitle(t=self.case+' Change in Temperature Study',
                          y = 1.05)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.fig.legend(bbox_to_anchor=(1.05, 1),
                        loc='upper left',
                        borderaxespad=0.)
        self.ax.grid(True)
        if SaveFig:
            self.fig.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/'+self.case+'_Temperature_Study_Delta.png',
                             bbox_inches='tight' )
    # -------------------------------------------------------------------------
    # Get max, min, max delta, and eclipse times for object
    # -------------------------------------------------------------------------                
    def get_info(self):
        print('Max temperature:', round(self.max.max(),2),'(K)      @ Time:', round(self.max.idxmax(),2),'(s)')
        print('Min temperature:', round(self.min.min(),2),'(K)      @ Time:', round(self.min.idxmin(),2),'(s)')
        print('Max temperature delta:', round(self.delta.max(),2),'(K)      @ Time:', round(self.delta.idxmax(),2),'(s)')
        print('Time of Eclipse:', round(self.eclipse_exit-self.eclipse_entry,2),'(s)')
        # print('Total Orbit time:', round(self.time.loc[self.shape[0]-2,'TIME'],2),'(s)')
        
# -----------------------------------------------------------------------------
# ESATAN DATA
# -----------------------------------------------------------------------------
print('ESATAN data loading...')
LEO = ESATAN('D:/ESATAN_TMS/Model/Reflector/esatan/LEO_T')
GEO = ESATAN('D:/ESATAN_TMS/Model/Reflector/esatan/GEO_T')
Summer = ESATAN('D:/ESATAN_TMS/Model/Reflector/esatan/Summer_T')
Winter = ESATAN('D:/ESATAN_TMS/Model/Reflector/esatan/Winter_T')
D10 = ESATAN('D:/ESATAN_TMS/Model/Reflector_10m/esatan/LEO_T')
F800 = ESATAN('D:/ESATAN_TMS/Model/Reflector_F800/esatan/LEO_T')
SiO = ESATAN('D:/ESATAN_TMS/Model/Reflector_SiO/esatan/LEO_T')
White = ESATAN('D:/ESATAN_TMS/Model/Reflector_White/esatan/LEO_T')
print('ESATAN data done loading')

# -----------------------------------------------------------------------------
# ABAQUS DATA
# -----------------------------------------------------------------------------
print('ABAQUS data loading...')
Ref_LEO = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/Reflector_LEO')
Ref_Summer = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/Reflector_Summer')
Ref_Winter = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/Reflector_Winter')
Ref_SiO_LEO = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/Reflector_SiO_LEO')

SDA_LEO = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_LEO')
SDA_Summer = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_Summer')
SDA_Winter = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_Winter')
SDA_CFRP = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_CFRP')
SDA_SW1 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_SW1')
SDA_SP1 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_SP1')
SDA_SP2 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_SP2')
SDA_CTE1 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_CTE1')
SDA_CTE2 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_CTE2')
SDA_TH1 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_TH1')
SDA_TH2 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_TH2')
SDA_K1 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_K1')
SDA_K2 = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_K2')
SDA_Improved = ABAQUS('D:/1-Master Thesis/ABAQUS/Reports/SDA_Optimal')
print('Abaqus data done loading')

# ------------------------------------------------------------------------------
# RMS Calculations for Different Frequency Bands
# ------------------------------------------------------------------------------
freq = np.array([1.5E9,6E9,10E9,16E9,33.25E9])
l = C/freq
Gl = np.array([[15],[20],[30],[50]])
r = l/Gl*1000
rmsf = pd.DataFrame(r,
                    columns=['1~2:L','4~8:S','8~12:C','12~18:Ku','26.5~40:Ka'],
                    index=['lambda/15','lambda/20','lambda/30','lambda/50'])
rmsf = rmsf.transpose()


# ------------------------------------------------------------------------------
# ANALYSIS PLOTS
# ------------------------------------------------------------------------------

fig_size = (6,4) # Parameter for the size of the figure

# ------------------------------------------------------------------------------
# Plot the displacent of SDA_CFRP and Reflector to look at the mean error / Plots
# ------------------------------------------------------------------------------
def A_mean_error():
    plt.figure(num='ABAQUS_mean_error',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_LEO.mean_u*1000,
             marker='.',
             label='Reflector mean')
    plt.plot(SDA_CFRP.time,
             SDA_CFRP.mean_u*1000,
             marker='.',
             label='CFRP SDA mean')
    
    plt.plot(Ref_LEO.time,
             Ref_LEO.max_u*1000,
             marker='.',
             label='Reflector max')
    plt.plot(SDA_CFRP.time,
             SDA_CFRP.max_u*1000,
             marker='.',
             label='CFRP SDA max')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.autoscale(axis='x',
                  tight=True)
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Mean Displacement Error Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Mean_Error.png',
                      bbox_inches='tight')
        
        
    plt.figure(num='ABAQUS_mean_error difference',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             abs(Ref_LEO.mean_u*1000-SDA_CFRP.mean_u*1000),
             marker='.',
             label='Mean absolute difference')
    plt.plot(Ref_LEO.time,
             abs(Ref_LEO.max_u*1000-SDA_CFRP.max_u*1000),
             marker='.',
             label='Max absolute difference')
   
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.autoscale(axis='x',
                  tight=True)
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Mean Displacement Error Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Mean_Error_Difference.png',
                      bbox_inches='tight')

# ------------------------------------------------------------------------------
# Plot the displacent of SDA_CFRP and Reflector to look at the RMS error / Plots
# ------------------------------------------------------------------------------
def A_rms_error():
    plt.figure(num='ABAQUS_RMS_error',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_LEO.rms,
             marker='.',
             label='Ref_LEO rms')
    plt.plot(SDA_CFRP.time,
             SDA_CFRP.rms,
             marker='.',
             label='SDA_CFRP rms')
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.autoscale(axis='x',
                  tight=True)
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('RMS Error Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/RMS_Error.png',
                      bbox_inches='tight')

# ------------------------------------------------------------------------------
# Plot the temperature of SDA_CFRP and Reflector to look at the temperature error / Plots
# ------------------------------------------------------------------------------
def A_temp_error():
    plt.figure(num='ABAQUS_temp_error',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_LEO.mean_t, 
             marker='.',
             color='red',
             label='Ref_LEO')
    plt.fill_between(Ref_LEO.time.loc[:,'time'],
                     Ref_LEO.max_t.loc[:,'max_t'],
                     Ref_LEO.min_t.loc[:,'min_t'],
                     color='red',alpha=0.3)
    
    plt.plot(SDA_CFRP.time,
             SDA_CFRP.mean_t,
             marker='.',
             color='blue',
             label='SDA_CFRP')
    plt.fill_between(SDA_CFRP.time.loc[:,'time'],
                     SDA_CFRP.max_t.loc[:,'max_t'],
                     SDA_CFRP.min_t.loc[:,'min_t'],
                     color='blue',
                     alpha=0.2)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.autoscale(axis='x',
                  tight=True)
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Temp Error Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Temp_Error.png',
                      bbox_inches='tight')
        
    plt.figure(num='ABAQUS_temp_error2',
               figsize=fig_size)
    plt.plot(SDA_LEO.time,
             SDA_LEO.mean_t,
             marker='.',
             color='red',
             label='SDA_LEO')
    plt.fill_between(SDA_LEO.time.loc[:,'time'],
                     SDA_LEO.max_t.loc[:,'max_t'],
                     SDA_LEO.min_t.loc[:,'min_t'],
                     color='red',
                     alpha=0.3)
    
    plt.plot(SDA_CFRP.time,
             SDA_CFRP.mean_t,
             marker='.',
             color='blue',
             label='SDA_CFRP')
    plt.fill_between(SDA_CFRP.time.loc[:,'time'],
                     SDA_CFRP.max_t.loc[:,'max_t'],
                     SDA_CFRP.min_t.loc[:,'min_t'],
                     color='blue',
                     alpha=0.2)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.autoscale(axis='x',
                  tight=True)
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Temp Error Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Temp_Error2.png',
                      bbox_inches='tight')

#------------------------------------------------------------------------------
# Diameter Study / Plots
#------------------------------------------------------------------------------
def diameter_study():
    plt.figure(num='Diameter_Study',
               figsize=fig_size)
    plt.plot(LEO.mean,
             marker='.',
             label='Avg.Temp D=3m')
    plt.plot(D10.mean,
             marker='.',
             label='Avg.Temp D=10m')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Diameter Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Diameter_Study.png',
                      bbox_inches='tight')

#------------------------------------------------------------------------------
# ESATAN to ABAQUS Study / Plots
#------------------------------------------------------------------------------
def EA_study():
    #--------------------------------------------------------------------------
    # LEO
    #--------------------------------------------------------------------------
    plt.subplots(nrows=1,
                 ncols=3,
                 figsize=(14,5),
                 num='ESATAN-ABAQUS LEO Study')
   
    plt.subplot(1, 3, 1)
    plt.plot(LEO.time,
              LEO.min,
              color='blue',
              marker='.')
    plt.plot(Ref_LEO.time,
              Ref_LEO.min_t,
              color='orange',
              marker='.')
    plt.plot(SDA_LEO.time,
              SDA_LEO.min_t,
              color='green',
              marker='.')
    plt.title('Min Temperature')
    plt.grid(True)
    plt.ylabel('Temperature (K)')
    plt.xlabel('Time (s)')
 
    plt.subplot(1, 3, 2)
    plt.plot(LEO.time,
              LEO.mean,
              color='blue',
              marker='.')
    plt.plot(Ref_LEO.time,
              Ref_LEO.mean_t,
              color='orange',
              marker='.')
    plt.plot(SDA_LEO.time,
              SDA_LEO.mean_t,
              color='green',
              marker='.')
    plt.title('Mean Temperature')
    plt.grid(True)
    plt.xlabel('Time (s)')
     
    plt.subplot(1, 3, 3)
    plt.plot(LEO.time,
              LEO.max,
              marker='.',
              color='blue',
              label='LEO\n(ESATAN)')
    plt.plot(Ref_LEO.time,
              Ref_LEO.max_t,
              marker='.',
              color='orange',
              label='Ref_LEO\n(ABAQUS)')
    plt.plot(SDA_LEO.time,
              SDA_LEO.max_t,
              marker='.',
              color='green',
              label='SDA_LEO\n(ABAQUS)')
    plt.title('Max Temperature')
    plt.grid(True)
    plt.xlabel('Time (s)')
    plt.legend(bbox_to_anchor=(1.05, 1, 0, 0),
                loc='upper left',
                borderaxespad=0.)
    
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/ESATAN_ABAQUS_LEO_Study.png',
                      bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Summer
    #--------------------------------------------------------------------------
    plt.subplots(nrows=1,
                 ncols=3,
                 figsize=(14,5),
                 num='ESATAN-ABAQUS Summer Study')
   
    plt.subplot(1, 3, 1)
    plt.plot(Summer.time,
              Summer.min,
              color='blue',
              marker='.')
    plt.plot(Ref_Summer.time,
              Ref_Summer.min_t,
              color='orange',
              marker='.')
    plt.plot(SDA_Summer.time,
              SDA_Summer.min_t,
              color='green',
              marker='.')
    plt.title('Min Temperature')
    plt.grid(True)
    plt.ylabel('Temperature (K)')
    plt.xlabel('Time (s)')
 
    plt.subplot(1, 3, 2)
    plt.plot(Summer.time,
             Summer.mean,
              color='blue',
              marker='.')
    plt.plot(Ref_Summer.time,
              Ref_Summer.mean_t,
              color='orange',
              marker='.')
    plt.plot(SDA_Summer.time,
              SDA_Summer.mean_t,
              color='green',
              marker='.')
    plt.title('Mean Temperature')
    plt.grid(True)
    plt.xlabel('Time (s)')
     
    plt.subplot(1, 3, 3)
    plt.plot(Summer.time,
             Summer.max,
              marker='.',
              color='blue',
              label='Summer\n(ESATAN)')
    plt.plot(Ref_Summer.time,
              Ref_Summer.max_t,
              marker='.',
              color='orange',
              label='Ref_Summer\n(ABAQUS)')
    plt.plot(SDA_Summer.time,
              SDA_Summer.max_t,
              marker='.',
              color='green',
              label='SDA_Summer\n(ABAQUS)')
    plt.title('Max Temperature')
    plt.grid(True)
    plt.xlabel('Time (s)')
    plt.legend(bbox_to_anchor=(1.05, 1, 0, 0),
                loc='upper left',
                borderaxespad=0.)
    
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/ESATAN_ABAQUS_Summer_Study.png',
                      bbox_inches='tight')

    #--------------------------------------------------------------------------
    # Winter
    #--------------------------------------------------------------------------
    plt.subplots(nrows=1,
                 ncols=3,
                 figsize=(14,4),
                 num='ESATAN-ABAQUS Winter Study')
   
    plt.subplot(1, 3, 1)
    plt.plot(Winter.time,
              Winter.min,
              color='blue',
              marker='.')
    plt.plot(Ref_Winter.time,
              Ref_Winter.min_t,
              color='orange',
              marker='.')
    plt.plot(SDA_Winter.time,
              SDA_Winter.min_t,
              color='green',
              marker='.')
    plt.title('Min Temperature')
    plt.grid(True)
    plt.ylabel('Temperature (K)')
    plt.xlabel('Time (s)')
 
    plt.subplot(1, 3, 2)
    plt.plot(Winter.time,
             Winter.mean,
              color='blue',
              marker='.')
    plt.plot(Ref_Winter.time,
              Ref_Winter.mean_t,
              color='orange',
              marker='.')
    plt.plot(SDA_Winter.time,
              SDA_Winter.mean_t,
              color='green',
              marker='.')
    plt.title('Mean Temperature')
    plt.grid(True)
    plt.xlabel('Time (s)')
     
    plt.subplot(1, 3, 3)
    plt.plot(Winter.time,
             Winter.max,
              marker='.',
              color='blue',
              label='Winter\n(ESATAN)')
    plt.plot(Ref_Winter.time,
              Ref_Winter.max_t,
              marker='.',
              color='orange',
              label='Ref_Winter\n(ABAQUS)')
    plt.plot(SDA_Winter.time,
              SDA_Winter.max_t,
              marker='.',
              color='green',
              label='SDA_Winter\n(ABAQUS)')
    plt.title('Max Temperature')
    plt.grid(True)
    plt.xlabel('Time (s)')
    plt.legend(bbox_to_anchor=(1.05, 1, 0, 0),
                loc='upper left',
                borderaxespad=0.)
    
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/ESATAN_ABAQUS_Winter_Study.png',
                      bbox_inches='tight')

        
# ------------------------------------------------------------------------------
# ESATAN to ABAQUS Differential Study / Plots
# ------------------------------------------------------------------------------
def EA_diff_study():
    plt.figure(num='ESATAN-ABAQUS LEO Diff Study',
               figsize=fig_size)
    plt.plot(LEO.time,
             LEO.mean.diff(),
             marker='.',
             label=' \u0394TT_AVG LEO\n(ESATAN)')
    plt.plot(Ref_LEO.time,
             Ref_LEO.mean_t.diff(),
             marker='.',
             label='\u0394T_AVG Ref_LEO\n(ABAQUS)')
    plt.plot(SDA_LEO.time,
             SDA_LEO.mean_t.diff(),
             marker='.',
             label='\u0394T_AVG SDA_LEO\n(ABAQUS)')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('ESATAN-ABAQUS LEO Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/ESATAN_ABAQUS_LEO_Diff_Study.png',
                      bbox_inches='tight')
    
    plt.figure(num='ESATAN-ABAQUS Summer Diff Study',
               figsize=fig_size)
    plt.plot(Summer.time,
             Summer.mean.diff(),
             marker='.',
             label='\u0394T_AVG Summer\n(ESATAN)')
    plt.plot(Ref_Summer.time,
             Ref_Summer.mean_t.diff(),
             marker='.',
             label='\u0394T_AVG Ref_Summer\n(ABAQUS)')
    plt.plot(SDA_Summer.time,
             SDA_Summer.mean_t.diff(),
             marker='.',
             label='\u0394T_AVG SDA_Summer\n(ABAQUS)')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('ESATAN-ABAQUS Summer Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/ESATAN_ABAQUS_Summer_Diff_Study.png',
                      bbox_inches='tight')
        
    plt.figure(num='ESATAN-ABAQUS Winter Diff Study',
               figsize=fig_size)
    plt.plot(Winter.time,
             Winter.mean.diff(),
             marker='.',
             label='\u0394T_AVG Winter\n(ESATAN)')
    plt.plot(Ref_Winter.time,
             Ref_Winter.mean_t.diff(),
             marker='.',
             label='\u0394T_AVG Ref_Winter\n(ABAQUS)')
    plt.plot(SDA_Winter.time,
             SDA_Winter.mean_t.diff(),
             marker='.',
             label='\u0394T_AVG SDA_Winter\n(ABAQUS)')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('ESATAN-ABAQUS Wnter Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/ESATAN_ABAQUS_Winter_Diff_Study.png',
                      bbox_inches='tight')

# ------------------------------------------------------------------------------
# Ref vs SDA Displacement Study / Plots
# ------------------------------------------------------------------------------
def RS_disp_study():
    #--------------------------------------------------------------------------
    # LEO displacement study / Plots
    #--------------------------------------------------------------------------
    plt.figure(num='Ref vs SDA LEO Disp Study',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_LEO.mean_u*1000,
             marker='.',
             label='U_Mean Ref_LEO')
    plt.plot(SDA_LEO.time,
             SDA_LEO.mean_u*1000,
             marker='.',
             label='U_Mean SDA_LEO')
    
    plt.plot(Ref_LEO.time,
             Ref_LEO.max_u*1000,
             marker='.',
             label='U_Max Ref_LEO')
    plt.plot(SDA_LEO.time,
             SDA_LEO.max_u*1000,
             marker='.',
             label='U_Max SDA_LEO')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Ref vs SDA LEO Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Ref_vs_SDA_LEO_Disp_Study.png',
                      bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Summer Displacement study / Plots
    #--------------------------------------------------------------------------
    plt.figure(num='Ref vs SDA Summer Disp Study',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_Summer.mean_u*1000,
             marker='.',
             label='U_Mean Ref_Summer')
    plt.plot(SDA_LEO.time,
             SDA_Summer.mean_u*1000,
             marker='.',
             label='U_Mean SDA_Summer')
    
    plt.plot(Ref_LEO.time,
             Ref_Summer.max_u*1000,
             marker='.',
             label='U_Max Ref_Summer')
    plt.plot(SDA_LEO.time,
             SDA_Summer.max_u*1000,
             marker='.',
             label='U_Max SDA_Summer')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Ref vs SDA Summer Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Ref_vs_SDA_Summer_Disp_Study.png',
                      bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Winter Displacement Study / Plots
    #--------------------------------------------------------------------------
    plt.figure(num='Ref vs SDA Winter Disp Study',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_Winter.mean_u*1000,
             marker='.',
             label='U_Mean Ref_Winter')
    plt.plot(SDA_LEO.time,
             SDA_Winter.mean_u*1000,
             marker='.',
             label='U_Mean SDA_Winter')
    
    plt.plot(Ref_LEO.time,
             Ref_Winter.max_u*1000,
             marker='.',
             label='U_Max Ref_Winter')
    plt.plot(SDA_LEO.time,
             SDA_Winter.max_u*1000,
             marker='.',
             label='U_Max SDA_Winter')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Ref vs SDA Winter Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Ref_vs_SDA_Winter_Disp_Study.png',
                      bbox_inches='tight')
        
# ------------------------------------------------------------------------------
# Ref vs SDA Mean Displacement Study / Plots
# ------------------------------------------------------------------------------
def RS_disp_mean_study():
    #--------------------------------------------------------------------------
    # LEO displacement study / Plots
    #--------------------------------------------------------------------------
    plt.figure(num='Ref vs SDA LEO Disp Study',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_LEO.mean_u*1000,
             marker='.',
             label='U_AVG Ref_LEO')
    plt.plot(SDA_LEO.time,
             SDA_LEO.mean_u*1000,
             marker='.',
             label='U_AVG SDA_LEO')
    # plt.plot(SDA_CFRP.time,
    #          SDA_CFRP.mean_u*1000,
    #          marker='.',
    #          label='U_AVG SDA_CFRP')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Ref vs SDA LEO Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Ref_vs_SDA_LEO_Mean_Disp_Study.png',
                      bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Summer Displacement study / Plots
    #--------------------------------------------------------------------------
    plt.figure(num='Ref vs SDA Summer Disp Study',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_Summer.mean_u*1000,
             marker='.',
             label='U_AVG Ref_Summer')
    plt.plot(SDA_LEO.time,
             SDA_Summer.mean_u*1000,
             marker='.',
             label='U_AVG SDA_Summer')
     
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Ref vs SDA Summer Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Ref_vs_SDA_Summer_Mean_Disp_Study.png',
                      bbox_inches='tight')
    
    #--------------------------------------------------------------------------
    # Winter Displacement Study / Plots
    #--------------------------------------------------------------------------
    plt.figure(num='Ref vs SDA Winter Disp Study',
               figsize=fig_size)
    plt.plot(Ref_LEO.time,
             Ref_Winter.mean_u*1000,
             marker='.',
             label='U_AVG Ref_Winter')
    plt.plot(SDA_LEO.time,
             SDA_Winter.mean_u*1000,
             marker='.',
             label='U_AVG SDA_Winter')
        
    plt.xlabel('Time (s)')
    plt.ylabel('Displacement (mm)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Ref vs SDA Winter Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Ref_vs_SDA_Winter_Mean_Disp_Study.png',
                      bbox_inches='tight')
    
#------------------------------------------------------------------------------
# Orbit Study / Plots
#------------------------------------------------------------------------------
def orbit_study():
    plt.figure(num='Orbit Study',
               figsize=fig_size)
    plt.fill_between(x=Summer.time,
                     y1=Summer.min,
                     y2=Summer.max,
                     alpha=0.5, label='Summer Orbit')
    plt.fill_between(x=Winter.time,
                     y1=Winter.min,
                     y2=Winter.max,
                     alpha=0.5,
                     label='Winter Orbit')
    plt.fill_between(x=LEO.time,
                     y1=LEO.min,
                     y2=LEO.max,
                     alpha=0.5,
                     label="LEO Orbit")
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (C)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Orbit Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Orbit_Study.png',
                    bbox_inches='tight')

#------------------------------------------------------------------------------
# Coating Study / Plots
#------------------------------------------------------------------------------
def coating_study():        
    plt.figure(num='Coating Study',
               figsize=fig_size)
    plt.fill_between(x=SiO.time,
                     y1=SiO.min,
                     y2=SiO.max,
                     color='purple',
                     alpha=0.5,
                     label='SiO2 Coating')
    plt.fill_between(x=White.time,
                     y1=White.min,
                     y2=White.max,
                     alpha=0.5,
                     color='cyan',
                     label='White Paint Coating')
    plt.fill_between(x=LEO.time,
                     y1=LEO.min,
                     y2=LEO.max,
                     alpha=0.5,
                     color='yellow',
                     label="CFRP Coating")
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Coating Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/Coating_Study.png',
                    bbox_inches='tight')
        
#------------------------------------------------------------------------------
# Spiral Width Study / Plots
#------------------------------------------------------------------------------
def SW_study():        
    #------------------------------------------------------------------------------
    # Spiral Width Temperature Study / Plots
    #------------------------------------------------------------------------------
    plt.figure(num='SW_t Study',
               figsize=fig_size)

    plt.plot([0.01, 0.05], [SDA_LEO.mean_t.loc[:,'mean_t'].mean(), SDA_SW1.mean_t.loc[:,'mean_t'].mean()],
             label='Mean_T')
    plt.plot([0.01, 0.05], [max(SDA_LEO.max_t.loc[:,'max_t']), max(SDA_SW1.max_t.loc[:,'max_t'])],
             label='Max_T')
    plt.plot([0.01, 0.05], [min(SDA_LEO.min_t.loc[:,'min_t']), min(SDA_SW1.min_t.loc[:,'min_t'])],
             label='Min_T')
    
    plt.xlabel('Spiral Width (m)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Spiral Width Temperature Study')
    
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/SW_t_Study.png',
                    bbox_inches='tight')

    #------------------------------------------------------------------------------
    # Spiral Width Displacement Study / Plots
    #------------------------------------------------------------------------------
    plt.figure(num='SW_u Study',
               figsize=fig_size)

    plt.plot([0.01, 0.05], [SDA_LEO.mean_u.loc[:,'mean_u'].mean(), SDA_SW1.mean_u.loc[:,'mean_u'].mean()],
             label='Mean_U')
    plt.plot([0.01, 0.05], [max(SDA_LEO.max_u.loc[:,'max_u']), max(SDA_SW1.max_u.loc[:,'max_u'])],
             label='Max_U')
    plt.plot([0.01, 0.05], [min(SDA_LEO.min_u.loc[:,'min_u']), min(SDA_SW1.min_u.loc[:,'min_u'])],
             label='Min_U')
    
    plt.xlabel('Spiral Width (m)')
    plt.ylabel('Displacement (m)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Spiral Width Displacement Study')
    
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/SW_u_Study.png',
                    bbox_inches='tight')

#------------------------------------------------------------------------------
# Spiral Parameter Study / Plots
#------------------------------------------------------------------------------
def SP_study():   
    #------------------------------------------------------------------------------
    # Spiral Parameter Temperature Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='SP_t Study',
               figsize=fig_size)
    
    plt.plot([0.04, 0.06, 0.08], [SDA_SP2.mean_t.loc[:,'mean_t'].mean(), SDA_LEO.mean_t.loc[:,'mean_t'].mean(), SDA_SP1.mean_t.loc[:,'mean_t'].mean()],
             label='Mean_T')
    plt.plot([0.04, 0.06, 0.08], [max(SDA_SP2.max_t.loc[:,'max_t']), max(SDA_LEO.max_t.loc[:,'max_t']), max(SDA_SP1.max_t.loc[:,'max_t'])],
             label='Max_T')
    plt.plot([0.04, 0.06, 0.08], [min(SDA_SP2.min_t.loc[:,'min_t']), min(SDA_LEO.min_t.loc[:,'min_t']), min(SDA_SP1.min_t.loc[:,'min_t'])],
             label='Min_T')
    
    plt.xlabel('Spiral Parameter')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Spiral Parameter Temperature Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/SP_t_Study.png',
                    bbox_inches='tight')
        
    #------------------------------------------------------------------------------
    # Spiral Parameter Displacement Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='SP_u Study',
               figsize=fig_size)
    
    plt.plot([0.04, 0.06, 0.08], [SDA_SP2.mean_u.loc[:,'mean_u'].mean(), SDA_LEO.mean_u.loc[:,'mean_u'].mean(), SDA_SP1.mean_u.loc[:,'mean_u'].mean()],
             label='Mean_U')
    plt.plot([0.04, 0.06, 0.08], [max(SDA_SP2.max_u.loc[:,'max_u']), max(SDA_LEO.max_u.loc[:,'max_u']), max(SDA_SP1.max_u.loc[:,'max_u'])],
             label='Max_U')
    plt.plot([0.04, 0.06, 0.08], [min(SDA_SP2.min_u.loc[:,'min_u']), min(SDA_LEO.min_u.loc[:,'min_u']), min(SDA_SP1.min_u.loc[:,'min_u'])], 
             label='Min_U')
    
    plt.xlabel('Spiral Parameter')
    plt.ylabel('Displacement (m)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Spiral Parameter Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/SP_u_Study.png',
                    bbox_inches='tight')
        
#------------------------------------------------------------------------------
# CTE Study / Plots
#------------------------------------------------------------------------------
def CTE_study():   
    #------------------------------------------------------------------------------
    # CTE Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='CTE_t Study',
               figsize=fig_size)
    
    plt.plot([-0.0000011, -0.0000005, 0.0], [SDA_CFRP.mean_t.loc[:,'mean_t'].mean(), SDA_CTE1.mean_t.loc[:,'mean_t'].mean(), SDA_CTE2.mean_t.loc[:,'mean_t'].mean()], 
             label='Mean_T')
    plt.plot([-0.0000011, -0.0000005, 0.0], [max(SDA_CFRP.max_t.loc[:,'max_t']), max(SDA_CTE1.max_t.loc[:,'max_t']), max(SDA_CTE2.max_t.loc[:,'max_t'])],
             label='Max_T')
    plt.plot([-0.0000011, -0.0000005, 0.0], [min(SDA_CFRP.min_t.loc[:,'min_t']), min(SDA_CTE1.min_t.loc[:,'min_t']), min(SDA_CTE2.min_t.loc[:,'min_t'])], 
             label='Min_T')
    
    plt.xlabel('CTE (/K)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('CTE Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/CTE_t_Study.png',
                    bbox_inches='tight')
        
    #------------------------------------------------------------------------------
    # CTE Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='CTE_u Study',
               figsize=fig_size)
    
    plt.plot([-0.0000011, -0.0000005, 0.0], [SDA_CFRP.mean_u.loc[:,'mean_u'].mean(), SDA_CTE1.mean_u.loc[:,'mean_u'].mean(), SDA_CTE2.mean_u.loc[:,'mean_u'].mean()],
             label='Mean_U')
    plt.plot([-0.0000011, -0.0000005, 0.0], [max(SDA_CFRP.max_u.loc[:,'max_u']), max(SDA_CTE1.max_u.loc[:,'max_u']), max(SDA_CTE2.max_u.loc[:,'max_u'])],
             label='Max_U')
    plt.plot([-0.0000011, -0.0000005, 0.0], [min(SDA_CFRP.min_u.loc[:,'min_u']), min(SDA_CTE1.min_u.loc[:,'min_u']), min(SDA_CTE2.min_u.loc[:,'min_u'])],
             label='Min_U')
    
    plt.xlabel('CTE (/K)')
    plt.ylabel('Displacement (m)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('CTE Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/CTE_u_Study.png',
                    bbox_inches='tight')

#------------------------------------------------------------------------------
# Thickness Study / Plots
#------------------------------------------------------------------------------
def TH_study():   
    #------------------------------------------------------------------------------
    # TH Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='TH_t Study',
               figsize=fig_size)
    
    plt.plot([0.0004, 0.001, 0.005], [SDA_LEO.mean_t.loc[:,'mean_t'].mean(), SDA_TH1.mean_t.loc[:,'mean_t'].mean(), SDA_TH2.mean_t.loc[:,'mean_t'].mean()],
             label='Mean_T')
    plt.plot([0.0004, 0.001, 0.005], [max(SDA_LEO.max_t.loc[:,'max_t']), max(SDA_TH1.max_t.loc[:,'max_t']), max(SDA_TH2.max_t.loc[:,'max_t'])],
             label='Max_T')
    plt.plot([0.0004, 0.001, 0.005], [min(SDA_LEO.min_t.loc[:,'min_t']), min(SDA_TH1.min_t.loc[:,'min_t']), min(SDA_TH2.min_t.loc[:,'min_t'])],
             label='Min_T')
    
    plt.xlabel('Thickness (m)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Thickness Temperature Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/TH_t_Study.png',
                    bbox_inches='tight')
        
    #------------------------------------------------------------------------------
    # TH Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='TH_u Study', figsize=fig_size)
    
    plt.plot([0.0004, 0.001, 0.005], [SDA_LEO.mean_u.loc[:,'mean_u'].mean(), SDA_TH1.mean_u.loc[:,'mean_u'].mean(), SDA_TH2.mean_u.loc[:,'mean_u'].mean()],
             label='Mean_U')
    plt.plot([0.0004, 0.001, 0.005], [max(SDA_LEO.max_u.loc[:,'max_u']), max(SDA_TH1.max_u.loc[:,'max_u']), max(SDA_TH2.max_u.loc[:,'max_u'])],
             label='Max_U')
    plt.plot([0.0004, 0.001, 0.005], [min(SDA_LEO.min_u.loc[:,'min_u']), min(SDA_TH1.min_u.loc[:,'min_u']), min(SDA_TH2.min_u.loc[:,'min_u'])],
             label='Min_U')
    
    plt.xlabel('Thickness (m)')
    plt.ylabel('Displacement (m)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Thickness Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/TH_u_Study.png',
                    bbox_inches='tight')

#------------------------------------------------------------------------------
# Conductivity Study / Plots
#------------------------------------------------------------------------------
def K_study():   
    #------------------------------------------------------------------------------
    # Conductivity Temperature Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='K_t Study',
               figsize=fig_size)
    
    plt.plot([1, 20, 46.97], [SDA_K1.mean_t.loc[:,'mean_t'].mean(), SDA_K2.mean_t.loc[:,'mean_t'].mean(), SDA_CFRP.mean_t.loc[:,'mean_t'].mean()],
             label='Mean_T')
    plt.plot([1, 20, 46.97], [max(SDA_K1.max_t.loc[:,'max_t']), max(SDA_K2.max_t.loc[:,'max_t']), max(SDA_CFRP.max_t.loc[:,'max_t'])],
             label='Max_T')
    plt.plot([1, 20, 46.97], [min(SDA_K1.min_t.loc[:,'min_t']), min(SDA_K2.min_t.loc[:,'min_t']), min(SDA_CFRP.min_t.loc[:,'min_t'])],
             label='Min_T')
    
    plt.xlabel('Conductivty (m.S.K)')
    plt.ylabel('Temperature (K)')
    plt.legend(bbox_to_anchor=(1.05, 1), 
               loc='upper left', 
               borderaxespad=0.)
    plt.title('Conductivity Temperature Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/K_t_Study.png',
                    bbox_inches='tight')
        
    #------------------------------------------------------------------------------
    # Spiral Parameter Displacement Study / Plots
    #------------------------------------------------------------------------------     
    plt.figure(num='K_u Study',
               figsize=fig_size)
    
    plt.plot([1, 20, 46.97], [SDA_K1.mean_u.loc[:,'mean_u'].mean(), SDA_K2.mean_u.loc[:,'mean_u'].mean(), SDA_CFRP.mean_u.loc[:,'mean_u'].mean()], 
             label='Mean_U')
    plt.plot([1, 20, 46.97], [max(SDA_K1.max_u.loc[:,'max_u']), max(SDA_K2.max_u.loc[:,'max_u']), max(SDA_CFRP.max_u.loc[:,'max_u'])],
             label='Max_U')
    plt.plot([1, 20, 46.97], [min(SDA_K1.min_u.loc[:,'min_u']), min(SDA_K2.min_u.loc[:,'min_u']), min(SDA_CFRP.min_u.loc[:,'min_u'])],
             label='Min_U')
    
    plt.xlabel('Conductivity (m.S.K)')
    plt.ylabel('Displacement (m)')
    plt.legend(bbox_to_anchor=(1.05, 1),
               loc='upper left',
               borderaxespad=0.)
    plt.title('Conductivity Displacement Study')
    plt.grid(True)
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/K_u_Study.png',
                    bbox_inches='tight')
    
#------------------------------------------------------------------------------
# RMS surface error for different Frequency Bands / Plots
#------------------------------------------------------------------------------
def RMS():
    fig, ax = plt.subplots()
    ax.plot(rmsf['lambda/15'],
            linewidth=0,
            marker='.',
            label='lambda/15')
    ax.plot(rmsf['lambda/20'],
            linewidth=0,
            marker='.',
            label='lambda/20')
    ax.plot(rmsf['lambda/30'],
            linewidth=0,
            marker='.',
            label='lambda/30')
    ax.plot(rmsf['lambda/50'],
            linewidth=0,
            marker='.', 
            label='lambda/50')
    ax.axhline(y = max(Ref_LEO.rms.loc[:,'RMS']),
               color='red',
               linestyle = '-',
               label='Max Ref_LEO RMS')
    ax.axhline(y = max(SDA_LEO.rms.loc[:,'RMS']),
               color='blue',
               linestyle = '-',
               label='Max SDA_LEO RMS')
    ax.axhline(y = max(SDA_Summer.rms.loc[:,'RMS']),
               color='purple',
               linestyle = '-',
               label='Max SDA_Summer RMS')
    ax.axhline(y = max(SDA_Winter.rms.loc[:,'RMS']),
               color='yellow',
               linestyle = '-',
               label='Max SDA_Winter RMS')
    
    ax.set_yscale(value='log')
    ax.set_ylim(ymin=0.05)
    ax.set_xlabel('Frequency Band (GHz)')
    ax.set_ylabel('RMS Surface Error (mm)')
    ax.set_title('RMS Study')
    fig.legend(bbox_to_anchor=(0.95, 0.85, 0.01, 0.0),
               loc='upper left',
               borderaxespad=0.)
    ax.grid(visible=True,
            which='both',
            axis='both')
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/RMS_Study.png',
                      bbox_inches='tight')
        
def RMS1():
    fig, ax = plt.subplots()
    ax.plot(rmsf['lambda/15'],
            linewidth=0,
            marker='.',
            label='lambda/15')
    ax.plot(rmsf['lambda/20'],
            linewidth=0,
            marker='.',
            label='lambda/20')
    ax.plot(rmsf['lambda/30'],
            linewidth=0,
            marker='.',
            label='lambda/30')
    ax.plot(rmsf['lambda/50'],
            linewidth=0,
            marker='.',
            label='lambda/50')
    ax.axhline(y = max(SDA_Improved.rms.loc[:,'RMS']),
               color='blue',
               linestyle = '-',
               label='Max SDA_Improved RMS')
    ax.axhline(y = 0.3,
               color='red',
               linestyle = '-',
               label='Park')
    
    ax.set_yscale(value='log')
    ax.set_ylim(ymin=0.05)
    ax.set_xlabel('Frequency Band (GHz)')
    ax.set_ylabel('RMS Surface Error (mm)')
    ax.set_title('RMS Study')
    fig.legend(bbox_to_anchor=(0.95, 0.85, 0.01, 0.0),
               loc='upper left',
               borderaxespad=0.)
    ax.grid(visible=True,
            which='both',
            axis='both')
    if SaveFig == True:
        plt.savefig('D:/1-Master Thesis/6-Master Thesis Text/Figures/Results/RMS1_Study.png',
                      bbox_inches='tight')



print('Program done running')
print('Elapsed time:', time.process_time() - t)


