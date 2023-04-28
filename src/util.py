# Standard library imports
import numpy as np
import pandas as pd    
import random
import time 
import scipy   
from scipy.interpolate import interp1d
from scipy import signal

def get_data(name):
    if name=='Benthic Stack':
        url = 'https://raw.githubusercontent.com/seonminahn/HMM-Stack/master/Prob_stack.txt'
        df = stack = pd.read_csv(url,sep='  ',engine='python',header=None)
        df.rename(columns={0:'age',1:'d18o',2:'std',3:'c95',4:'c5'},inplace=True)
        df = df.set_index('age')
        t_new = np.arange(0,3000,2.5)
        df_ = pd.DataFrame(index=t_new)
        df_['data'] = interp1d(df.index.values,df['d18o'])(t_new)
        df_['std'] = interp1d(df.index.values,df['std'])(t_new)
    
    elif name=='Ice Core CO2':
        url = 'https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/antarctica2015co2composite.txt'
        df = pd.read_csv(url,header=137,sep='\t')
        df['age_gas_calBP'] = df['age_gas_calBP']/1000
        t_new = np.arange(0,800,0.5)
        df_ = pd.DataFrame(index=t_new)
        df_['data'] = interp1d(df['age_gas_calBP'],df['co2_ppm'])(t_new)
        df_['std'] = interp1d(df['age_gas_calBP'],df['co2_1s_ppm'])(t_new)
        
    return df_

def next_power_of_2(x):  
    return 1 if x == 0 else 2**(x - 1).bit_length()

def periodogram(x,y):
    fs = 1/(np.mean(np.diff(x)))
    f,Pxx = signal.periodogram(y,fs=fs,nfft=next_power_of_2(len(y)))
    return f,Pxx

    
class SSA(object):
    """ 
    SOURCE: 
        - This code is taken directly from Braun et al., 2023: https://www.nature.com/articles/s43247-023-00717-5 
        - All credit to Braun et al., for the excellent paper and great code! 
        - The code itself was accessed via https://zenodo.org/record/7104976#.ZEvKvuxuerc
        - Braun cites this blog as the original source: https://www.kaggle.com/code/jdarcy/introducing-ssa-for-time-series-decomposition/notebook 
        
    """

    
    __supported_types = (pd.Series, np.ndarray, list)
    
    def __init__(self, tseries, L, save_mem=True):
        """
        Decomposes the given time series with a singular-spectrum analysis. Assumes the values of the time series are
        recorded at equal intervals.

        """
        
        # Tedious type-checking for the initial time series
        if not isinstance(tseries, self.__supported_types):
            raise TypeError("Unsupported time series object. Try Pandas Series, NumPy array or list.")
        
        # Checks to save us from ourselves
        self.N = len(tseries)
        if not 2 <= L <= self.N/2:
            raise ValueError("The window length must be in the interval [2, N/2].")
        
        self.L = L
        self.orig_TS = pd.Series(tseries)
        self.K = self.N - self.L + 1
        
        # Embed the time series in a trajectory matrix
        self.X = np.array([self.orig_TS.values[i:L+i] for i in range(0, self.K)]).T
        
        # Decompose the trajectory matrix
        self.U, self.Sigma, VT = np.linalg.svd(self.X)
        self.d = np.linalg.matrix_rank(self.X)
        
        self.TS_comps = np.zeros((self.N, self.d))
        
        if not save_mem:
            # Construct and save all the elementary matrices
            self.X_elem = np.array([ self.Sigma[i]*np.outer(self.U[:,i], VT[i,:]) for i in range(self.d) ])

            # Diagonally average the elementary matrices, store them as columns in array.           
            for i in range(self.d):
                X_rev = self.X_elem[i, ::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]
            
            self.V = VT.T
        else:
            # Reconstruct the elementary matrices without storing them
            for i in range(self.d):
                X_elem = self.Sigma[i]*np.outer(self.U[:,i], VT[i,:])
                X_rev = X_elem[::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]
            
            self.X_elem = "Re-run with save_mem=False to retain the elementary matrices."
            
            # The V array may also be very large under these circumstances, so we won't keep it.
            self.V = "Re-run with save_mem=False to retain the V matrix."
        
        # Calculate the w-correlation matrix.
        self.calc_wcorr()
            
    def components_to_df(self, n=0):
        """
        Returns all the time series components in a single Pandas DataFrame object.
        """
        if n > 0:
            n = min(n, self.d)
        else:
            n = self.d
        
        # Create list of columns - call them F0, F1, F2, ...
        cols = ["F{}".format(i) for i in range(n)]
        return pd.DataFrame(self.TS_comps[:, :n], columns=cols, index=self.orig_TS.index)
            
    
    def reconstruct(self, indices):
        """
        Reconstructs the time series from its elementary components, using the given indices. 
        """
        if isinstance(indices, int): indices = [indices]
        
        ts_vals = self.TS_comps[:,indices].sum(axis=1)
        return pd.Series(ts_vals, index=self.orig_TS.index)
    
    def calc_wcorr(self):
        """
        Calculates the w-correlation matrix
        """
             
        # Calculate the weights
        w = np.array(list(np.arange(self.L)+1) + [self.L]*(self.K-self.L-1) + list(np.arange(self.L)+1)[::-1])
        
        def w_inner(F_i, F_j):
            return w.dot(F_i*F_j)
        
        # Calculated weighted norms, ||F_i||_w, then invert.
        F_wnorms = np.array([w_inner(self.TS_comps[:,i], self.TS_comps[:,i]) for i in range(self.d)])
        F_wnorms = F_wnorms**-0.5
        
        # Calculate Wcorr.
        self.Wcorr = np.identity(self.d)
        for i in range(self.d):
            for j in range(i+1,self.d):
                self.Wcorr[i,j] = abs(w_inner(self.TS_comps[:,i], self.TS_comps[:,j]) * F_wnorms[i] * F_wnorms[j])
                self.Wcorr[j,i] = self.Wcorr[i,j]
    
    