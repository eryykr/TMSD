

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from hampel import hampel
import pandas as pd 

#loading fluorescence data
data_file = 'example.csv'
df = pd.read_csv(f'{data_file}')

#getting time, converting to seconds
time = df.iloc[:,0][1:-1] #<--- first number refers to index of column, second and third are indices of data within the column
time = np.array([float(i.split()[0]) for i in time])
time = time*60

#getting fluorescence
flrs = df.iloc[:,1][1:-1] #<--- first number refers to index of column, second and third are indices of data within the column
flrs = np.array([float(i.split()[0]) for i in flrs])



#filtering to get rid of spike at t = 0 (change parameters to get desired result)
result = hampel(flrs, window_size=100, n_sigma=1.0)
flrs = result.filtered_data




#concentration as a function of time (second order kinetics)
def conc_of_t(t, t_0, A_0, A_inf, k):
	out = []
	for t_i in t:
		if t_i <= t_0:
			out.append(A_0)
		else:
			out.append( ( A_0 +  (A_inf-A_0)*(1-1/(1 + abs(k)*(t_i-t_0)) ) ) )
	return out


#fitting (p0 are initial guesses of parameters)
popt, pcov = curve_fit(conc_of_t, time, flrs, p0 = [10,10,10,10])

#fitted curve
flrs_fit = conc_of_t(time,popt[0], popt[1], popt[2], popt[3])




#plotting experimental data and fit
plt.plot(time, np.array(flrs), label = 'Data', linewidth=1)
plt.plot(time, np.array(flrs_fit), label = 'Fit', color='black', linestyle='dashed', linewidth=2)
plt.show()


