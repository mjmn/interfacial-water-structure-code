#!/home/mniesen/Software/Anaconda/anaconda3/bin/python

# Aggregate local water structure data into a single file

import numpy as np

for i in range(1,100001,100):
  try:
    fn = "t%i/local_wstructure_histogram.txt" % i
    data = np.loadtxt(fn,dtype=int,usecols=(3))
    try:
      sum_data
    except NameError: # sum_data doesn't excist yet
      sum_data = np.zeros(np.shape(data),dtype=int)
    sum_data += data
  except:
    print('Failed to load data from t',i)

np.savetxt("Aggregate_local_wstructure_histogram.txt",sum_data)
