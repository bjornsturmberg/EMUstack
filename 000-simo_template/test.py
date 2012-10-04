import numpy as np

file_name_1 = "wl0001_Tnet.txt" 
data   = np.loadtxt(file_name_1)
num_1  = max(data[:,0])
num_2  = max(data[:,1])
matrix = np.mat(data[:,2] + data[:,3]*(0+1j))
matrix1 = np.reshape(matrix, (num_2, num_1))
file_name_2 = "TLambda.txt" 
data   = np.loadtxt(file_name_2)
num_1  = max(data[:,0])
num_2  = max(data[:,1])
matrix = np.mat(data[:,2] + data[:,3]*(0+1j))
matrix2 = np.reshape(matrix, (num_2, num_1))

dif = matrix1 - matrix2
print dif.max()