# -*- coding: UTF-8 -*-
import numpy as np
#from pysnptools.snpreader import Bed
from scipy import linalg
import pandas as pd
import sys

grm=open(sys.argv[1])
ginv = open(sys.argv[2], 'w')
#id=open('ID.txt','w')

matrix=[]
for i in grm.readlines()[0:]:
	f=i.split(',')
	matrix.append('\t'.join(f[0:])+'\n')

       
lines = matrix          
n=len(lines)
A = np.zeros((n,n),dtype=float)    
A_row = 0                       
for line in lines:              
    list = line.strip('\n').split()      
    A[A_row:] = list[0:n]                   
    A_row+=1                                



d = np.diag(A) + 0.001
#替换原来的G矩阵对角线
np.fill_diagonal(A, d)
#求逆
Ginv = linalg.inv(A)
np.savetxt(sys.argv[2],Ginv,delimiter=" ")