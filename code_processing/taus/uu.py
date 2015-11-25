from pylab import *
from numpy import *

def acc_i(v, i):
	sum = 0.0
	for j in range(i+1):
		sum += v[j]
	return sum

def acc(v):
	n = size(v)
	acc_v = zeros(n)
	for i in range(n):
		acc_v[i] = acc_i(v, i)
	return acc_v
