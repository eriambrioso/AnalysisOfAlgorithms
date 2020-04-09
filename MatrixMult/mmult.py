# An implementation of the Strassen methood for matrix multiplication
# https://en.wikipedia.org/wiki/Strassen_algorithm
"""
Strassen method for matrix multiplicaton a square matrix 2^n x 2^n matrix into 2^(n-1) square submatices
"""
import random as rd

# This functions prints out a matrix line by line so that it is easy to read.
def printMat(M):
	s = [[str(e) for e in row] for row in M]
	lens = [max(map(len, col)) for col in zip(*s)]
	fmt = ' '.join('{{:{}}}'.format(x) for x in lens)
	table = [fmt.format(*row) for row in s]
	print('\n'.join(table))

def make_zeros(n_rows: int, n_columns: int):
    matrix = []
    for i in range(n_rows):
        matrix.append([0] * n_columns)
    return matrix

# Matrix multiplication using the standard algorithm and list comprehension.
def my_matmul(X, Y):
	result = [[sum(a*b for a,b in zip(X_row, Y_col)) for Y_col in zip(*Y)] for X_row in X]
	return result

# Matrix multiplication using the standard algorithm and for loops
def standard(A, B, solution):
	size = len(A)
	for i in range(size):
		for j in range(size):
			for k in range(size):
				pass
				#solution[i][j] += A[i][k] * B[k][j]
	return solution

"""
We define the functions we need for the Strassen algorithm.
"""
# Build up a matrix to the nearest power of 2
def buildPowerTwo(M):
	result = M

	powersOfTwo = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
	size = len(M)
	if len(M) not in powersOfTwo:
		result = M
		list = [i-size for i in powersOfTwo if i-size > 0]
		diff = abs(min(list))
		size = size + diff
		# build matrix up to new size
		for j in range(diff):
			for row in result:
				row.append(0)
		for k in range(diff):
			result.append([0] * size)

	return result
# Adds two matrices
def add(A,B):
	n = len(A)
	return [[A[i][j] + B[i][j] for j in range(n)] for i in range(n)]

# Subtracts two matrices
def subtract(A,B):
	n=len(A)
	return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]

#A matrix chopping function.  Cuts a 2^n x 2^n matrix into four square submatrices.
def chop_mat(M):
	n = len(M)
	A00 = [[M[i][j] for j in range(n//2)] for i in range(n//2)]
	A01 = [[M[i][j] for j in range(n//2)] for i in range(n//2,n)]
	A10 = [[M[i][j] for j in range(n//2,n)] for i in range(n//2)]
	A11 = [[M[i][j] for j in range(n//2,n)] for i in range(n//2,n)]
	if len(A11) == 1:
		return [A00, A10, A01, A11]
	else:
		return [A00, A10, A01, A11]

# This the main part of the algorithm
def strassen_mul(A,B):
	size = len(A)
	if size == 2:
		M1 = (A[0][0] + A[1][1])*(B[0][0] + B[1][1])
		M2 = (A[1][0] + A[1][1])*B[0][0]
		M3 = A[0][0]*(B[0][1] - B[1][1])
		M4 = A[1][1]*(B[1][0] - B[0][0])
		M5 = (A[0][0] + A[0][1])*B[1][1]
		M6 = (A[1][0] - A[0][0])*(B[0][0] + B[0][1])
		M7 = (A[0][1] - A[1][1])*(B[1][0] + B[1][1])
		C00 = M1 + M4 - M5 + M7
		C01 = M3 + M5
		C10 = M2 + M4
		C11 = M1 - M2 + M3 + M6
		return [[C00, C01], [C10, C11]]
	if size > 2:
		P1 = chop_mat(A)
		P2 = chop_mat(B)
		A00 = P1[0]
		A01 = P1[1]
		A10 = P1[2]
		A11 = P1[3]
		B00 = P2[0]
		B01 = P2[1]
		B10 = P2[2]
		B11 = P2[3]
		p1_1 = add(A00, A11)
		p1_2 = add(B00, B11)
		p2 = add(A10, A11)
		p3 = subtract(B01, B11)
		p4 = subtract(B10, B00)
		p5 = add(A00, A01)
		p6_1 = subtract(A10,A00)
		p6_2 = add(B00, B01)
		p7_1 = subtract(A01, A11)
		p7_2 = add(B10, B11)
		M1 = strassen_mul(p1_1, p1_2)
		M2 = strassen_mul(p2, B00)
		M3 = strassen_mul(A00, p3)
		M4 = strassen_mul(A11, p4)
		M5 = strassen_mul(p5, B11)
		M6 = strassen_mul(p6_1, p6_2)
		M7 = strassen_mul(p7_1, p7_2)
		C00 = subtract(add(add(M1, M4), M7), M5)
		C01 = add(M3, M5)
		C10 = add(M2, M4)
		C11 = subtract(add(add(M1,M3), M6), M2)
		C = [[0 for j in range(size)] for i in range(size)]
		new_size = size//2
		for i in range(new_size):
			for j in range(new_size):
				C[i][j] = C00[i][j]
				C[i][j + new_size] = C01[i][j]
				C[i + new_size][j] = C10[i][j]
				C[i + new_size][j + new_size] = C11[i][j]
		return C
	else:
		print('Error in strassen_mul')

"""
The shows how to use a context manager to time a block of code.
The code idea comes from Beazley and Jones's Python Cookbook (pp. 588-599)
"""
import time
import numpy as np
import matplotlib.pylab as plt


# Timing the algorithms
strassen_timing = {}
standard_timing = {}
comp_timing = {}
numpy_timing = {}

for n in range(2, 250, 10):
	print(n)
	E = [[rd.randint(1, 3) for j in range(n)] for i in range(n)]
	F = [[rd.randint(1, 3) for j in range(n)] for i in range(n)]
	E = buildPowerTwo(E)
	F = buildPowerTwo(F)
	
	
	# time strassen
	start = time.perf_counter()
	strassen_mul(E, F)
	end = time.perf_counter()
	strassen_timing[n] = end-start
	
	# time standard
	start = time.perf_counter()
	standard(E, F, make_zeros(len(E), len(E)))
	end = time.perf_counter()
	standard_timing[n] = end-start
	"""
	# time list comp
	start = time.perf_counter()
	my_matmul(E, F)
	end = time.perf_counter()
	comp_timing[n] = end-start
	"""

	# time numpy
	start = time.perf_counter()
	E = np.array(E)
	F = np.array(F)
	solution = E.dot(F)
	end = time.perf_counter()
	numpy_timing[n] = end-start
	

# Plotting the time complexity
ax = plt.subplots()
#plt.plot(list(comp_timing.keys()), list(comp_timing.values()), label = "Comp")
plt.plot(list(standard_timing.keys()), list(standard_timing.values()), label = 'Standard')
plt.plot(list(numpy_timing.keys()), list(numpy_timing.values()), label = 'Numpy')
plt.plot(list(strassen_timing.keys()), list(strassen_timing.values()), label = 'Strassen')
plt.title("Matrix Multiplication time complexity")
plt.ylabel("Time")
plt.xlabel("n, size of input")
plt.legend()
plt.show()