import numpy as np
import factorization as fc
import solution as slt

np.set_printoptions(suppress = True,precision = 3)

#Householder reduction solution
A = np.array([[4,-3,4],
              [2,-14,-3],
              [-2,14,0],
              [1,-7,15]])
b = np.array([[5],
              [-15],
              [0],
              [30]])
slt.solution(A,'Householder',b)

#partial pivot LU reduction solution
A = np.array([[1,2,-3,4],
              [4,8,12,-8],
              [2,3,2,1],
              [-3,-1,1,-4]])
b = np.array([[3],
              [60],
              [1],
              [5]])
slt.solution(A,'LU',b)

#Gram-Schmidt reduction solution
A = np.array([[2,4],
              [3,-5],
              [1,2],
              [2,1]])
b = np.array([[11],
              [3],
              [6],
              [7]])
slt.solution(A,'GS',b)

#URV reduction
A = np.array([[-4,-2,-4,-2],
              [2,-2,2,1],
              [-4,1,-4,-2]])
fc.factorization(A,'URV')

#Givens reduction solution
A = np.array([[4,-3,4],
              [2,-14,-3],
              [-2,14,0],
              [1,-7,15]])
b = np.array([[5],
              [-15],
              [0],
              [30]])
slt.solution(A,'Givens',b)