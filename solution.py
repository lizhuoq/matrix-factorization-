import numpy as np
import factorization as fc

np.set_printoptions(suppress = True,precision = 3)

def sub_matrix(A,r,c,verbose = False):
    m,n = A.shape
    SM = [[A[x,y] for y in range(n) if y != c] for x in range(m) if x != r]
    SM = np.array(SM)
    if verbose == True:
        print(SM)

    return SM


def det(A,verbose = False):
    A = np.array(A,dtype = np.float64)
    m,n = A.shape
    if m == 1 & n == 1:
        det_value = A[0,0]
    else:
        det_value = 0
        for i in range(n):
            det_value += ((-1) ** (i + 2)) * A[0,i] * det(sub_matrix(A,0,i))
    if verbose == True:
        print('det(A): \n',det_value)

    return det_value


def ly_pb(L,b,P = None,verbose = False):
    m = L.shape[0]
    if P is None:
        P = np.identity(m)
    b = np.array(b,dtype = np.float64)
    b = np.dot(P,b)
    y = np.zeros_like(b)
    for i in range(m):
        t = 0
        for j in range(i):
            t += L[i,j] * y[j,0]
        y[i,0] = (b[i,0] - t) / L[i,i]
    if verbose == True:
        print('y: \n',y)

    return y


def ux_y(U,y):
    m = U.shape[0]
    x = np.zeros_like(y)
    for i in range(m - 1,-1,-1):
        t = 0
        for j in range(i + 1,m):
            t += U[i,j] * x[j,0]
        t = y[i,0] - t
        x[i,0] = t / U[i,i]
    print('x: \n',x)

    return x


def qr_resolution(Q,R,b):
    R = R[:np.linalg.matrix_rank(R),:]
    c = np.dot(Q.T,b)
    c = c[:np.linalg.matrix_rank(R),:]
    x = ux_y(R,c)

    return x


def householder_givens_resolution(P,T,b):
    R = T[:np.linalg.matrix_rank(T),:]
    c = np.dot(P,b)
    c = c[:np.linalg.matrix_rank(T),:]
    x = ux_y(R,c)

    return x


def solution(A,method,b):
    np.set_printoptions(precision=3, suppress=True)
    if A.shape[0] == A.shape[1]:
        det_a = det(A,verbose = True)
    else:
        print('A不是方阵，没有特征值')
    if method == 'LU':
        P,L,U = fc.factorization(A,method)
        y = ly_pb(L,b,P)
        x = ux_y(U,y)

        return P,L,U,x
    elif method == 'GS':
        Q,R = fc.factorization(A,method)
        x = qr_resolution(Q,R,b)

        return Q,R,x
    elif method == 'Givens' or method =='Householder':
        P,T = fc.factorization(A,method)
        x = householder_givens_resolution(P,T,b)

        return P,T,x