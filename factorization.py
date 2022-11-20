import numpy as np

np.set_printoptions(suppress = True,precision = 3)

def pp_lu(A):
    A = np.array(A,dtype = np.float64)
    L = np.zeros_like(A)
    n = A.shape[0]
    P = np.identity(n)
    for i in range(n - 1):
        xrow = i + np.argmax(abs(A[i:,i]))
        if xrow != i:
            A[[i,xrow],:] = A[[xrow,i],:]
            L[[i,xrow],:] = L[[xrow,i],:]
            P[[i,xrow],:] = P[[xrow,i],:]
        L[i,i] = 1
        for j in range(i + 1,n):
            L[j,i] = A[j,i] / A[i,i]
            A[j,:] = A[j,:] - L[j,i] * A[i,:]
    L[i + 1,i + 1] = 1
    U = A
    print('P: \n',P)
    print('L: \n',L)
    print('U: \n',U)
    return P,L,U


def qr_gs(A):
    A = np.array(A,dtype = np.float64)
    Q = np.zeros_like(A)
    m,n = A.shape
    R = np.zeros((n,n))
    for i in range(n):
        a = A.T[i]
        for j in range(0,i):
            r = np.dot(Q[:,j],a)
            R[j,i] = r
            a -= np.dot(r,Q[:,j])
        anorm = np.linalg.norm(a)
        q = a / anorm
        R[i,i] = anorm
        Q[:,i] = q
    print('Q: \n',Q)
    print('R: \n',R)

    return Q,R


def householder(A,verbose = False):
    A = np.array(A,dtype = np.float64)
    T = np.copy(A)
    m,n = A.shape
    P = np.identity(m)
    for i in range(m - 1):
        if i > n - 1:
            break
        a = T[i:,i]
        e = np.zeros_like(a)
        e[0] = np.linalg.norm(a)
        u = a - e
        R = np.identity(m)
        R[i:,i:] -= 2.0 * np.outer(u,u) / np.inner(u,u)
        P = np.dot(R,P)
        T = np.dot(R,T)
    if verbose == True:
        print('P: \n',P)
        print('T: \n',T)

    return P,T


def givens(A):
    A = np.array(A,dtype = np.float64)
    m,n = A.shape
    P = np.identity(m)
    T = np.copy(A)
    r,c = np.tril_indices(m,-1,n)
    for i,j in zip(r,c):
        if T[i,j] != 0:
            s = np.hypot(T[j,j],T[i,j])
            cos = T[j,j] / s
            sin = T[i,j] / s
            R = np.identity(m)
            R[j,j] = cos
            R[i,i] = cos
            R[i,j] = -sin
            R[j,i] = sin
            T = np.dot(R,T)
            P = np.dot(R,P)
    print('P: \n',P)
    print('T: \n',T)

    return P,T


def urv(A):
    A = np.array(A,dtype = np.float64)
    P,B = householder(A)
    B = B[:np.linalg.matrix_rank(B),:].T
    Q,T = householder(B)
    T = T[:np.linalg.matrix_rank(T),:]
    V = Q.T
    U = P.T
    R = np.zeros_like(A,dtype = np.float64)
    R[:np.linalg.matrix_rank(T),:np.linalg.matrix_rank(T)] = T.T
    print('U: \n',U)
    print('R: \n',R)
    print('V: \n',V)

    return U,R,V


def factorization(A,method):
    np.set_printoptions(suppress = True,precision = 3)
    verbose = True
    if method == 'LU':
        P,L,U = pp_lu(A)

        return P,L,U
    elif method == 'Givens':
        P,T = givens(A)

        return P,T
    elif method == 'Householder':
        P,T = householder(A,verbose)

        return P,T
    elif method == 'GS':
        Q,R = qr_gs(A)

        return Q,R
    elif method == 'URV':
        U,R,V = urv(A)

        return U,R,V