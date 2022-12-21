import pprint
import numpy as np
import scipy.linalg as la
import matplotlib.pylab as plt

def print_sparsity(A):
    n,m = A.shape
    for i in range(n):
        for j in range(m):
            if abs(A[i][j]) != 0:
                print("x",end=' ')
            else:
                print(" ",end=' ')
        print()

def run(n, off):
    # Matrix
    d = 2.*np.ones(n)
    s = -1.*np.ones(n-off)
    A= np.diag(d) + np.diag(s,off) + np.diag(s,-off)

    # Decomposition and inverses
    P,L,U = la.lu(A)
    Uinv = np.linalg.inv(U)
    Linv = np.linalg.inv(L)
    Ainv = np.linalg.inv(A)


    # Plot sparsity
    plt.subplot(2,3,1)
    plt.spy(A, markersize=1)
    plt.xticks([])
    plt.yticks([])
    plt.title('A')
    plt.subplot(2,3,2)
    plt.spy(L, markersize=1)
    plt.xticks([])
    plt.yticks([])
    plt.title('L')
    plt.subplot(2,3,3)
    plt.spy(U, markersize=1)
    plt.xticks([])
    plt.yticks([])
    plt.title('U')

    plt.subplot(2,3,4)
    plt.spy(Ainv, markersize=1)
    plt.xticks([])
    plt.yticks([])
    plt.title('Ainv')
    plt.subplot(2,3,5)
    plt.spy(Linv, markersize=1)
    plt.xticks([])
    plt.yticks([])
    plt.title('Linv')
    plt.subplot(2,3,6)
    plt.spy(Uinv, markersize=1)
    plt.xticks([])
    plt.yticks([])
    plt.title('Uinv')

    plt.show()

    # # Output
    # print("A:")
    # print_sparsity(A)
    # # print("P:")
    # # pprint.pprint(P)
    # print("L:")
    # print_sparsity(L)
    # print("U:")
    # print_sparsity(U)
    # print("Uinv:")
    # print_sparsity(Uinv)
    # print("Linv:")
    # print_sparsity(Linv)
    # print("Ainv:")
    # print_sparsity(Ainv)

# Run Program
n = 20
off = 2
run(n, off)

n = 20
off = 4
run(n, off)
