import numpy as np
from random import random
import sympy
import math
def calculation(alpha,beta,m):
    # m = 10 # 精度
    a = np.zeros((m + 2,m + 2))
    b = np.zeros((m + 2,m + 2))
    i = j = 1
    shift = True
    a[i,i] = beta[i]
    # print('a',i,i,a[i,i])
    while True:
        if i == 1 and shift:
            j += 1
            if j > m:
                break
            i = j
            a[i,i] = beta[i]
            shift = True
            # print('a',i,i,a[i,i])
        else:
            if shift:
                i -= 1
                b[i,j] = np.dot(alpha[:(j - i)],a[(i + 1):(j + 1),j])
                shift = False
                # print('b',i,j,b[i,j])
            else:
                temp = [b[x,j - i + x] for x in range(1,i + 1)]
                a[i,j] = np.dot(beta[(i - 1)::-1],temp)
                shift = True
                # print('a',i,j,a[i,j])
    a = a[1::,1::]
    b = b[1::,1::]
    # np.delete(a,0,axis=0)
    # np.delete(a,0,axis=1)
    # np.delete(b,0,axis=0)
    # np.delete(b,0,axis=1)
    return a,b

def get_alpha(m):
    x = sympy.symbols('x')
    f = (sympy.exp(-x) + 0.5 * sympy.exp(-0.5 * x)) / 2
    alpha = [f]
    while len(alpha) < m:
        alpha.append(sympy.diff(alpha[-1]))
    alpha = [i.subs(x,0) for i in alpha]
    return alpha

def get_beta(m):
    t = sympy.symbols('t')
    f = 1 / (1 - t)
    beta = [f]
    while len(beta) < m:
        beta.append(sympy.diff(beta[-1]))
    beta = [j.subs(t,0) / math.factorial(i) for i,j in enumerate(beta)]
    return beta

def get_gamma(m):
    return get_beta(m)
m = 40
alpha = get_alpha(m + 1)
beta = get_beta(m + 1)
a,b = calculation(alpha,beta,m)
delay = b[0,:].sum()
print(delay)