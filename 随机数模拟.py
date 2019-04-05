import random
class Generate():
    def __init__(self):
        print()
    
    def uniform(self,n=1,seed=3):
        a = 7 ** 5
        m = 2 ** 31 - 1
        c = 4
        x = [seed]
        for i in range(n):
            x.append(a * x[-1] + c - (a * x[-1] + c) // m * m)
        x = [i / m for i in x]
        if len(x) == 2:
            return x[1]
        else:
            return x[1:]
    
    def discrete(self,prob,value):
        # U = self.uniform(100)[-1]
        if sum(prob) != 1:
            return print('错误的概率分布')
        U = random.uniform(0,1)
        sum_p = prob[0]
        for i in range(len(prob)):
            if U < sum_p:
                return value[i]
            else:
                sum_p += prob[i + 1]
    
    def permutation(self,n,r=None):
        arangement = [i for i in range(1,n+1)]
        if not r:
            r = n
        k = n
        # flag = False
        # if r > n // 2:
        #     r = n - r
        #     flag = True
        while k > n - r:
            index = int(k * random.uniform(0,1))
            arangement[k-1],arangement[index] = arangement[index],arangement[k-1]
            k -= 1
            print('$')
        # if flag:
        #     return arangement[:n-r+1]
        # else:
        return arangement[n-r::]
    
    def permutation2(self,n,r=None):
        arangement = [i for i in range(1,n+1)]
        if not r:
            r = n
        d = {random.uniform(0,1):arangement[i] for i in range(n)}
        return [d[i] for i in sorted(d)][:r]
    
    def permutation3(self,n):
        if n == 1:
            return [1]
        else:
            arangement = self.permutation3(n-1)
            arangement += [n]
            a = int(random.uniform(0,1) * n)
            arangement[-1],arangement[a] = arangement[a],arangement[-1]
            return arangement

    
    def bernoulli(self,p):
        if random.uniform(0,1) < p:
            return 1
        else:
            return 0
    
    def binomial(self,n,p):
        c = p / (1 - p)
        pr = (1 - p) ** n
        F = pr
        U = random.uniform(0,1)
        for i in range(n + 1):
            if U < F:
                return i
            else:
                pr *= (c * (n - i) / (i + 1))
                F += pr
    def reject(self):
        c = 11 / 6
        p = [i / 36 for i in list(range(1,6)) + list(range(6,0,-1))]
        q = [i / 11 for i in [1]*11]
        while 1:
            Y = int(11 * random.uniform(0,1)) + 2
            U = random.uniform(0,1)
            if U < p[Y-2] / q[Y-2] / c:
                # print(Y)
                return Y

def pr():
    import math
    U = 1
    for i in range(100):
        U *= random.uniform(0,1)
        if U < math.e ** (-3):
            return i
if __name__=='__main__':
    s = 0
    n = 1000
    generate = Generate()
    for i in range(n):
        s += pr()
    s /= n
    print(s)




