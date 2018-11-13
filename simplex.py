import numpy as np
import scipy as sp
import itertools
from scipy.optimize import linprog
from time import time

class Simplex():

    def __init__(self,A=[],b=[],c=[]):
        self.A = A
        self.b = b
        self.c = c
        # self.input()
        self.show()
        self.m = len(self.A)
        self.n = np.size(self.A,1)
        self.reduce_redundant()
        self.feasible = True
        self.iteration = 0
        self.unbounded = False
        self.error = 0

    def input(self):
        print('-'*9,'\t请输入消耗系数矩阵 A\t','-'*9)
        while True:
            temp = input()
            if temp:
                if ',' in temp:
                    self.A.append(list(map(eval,temp.split(','))))
                elif ' ' in temp:
                    self.A.append(list(map(eval,temp.split(' '))))
            else:
                break

        print('-'*9,'\t请输入资源约束向量 b\t','-'*9)
        temp = input()
        if temp:
            if ',' in temp:
                self.b = list(map(eval,temp.split(',')))
            elif ' ' in temp:
                self.b = list(map(eval,temp.split(' ')))

        print('-'*9,'\t请输入目标函数向量 c\t','-'*9)
        temp = input()
        if temp:
            if ',' in temp:
                self.c = list(map(eval,temp.split(',')))
            elif ' ' in temp:
                self.c = list(map(eval,temp.split(' ')))      

        self.A = np.array(self.A)
        self.b = np.array(self.b)
        self.c = np.array(self.c)

    def solve(self):
        self.initialize()
        if self.feasible:
            while self.feasible:
                if self.iteration > 20:
                    break
                if self.optimal_test():
                    break
                else:
                    self.get_basic_direction()
                    self.get_stepsize()
                    if self.unbounded:
                        print('无界解')
                        break
                    else:
                        self.update()
            else:
                print(self.x)
        else:
            print('无解')

    def reduce_redundant(self):
        if np.linalg.matrix_rank(np.c_[self.A,self.b]) < self.m:
            for i in itertools.combinations(range(self.m),np.linalg.matrix_rank(np.c_[self.A,self.b])):
                if np.linalg.matrix_rank(np.c_[self.A,self.b][i,:]) == np.linalg.matrix_rank(np.c_[self.A,self.b]):
                    self.A = self.A[i,:]
                    self.b = np.array([self.b[j] for j in i])
                    self.m = len(self.A)
                    break

    def initialize(self,method='two phase'):
        # try:
        #     if np.linalg.det(self.A[-self.m:]) == 0:
        #         A_two_phase = np.c_[self.A,np.eye(self.m)]
        #         c_two_phase = np.r_(self.c,np.ones(self.m))
        #         LP = Simplex(A_two_phase,self.b,c_two_phase)
        #         LP.solve()
        #         if all(LP.x[-LP.m:-1]==0):
        #             self.x = LP.x[:LP.n]
        #         else:
        #             self.feasible = False
        #     else:
        #         self.B = self.A[-self.m:]
        #         self.B_inverse = np.linalg.inv(self.B)
        #         self.basic_indices = range(self.n - self.m,self.n)
        #         self.nonbasic_indices = range(self.n - self.m)
        #         self.x = np.zeros(self.n)
        #         self.x[self.basic_indices] = np.dot(self.B_inverse,self.b)
        #         self.objective = np.dot(self.c,self.x)
        # except:
        for i in itertools.combinations(range(self.n),self.m):
            if np.linalg.det(self.A[:,i]) != 0:
                self.basic_indices = list(i)
                self.nonbasic_indices = [i for i in range(self.n) if i not in self.basic_indices]
                self.B = self.A[:,self.basic_indices]
                if np.all(np.linalg.solve(self.B,self.b) >= [0]):
                    self.x = np.zeros(self.n)
                    self.x[self.basic_indices] = np.linalg.solve(self.B,self.b)
                    self.objective = np.dot(self.c,self.x)
                    print('='*9,'\t\t初始化\t\t','='*9)
                    print('基\t\t：',end='')
                    [print(i,'\t',end='') for i in self.basic_indices]
                    print('\n决策变量\t：',end='')
                    [print(i,'\t',end='') for i in np.round(self.x,2)]
                    print('\n目标函数\t：',np.round(self.objective,2))
                    break
        else:
            self.feasible = False
    
    def optimal_test(self):
        self.reduced_cost = self.c - np.dot(np.dot(self.c[self.basic_indices],np.linalg.inv(self.B)),self.A)
        if min(self.reduced_cost) < -self.error:
            self.into_basis = np.argmin(self.reduced_cost)
            return False
        elif 0 in self.reduced_cost[self.nonbasic_indices]:
            print('*'*9,'\t\t无穷解\t\t','*'*9)
            self.find_more_solutions()
            return True
        else:
            print('*'*9,'\t\t最优解\t\t','*'*9)
            return True
    
    def find_more_solutions(self):
        test = np.where(self.reduced_cost[self.nonbasic_indices]==0)[0]
        test = np.array(self.nonbasic_indices)[test]
        base = [self.basic_indices]
        for i in test:
            basis = self.basic_indices[:]
            self.into_basis = i
            self.get_basic_direction()
            self.get_stepsize()
            basis.remove(self.out_basis)
            basis.append(self.into_basis)
            base.append(basis)
        solutions = []
        for i in base:
            x = np.zeros(self.n)
            x[i] = np.dot(np.linalg.inv(self.A[:,i]),self.b)
            solutions.append(x)
        print('最优解是下列向量的凸组合')
        for i in solutions:
            print('x = ',end='')
            [print(round(j,2),'\t',end='') for j in i]
            print('')

    def get_basic_direction(self):
        self.basic_direction = np.zeros(self.n)
        self.basic_direction[self.basic_indices] = - np.dot(np.linalg.inv(self.B),self.A[:,self.into_basis])
        self.basic_direction[self.into_basis] = 1

    def get_stepsize(self):
        self.filter_indices = [i for i,j in enumerate(self.basic_direction) if j < self.error]
        if self.filter_indices:
            self.stepsize = min(- self.x[self.filter_indices] / self.basic_direction[self.filter_indices])
            self.out_basis = self.filter_indices[np.argmin(- self.x[self.filter_indices] / self.basic_direction[self.filter_indices])]
        else:
            self.unbounded = True
    
    def update(self):
        self.iteration += 1
        print('='*9,'\t\t迭代',self.iteration,'\t\t','='*9)
        self.basic_indices.remove(self.out_basis)
        self.basic_indices.append(self.into_basis)
        self.nonbasic_indices = [i for i in range(self.n) if i not in self.basic_indices]
        self.B = self.A[:,self.basic_indices]
        self.x += self.stepsize * self.basic_direction
        self.objective = np.dot(self.c,self.x)
        print('入基变量：变量\t',self.into_basis)
        print('出基变量：变量\t',self.out_basis)
        print('基\t\t：',end='')
        [print(i,'\t',end='') for i in self.basic_indices]
        print('\n前进方向：',end='')
        [print(i,'\t',end='') for i in np.round(self.basic_direction,2)]
        print('\n前进值\t：',round(self.stepsize,2))
        print('决策变量：',end='')
        [print(i,'\t',end='') for i in np.round(self.x,2)]
        print('\n目标函数：',np.round(self.objective,2))
    
    def show(self):
        print('~'*9,'\t\t标准型\t\t','~'*9)
        for i in self.c:
            print(i,end='\t')
        print('')
        for i,k in enumerate(self.A):
            for j in k:
                print(j,end='\t')
            print(self.b[i])

if __name__=='__main__':
    # problem 1 报错
    # A = np.array([[-1,1,1,0,0],[-2,1,0,1,0],[4,1,0,0,1]])
    # b = np.array([3,2,16])
    # c = np.array([3,-2,0,0,0])
    # problem 2 无
    # A = np.array([[1,2,2,1,0,0],[2,1,2,0,1,0],[2,2,1,0,0,1]])
    # b = np.array([20,20,20])
    # c = np.array([-10,-12,-12,0,0,0])
    # problem 3 无
    # A = np.array([[1,2,3,0],[-1,2,6,0,],[0,4,9,0],[0,0,3,1]])
    # b = np.array([3,2,5,1])
    # c = np.array([1,1,1,0])
    # problem 4 无
    # A = np.array([[1,1,1,1,0,0,0],[1,0,0,0,1,0,0],[0,0,1,0,0,1,0],[0,3,1,0,0,0,1]])
    # b = np.array([4,2,3,6])
    # c = np.array([-1,-14,-6,0,0,0,0])
    # problem 5 冗余约束
    A = np.array([[3,5,1,0],[3,1,0,1],[6,6,1,1]])
    b = np.array([15,12,27])
    c = np.array([-2,-1,0,0])
    # problem 6 无界解
    # A = np.array([[1,1,1,-1,0,0],[-2,0,1,0,-1,0],[0,2,-1,0,0,-1]])
    # b = np.array([6,2,0])
    # c = np.array([-2,1,-2,0,0,0])
    # problem 7 循环？
    # A = np.array([[1,0,0,1/4,-8,-1,9],[0,1,0,1/2,-12,-1/2,3],[0,0,1,0,1,1,0]])
    # b = np.array([0,0,1])
    # c = np.array([0,0,0,-3/4,20,-1/2,6])
    # problem 8 无穷解
    # A = np.array([[1,3,1,-1,0,0],[1,2,-1,0,-1,0],[1,0,1,0,0,1]])
    # b = np.array([4,6,12])
    # c = np.array([1,2,-1,0,0,0])
    # problem 9 无
    # A = np.array([[1,1,0,0,1,0,0,0],[0,0,1,1,0,1,0,0],[1,0,1,0,0,0,-1,0],[0,1,0,1,0,0,0,-1]])
    # b = np.array([48,60,36,72])
    # c = np.array([6,8,4,3,0,0,0,0])
    # problem 10 无
    # A = np.array([[1,1,1,0,0],[1,2,0,1,0],[1,3,0,0,1]])
    # b = np.array([1,1,1])
    # c = np.array([-1,-4,0,0,0])
    # problem 11 冗余约束
    # A = np.array([[1,0,1],[1,1,0],[1,1,0]])
    # b = np.array([2,2,2])
    # c = np.array([4,2,1])
    begin = time()
    LP = Simplex(A,b,c)
    LP.solve()
    end = time()
    duration1 = end - begin
    begin = time()
    x = linprog(c,A_eq=A,b_eq=b,method='simplex')
    print(x)
    end = time()
    duration2 = end - begin
    print(duration1//duration2)