import numpy as np
import scipy as sp
import itertools

class simplex():

    def __init__(self):
        self.input()
        self.m = len(self.A)
        self.n = np.size(self.A,1)
        self.question_simplify()
        self.x = np.zeros(self.n)
        self.objective = np.dot(self.c,self.x)
        self.B = []
        self.feasible = True
        self.reduced_cost = []
        self.basic_indices = []
        self.into_basis = None
        self.out_basis = None
        self.iteration = 0
        self.unbounded = False

    def input(self):
        self.A = []
        print('-'*9,'\t请输入消耗系数矩阵 A\t','-'*9)
        while True:
            temp = input()
            if temp:
                if ',' in temp:
                    self.A.append(list(map(eval,temp.split(','))))
                elif ' ' in temp:
                    self.A.append(list(map(eval,temp.split(','))))
            else:
                break

        print('-'*9,'\t请输入资源约束向量 b\t','-'*9)
        temp = input()
        if temp:
            if ',' in temp:
                self.b = list(map(eval,temp.split(',')))
            elif ' ' in temp:
                self.b = list(map(eval,temp.split(',')))

        print('-'*9,'\t请输入目标函数向量 c\t','-'*9)
        temp = input()
        if temp:
            if ',' in temp:
                self.c = list(map(eval,temp.split(',')))
            elif ' ' in temp:
                self.c = list(map(eval,temp.split(',')))      

        self.A = np.array(self.A)
        self.b = np.array(self.b)
        self.c = np.array(self.c)
        self.show() 

    def solve(self):
        self.iteration_initial()
        if self.feasible:
            while self.feasible:
                self.get_reduced_cost()
                if self.optimal_test():
                    print('*'*9,'\t\t当前为最优解\t\t','*'*9)
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
            print('可行集为空!!!')

    def question_simplify(self):
        if np.linalg.matrix_rank(np.c_[self.A,self.b]) < self.m:
            for i in itertools.combinations(range(self.m),np.linalg.matrix_rank(np.c_[self.A,self.b])):
                if np.linalg.matrix_rank(np.c_[self.A,b][i,:]) == np.linalg.matrix_rank(np.c_[self.A,self.b]):
                    self.A = self.A[i,:]
                    self.b = np.array([self.b[j] for j in i])
                    self.m = len(self.A)
                    break

    def iteration_initial(self):
        for i in itertools.combinations(range(self.n),self.m):
            if np.linalg.det(self.A[:,i]) != 0:
                self.basic_indices = list(i)
                self.B = self.A[:,self.basic_indices]
                if np.all(np.linalg.solve(self.B,self.b) > [0]):
                    self.x[self.basic_indices] = np.linalg.solve(self.B,self.b)
                    self.objective = np.dot(self.c,self.x)
                    print('='*9,'\t\t初始化\t\t','='*9)
                    print('基：\t\t',self.basic_indices)
                    print('决策变量值：',np.round(self.x,2))
                    print('目标函数值：',np.round(self.objective,2))
                    break
        else:
            self.feasible = False
    
    def optimal_test(self):
        if min(self.reduced_cost) < -1e-12:
            self.into_basis = np.argmin(self.reduced_cost)
            return False
        else:
            return True
    
    def get_reduced_cost(self):
        self.reduced_cost = self.c - np.dot(np.dot(self.c[self.basic_indices],np.linalg.inv(self.B)),self.A)
    
    def get_basic_direction(self):
        self.basic_direction = np.zeros(self.n)
        self.basic_direction[self.basic_indices] = - np.dot(np.linalg.inv(self.B),self.A[:,self.into_basis])
        self.basic_direction[self.into_basis] = 1

    def get_stepsize(self):
        self.filter_indices = [i for i,j in enumerate(self.basic_direction) if j < 0]
        if self.filter_indices:
            self.stepsize = min(- self.x[self.filter_indices] / self.basic_direction[self.filter_indices])
            self.out_basis = self.filter_indices[np.argmin(- self.x[self.filter_indices] / self.basic_direction[self.filter_indices])]
        else:
            self.unbounded = True
    
    def update(self):
        self.iteration += 1
        print('='*9,'\t\t第',self.iteration,'次迭代\t\t','='*9)
        print('原基：\t\t',self.basic_indices)
        self.basic_indices.remove(self.out_basis)
        self.basic_indices.append(self.into_basis)
        self.B = self.A[:,self.basic_indices]
        print('本次入基变量：\t\t变量',self.into_basis)
        print('本次出基变量：\t\t变量',self.out_basis)
        print('新基：\t\t',self.basic_indices)
        print('原决策变量值：\t\t',np.round(self.x,2))
        print('本次基本方向：\t\t',np.round(self.basic_direction,2))
        print('本次前进值：\t\t',round(self.stepsize,2))
        self.x += self.stepsize * self.basic_direction
        print('新决策变量值：\t\t',np.round(self.x,2))
        print('原目标函数值：\t\t',np.round(self.objective,2))
        self.objective = np.dot(self.c,self.x)
        print('新目标函数值：',np.round(self.objective,2))
    
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
    # problem 1
    # A = np.array([[-1,1,1,0,0],[-2,1,0,1,0],[4,1,0,0,1]])
    # b = np.array([3,2,16])
    # c = np.array([3,-2,0,0,0])
    # problem 2
    # A = np.array([[1,2,2,1,0,0],[2,1,2,0,1,0],[2,2,1,0,0,1]])
    # b = np.array([20,20,20])
    # c = np.array([-10,-12,-12,0,0,0])
    # problem 3
    # A = np.array([[1,2,3,0],[-1,2,6,0,],[0,4,9,0],[0,0,3,1]])
    # b = np.array([3,2,5,1])
    # c = np.array([1,1,1,0])
    # problem 4
    # A = np.array([[1,1,1,1,0,0,0],[1,0,0,0,1,0,0],[0,0,1,0,0,1,0],[0,3,1,0,0,0,1]])
    # b = np.array([4,2,3,6])
    # c = np.array([-1,-14,-6,0,0,0,0])
    # problem 5
    # A = np.array([[3,5,1,0],[3,1,0,1]])
    # b = np.array([15,12])
    # c = np.array([-2,-1,0,0])
    problem = simplex()
    problem.solve()