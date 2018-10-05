import numpy as np
import scipy as sp
import itertools

class simplex():

    def __init__(self,A,b,c):
        self.A = A
        self.m = len(self.A)
        self.n = np.size(self.A,1)
        self.b = b
        self.question_simplify()
        self.x = np.zeros(self.n)
        self.c = c
        self.objective = np.dot(self.c,self.x)
        self.B = []
        self.feasible = True
        self.reduced_cost = []
        self.basic_indices = []
        self.into_basis = None
        self.out_basis = None
        self.iteration = 0

    def solve(self):
        self.iteration_initial()
        if self.feasible:
            while self.feasible:
                self.get_reduced_cost()
                if self.optimal_test():
                    print('*'*9,'当前为最优解','*'*9)
                    break
                else:
                    self.get_basic_direction()
                    self.get_theta()
                    self.update()
            else:
                print(self.x)
        else:
            print('可行集为空')

    def question_simplify(self):
        if np.linalg.matrix_rank(np.c_[self.A,b]) < self.m:
            for i in itertools.combinations(range(self.m),np.linalg.matrix_rank(np.c_[self.A,b])):
                if np.linalg.matrix_rank(np.c_[self.A,b][i,:]) == np.linalg.matrix_rank(np.c_[self.A,b]):
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
                    print('='*9,'初始化','='*9)
                    print('基：',self.basic_indices)
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

    def get_theta(self):
        self.filter_indices = [i for i,j in enumerate(self.basic_direction) if j < 0]
        self.theta = min(- self.x[self.filter_indices] / self.basic_direction[self.filter_indices])
        self.out_basis = self.filter_indices[np.argmin(- self.x[self.filter_indices] / self.basic_direction[self.filter_indices])]
    
    def update(self):
        self.iteration += 1
        print('='*9,'第',self.iteration,'次迭代','='*9)
        print('原基：',self.basic_indices)
        self.basic_indices.remove(self.out_basis)
        self.basic_indices.append(self.into_basis)
        self.B = self.A[:,self.basic_indices]
        print('本次入基变量：变量',self.into_basis)
        print('本次出基变量：变量',self.out_basis)
        print('新基：',self.basic_indices)
        print('原决策变量值：',np.round(self.x,2))
        print('本次基本方向：',np.round(self.basic_direction,2))
        print('本次沿基本方向前进值：',round(self.theta,2))
        self.x += self.theta * self.basic_direction
        print('新决策变量值：',np.round(self.x,2))
        print('原目标函数值：',np.round(self.objective,2))
        self.objective = np.dot(self.c,self.x)
        print('新目标函数值：',np.round(self.objective,2))
    
    def show(self):
        pass

if __name__=='__main__':
    import pprint
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
    A = np.array([[1,1,1,1,0,0,0],[1,0,0,0,1,0,0],[0,0,1,0,0,1,0],[0,3,1,0,0,0,1]])
    b = np.array([4,2,3,6])
    c = np.array([-1,-14,-6,0,0,0,0])
    problem = simplex(A,b,c)
    problem.show()
    problem.solve()