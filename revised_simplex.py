import numpy as np
import scipy as sp
import itertools

class revised_simplex():

    def __init__(self,A,b,c):
        self.A = A
        self.b = b
        self.c = c
        # self.input()
        self.show()
        self.m = len(self.A)
        self.n = np.size(self.A,1)
        self.question_simplify()
        self.feasible = True
        self.iteration = 0
        self.unbounded = False
        self.error = 1e-12

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
                if self.is_optimal():
                    print('*'*9,'\t\t最优解\t\t','*'*9)
                    break
                else:
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
                if np.linalg.matrix_rank(np.c_[self.A,self.b][i,:]) == np.linalg.matrix_rank(np.c_[self.A,self.b]):
                    self.A = self.A[i,:]
                    self.b = np.array([self.b[j] for j in i])
                    self.m = len(self.A)
                    break

    def iteration_initial(self):
        for i in itertools.combinations(range(self.n),self.m):
            if np.linalg.det(self.A[:,i]) != 0:
                self.basic_indices = list(i)
                self.nonbasic_indices = [i for i in range(self.n) if i not in self.basic_indices]
                self.B = self.A[:,self.basic_indices]
                self.B_inverse = np.linalg.inv(self.B)
                if all(np.linalg.solve(self.B,self.b) > 0):
                    self.x = np.zeros(self.n)
                    self.x[self.basic_indices] = np.dot(self.B_inverse,self.b)
                    self.objective = np.dot(self.c,self.x)
                    print('='*9,'\t\t初始化\t\t','='*9)
                    print('基：\t',self.basic_indices)
                    print('决策变量：',np.round(self.x,2))
                    print('目标函数：',np.round(self.objective,2))
                    break
        else:
            self.feasible = False
    
    def is_optimal(self):
        self.reduced_cost = np.zeros(self.n)
        self.reduced_cost[self.nonbasic_indices] = self.c[self.nonbasic_indices] - np.dot(np.dot(self.c[self.basic_indices],self.B_inverse),self.A[:,self.nonbasic_indices])
        if min(self.reduced_cost) < -self.error:
            self.into_basis = np.argmin(self.reduced_cost)
            return False
        else:
            return True

    def get_stepsize(self):
        self.u = np.dot(self.B_inverse,self.A[:,self.into_basis])
        if any(self.u > self.error):
            temp = self.x[np.array(self.basic_indices)[self.u > self.error]] / self.u[self.u > self.error]
            self.stepsize = min(temp)
            self.out_basis = np.array(self.basic_indices)[self.u==self.u[self.u > self.error][np.argmin(temp)]]
            self.out_basis = self.out_basis[0]
        else:
            self.unbounded = True

    def update(self):
        self.x[self.basic_indices] -= self.stepsize * self.u
        self.x[self.into_basis] = self.stepsize
        self.iteration += 1
        print('='*9,'\t\t迭代',self.iteration,'\t\t','='*9)
        temp = self.basic_indices.index(self.out_basis)
        self.basic_indices.remove(self.out_basis)
        self.basic_indices.append(self.into_basis)
        self.nonbasic_indices = [i for i in range(self.n) if i not in self.basic_indices]
        self.B = self.A[:,self.basic_indices]
        for i,j in enumerate(self.u):
            if i != temp:
                if j != 0:
                    self.B_inverse[i] -= self.B_inverse[temp] * (j / self.u[temp])
            else:
                self.B_inverse[i] /= self.u[temp]
        print('入基变量：变量',self.into_basis)
        print('出基变量：变量',self.out_basis)
        print('基：\t',self.basic_indices)
        print('前进方向：',np.round(- self.u,2))
        print('前进值：',round(self.stepsize,2))
        print('决策变量：',np.round(self.x,2))
        self.objective = np.dot(self.c,self.x)
        print('目标函数：',np.round(self.objective,2))
    
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
    # problem 6
    A = np.array([[1,1,1,-1,0,0],[-2,0,1,0,-1,0],[0,2,-1,0,0,-1]])
    b = np.array([6,2,0])
    c = np.array([-2,1,-2,0,0,0])
    B = np.array([[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]])
    print(np.linalg.inv(B))
    problem = revised_simplex(A,b,c)
    problem.solve()