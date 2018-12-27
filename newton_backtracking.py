import sympy

class Problem:
    def __init__(self):
        '''创建目标函数'''

        self.x = sympy.symbols('x')
        self.create_y()
    
    def create_y(self):
        '''创建目标函数及其一、二阶导数'''

        self.y = sympy.log(sympy.exp(self.x) + sympy.exp(- self.x)) # 目标函数
        self.dy = sympy.diff(self.y,self.x) # 一阶导
        self.dy2 = sympy.diff(self.dy,self.x) # 二阶导
    
    def newton(self,x0,x_tol=1e-8,fun_tol=1e-8,max_iter=20):
        """使用牛顿法求解
        
        Arguments:
            x0 {[float]} -- [初始值]
        
        Keyword Arguments:
            x_tol {[float]} -- [variable tolerance] (default: {1e-8})
            fun_tol {[float]} -- [function value tolerance] (default: {1e-8})
            max_iter {[float]} -- [最大迭代次数] (default: {20})
        
        Returns:
            [list] -- [返回每次迭代的变量值和函数值]
        """

        alpha = 0.2
        beta = 0.5
        x = x0 
        num = 0 # 记录迭代次数
        iter_x = [x0] # 记录每次迭代的 x 
        iter_y = [float(self.y.subs(self.x,iter_x[-1]))] # 记录每次迭代的 y
        while 1:
            '''牛顿迭代'''
            delta_x = float((- self.dy / self.dy2).subs(self.x,x))
            t = 1
            while float(self.y.subs(self.x,x + t * delta_x)) > float(self.y.subs(self.x,x)) + alpha * t * float(delta_x * self.dy.subs(self.x,x)): # backtracking stop condition
                '''backtracking'''    
                t = beta * t
            if abs(t * delta_x) < x_tol: 
                break # stop criterior 
            if abs(alpha * t * float(delta_x * self.dy.subs(self.x,x))) < fun_tol: 
                break # stop criterior 
            x = x + t * delta_x
            iter_x.append(x)
            iter_y.append(float(self.y.subs(self.x,iter_x[-1])))
            num += 1
            if num > max_iter:
                break
        return iter_x,iter_y

if __name__=='__main__':
    fun_tol = 1e-8
    x_tol = 1e-8
    max_iter = 5
    problem = Problem()
    iter_x,iter_y = problem.newton(1,x_tol,fun_tol,max_iter)
    print(iter_x)
    print(iter_y)