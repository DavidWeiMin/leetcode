def dichotomous(f,a,b,d,l):
    '''使用二分法搜索严格拟凸函数的最小值
    
    Arguments:
        f {fun} -- 目标函数
        a {float} -- 区间下界
        b {float} -- 区间上界
        d {float} -- distinguish constant:非常小的正数
        l {float} -- 不确定性区间长度
        l > 2 * d

    Returns:
        [list] -- 最小值所在区间
    '''
    while l < b - a:
        lamb = (a + b) / 2 - d
        mu = (a + b) / 2 + d
        if f(lamb) < f(mu):
            b = mu
        else:
            a = lamb
        print([a,b])
    return [a,b]

def golden_section(f,a,b,l):
    '''使用黄金分割法搜索严格拟凸函数的最小值
    
    Arguments:
        f {fun} -- 目标函数
        a {float} -- 区间下界
        b {float} -- 区间上界
        l {float} -- 不确定性区间长度
    
    Returns:
        [list] -- 最小值所在区间
    '''

    alpha = 0.61803399
    lamb = a + (1 - alpha) * (b - a)
    mu = a + alpha * (b - a)
    while l < b - a:
        if f(lamb) < f(mu):
            b = mu
            mu = lamb
            lamb = a + (1 - alpha) * (b - a)
        else:
            a = lamb
            lamb = mu
            mu = a + alpha * (b - a)
        print([a,b])
    return [a,b]

def fibonacci(n):
    if n==0 or n==1:
        return 1
    else:
        return fibonacci(n - 1) + fibonacci(n - 2)

def fibonacci_search(f,a,b,d,l):
    '''使用斐波那契法搜索严格拟凸函数的最小值
    
    Arguments:
        f {fun} -- 目标函数
        a {float} -- 区间下界
        b {float} -- 区间上界
        d {float} -- distinguish constant:非常小的正数
        l {float} -- 不确定性区间长度
    
    Returns:
        [list] -- 最小值所在区间
    '''
    F = []
    for i in range(1000):
        F.append(fibonacci(i))
        if fibonacci(i) > (b -a) / l:
            n = i
            break
    k = 1
    lamb = a + F[n - 3] / F[n - 1] * (b - a)
    mu = a + F[n - 2] / F[n - 1] * (b - a)
    while l < b - a:
        k += 1
        if k < n - 1:
            if f(lamb) < f(mu):
                b = mu
                mu = lamb
                lamb = a + F[n - k - 2] / F[n - k] * (b - a)
            else:
                a = lamb
                lamb = mu
                mu = a + F[n - k - 1] / F[n - k] * (b - a)
        else:
            if f((a + b) / 2 + d) < f((a + b) / 2):
                a = (a + b) / 2
            else:
                b = (a + b) / 2
        print([a,b])
    return [a,b]

f = lambda x: x ** 2 + 2 * x
# solution1 = dichotomous(f,-3,5,0.0001,0.001)
# solution2 = golden_section(f,-3,5,0.001)
solution3 = fibonacci_search(f,-3,5,0.000001,0.0001)
# print(solution1,solution2)
