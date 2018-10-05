class Solution():
    '''return the index of two numbers
    
    Given an array of integers, return indices of the two numbers such that they add up to a specific target.
    You may assume that each input would have exactly one solution, and you may not use the same element twice.
    '''
    def two_sum_1(self,nums,target):
        if len(nums) <= 1:
            return False
        buff_dict = {}
        for i in range(len(nums)):
            if nums[i] in buff_dict:
                return [buff_dict[nums[i]], i]
            else:
                buff_dict[target - nums[i]] = i

if __name__=='__main__':
    nums = [2 ** i for i in range(50)]
    s = Solution()
    solution = s.two_sum_1(nums,100)
    print('solution:',solution)

