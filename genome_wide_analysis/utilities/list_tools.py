'''
Created on 2013-05-19

@author: jyeung
'''

def chunks(my_list, n):
    '''
    Takes a long list and puts n into each bin until you reach the end
    of your list. 
    '''
    return [my_list[i:i+n] for i in range(0, len(my_list), n)]