import sys
import numpy as np
import symnmf

np.random.seed(1234)

'''
============================================PRINT_MATRIX=========================================
this is a simple helper method to print a matrix in 4 decimal places format, seperated by commas
'''
def print_matrix(mat):
    for row in mat:
        print(','.join(['%.4f' % num for num in row]))

'''
=================================================================================================
=============================================MAIN================================================
=================================================================================================
'''
def main():
    # read and parse the arguments
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    
    # parse the points into a python matrix using numpy's loadtxt
    np_array = np.loadtxt(file_name, delimiter = ',')
    data_points = np_array.tolist()
    
    # call the c module
    if goal == 'sym':
        res = symnmf.sym(data_points)
        print_matrix(res)
    elif goal == 'ddg':
        res = symnmf.ddg(data_points)
        print_matrix(res)
    elif goal == 'norm':
        res = symnmf.norm(data_points)
        print_matrix(res)
    #if we want the full run, we need the matrix W, and to calculate the average over all its entries
    # and then we can get the initial H and call symnmf
    elif goal == 'symnmf':
        W = symnmf.norm(data_points)
        m = np.array(W).mean() #convert to numpy array and get average cleanly
        upper_bound = 2 * np.sqrt(m/k)
        N = len(data_points)
        init_H = np.random.uniform(low=0, high=upper_bound, size=(N,k)) # a random matrix of size Nxk with values in the legal interval
        optimized_H = symnmf.symnmf(data_points, init_H.tolist()) #pass the data and the initial H in the form of a python matrix
        print_matrix(optimized_H)
        
    
if __name__ == "__main__":
    main()
