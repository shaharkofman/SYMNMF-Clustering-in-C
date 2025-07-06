import sys
import numpy as np
import symnmf
import kmeans #modified HW1 kmeans implementation
from sklearn.metrics import silhouette_score

np.random.seed(1234)

def main():
    # read and parse the arguments
    k = int(sys.argv[1])
    file_name = sys.argv[2]
    
    # parse the points into a python matrix using numpy's loadtxt
    np_array = np.loadtxt(file_name, delimiter = ',')
    data_points = np_array.tolist()
    N = len(data_points)
    
    # call c module to compute the optimized H matrix
    W = symnmf.norm(data_points)
    m = np.array(W).mean() #convert to numpy array and get average cleanly
    upper_bound = 2 * np.sqrt(m/k)
    N = len(data_points)
    init_H = np.random.uniform(low=0, high=upper_bound, size=(N,k)) # a random matrix of size Nxk with values in the legal interval
    optimized_H = symnmf.symnmf(data_points, init_H.tolist()) #pass the data and the initial H in the form of a python matrix

    #convert the H matrix to a numpy array and extract the custer assignment for each point using argmax()
    H_np = np.array(optimized_H)
    symnmf_clusters = np.argmax(H_np, axis=1)
    
    #get the cluster assignment from kmeans.py
    kmeans_clusters = kmeans.find_kmeans_labels(data_points, k)
    
    #get silhouette scores
    symnmf_score = silhouette_score(data_points, symnmf_clusters)
    kmeans_score = silhouette_score(data_points, kmeans_clusters)
    
    print(f"symnmf: {symnmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")

if __name__ == "__main__":
    main()
    
