import sys
import numpy as np #for labeling

#read from .txt file, each line becomes a list of coordinates
#output - list of points
def read_input():
    points = []
    for line in sys.stdin:
        if line.strip() != '':
            point = list(map(float, line.strip().split(',')))
            points.append(point)
    return points

def euclidean(p1, p2):
    total = 0
    for i in range(len(p1)):
        total += (p1[i] - p2[i]) ** 2
    return total ** 0.5
    
def assign_clusters(data, centroids):
    clusters = [[] for _ in centroids]
    for point in data:
        dist_lst = [euclidean(point, centroid) for centroid in centroids] #gather distances from each centroid for each point
        min_ind = dist_lst.index(min(dist_lst)) 
        clusters[min_ind].append(point) #assign to closest centroid
    return clusters
    
def update_centroids(clusters, dim):
    new = []
    for cluster in clusters:
        #if cluster isn't empty, take the mean of each coordinate to be the new centroid
        if cluster:
            centroid = []
            denominator = 1 / abs(len(cluster))
            for i in range(dim):
                centroid.append(sum([coord[i] for coord in cluster]) * denominator)
        else:
            centroid = [0.0000 for _ in range(dim)]
        new.append(centroid)
    return new
    
def has_converged(old, new, epsilon):
    for i in range(len(old)):
        if euclidean(old[i], new[i]) >= epsilon:
            return False
    return True
    
def format_point(point):
    return ','.join(f"{coord:.4f}" for coord in point)

#==========================ADDED METHOD==========================================
# this is a new method that returns a list with the indices of the closest cluster 
# of each point
def get_labels(data, centroids):
    labels = []
    for point in data:
        distances = [euclidean(point, c) for c in centroids]
        closest_index = np.argmin(distances)
        labels.append(closest_index)
    return labels

#==========================MODDED METHOD=========================================
#this method was modded to return the cluster assignment instead of a list of clusters   
def find_kmeans_labels(data, k, iter=300, epsilon=1e-4):
        
    #initialize first k points as initial centroids
    centroids = [point.copy() for point in data[:k]]
    dim = len(data[0])
    
    #perform clustering until reaching one of the stop conditions
    for _ in range(iter):
        clusters = assign_clusters(data, centroids)
        new_cen = update_centroids(clusters, dim)
        if has_converged(centroids, new_cen, epsilon):
            break
        centroids = new_cen
    
    final_labels = get_labels(data, centroids)
    return final_labels
            
        
if __name__ == "__main__":

    def main():
            
        try:
            k_raw = float(sys.argv[1])
            if not k_raw.is_integer() or k_raw <= 0:
                raise ValueError
            k = int(k_raw)
        except:
            print("Incorrect number of clusters!")
            sys.exit(1)

        if len(sys.argv) == 3:
            try:
                iter_raw = float(sys.argv[2])
                if not iter_raw.is_integer() or iter_raw <= 0:
                    raise ValueError
                iter = int(iter_raw)
            except:
                print("Incorrect maximum iteration!")
                sys.exit(1)
        else:
            iter = 400
            
        data = read_input()
        
        #assert number of clusters is not greater than number of vectors 
        if k > len(data):
            print("Incorrect number of clusters!")
            sys.exit(1)
        
        final_centroids = kmeans(data, k, iter)
        
        for centroid in final_centroids:
            print(format_point(centroid))
        
    main()
