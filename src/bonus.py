from sklearn.cluster import KMeans
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt

# load data
data = load_iris().data

#calculate inertias
inertias = [KMeans(n_clusters=k, init='k-means++',
                random_state=0).fit(data).inertia_ for k in range(1,11)]
#define indexes
k_indexes = range(1,11)

#plot the graph
plt.plot(k_indexes, inertias)
plt.xlabel("K")
plt.ylabel("Inertia")
plt.title('Elbow Method for selection of optimal "K" clusters')
plt.annotate("Elbow",xy=(2,inertias[1]),xytext=(4, 300),arrowprops=dict(arrowstyle="->"))

#save the graph into a png file
plt.savefig('elbow.png')