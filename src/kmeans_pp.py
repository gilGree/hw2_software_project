import pandas as pd
import numpy as np
import sys
import mykmeanssp as mk
def input_user():
    assert((len(sys.argv)==6)or(len(sys.argv)==5)), "An Error Has Occurred"
    if len(sys.argv) == 6:
        temp_k = sys.argv[1]
        max_iter = sys.argv[2]
        ep = sys.argv[3]
        file_name1 = sys.argv[4]
        file_name2 = sys.argv[5]
    elif len(sys.argv) == 5:
        temp_k = sys.argv[1]
        ep = sys.argv[2]
        max_iter = "300"
        file_name1 = sys.argv[3]
        file_name2 = sys.argv[4]
    assert((temp_k.isdigit())and(int(temp_k)>1)), "invalid number of clusters!"
    assert((max_iter.isdigit())and(1000>int(max_iter)>1)), "Invalid maximum iteration!"
    assert(((file_name1[len(file_name1)-1]=="t")
    and(file_name1[len(file_name1)-2]=="x")
    and(file_name1[len(file_name1)-3]=="t")
    and(file_name1[len(file_name1)-4]=="."))or
    ((file_name1[len(file_name1)-1]=="v")
    and(file_name1[len(file_name1)-2]=="s")
    and(file_name1[len(file_name1)-3]=="c")
    and(file_name1[len(file_name1)-4]=="."))),"NA"
    assert(((file_name2[len(file_name2)-1]=="t")
    and(file_name2[len(file_name2)-2]=="x")
    and(file_name2[len(file_name2)-3]=="t")
    and(file_name2[len(file_name2)-4]=="."))or
    ((file_name2[len(file_name2)-1]=="v")
    and(file_name2[len(file_name2)-2]=="s")
    and(file_name2[len(file_name2)-3]=="c")
    and(file_name2[len(file_name2)-4]=="."))),"NA"
    assert(float(ep)>=0),"Invalid epsilon!"

    return int(temp_k), int(max_iter), float(ep), file_name1, file_name2
  
def main():
    K,max_iter,ep,file_name_1,file_name_2 = input_user()
    # build the dataframes:
    df1 = pd.read_csv(file_name_1, header=None)# none header to not override the first row.
    df2 = pd.read_csv(file_name_2, header=None)# none header to not override the first row.
    # relabel the first column (key):
    df1.rename(columns={df1.columns[0]: 'id'}, inplace=True)
    df2.rename(columns={df2.columns[0]: 'id'}, inplace=True)

    # merging:
    df = pd.merge(df1, df2, on="id")

    # Step 3:
    df = df.sort_values(by="id", ascending=True)

    # Step 4:
    df.drop(columns='id', inplace=True)# unnecessary data column
    n = df.shape[0]
    d = df.shape[1]
    assert(K<n), "invalid number of clusters!"

    chosen_indices = []# the centroids
    np.random.seed(0)# the required seed
    center_index = np.random.choice(n)# uniform selection for 1st center.

    center = df.iloc[center_index].values# row as np array
    chosen_indices.append(center_index)

    while(len(chosen_indices) < K):
        distances = [np.linalg.norm(center - df.iloc[i].values) if i not in chosen_indices else 0 for i in range(n)]# dist to newest center
        for j in chosen_indices:
            maybe_better_distances = [np.linalg.norm(df.iloc[j].values - df.iloc[i].values) if i not in chosen_indices else 0 for i in range(n)]
            distances = [min(maybe_better_distances[i], distances[i]) for i in range(len(distances))]
        D = sum(distances)
        probabilities = np.array([float(x)/D for x in distances])# as required
        center_index  = np.random.choice(a=np.array([x for x in range(n)]), p=probabilities)# update as required 
        center = df.iloc[center_index].values# update
        chosen_indices.append(center_index)
    for i in range(len(chosen_indices)-1):
        print(str(chosen_indices[i])+",",end='')
    print(format(chosen_indices[len(chosen_indices)-1]))
    DB = df.values.tolist()
    mu = []
    for i in range(len(chosen_indices)):
        mu.append(DB[chosen_indices[i]])
    mu = mk.fit(DB, mu, max_iter, n, d, K, ep)
    for i in range(len(mu)):
        for j in range(len(mu[i])-1):
            print('{:.4f}'.format(mu[i][j]),end=",")
        print('{:.4f}'.format(mu[i][len(mu[i])-1]))
if __name__ == "__main__":
    main()
