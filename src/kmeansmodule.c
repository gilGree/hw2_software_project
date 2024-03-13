
# define PY_SSIZE_T_CLEAN
# include <Python.h>
#include <math.h>

/*distance between 2 points*/
double dist(double *p, double *q, int d){
    double sumofsquares = 0;
    int i;
    for(i = 0 ; i < d ; ++i){
        sumofsquares += pow(*(p+i)-*(q+i), 2);
    }
    return sqrt(sumofsquares);
}

/*assign to xi the cluster. returns the index.*/
int clusterSelection(double *x, double *mu, int K, int d){
    double minDist = dist(x, mu, d);/*the first mu*/
    int index = 0;
    int i;
    for(i = 1; i < K ; ++i){
        double tmp=dist(x, mu+(i*d), d);
        if(tmp < minDist) {
            minDist = tmp;
            index = i;
        }
    }
    return index;
}

/*assign to all x the clusters*/
void assign(double *mu, double *DB, int *association, int d, int K, int n){
    int i;
    for(i = 0 ; i < n ; i++){
        *(association+i) = clusterSelection(DB+(i*d),mu,K,d);/*the ith point.*/
    }
}

/*update mu and returns if the delta condition is violated
// 1 - keep going
// 0 - stop:*/
int updateMu(double *mu, double *DB, int *association, int d, int K, int n,int ep){
    int *clusterSize = (int *)malloc(K*sizeof(int));/*stores sizes of clusters.*/
    double *muTmp = (double *)malloc(K*d*sizeof(double));/*stores the acccumelated positions.*/
    int i,j;
    int result = 0;
    for(i = 0 ; i < K*d ; i++) *((double *)muTmp+i) = 0;/* initiate muTmp to zeros.*/
    for(i = 0 ; i < K ; i++) clusterSize[i]=0;
    for(i = 0 ; i < n; ++i){
        clusterSize[association[i]]++;
        for(j = 0;j<d;j++) muTmp[association[i]*d+j] += DB[i*d+j];
    }
    for(i = 0 ; i < K ; i++){
        for(j = 0 ; j < d ; j++){
            muTmp[i*d+j] /= clusterSize[i];
        }
    }

    for(i = 0 ; i < K ; i++){
        /*If even 1 of the centroids is far we keep going.*/
        if(dist(mu+(i*d),(muTmp+i*d),d)>=ep){
            result = 1;
        }
    }
    for(i = 0 ; i < K ; i++){
        for(j = 0 ; j < d ; j++){
            *(mu+(i*d+j)) = muTmp[i*d+j];
        }
    }
    free(clusterSize);
    free(muTmp);
    return result;
}

/*kmeans*/
void kmeans(double *DB, int d, int K, int n, int iter,int ep,double *mu){ 
    int count = 0;
    int deltaCondition;
    int *association = (int *)malloc(n*sizeof(int));;

    do{
        assign((double *)mu, DB, association, d, K, n);
        deltaCondition = updateMu((double *)mu, DB, (int *) association, d, K, n,ep);
    }while(deltaCondition && (++count)<iter);
    free(association);
}

static PyObject *fit(PyObject *self, PyObject *args) {
    PyObject *py_points,*py_centroids,*point,*item,*py_centroinds2return;
    int d,K,n,iter,i,j;
    double ep;
    /*allocate space*/
    double *DB,*mu;
    /* check if we managed to pass them well*/
    if (!PyArg_ParseTuple(args, "OOiiiid", &py_points, &py_centroids, &iter, &n, &d, &K, &ep)){
        printf("An Error Has Occurred1\n");
        return NULL;
    }
    if (!PyList_Check(py_centroids)){ /*check if its an instance of a subtype of list*/
        printf("An Error Has Occurred2\n");
        return NULL;
    }
    if (!PyList_Check(py_points)){ /*check if its an instance of a subtype of list*/
        printf("An Error Has Occurred3\n");
        return NULL;
    }

    DB = (double *)malloc(n*d * sizeof(double));
    mu = (double *)malloc(K*d * sizeof(double));

    if (DB == NULL||mu==NULL) {
        printf("An Error Has Occurred\n");
        return NULL;
    }
    /* get input points from python */
    for (i = 0; i < n; i++) {
        point = PyList_GetItem(py_points, i); /* get the row */
        if (!PyList_Check(point)){
            continue;
        }
        /* save point in points */
        for (j = 0; j < d; j++) {
            item = PyList_GetItem(point, j);
            DB[i*d + j] = PyFloat_AsDouble(item);
        }
    }
    /* get input centroids from python */
    for (i = 0; i < K; i++) {
        point = PyList_GetItem(py_centroids, i); /* get the row */
        if (!PyList_Check(point)){
            continue;
        }
        /* save point in mu */
        for (j = 0; j < d; j++) {
            item = PyList_GetItem(point, j);
            mu[i*d + j] = PyFloat_AsDouble(item);
        }
    }
    kmeans(DB,d,K,n,iter,ep,mu);
    free(DB);
    py_centroinds2return = PyList_New(K); /* create a list */
    if (py_centroinds2return == NULL){
        printf("An Error Has Occurred\n");
        return NULL;
    }
    for (i = 0; i < K; i++)
    {
        point = PyList_New(d); /* create a row */
        if (point == NULL){
       	    printf("An Error Has Occurred\n");
            return NULL;
        }
        for (j = 0; j < d; j++)/*copying the row*/
        {
            PyList_SetItem(point, j, Py_BuildValue("d", mu[i*d+j]));
        }
        PyList_SetItem(py_centroinds2return, i, Py_BuildValue("O", point)); /* adding the row */
    }
    free(mu);
    return py_centroinds2return;
}


static PyMethodDef capiMethods[] = {
    {"fit",                   
      (PyCFunction) fit,
      METH_VARARGS,         
      PyDoc_STR("expected input (in this order): points, centroids, max_iter, number of point, dimension of points,number of cluster points, epsilon")},
    {NULL, NULL, 0, NULL}     
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,  
    capiMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
