/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Build_cBiK.h
 * Author: carlos
 *
 * Created on January 16, 2018, 3:28 PM
 */

#ifndef BUILD_CBIK_H
#define BUILD_CBIK_H

#include <stdbool.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <list>
#include <iostream>
#include <queue>
#include "NodeSKQ.h"
#include "Node.h"
#include "PriorityQueue.h"
#include "MinHeap.h"
#include <set>
//Para los tiempos
#include <time.h>
#include <sys/time.h>
//Para leer archivos
#include <sstream>
#include <fstream>
//Para la funcion hash
#include <unordered_map>
#include <string>
#include <locale>
#include "utils.h"
using namespace std;

class Build_cBiK {
public:
    Build_cBiK(string datasetName);
    Build_cBiK();
    
    void load_cBiK(string datasetNameA, string datasetNameB);
    
    double hdKD1();
    double hdKD2();
    double hdKD2v2();

    double bruteForce();
    double getStorageSize();
    std::vector<Point> extractPointKDTree();
    int getHeapMaxElements();
    int getHeapExitElements();

private:
     int heapMaxElements;
     int exitHeapElements;

    //Atributos generales de construccion
    string datasetNameBuild;
    string datasetNameShort;
    string datasetNameA;
    string datasetNameB;
    //En esta variable se guarda los puntos recibidos
    vector<NodeSKQ> nodosSKQ;
    //En esta variable se guardan las referencias finales del KD-Tree
    vector<NodeSKQ> KdTreeBuild;
    vector<NodeSKQ> KdTreeA;
    vector<NodeSKQ> KdTreeB;
    
    //para convertir las palabras
    locale loc;
    
    int* pNN = new int[2];
    
    //region completa de conjuntoA, se inicializan como infinitos
    int xMenorA = numeric_limits<int>::max();
    int yMenorA = numeric_limits<int>::max();
    int xMayorA = numeric_limits<int>::min();
    int yMayorA = numeric_limits<int>::min();
    
    //region completa, se inicializan como infinitos
    int xMenorB = numeric_limits<int>::max();
    int yMenorB = numeric_limits<int>::max();
    int xMayorB = numeric_limits<int>::min();
    int yMayorB = numeric_limits<int>::min();
    
    
    void loadDataset();
    void create_iKDTree();
    void export_iKDTree();
    
    void initializeReference(vector<NodeSKQ>& coordinates, vector<NodeSKQ>& reference);
    double superKeyCompare(const int *a, const int *b, const long p, const long dim);
    void mergeSort(vector<NodeSKQ> &reference, vector<NodeSKQ>& temporary, const long low, const long high, const long p, const long dim);
    long removeDuplicates(vector<NodeSKQ>& reference, const long i, const long dim);
    void buildKdTree(vector< vector<NodeSKQ> >& references, vector<NodeSKQ>& temporary, const long start, const long end, const long dim, const long depth);
    
    double getEuclideanDistance(int p1[], int p2[]);
    double minDist(const int* q, int* r1, int* r2, int p);
    double maxDist(int* q, int* r1, int* r2, int p);
    
    double timeval_diff(struct timeval *a, struct timeval *b);
    void  saveMaxPqueue(unsigned long num) ;
    void  saveExitHeap(unsigned long num) ;
    //HDKP2
    double nearestNeighbor(int query[], int start, int end, const int profundidad, double minNN);
    double HDKPruning(double cmax);
    double HDKPruningv2(double cmax);

    double isCandidate(int region1[], int region2[], double cmax);
    double MaxMaxDist(int* region1A, int* region2A, int* region1B, int* region2B);
    
    //Metodos compartidos
    void getRegion(int region[], Node nodoQ, bool iod, bool aob);
    
    double NNMax(int pointP[], double cmax, double minNN);
    bool update(MinHeap &pQ, Node nodoQ, double cmax, double minNN, int* pointP, bool iod);
    
    void load_KDTree(string filename, bool first);
    
    void printTree();
    void printKdTree(vector<NodeSKQ>& tree, int start, int end, int depth);
    void printTuple(int* tuple);
    string getPoint(int* tuple);


};

#endif /* BUILD_CBIK_H */

