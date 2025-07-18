/* 
 * File:   Build_cBiK.cpp
 * Author: carlos
 * 
 * Created on January 16, 2018, 3:28 PM
 */

#include "Build_cBiK.h"
#include "MinHeap.h"
#include "MaxHeap.h"

std::string extractFilename2(const std::string &path) {
    size_t pos = path.find_last_of('/');
    if (pos == std::string::npos) {
        return path; // Path doesn't contain a '/' separator
    } else {
        return path.substr(pos + 1); // Extract filename from position after '/'
    }
}
Build_cBiK::Build_cBiK(string datasetName) {
    this->datasetNameBuild = datasetName;
    this->datasetNameShort = extractFilename2(datasetName);
    loadDataset();
    create_iKDTree();
    export_iKDTree();
}

Build_cBiK::Build_cBiK() {
    
}

void Build_cBiK::load_cBiK(string datasetNameA, string datasetNameB) {
    this->datasetNameA = datasetNameA;
    this->datasetNameB = datasetNameB;
    
    load_KDTree(datasetNameA, true);//Carga datos de A
    load_KDTree(datasetNameB, false);//Carga datos de B
    
    //cout << "(" << xMenor << "," << yMenor << ")" << " y " << "(" << xMayor << "," << yMayor << ")" << endl;
}

void Build_cBiK::load_KDTree(string filename, bool first) {
    /**************************************************************************/
    /************************** CARGAR PUNTOS *******************************/
    /**************************************************************************/
  //  cout << "CARGANDO PUNTOS >> " << "iKDTree_"+filename+"_points" << endl;
    int totalObjects = 0;
    ifstream lecturaPalabras;
    lecturaPalabras.open("iKDTree_"+filename+"_points",ios::in);
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        totalObjects++;
    }
    lecturaPalabras.close();
 //   cout << "Total de objetos >> " << totalObjects << endl;
    lecturaPalabras.open("iKDTree_"+filename+"_points",ios::in);

    //arreglo para guardar los puntos
    int **coordinates;
    coordinates = new int *[totalObjects];
    for(int i = 0; i<totalObjects; i++) {
        coordinates[i] = new int[2];//dos dimensiones
    }
    
    //Se lee cada linea del archivo de texto
    int contToken=0;
    int i=0;
    for(string linea; getline(lecturaPalabras, linea); ) {
        stringstream registro(linea);
        
        contToken = 0;
        //Lee cada elemento de la linea
        for (string dato; getline(registro, dato, ' '); ) {
            if(contToken == 0) {
                //LATITUD
                coordinates[i][0] = stod(dato);
            }else if(contToken == 1) {
                //LONGITUD
                coordinates[i][1] = stod(dato);
            }
            contToken++;
        }
        i++;
    }
    lecturaPalabras.close();


    if(first) {
        KdTreeA.resize(totalObjects);
        for (int i = 0; i < KdTreeA.size(); i++) {
            KdTreeA[i] = NodeSKQ( &(coordinates[i][0]) );
            
            //Obtiene la region del espacio evaluado A
            if(coordinates[i][0] < xMenorA) {
                xMenorA = coordinates[i][0];
            }
            
            if(coordinates[i][0] > xMayorA) {
                xMayorA = coordinates[i][0];
            }
            
            if(coordinates[i][1] < yMenorA) {
                yMenorA = coordinates[i][1];
            }
            
            if(coordinates[i][1] > yMayorA) {
                yMayorA = coordinates[i][1];
            }
        }
    }else {
        KdTreeB.resize(totalObjects);
        for (int i = 0; i < KdTreeB.size(); i++) {
            KdTreeB[i] = NodeSKQ( &(coordinates[i][0]) );
            
            //Obtiene la region del espacio evaluado
            if(coordinates[i][0] < xMenorB) {
                xMenorB = coordinates[i][0];
            }
            
            if(coordinates[i][0] > xMayorB) {
                xMayorB = coordinates[i][0];
            }
            
            if(coordinates[i][1] < yMenorB) {
                yMenorB = coordinates[i][1];
            }
            
            if(coordinates[i][1] > yMayorB) {
                yMayorB = coordinates[i][1];
            }
        }
    }
}

void Build_cBiK::export_iKDTree() {
    /**************************************************************************/
    /************************** EXPORTAR PUNTOS *******************************/
    /**************************************************************************/
 //   cout << "EXPORTANDO PUNTOS >> " << "iKDTree_"+datasetNameBuild+"_points" << endl;
    ofstream archivoP("iKDTree_"+datasetNameShort+"_points");
    for(int i=0; i<KdTreeBuild.size(); i++) {
        archivoP << KdTreeBuild[i].coordenadas[0] << " " << KdTreeBuild[i].coordenadas[1] << endl;
    }
    archivoP.close();
}
std::vector<Point> Build_cBiK::extractPointKDTree(){
    std::vector<Point> puntillos;

    for (int i = 0; i < KdTreeA.size(); i++) {
        Point aux;
        aux.x = KdTreeA[i].coordenadas[0];
        aux.y = KdTreeA[i].coordenadas[1];
        puntillos.push_back(aux);
    }
    return puntillos;
}

/*agregado por fernando */
double  Build_cBiK::getStorageSize(){
    // Calcular el tamaño en memoria de un NodeSKQ (solo el puntero)
    size_t nodeSize = sizeof(int*); // Tamaño del puntero, no del array al que apunta

    // Calcular el tamaño total de los punteros en el vector
    size_t totalNodeSize = KdTreeA.size() * nodeSize;

    // Calcular el tamaño total de los arrays de coordenadas
    size_t totalCoordSize = 0;
    for (const auto& node : KdTreeBuild) {
        totalCoordSize += sizeof(int) * 2; // Asumimos que cada NodeSKQ tiene un array de 2 int
    }
    return (totalNodeSize + totalCoordSize);
}

void Build_cBiK::loadDataset() {
    long long int totalObjects = 0;
    
  //  cout << "LEYENDO DATASET >> " << datasetNameBuild << endl;
    
    /**************************************************************************/
    /*************************** CUENTA OBJETOS *******************************/
    /**************************************************************************/
    ifstream lecturaPalabras;
    lecturaPalabras.open(datasetNameBuild,ios::in);
    
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        totalObjects++;
    }
    lecturaPalabras.close();
    //cout << "Total de objetos >> " << totalObjects << endl;
    
    /**************************************************************************/
    /*************************** CARGAR OBJETOS *******************************/
    /**************************************************************************/
    lecturaPalabras.open(datasetNameBuild,ios::in);
    
    //arreglo para guardar los puntos
    int **coordinates;
    coordinates = new int *[totalObjects];
    for(int i = 0; i<totalObjects; i++) {
        coordinates[i] = new int[2];//dos dimensiones
    }
    
    //HASHING PALABRA-ID
    int contToken=0;
    int i=0;
    
    //Se lee cada linea del archivo de texto
    for(string linea; getline(lecturaPalabras, linea); ) {
        stringstream registro(linea);
        contToken = 0;
        
        //Lee cada elemento de la linea
        for (string dato; getline(registro, dato, ','); ) {
            if(contToken == 0) {
                //LATITUD
                coordinates[i][0] = stod(dato);
            }else if(contToken == 1) {
                //LONGITUD
                coordinates[i][1] = stod(dato);
            }
            contToken++;
        }
        i++;
    }
    lecturaPalabras.close();
    
    /**************************************************************************/
    /*************************** VECTOR nodosSKQ*******************************/
    /**************************************************************************/
    
    //Se guardan las coordenadas y las palabras en un objeto
    this->nodosSKQ.clear();
    this->nodosSKQ.resize(totalObjects);
    for (int i = 0; i < this->nodosSKQ.size(); i++) {
        this->nodosSKQ[i] = NodeSKQ( &(coordinates[i][0]) );
    }
}

void Build_cBiK::create_iKDTree(){
   // cout << "\nCREANDO iKD-Tree >> " << endl;
    
    struct timeval iniTime2, endTime2;
    double secs;
    
    //Marca el tiempo de inicio
  //  gettimeofday(&iniTime2, NULL);

    vector<NodeSKQ>& coordinates = this->nodosSKQ;
    const long numDimensions = 2;//numero de dimensiones
    
    //Inicializa y ordena la referencia del arreglo
    vector< vector<NodeSKQ> > references(numDimensions, vector<NodeSKQ>(coordinates.size()));
    vector<NodeSKQ> temporary(coordinates.size());

    for (long i = 0; i < references.size(); i++) {
        initializeReference(coordinates, references.at(i));
        mergeSort(references[i], temporary, 0, references[i].size() - 1, i, numDimensions);
    }
    
    //Elimine las referencias a coordenadas duplicadas mediante una pasada a través de cada matriz de referencia.
    vector<long> end(references.size());
    for (long i = 0; i < end.size(); i++) {
        end[i] = removeDuplicates(references[i], i, numDimensions);
    }

    //Verifica que se elimino la misma cantidad de referencias de cada matriz de referencia.
    for (long i = 0; i < end.size() - 1; i++) {
        for (long j = i + 1; j < end.size(); j++) {
            if (end.at(i) != end.at(j)) {
                cout << "reference removal error" << endl;
                exit(1);
            }
        }
    }
    
    //Se crea el arreglo final
    long long int totalSize = (end.at(0) + 1);
    KdTreeBuild.resize(totalSize);

    //Crea el KD-Tree por Medianas
    buildKdTree(references, temporary, 0, end.at(0), numDimensions, 0);

    //marco el tiempo final
 //   gettimeofday(&endTime2, NULL);
  //  secs = timeval_diff(&endTime2, &iniTime2);
    
   // printf("Tiempo: %.16g microseg -- %.16g miliseg -- %.16g seg\n", secs * 1000000.0, secs * 1000.0, secs);
  //  cout << endl;
}

/*
 * Initialize a reference array by creating references into the coordinates array.
 *
 * calling parameters:
 *
 * coordinates - a vector<NodeSKQ> of pointers to each of the (x, y, z, w...) tuples
 * reference - a vector<NodeSKQ> that represents one of the reference arrays
 */
void Build_cBiK::initializeReference(vector<NodeSKQ>& coordinates, vector<NodeSKQ>& reference) {
    for (long j = 0; j < coordinates.size(); j++) {
        reference.at(j) = coordinates.at(j);
    }
}

/*
 * The superKeyCompare method compares two double arrays in all k dimensions,
 * and uses the sorting or partition coordinate as the most significant dimension.
 *
 * calling parameters:
 *
 * a - a double*
 * b - a double*
 * p - the most significant dimension
 * dim - the number of dimensions
 *
 * returns: a double result of comparing two double arrays
 */
double Build_cBiK::superKeyCompare(const int *a, const int *b, const long p, const long dim) {
    double diff = 0;
    for (long i = 0; i < dim; i++) {
        long r = i + p;
        // A fast alternative to the modulus operator for (i + p) < 2 * dim.
        r = (r < dim) ? r : r - dim;
        diff = a[r] - b[r];
        if (diff != 0) {
            break;
        }
    }
    return diff;
}

/*
 * The mergeSort function recursively subdivides the array to be sorted
 * then merges the elements. Adapted from Robert Sedgewick's "Algorithms
 * in C++" p. 166. Addison-Wesley, Reading, MA, 1992.
 *
 * calling parameters:
 *
 * reference - a vector<NodeSKQ> that represents the reference array to sort
 * temporary - a temporary array into which to copy intermediate results;
 *             this array must be as large as the reference array
 * low - the start index of the region of the reference array to sort
 * high - the high index of the region of the reference array to sort
 * p - the sorting partition (x, y, z, w...)
 * dim - the number of dimensions
 * depth - the depth of subdivision
 */
void Build_cBiK::mergeSort(vector<NodeSKQ> &reference, vector<NodeSKQ>& temporary, const long low, const long high, const long p, const long dim) {
    long i, j, k;

    if (high > low) {
        // Evite el desbordamiento al calcular la mediana.
        const long mid = low + ((high - low) >> 1);

        // Recursivamente subdivide las mitades inferior y superior de la matriz.
        mergeSort(reference, temporary, low, mid, p, dim);
        mergeSort(reference, temporary, mid + 1, high, p, dim);

        // Combina los resultados para este nivel de subdivisión.
        for (i = mid + 1; i > low; i--) {
            temporary[i - 1] = reference[i - 1];
        }
        for (j = mid; j < high; j++) {
            temporary[mid + (high - j)] = reference[j + 1]; // Evite el desbordamiento de direcciones.
        }
        for (k = low; k <= high; k++) {
            if(superKeyCompare(temporary[i].coordenadas, temporary[j].coordenadas, p, dim) < 0) {
                reference[k] = temporary[i++];
            }else {
                reference[k] = temporary[j--];
            }
        }
    }
}

/*
 * Check the validity of the merge sort and remove duplicates from a reference array.
 *
 * calling parameters:
 *
 * reference - a vector<NodeSKQ> that represents one of the reference arrays
 * i - the leading dimension for the super key
 * dim - the number of dimensions
 *
 * returns: the end index of the reference array following removal of duplicate elements
 */
long Build_cBiK::removeDuplicates(vector<NodeSKQ>& reference, const long i, const long dim) {
    long end = 0;
    for (long j = 1; j < reference.size(); j++) {
        double compare = superKeyCompare(reference[j].coordenadas, reference[j - 1].coordenadas, i, dim);
        if (compare < 0) {
            cout << "merge sort failure: superKeyCompare(ref[" << j << "], ref["
                    << j - 1 << "], (" << i << ") = " << compare << endl;
            exit(1);
        } else if (compare > 0) {
            reference[++end] = reference[j];
        }
    }
    return end;
}

/*
 * This function builds a k-d tree by recursively partitioning the
 * reference arrays and adding kdNodes to the tree.  These arrays
 * are permuted cyclically for successive levels of the tree in
 * order that sorting occur on x, y, z, w...
 *
 * calling parameters:
 *
 * references - a vector< vector<NodeSKQ> > of pointers to each of the (x, y, z, w...) tuples
 * temporary - a vector<NodeSKQ> that is used as a temporary array
 * start - start element of the reference arrays
 * end - end element of the reference arrays
 * dim - the number of dimensions
 * depth - the depth in the tree
 *
 * returns: a KdNode pointer to the root of the k-d tree
 */
void Build_cBiK::buildKdTree(vector< vector<NodeSKQ> >& references, vector<NodeSKQ>& temporary, const long start, const long end, const long dim, const long depth) {
    // The axis permutes as x, y, z, w... and addresses the referenced data.
    long axis = depth % dim;
    
    if (end == start) {
        //Se agrega la key al punto, siempre es una hoja
        KdTreeBuild[end] = references[0][end];
        //cout << "Start: " << start << endl;
    } else if (end == start + 1) {
        // Two references were passed to this function in sorted order, so store the start
        // element at this level of the tree and store the end element as the > child.
        KdTreeBuild[start] = references[0][start];
        KdTreeBuild[end] = references[0][end];

        //cout << "Start: " << start << " End: " << end << endl;
    } else if (end == start + 2) {
        // Three references were passed to this function in sorted order, so
        // store the median element at this level of the tree, store the start
        // element as the < child and store the end element as the > child.
        KdTreeBuild[start + 1] = references[0][start + 1];
        KdTreeBuild[start] = references[0][start];
        KdTreeBuild[end] = references[0][end];

        //cout << "Start: " << start << " Start + 1: " << (start+1) << " End: " << end << endl;
    } else if (end > start + 2) {
        // More than three references were passed to this function, so
        // the median element of references[0] is chosen as the tuple about
        // which the other reference arrays will be partitioned.  Avoid
        // overflow when computing the median.
        const long median = start + ((end - start) / 2);
        
        //Se guardan los puntos de la mediana del arreglo
        KdTreeBuild[median] = references[0][median];

        //cout << "Median: " << median << endl;
        
        // Copy references[0] to the temporary array before partitioning.
        for (long i = start; i <= end; i++) {
            temporary[i] = references[0][i];
        }

        // Process each of the other reference arrays in a priori sorted order
        // and partition it by comparing super keys.  Store the result from
        // references[i] in references[i-1], thus permuting the reference
        // arrays.  Skip the element of references[i] that that references
        // a point that equals the point that is stored in the new k-d node.
        long lower, upper, lowerSave, upperSave;
        for (long i = 1; i < dim; i++) {

            // Process one reference array.  Compare once only.
            lower = start - 1;
            upper = median;
            for (long j = start; j <= end; j++) {
                double compare = superKeyCompare(references[i][j].coordenadas, KdTreeBuild[median].coordenadas, axis, dim);
                if (compare < 0) {
                    references[i - 1][++lower] = references[i][j];
                } else if (compare > 0) {
                    references[i - 1][++upper] = references[i][j];
                }
            }

            // Check the new indices for the reference array.
            if (lower < start || lower >= median) {
                cout << "incorrect range for lower at depth = " << depth << " : start = "
                        << start << "  lower = " << lower << "  median = " << median << endl;
                exit(1);
            }

            if (upper <= median || upper > end) {
                cout << "incorrect range for upper at depth = " << depth << " : median = "
                        << median << "  upper = " << upper << "  end = " << end << endl;
                exit(1);
            }

            if (i > 1 && lower != lowerSave) {
                cout << " lower = " << lower << "  !=  lowerSave = " << lowerSave << endl;
                exit(1);
            }

            if (i > 1 && upper != upperSave) {
                cout << " upper = " << upper << "  !=  upperSave = " << upperSave << endl;
                exit(1);
            }

            lowerSave = lower;
            upperSave = upper;
        }

        // Copy the temporary array to references[dim-1] to finish permutation.
        for (long i = start; i <= end; i++) {
            references[dim - 1][i] = temporary[i];
        }

        // Recursively build the < branch of the tree.
        buildKdTree(references, temporary, start, lower, dim, depth + 1);
        
        // Recursively build the > branch of the tree.
        buildKdTree(references, temporary, median + 1, upper, dim, depth + 1);
    } else if (end < start) {

        // This is an illegal condition that should never occur, so test for it last.
        cout << "error has occurred at depth = " << depth << " : end = " << end
                << "  <  start = " << start << endl;
        exit(1);
    }
}

double Build_cBiK::bruteForce() {
    struct timeval iniTime2, endTime2;
    double secs;
    
    //Inicio tiempo
    gettimeofday(&iniTime2, NULL);
    int bigI = 0;
    int bigJ = 0;
    
    double hd = 0;
    for (int i = 0; i < KdTreeA.size(); i++) {
        //cout << "Punto " << (i+1) << endl;
        double distancia = numeric_limits<double>::max();
        int localJ=0;
        for (int j = 0; j < KdTreeB.size(); j++) {
            double distEu = getEuclideanDistance(KdTreeA[i].coordenadas, KdTreeB[j].coordenadas);
            //cout << "P: " << getPoint(KdTreeA[i].coordenadas) << " - Q:" << getPoint(KdTreeA[j].coordenadas) << " - D:" << distEu << endl;
            if(distEu < distancia) {
                distancia = distEu;
                localJ = j;
            }
        }
        if(distancia > hd) {
            bigI = i;
            bigJ = localJ;
            hd = distancia;
        }
    }
    
    gettimeofday(&endTime2, NULL);
    secs = timeval_diff(&endTime2, &iniTime2);
    
    cout << "\t\t\t*** Fuerza bruta: " << hd << " ***" << " Tiempo: " << secs * 1000.0 << " ms" << " i: " << bigI << " j: " << bigJ <<endl;
}

void  Build_cBiK::saveMaxPqueue(unsigned long num) {
    if (num > heapMaxElements){
        heapMaxElements = num;
    }
}
void  Build_cBiK::saveExitHeap(unsigned long num) {
    exitHeapElements = num;
}


double Build_cBiK::hdKD2() {
  //  struct timeval iniTime2, endTime2;
   // double secs;
    
    //Inicio tiempo
    //gettimeofday(&iniTime2, NULL);
    
    double minNN = numeric_limits<double>::max();
    
    //Distancia de Hausdorff
    double cmax = 0;
    int anyPointA[2];
    anyPointA[0] = KdTreeA[0].coordenadas[0];
    anyPointA[1] = KdTreeA[0].coordenadas[1];
    cmax = nearestNeighbor(anyPointA, 0, (KdTreeB.size() - 1), 0, minNN);
    
    cmax = HDKPruning(cmax);
   /* gettimeofday(&endTime2, NULL);
    secs = timeval_diff(&endTime2, &iniTime2);
    
    
   // cout << "********************************************************************" << endl;
   // cout << "Distancia H: " << cmax << endl;
    
    //Muestra el tiempo por cada consulta
    //cout << "Tiempo: " << secs * 1000.0 << " ms" << endl;
    */
    return cmax;
}

double Build_cBiK::hdKD2v2() {
    //  struct timeval iniTime2, endTime2;
    // double secs;

    //Inicio tiempo
    //gettimeofday(&iniTime2, NULL);

    double minNN = numeric_limits<double>::max();

    //Distancia de Hausdorff
    double cmax = 0;
    int anyPointA[2];
    anyPointA[0] = KdTreeA[0].coordenadas[0];
    anyPointA[1] = KdTreeA[0].coordenadas[1];
    cmax = nearestNeighbor(anyPointA, 0, (KdTreeB.size() - 1), 0, minNN);

    cmax = HDKPruningv2(cmax);
    /* gettimeofday(&endTime2, NULL);
     secs = timeval_diff(&endTime2, &iniTime2);


    // cout << "********************************************************************" << endl;
    // cout << "Distancia H: " << cmax << endl;

     //Muestra el tiempo por cada consulta
     //cout << "Tiempo: " << secs * 1000.0 << " ms" << endl;
     */
    return cmax;
}

double Build_cBiK::HDKPruning(double cmax) {
    //Creacion MaxHeap
    MaxHeap pQ;
    
    //Valores KDTree
    int start = 0;
    int end = KdTreeA.size() - 1;
    int depth = 0;
    
    //Region inicial
    int r1[2];
    r1[0] = xMenorA;
    r1[1] = yMenorA;
    int r2[2];
    r2[0] = xMayorA;
    r2[1] = yMayorA;
    
    //Inserta la region inicial con distancia cero
    double d = 0;
  //  double d =  sqrt(pow(xMayorA - xMenorA,2)+ pow(yMayorA-yMenorB,2)); //DISTANCIA DE LA DIAGONAL DEL RECTANGULO.
    pQ.push( Node(start, end, depth, d, r1, r2) );
    
    int mid;
    Node nodoQ;
    int region[4];
    int regionH1[2];
    int regionH2[2];
    double nn;
    while(!pQ.empty()) {
      //  cout << "HEAP SIZE: " << sizeof(pQ) << endl;
        saveMaxPqueue(pQ.size()); // para saber el maximo ocupado por pqueue y reportar el espacio ocupado..
        nodoQ = pQ.top();
        pQ.pop();

       /* if(nodoQ.distance <= cmax){
            return cmax;
        }*/
        mid = (nodoQ.start + nodoQ.end) * 0.5;
        //cout <<"mid is:"<< mid <<endl;
        //cout << "Punto:" << getPoint(KdTreeA[mid].coordenadas) << " INDICES: [" << nodoQ.start << "," << nodoQ.end << "] => " << mid << endl;
        
        //Prepara el punto para consultarlo en el conjunto B
        nn = NNMax(KdTreeA[mid].coordenadas, cmax, numeric_limits<double>::max());
        //cout << "nn: " << nn << " > " << "cmax: " << cmax << endl;
        if(nn > cmax) {
            cmax = nn;
        }
        //Es un nodo interno
        if(nodoQ.start != nodoQ.end) {
            //Si existe un hijo izquierdo
            //cout << "ES UN NODO INTERNO" << endl;
            if( (nodoQ.start <= mid-1) ){
                //cout << "buscando region HI" << endl;
                getRegion(region, nodoQ, true, true);
                regionH1[0] = region[0];
                regionH1[1] = region[1];
                regionH2[0] = region[2];
                regionH2[1] = region[3];
                
                d = isCandidate(regionH1, regionH2, cmax);
                if( d != -1 ) {
                    //cout << "Ingresando HI start:" << nodoQ.start << " end: " << mid-1 << endl;
                    pQ.push( Node(nodoQ.start, mid-1, depth+1, d, regionH1, regionH2) );
                }
            }
            
            //Si existe un hijo derecho
            if( (mid+1 <= nodoQ.end) ){
                //cout << "buscando region HD" << endl;
                getRegion(region, nodoQ, false, true);
                regionH1[0] = region[0];
                regionH1[1] = region[1];
                regionH2[0] = region[2];
                regionH2[1] = region[3];
                
                d = isCandidate(regionH1, regionH2, cmax);
                //cout << "DISTANCIA ENTRE REGIONES: " << d << endl;
                if( d != -1 ) {
                    //cout << "Ingresando HD start:" << mid+1 << " end: " << nodoQ.end << endl;
                    pQ.push( Node(mid+1, nodoQ.end, depth+1, d, regionH1, regionH2) );
                }
            }
        }
    }
    saveExitHeap(pQ.size()); //agregado por fernando.
    return cmax;
}

double Build_cBiK::HDKPruningv2(double cmax) {
    //Creacion MaxHeap
    MaxHeap pQ;

    //Valores KDTree
    int start = 0;
    int end = KdTreeA.size() - 1;
    int depth = 0;

    //Region inicial
    int r1[2];
    r1[0] = xMenorA;
    r1[1] = yMenorA;
    int r2[2];
    r2[0] = xMayorA;
    r2[1] = yMayorA;

    //Inserta la region inicial con distancia cero
    // double d = 0;
    double d =  sqrt(pow(xMayorA - xMenorA,2)+ pow(yMayorA-yMenorB,2)); //DISTANCIA DE LA DIAGONAL DEL RECTANGULO.
    pQ.push( Node(start, end, depth, d, r1, r2) );

    int mid;
    Node nodoQ;
    int region[4];
    int regionH1[2];
    int regionH2[2];
    double nn;
    while(!pQ.empty()) {
        //  cout << "HEAP SIZE: " << sizeof(pQ) << endl;
        saveMaxPqueue(pQ.size()); // para saber el maximo ocupado por pqueue y reportar el espacio ocupado..
        nodoQ = pQ.top();
        pQ.pop();

        if(nodoQ.distance <= cmax){
            saveExitHeap(pQ.size()); //agregado por fernando.
            return cmax;
        }
        mid = (nodoQ.start + nodoQ.end) * 0.5;
        //cout <<"mid is:"<< mid <<endl;
        //cout << "Punto:" << getPoint(KdTreeA[mid].coordenadas) << " INDICES: [" << nodoQ.start << "," << nodoQ.end << "] => " << mid << endl;

        //Prepara el punto para consultarlo en el conjunto B
        nn = NNMax(KdTreeA[mid].coordenadas, cmax, numeric_limits<double>::max());
        //cout << "nn: " << nn << " > " << "cmax: " << cmax << endl;
        if(nn > cmax) {
            cmax = nn;
        }
        //Es un nodo interno
        if(nodoQ.start != nodoQ.end) {
            //Si existe un hijo izquierdo
            //cout << "ES UN NODO INTERNO" << endl;
            if( (nodoQ.start <= mid-1) ){
                //cout << "buscando region HI" << endl;
                getRegion(region, nodoQ, true, true);
                regionH1[0] = region[0];
                regionH1[1] = region[1];
                regionH2[0] = region[2];
                regionH2[1] = region[3];

                d = isCandidate(regionH1, regionH2, cmax);
                if( d != -1 ) {
                    //cout << "Ingresando HI start:" << nodoQ.start << " end: " << mid-1 << endl;
                    pQ.push( Node(nodoQ.start, mid-1, depth+1, d, regionH1, regionH2) );
                }
            }

            //Si existe un hijo derecho
            if( (mid+1 <= nodoQ.end) ){
                //cout << "buscando region HD" << endl;
                getRegion(region, nodoQ, false, true);
                regionH1[0] = region[0];
                regionH1[1] = region[1];
                regionH2[0] = region[2];
                regionH2[1] = region[3];

                d = isCandidate(regionH1, regionH2, cmax);
                //cout << "DISTANCIA ENTRE REGIONES: " << d << endl;
                if( d != -1 ) {
                    //cout << "Ingresando HD start:" << mid+1 << " end: " << nodoQ.end << endl;
                    pQ.push( Node(mid+1, nodoQ.end, depth+1, d, regionH1, regionH2) );
                }
            }
        }
    }
    saveExitHeap(pQ.size()); //agregado por fernando.
    return cmax;
}

double Build_cBiK::isCandidate(int region1[], int region2[], double cmax) {
    MinHeap pR;
    
    int start = 0;
    int end = KdTreeB.size() - 1;
    int depth = 0;
    
    int r1[2];
    r1[0] = xMenorB;
    r1[1] = yMenorB;
    int r2[2];
    r2[0] = xMayorB;
    r2[1] = yMayorB;
    
    //cout << "A: [" << region1[0] << "," << region1[1] << "] [" << region2[0] << "," << region2[1] << "] == ";
    //cout << "B: [" << r1[0] << "," << r1[1] << "] [" << r2[0] << "," << r2[1] << "]" << endl;
    
    double d = MaxMaxDist(region1, region2, r1, r2);
    //cout << "(isCANDIDATE) d: " << d << endl;
    pR.push( Node(start, end, depth, d, r1, r2) );
    
    int mid;
    int region[4];
    int regionH1[2];
    int regionH2[2];
    double maxDist = 0;
    Node nodoQ;
    while(!pR.empty()) {
        nodoQ = pR.top();
        pR.pop();
        mid = (nodoQ.start + nodoQ.end) * 0.5;
        
        if(nodoQ.start == nodoQ.end) {//Es una hoja
          //  cout <<"distance is:"<<nodoQ.distance<<endl;
            return nodoQ.distance;
        }else {//Es un nodo interno
            //Si existe un hijo izquierdo
            //cout << "ES UN NODO INTERNO" << endl;
            if( (nodoQ.start <= mid-1) ){
                getRegion(region, nodoQ, true, false);
                regionH1[0] = region[0];
                regionH1[1] = region[1];
                regionH2[0] = region[2];
                regionH2[1] = region[3];
                
                maxDist = MaxMaxDist(region1, region2, regionH1, regionH2);
                
                //cout << "HI:: MaxDist: " << maxDist << " - " << "cmax: " << cmax << endl;
                
                if(maxDist <= cmax) {
                    return -1;//Regla 4
                }
             //   cout <<"insertando:"<< maxDist<<endl;
                pR.push( Node(nodoQ.start, mid-1, depth+1, maxDist, regionH1, regionH2) );
            }
            
            //Si existe un hijo derecho
            if( (mid+1 <= nodoQ.end) ){
                getRegion(region, nodoQ, false, false);
                regionH1[0] = region[0];
                regionH1[1] = region[1];
                regionH2[0] = region[2];
                regionH2[1] = region[3];
                
                maxDist = MaxMaxDist(region1, region2, regionH1, regionH2);
                
                //cout << "HD:: MaxDist: " << maxDist << " - " << "cmax: " << cmax << endl;
                
                if(maxDist <= cmax) {
                    return -1;//Regla 4
                }
                //cout <<"insertando:"<< maxDist<<endl;

                pR.push( Node(mid+1, nodoQ.end, depth+1, maxDist, regionH1, regionH2) );
            }
        }
    }
    
    return d;
}

double Build_cBiK::MaxMaxDist(int* region1A, int* region2A, int* region1B, int* region2B) {
    //cout << "A: [" << region1A[0] << "," << region1A[1] << "] [" << region2A[0] << "," << region2A[1] << "] == ";
    //cout << "B: [" << region1B[0] << "," << region1B[1] << "] [" << region2B[0] << "," << region2B[1] << "]" << endl;
    
    int q[2];
    int r[2];
    double mayor = 0;
    double candidate = 0;
    
    //INF-IZQ de A
    r[0] = region1A[0];
    r[1] = region1A[1];
    q[0] = ( abs(r[0]-region1B[0]) > abs(r[0]-region2B[0]) ) ? region1B[0] : region2B[0];
    q[1] = ( abs(r[1]-region1B[1]) > abs(r[1]-region2B[1]) ) ? region1B[1] : region2B[1];
    
    mayor = getEuclideanDistance(r, q);
    //cout << "MAYOR 1: " << to_string(mayor) << " region1A[" << region1A[0] << "," << region1A[1] << "]" << " Q[" << q[0] << "," << q[1] << "]" << endl;
    
    //SUP-IZQ de A
    r[0] = region1A[0];
    r[1] = region2A[1];
    q[0] = ( abs(r[0]-region1B[0]) > abs(r[0]-region2B[0]) ) ? region1B[0] : region2B[0];
    q[1] = ( abs(r[1]-region1B[1]) > abs(r[1]-region2B[1]) ) ? region1B[1] : region2B[1];
    
    candidate = getEuclideanDistance(r, q);
    if(candidate > mayor) {
        mayor = candidate;
    }
    
    //INF-DER de A
    r[0] = region2A[0];
    r[1] = region1A[1];
    q[0] = ( abs(r[0]-region1B[0]) > abs(r[0]-region2B[0]) ) ? region1B[0] : region2B[0];
    q[1] = ( abs(r[1]-region1B[1]) > abs(r[1]-region2B[1]) ) ? region1B[1] : region2B[1];
    
    candidate = getEuclideanDistance(r, q);
    if(candidate > mayor) {
        mayor = candidate;
    }
    
    //SUP-DER de A
    r[0] = region2A[0];
    r[1] = region2A[1];
    q[0] = ( abs(r[0]-region1B[0]) > abs(r[0]-region2B[0]) ) ? region1B[0] : region2B[0];
    q[1] = ( abs(r[1]-region1B[1]) > abs(r[1]-region2B[1]) ) ? region1B[1] : region2B[1];
    
    candidate = getEuclideanDistance(r, q);
    if(candidate > mayor) {
        mayor = candidate;
    }
    
    return mayor;
}

//(Izquierda o Derecha) iod => true -> Izquierda y false -> Derecha
//(Conjunto A o B) aob => true -> conjuntoA y false -> conjuntoB
void Build_cBiK::getRegion(int region[], Node nodoQ, bool iod, bool aob) {
    int r1Child[2];
    int r2Child[2];
    
    int mid = ((nodoQ.start + nodoQ.end) * 0.5);
    int p = nodoQ.depth % 2;
    
    //cout << "\t\tRegion recibida MID:" << mid << " R1:" << getPoint(nodoQ.r1) << " R2: " << getPoint(nodoQ.r2) << endl;
    
    if(p == 0) {//Eje X
        if(iod) {//subespacio izquierdo
            //r1 se mantiene
            r1Child[0] = nodoQ.r1[0];
            r1Child[1] = nodoQ.r1[1];
            //cout << "R1: " << nodoQ.r1[0] << ", " << nodoQ.r1[1] << endl;

            //r2 debe actualizarce como una combinacion del eje de particion junto con r2
            r2Child[0] = aob ? KdTreeA[mid].coordenadas[0] : KdTreeB[mid].coordenadas[0];
            r2Child[1] = nodoQ.r2[1];
        }else {//subespacio derecho
            r1Child[0] = aob ? KdTreeA[mid].coordenadas[0] : KdTreeB[mid].coordenadas[0];
            r1Child[1] = nodoQ.r1[1];

            //se mantiene
            r2Child[0] = nodoQ.r2[0];
            r2Child[1] = nodoQ.r2[1];
        }
    }else {//Eje Y
        if(iod) {//subespacio inferior
            //r1 se mantiene
            r1Child[0] = nodoQ.r1[0];
            r1Child[1] = nodoQ.r1[1];
            
            r2Child[0] = nodoQ.r2[0];
            r2Child[1] = aob ? KdTreeA[mid].coordenadas[1] : KdTreeB[mid].coordenadas[1];
        }else {//subespacio superior
            //r1 se mantiene
            r1Child[0] = nodoQ.r1[0];
            r1Child[1] = aob ? KdTreeA[mid].coordenadas[1] : KdTreeB[mid].coordenadas[1];
            
            r2Child[0] = nodoQ.r2[0];
            r2Child[1] = nodoQ.r2[1];
        }
    }
    
    region[0] = r1Child[0];
    region[1] = r1Child[1];
    region[2] = r2Child[0];
    region[3] = r2Child[1];
}

/*
double Build_cBiK::nearestNeighbor(int query[], int start, int end, const int depth, double minNN) {
    //Determina el eje x - y dada la profundidad
    const int p = depth % 2;//Debido a las dos dimensiones
    int mid = (start + end) * 0.5;
    
//    cout << "Punto visitado: ";
//    printTuple(KdTreeB[mid].coordenadas);
//    cout << endl;
    
    //Es un nodo interno
    //cout << "Start: " << start << " -- End: " << end << " -- Mid: " << mid << endl;
    if(start != end) {
        //Busca por la izquierda
        //cout << "BI:: " << endl;
        if((start <= mid-1) && query[p] <= KdTreeB[mid].coordenadas[p]) {
            minNN = nearestNeighbor(query, start, mid-1, depth + 1, minNN);

            //Busca si intersecta un subespacio mediante la circunferencia
            if ( abs(KdTreeB[mid].coordenadas[p] - query[p])*abs(KdTreeB[mid].coordenadas[p] - query[p]) < minNN*minNN ) {
                //Busca por la derecha
                minNN = nearestNeighbor(query, mid+1, end, depth + 1, minNN);
            }
        }
        
        //Busca por la derecha
        if((mid+1 <= end) && query[p] >= KdTreeB[mid].coordenadas[p]) {
            minNN = nearestNeighbor(query, mid+1, end, depth + 1, minNN);

            //Busca si intersecta un subespacio mediante la circunferencia
            if ( abs(KdTreeB[mid].coordenadas[p] - query[p])*abs(KdTreeB[mid].coordenadas[p] - query[p]) < minNN*minNN ) {
                //Busca por la izquierda
                minNN = nearestNeighbor(query, start, mid-1, depth + 1, minNN);
            }
        }
    }
    
    //Revisa el punto
    if ( getEuclideanDistance(query, KdTreeB[mid].coordenadas) < minNN ) {
        minNN = getEuclideanDistance(query, KdTreeB[mid].coordenadas);
    }

    return minNN;
}
*/
double Build_cBiK::nearestNeighbor(int query[], int start, int end, const int depth, double minNN) {
    // Determina el eje x - y dada la profundidad
    const int p = depth % 2; // Debido a las dos dimensiones

    // Cálculo del índice medio usando división entera
    int mid = (start + end) / 2;

    // Es un nodo interno
    if (start <= end) {
        // Busca por la izquierda
        if (start <= mid - 1 && query[p] <= KdTreeB[mid].coordenadas[p]) {
            minNN = nearestNeighbor(query, start, mid - 1, depth + 1, minNN);

            // Busca si intersecta un subespacio mediante la circunferencia
            if (std::abs(KdTreeB[mid].coordenadas[p] - query[p]) < minNN) {
                // Busca por la derecha
                minNN = nearestNeighbor(query, mid + 1, end, depth + 1, minNN);
            }
        }

        // Busca por la derecha
        if (mid + 1 <= end && query[p] >= KdTreeB[mid].coordenadas[p]) {
            minNN = nearestNeighbor(query, mid + 1, end, depth + 1, minNN);

            // Busca si intersecta un subespacio mediante la circunferencia
            if (std::abs(KdTreeB[mid].coordenadas[p] - query[p]) < minNN) {
                // Busca por la izquierda
                minNN = nearestNeighbor(query, start, mid - 1, depth + 1, minNN);
            }
        }
    }

    // Revisa el punto
    double distance = getEuclideanDistance(query, KdTreeB[mid].coordenadas);
    if (distance < minNN) {
        minNN = distance;
    }

    return minNN;
}

double Build_cBiK::hdKD1() {
    struct timeval iniTime2, endTime2;
    double secs;
    
    //Inicio tiempo
    //gettimeofday(&iniTime2, NULL);
    
    //pNN es el punto mas cercano
    pNN[0] = numeric_limits<int>::max();
    pNN[1] = numeric_limits<int>::max();
    
    //Distancia de Hausdorff
    double cmax = 0;
    double minNN;
    //Por cada punto de A
    //cout << "CALCULANDO HD" << endl;
    for (int i = 0; i < KdTreeA.size(); i++) {
        minNN = getEuclideanDistance(KdTreeA[i].coordenadas, pNN);
        
//        cout << "\n############## {"<< i <<"} Dist de p: [";
//        printTuple(KdTreeA[i].coordenadas);
//        cout << "] a pNN: [";
//        printTuple(pNN);
//        cout << "] == minNN: " <<  minNN << " > cmax actual:" << cmax << " ##############" << endl;
        
        //Regla 1 
        if(minNN > cmax) {
            //Busca en B el vecino mas cercano a p en A
            cmax = NNMax(KdTreeA[i].coordenadas, cmax, minNN);
        }
    }
  /*  gettimeofday(&endTime2, NULL);
    secs = timeval_diff(&endTime2, &iniTime2);
    
    cout << "********************************************************************" << endl;
    cout << "Distancia H: " << cmax << endl;
    
    //Muestra el tiempo por cada consulta
    cout << "Tiempo: " << secs * 1000.0 << " ms" << endl;
    */
  //  return secs * 1000.0;//retorna las milesimas
    return cmax;  //puesto por fernando
}

double Build_cBiK::NNMax(int pointP[], double cmax, double minNN) {
    //PriorityQueue pQ;
    MinHeap pQ;
    
    int start = 0;
    int end = KdTreeB.size() - 1;
    int depth = 0;
    
    int r1[2];
    r1[0] = xMenorB;
    r1[1] = yMenorB;
    int r2[2];
    r2[0] = xMayorB;
    r2[1] = yMayorB;
    
    double d = maxDist(pointP, r1, r2, (depth%2));
    pQ.push( Node(start, end, depth, d, r1, r2) );
    
    int pNN1[2];
    pNN1[0] = pNN[0];
    pNN1[1] = pNN[1];
    int mid;
    double dist;
    Node nodoQ;
    while(!pQ.empty()) {
        nodoQ = pQ.top();
        pQ.pop();
        mid = (nodoQ.start + nodoQ.end) * 0.5;
        
        dist = getEuclideanDistance(KdTreeB[mid].coordenadas, pointP);
        
        //Regla 1
        //cout << "<<Rule 1 {" << mid << "}>> Dist (Q" << getPoint(KdTreeB[mid].coordenadas) << " to P" << getPoint(pointP) << "): " << dist << " <= cmax: " << cmax << endl;
        if(dist <= cmax) {
            return cmax;
        }else {
            //cout << "\tEvalua dist: " << dist << " < minNN: " << minNN << endl;
            if(dist < minNN) {
                //cout << "\t\tnuevo MinNN: " << dist << " == Punto " << "{" << mid << "}: " << getPoint(KdTreeB[mid].coordenadas) << endl;
                minNN = dist;
                pNN1[0] = KdTreeB[mid].coordenadas[0];
                pNN1[1] = KdTreeB[mid].coordenadas[1];
            }
        }
        
        //Si es un nodo interno, se agregan los subespacios
        if(nodoQ.start != nodoQ.end) {
            //Si existe un hijo izquierdo
            if( (nodoQ.start <= mid-1) && !update(pQ, nodoQ, cmax, minNN, pointP, true) ){
                //cout << "RETURN LEFT - cmax: " << cmax << endl;
                return cmax;
            }
            
            //Si existe un hijo derecho
            if( (mid+1 <= nodoQ.end) && !update(pQ, nodoQ, cmax, minNN, pointP, false) ){
                //cout << "RETURN RIGHT - cmax: " << cmax << endl;
                return cmax;
            }
        }
    }
    pNN[0] = pNN1[0];
    pNN[1] = pNN1[1];
    
    //cout << "\t\t\t\tRETURN MIN-NN: " << minNN << endl;
    return minNN;
}

bool Build_cBiK::update(MinHeap &pQ, Node nodoQ, double cmax, double minNN, int* pointP, bool iod) {
    int r1Child[2];
    int r2Child[2];
    
    int mid = ((nodoQ.start + nodoQ.end) * 0.5);
    int p = nodoQ.depth%2;
    
    //cout << "\t\tRegion recibida MID:" << mid << " R1:" << getPoint(nodoQ.r1) << " R2: " << getPoint(nodoQ.r2) << endl;
    
    if(p == 0) {//Eje X
        if(iod) {//subespacio izquierdo
            //r1 se mantiene
            r1Child[0] = nodoQ.r1[0];
            r1Child[1] = nodoQ.r1[1];
            //cout << "R1: " << nodoQ.r1[0] << ", " << nodoQ.r1[1] << endl;

            //r2 debe actualizarce como una combinacion del eje de particion junto con r2
            r2Child[0] = KdTreeB[mid].coordenadas[0];
            r2Child[1] = nodoQ.r2[1];
        }else {//subespacio derecho
            r1Child[0] = KdTreeB[mid].coordenadas[0];
            r1Child[1] = nodoQ.r1[1];

            //se mantiene
            r2Child[0] = nodoQ.r2[0];
            r2Child[1] = nodoQ.r2[1];
        }
    }else {//Eje Y
        if(iod) {//subespacio inferior
            //r1 se mantiene
            r1Child[0] = nodoQ.r1[0];
            r1Child[1] = nodoQ.r1[1];
            
            r2Child[0] = nodoQ.r2[0];
            r2Child[1] = KdTreeB[mid].coordenadas[1];
        }else {//subespacio superior
            //r1 se mantiene
            r1Child[0] = nodoQ.r1[0];
            r1Child[1] = KdTreeB[mid].coordenadas[1];
            
            r2Child[0] = nodoQ.r2[0];
            r2Child[1] = nodoQ.r2[1];
        }
    }
    //cout << "\t\tRegion evaluada MID:" << mid << " R1:" << getPoint(r1Child) << " R2: " << getPoint(r2Child) << endl;
    
    double distMax = maxDist(pointP, r1Child, r2Child, p);
    
    //Regla 2
    //cout << "\t<<Rule 2>>" << " maxDist: " << distMax << " <= cmax: " << cmax << " Axis[" << p << "]" << endl;
    if(distMax <= cmax) {
        return false;
    }
    
    //Regla 3
    double distMin = minDist(pointP, r1Child, r2Child, p);
    //cout << "\t<<Rule 3>>" << " minDist: " << distMin << " <= minNN: " << minNN << endl;
    if(distMin <= minNN) {
        //cout << "\t\t" << (iod ? "[Izq]" : "[Der]") << " d: " << distMax << endl;
        
        //cout << "\t\tINSERTANDO NODO R1:" << getPoint(r1Child) << " - R2: " << getPoint(r2Child) << endl;
        if(iod) {//subespacio izquierdo
            pQ.push( Node(nodoQ.start, mid-1, nodoQ.depth+1, distMax, r1Child, r2Child) );
        }else {
            pQ.push( Node(mid+1, nodoQ.end, nodoQ.depth+1, distMax, r1Child, r2Child) );
        }
    }
    
    return true;
}

double Build_cBiK::minDist(const int* q, int* r1, int* r2, int p) {
    if(q[p] < r1[p]){
        return abs(q[p] - r1[p]);//*abs(q[p] - menor);
    }else if(q[p] > r2[p]) {
        return abs(q[p] - r2[p]);//*abs(q[p] - mayor);
    }else {
        return q[p] - q[p];
    }
}

double Build_cBiK::maxDist(int* q, int* r1, int* r2, int p) {
    int farthest[2];
    int otherAxis = p == 0 ? 1 : 0;
    
    if(q[p] <= r1[p]) {//Descarta eje de R1 ya que está más cerca
        farthest[p] = r2[p];
        
        if(q[otherAxis] <= r1[otherAxis]) {
            farthest[otherAxis] = r2[otherAxis];
        }else if(q[otherAxis] >= r2[otherAxis]) {
            farthest[otherAxis] = r1[otherAxis];
        }else {//Busca el más lejano entre ambos
            farthest[otherAxis] = (abs(q[otherAxis]-r1[otherAxis]) > abs(q[otherAxis]-r2[otherAxis])) ? r1[otherAxis] : r2[otherAxis];
        }
        
        return getEuclideanDistance(q, farthest);
    }else if(q[p] >= r2[p]) {//Descarta eje de R2 ya que está más cerca
        farthest[p] = r1[p];
        
        if(q[otherAxis] <= r1[otherAxis]) {
            farthest[otherAxis] = r2[otherAxis];
        }else if(q[otherAxis] >= r2[otherAxis]) {
            farthest[otherAxis] = r1[otherAxis];
        }else {//Busca el más lejano entre ambos
            farthest[otherAxis] = (abs(q[otherAxis]-r1[otherAxis]) > abs(q[otherAxis]-r2[otherAxis])) ? r1[otherAxis] : r2[otherAxis];
        }
        
        return getEuclideanDistance(q, farthest);
    }else {//Está entre R1 Y R2, debe calcularlos todos
        farthest[p] = (abs(q[p]-r1[p]) > abs(q[p]-r2[p])) ? r1[p] : r2[p];
        farthest[otherAxis] = (abs(q[otherAxis]-r1[otherAxis]) > abs(q[otherAxis]-r2[otherAxis])) ? r1[otherAxis] : r2[otherAxis];
        
        return getEuclideanDistance(q, farthest);
    }
    
    /*int supIzq[2];
    int supDer[2];
    int infIzq[2];
    int infDer[2];
    
    if(r1[0] < r2[0]) {
        infIzq[0] = r1[0];
        supIzq[0] = r1[0];
        
        infDer[0] = r2[0];
        supDer[0] = r2[0];
    }else {
        infIzq[0] = r2[0];
        supIzq[0] = r2[0];
        
        infDer[0] = r1[0];
        supDer[0] = r1[0];
    }
    
    if(r1[1] < r2[1]) {
        infIzq[1] = r1[1];
        supIzq[1] = r2[1];
        
        infDer[1] = r1[1];
        supDer[1] = r2[1];
    }else {
        infIzq[1] = r2[1];
        supIzq[1] = r1[1];
        
        infDer[1] = r2[1];
        supDer[1] = r1[1];
    }
    
    //cout << "\t\t\tREGIONES: " << getPoint(infIzq) << getPoint(supDer) << getPoint(supIzq) << getPoint(infDer) << endl;
    //cout << "Q: " << getPoint(q) << "  EVALUANDO: " << getPoint(r1) << " y " << getPoint(r2) << " MAXDIST: ";
    
    double d1 = getEuclideanDistance(q, infIzq);
    double d2 = getEuclideanDistance(q, supDer);
    double d3 = getEuclideanDistance(q, supIzq);
    double d4 = getEuclideanDistance(q, infDer);
    double elMayor = 0;
    
    elMayor = (d1 > d2) ? d1 : d2;
    elMayor = (elMayor > d3) ? elMayor : d3;
    elMayor = (elMayor > d4) ? elMayor : d4;
    
    return elMayor;*/
}

double Build_cBiK::getEuclideanDistance(int p1[], int p2[]) {
    //calcula la distancia para dos dimensiones
    //cout << "(ED) " << getPoint(p1) << " == " << getPoint(p2) << endl;
    //cout << abs(p2[0] - p1[0]) << "," << abs(p2[1] - p1[1]) << endl;
    //cout << "GED: " << sqrt( fabs(p2[0] - p1[0])*fabs(p2[0] - p1[0]) + fabs(p2[1] - p1[1])*fabs(p2[1] - p1[1]) ) << endl;
    return sqrt( fabs(p2[0] - p1[0])*fabs(p2[0] - p1[0]) + fabs(p2[1] - p1[1])*fabs(p2[1] - p1[1]) );
}

double Build_cBiK::timeval_diff(struct timeval *a, struct timeval *b){
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

void Build_cBiK::printTree() {
    cout << endl;
    cout << "************************** ARBOL KD-TREE A ***************************" << endl;
    printKdTree(KdTreeA, 0, (KdTreeA.size() - 1), 0);
    cout << "********************************************************************" << endl;
    cout << endl;
    
    cout << endl;
    cout << "************************** ARBOL KD-TREE B ***************************" << endl;
    printKdTree(KdTreeB, 0, (KdTreeB.size() - 1), 0);
    cout << "********************************************************************" << endl;
    cout << endl;
}

void Build_cBiK::printKdTree(vector<NodeSKQ>& tree, int start, int end, int depth) {
    int mid = (start + end) / 2;

    if (end > mid) {
        printKdTree(tree, mid + 1, end, depth + 1);
    }

    for (int i = 0; i < depth; i++) {
        cout << "         ";
    }
    cout << "[" << start << "," << end << "] IDX:" << mid;
    cout << " (";
    printTuple(tree[mid].coordenadas);
    cout << ") " << endl;
    
    if (start < mid) {
        printKdTree(tree, start, mid - 1, depth + 1);
    }
}

void Build_cBiK::printTuple(int* tuple) {
    //printf("%.16g, %.16g", tuple[0], tuple[1]); 
    cout << "(" << tuple[0] << "," << tuple[1] << ")";
}

string Build_cBiK::getPoint(int* tuple) {
    
    return "[" + to_string(tuple[0]) + ", " + to_string(tuple[1]) + "]" ;
}

int Build_cBiK::getHeapMaxElements(){
    return heapMaxElements;
}
int Build_cBiK::getHeapExitElements(){
    return exitHeapElements;
}