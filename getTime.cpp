//
// Created by fernando-lenovo on 30-07-24.
//
#include <iostream>
#include <vector>
#include <sstream>   // For std::stringstream
#include <fstream>
#include <cstring>
#include "TimeMesure.h"
#include "utils.h"
#include "Build_cBiK.h"
#include "NodeSKQ.h"
#include <dirent.h> //Agregado por Fernando
#include "ryu-kamata.h" //Agregado por Fernando
#include "MemoryMeasure.h"

int cuentaPuntos(string ruta1) {
    ifstream set1;
    set1.open(ruta1);
    string linea1;
    int nelementos1 = 0;
    //Necesitamos saber cuantos elementos tiene el
    while (getline(set1, linea1)) {
        nelementos1++;
    }
    //cout << nelementos1 << endl;
    set1.close();
    return nelementos1;
}
void printHeader() {
    printf("%-19s|%-8s|%-12s|%-12s|"
           "%-15s|%-15s|%-15s|%-15s|"
           "%-12s|%-7s|%-15s|%-15s|%-15s|%-15s|%-15s|\n",
           "SET NAME", "CATEGORY", "ALGORITHM", "DISTRIBUTION",
           "BUILD TIME", "QUERY TIME", "QUERY WALL TIME", "EXTRACT TIME",
           "Haus Dist", "LAMBDA", "STORAGE", "PEAK REAL MEM", "PEAK VIRT MEM",
           "MXHEAP", "EXIT HEAP"); // Añadir nuevos encabezados
}

void printStats(const Stats &stats) {
    printf("%-19s|%-8s|%-12s|%-12s|%15.01f|%15.01f|%15.01f|%15.01f|%12.01f|%7d|%15.01f|%15lu|%15lu|%15d|%15d|\n",
           stats.dataset.c_str(), stats.datasetCategory.c_str(),
           stats.algoritmo.c_str(), stats.distribution.c_str(),
           stats.buildTime, stats.queryTime, stats.queryTimeWall, stats.extractTime,
           stats.hauss_distance, stats.lambda, stats.storage,
           stats.peakRealMem, stats.peakVirtMem, stats.mxheap, stats.exit_heap); // Imprimir nuevos parámetros
}
std::vector<Point> llenaQueryPoints(char *file) {
    std::ifstream input;
    uint32_t a, b;
    std::vector<Point> puntos;

    input.open(file);
    if (!input.is_open()) {
        std::cout << "ERROR: THE FILE " << file << " DOES NOT EXIST" << std::endl;
        exit(211);
    }
    int i = 0;
    while (input >> a >> b) {
        Point nuevoPunto;
        nuevoPunto.x = a;
        nuevoPunto.y = b;
        puntos.push_back(nuevoPunto);
    }
    input.close();
    return puntos;
}
std::string extractFilename(const std::string &path) {
    size_t pos = path.find_last_of('/');
    if (pos == std::string::npos) {
        return path; // Path doesn't contain a '/' separator
    } else {
        return path.substr(pos + 1); // Extract filename from position after '/'
    }
}
std::string extractCategoryName(const std::string &path) {
    // Encuentra la última ocurrencia de "sets"
    size_t startPos = path.find("sets");
    if (startPos != std::string::npos) {
        startPos += 0; // Mueve el puntero después de "sets"

        // Encuentra la siguiente barra '/' después de "sets"
        size_t endPos = path.find('/', startPos);
        if (endPos != std::string::npos) {
            return path.substr(startPos, endPos - startPos);
        } else {
            // Si no se encuentra '/', toma el resto de la cadena
            return path.substr(startPos);
        }
    }
    return ""; // Retorna una cadena vacía si "sets" no se encuentra
}
std::string extractDistribution(const std::string &path) {
    size_t pos = path.find("Conjunto");
    if (pos != std::string::npos) {
        // Encuentra la posición del último número después de "Conjunto"
        pos = path.find_first_of("0123456789", pos);
        if (pos != std::string::npos && pos + 1 < path.length()) {
            char distChar = path[pos + 1]; // Captura el carácter justo después del número

            // Retorna el nombre completo de la distribución según la letra
            switch (distChar) {
                case 'R':
                    return "RANDOM";
                case 'M':
                    return "MIXTA";
                case 'G':
                    return "GAUSS";
                case 'T':
                    return "REAL";
                default:
                    return "UNKNOWN"; // En caso de que no sea R, M o G
            }
        }
    }
    return "UNKNOWN"; // Retorna "UNKNOWN" si no se encuentra la distribución
}
void imprimeResultados(Stats & resultados){
    char *homedir = getenv("HOME");
    strcat(homedir, "/RESULTADOS/resultHauss.txt");
    if (checkFileExist(homedir) == 1 && fileIsEmpty(homedir) == 1) {
        printHeader();
        fflush(stdout);
    }
    printStats(resultados);
    fflush(stdout);

}
void printVector(const std::vector<Point>& vec) {
    std::cout << "Contenido del vector:" << std::endl;
    for (const auto& p : vec) {
        std::cout << "(" << p.x << ", " << p.y << ")"<<endl;
    }
    std::cout << std::endl;
}
void HDKD1(const std::string &dataSetA, const std::string &dataSetB, Stats &resultados){
    Cronometer *crono = cCronometer();

    Build_cBiK cBiK1 = Build_cBiK();
    Build_cBiK cBiK2 = Build_cBiK();

    cBiK1.load_cBiK(extractFilename(dataSetA), extractFilename(dataSetB));
    cBiK2.load_cBiK(extractFilename(dataSetB), extractFilename(dataSetA));

    double time = 0.0;
    start_clock(crono);
        auto haussCompleto = max(cBiK1.hdKD1(),cBiK2.hdKD1());
    time = (double) (stop_clock(crono) * 1000000.0);
    resultados.algoritmo = "KD-HDK1";
    resultados.hauss_distance = haussCompleto;
    resultados.queryTime = time;
    resultados.queryTimeWall =  wallTime(crono)* 1000000.0;
    resultados.storage = cBiK1.getStorageSize() + cBiK2.getStorageSize();
    resultados.lambda = 0; //no ocupa.
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    imprimeResultados(resultados);
}
void HDKD2(const std::string &dataSetA, const std::string &dataSetB,Stats &resultados){
    Cronometer *crono = cCronometer();

    Build_cBiK cBiK1 = Build_cBiK();
    Build_cBiK cBiK2 = Build_cBiK();

    cBiK1.load_cBiK(extractFilename(dataSetA), extractFilename(dataSetB));
    cBiK2.load_cBiK(extractFilename(dataSetB), extractFilename(dataSetA));

    double time = 0.0;
    start_clock(crono);
        auto haussCompleto = max(cBiK1.hdKD2(),cBiK2.hdKD2());
    time = (double) (stop_clock(crono) * 1000000.0);
    resultados.algoritmo = "KD-HDK2";
    resultados.hauss_distance = haussCompleto;
    resultados.queryTime = time;
    resultados.queryTimeWall =  wallTime(crono)* 1000000.0;
    resultados.storage = cBiK1.getStorageSize() + cBiK2.getStorageSize();
    resultados.lambda = 0;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    resultados.exit_heap = cBiK1.getHeapExitElements()+cBiK2.getHeapExitElements();
    resultados.mxheap = cBiK1.getHeapMaxElements()+cBiK2.getHeapMaxElements();
    imprimeResultados(resultados);
}
void HDKD2v2(const std::string &dataSetA, const std::string &dataSetB,Stats &resultados){
    Cronometer *crono = cCronometer();

    Build_cBiK cBiK1 = Build_cBiK();
    Build_cBiK cBiK2 = Build_cBiK();

    cBiK1.load_cBiK(extractFilename(dataSetA), extractFilename(dataSetB));
    cBiK2.load_cBiK(extractFilename(dataSetB), extractFilename(dataSetA));

    double time = 0.0;
    start_clock(crono);
    auto haussCompleto = max(cBiK1.hdKD2v2(),cBiK2.hdKD2v2());
    time = (double) (stop_clock(crono) * 1000000.0);
    resultados.algoritmo = "KD-HDK2v2";
    resultados.hauss_distance = haussCompleto;
    resultados.queryTime = time;
    resultados.queryTimeWall =  wallTime(crono)* 1000000.0;
    resultados.storage = cBiK1.getStorageSize() + cBiK2.getStorageSize();
    resultados.lambda = 0;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    resultados.exit_heap = cBiK1.getHeapExitElements()+cBiK2.getHeapExitElements();
    resultados.mxheap = cBiK1.getHeapMaxElements()+cBiK2.getHeapMaxElements();
    imprimeResultados(resultados);
}

void KAMATA(const std::string &dataSetA, const std::string &dataSetB, int cantPuntosSet1, int cantPuntosSet2,Stats &resultados){
    Cronometer *crono = cCronometer();
    Build_cBiK KDTreeA = Build_cBiK();
    Build_cBiK KDTreeB = Build_cBiK();

    KDTreeA.load_cBiK(extractFilename(dataSetA), extractFilename(dataSetB));
    KDTreeB.load_cBiK(extractFilename(dataSetB), extractFilename(dataSetA));

    double extrTime = 0.0;
    start_clock(crono);
    auto dataset1 = KDTreeA.extractPointKDTree();
    auto dataset2 = KDTreeB.extractPointKDTree();
    extrTime = (double) (stop_clock(crono) * 1000000.0);
    double time = 0.0;
    start_clock(crono);
    auto haussCompleto = hausKamata(dataset1,dataset2, 3);
    time = (double) (stop_clock(crono) * 1000000.0);

    resultados.algoritmo = "KD-KAMATA";
    resultados.queryTime = time;
    resultados.extractTime = extrTime;
    resultados.queryTimeWall =  wallTime(crono)* 1000000.0;
    resultados.storage= (((sizeof(int)*2)*cantPuntosSet1) + ((sizeof(int)*2)*cantPuntosSet2)) + (sizeof(double)*cantPuntosSet1 + sizeof(double)*cantPuntosSet2) ;
    resultados.hauss_distance = haussCompleto;
    resultados.lambda = 3;
    resultados.peakRealMem = getResidentMemory();
    resultados.peakVirtMem = getVirtualMemory();
    imprimeResultados(resultados);
}
void ejecutaExperimentos(const std::string &dataset1, const std::string &dataset2, int experimento) {
    Cronometer *crono = cCronometer();
    Stats resultados; // Para guardar el resultado de los experimentos.

    int cantPuntosSet1 = cuentaPuntos(dataset1);
    int cantPuntosSet2 = cuentaPuntos(dataset2);
/*
    Build_cBiK KDTreeA, KDTreeB;
    double buildTime = 0.0;


    // Construir solo si se indica
    if (construir == 2) {
        start_clock(crono);
        KDTreeA = Build_cBiK(dataset1);
        KDTreeB = Build_cBiK(dataset2);
        buildTime = (double)(stop_clock(crono) * 1000000.0);

        resultados.dataset = extractFilename(dataset1);
        resultados.datasetCategory = extractCategoryName(dataset1);
        resultados.distribution = extractDistribution(resultados.dataset);
        resultados.buildTime = buildTime;
    }
*/
    resultados.dataset = extractFilename(dataset1);
    resultados.datasetCategory = extractCategoryName(dataset1);
    resultados.distribution = extractDistribution(resultados.dataset);

    // Ejecución del experimento según el tipo especificado
    switch (experimento) {
        case 0: // HDKD1
            HDKD1(dataset1, dataset2,resultados);
            break;
        case 1: // HDKD2
            HDKD2(dataset1, dataset2, resultados);
            break;
        case 2: // KAMATA
            KAMATA(dataset1, dataset2, cantPuntosSet1, cantPuntosSet2, resultados);
            break;
        case 3: // HDKD2 version 2
            HDKD2v2(dataset1, dataset2, resultados);
            break;
        default:
            std::cerr << "Experimento desconocido." << std::endl;
            return;
    }

    resultados.reset();
}
/*
void ejecutaExperimentos(const std::string &dataset1, const std::string &dataset2){

    Cronometer *crono = cCronometer();
    Stats resultados; //para guardar el resultado de los experimentos.

    int cantPuntosSet1 = cuentaPuntos(dataset1);
    int cantPuntosSet2 = cuentaPuntos(dataset2);

    */
/* para crear kdtree *//*

    double buildTime = 0.0;
    start_clock(crono);
        Build_cBiK KDTreeA = Build_cBiK(dataset1);
        Build_cBiK KDTreeB = Build_cBiK(dataset2);
    buildTime = (double) (stop_clock(crono) * 1000000.0);
    resultados.dataset = extractFilename(dataset1);
    resultados.datasetCategory = extractCategoryName(dataset1);
    resultados.distribution = extractDistribution(resultados.dataset);
    resultados.buildTime = buildTime;

    //EJECUCION DE LAS CONSULTAS...

    HDKD1(dataset1,dataset2,KDTreeA,KDTreeB, resultados);
    resultados.reset();
    HDKD2(dataset1,dataset2,KDTreeA,KDTreeB, resultados);
    resultados.reset();
    KAMATA(KDTreeA,KDTreeB, cantPuntosSet1, cantPuntosSet2,resultados);
    resultados.reset();
}
*/
/*
int main(int argc, char **argv) {

    if (argc < 1) {
        printf("%s <PATH DATASETS> \n", argv[0]);
        return -1;
    }

    string path = argv[1];

    std::vector<std::string> distributions = {"R", "G"};
    std::vector<std::string> distributionMix = {"M"};
    int start = 4;
    int end = 7;
    DIR *dir;

    if ((dir = opendir(path.c_str())) != NULL) {

        //para distribucion gauss y random
        for (const auto &dist: distributions) {
            for (int i = start; i <= end; ++i) {
                std::string dataset1 = path + "set1Conjunto" + std::to_string(i) + dist + ".csv";
                std::string dataset2 = path + "set2Conjunto" + std::to_string(i) + dist + ".csv";
                ejecutaExperimentos(dataset1, dataset2);
            }
        }

        //solo para mix distribution
        for (const auto &dist: distributionMix) {
            for (int i = start; i <= end; ++i) {
                std::string dataset1 = path + "set1Conjunto" + std::to_string(i) + dist + ".csv";
                std::string dataset2 = path + "set2Conjunto" + std::to_string(i) + dist + ".csv";
                ejecutaExperimentos(dataset1, dataset2);
            }
        }
        closedir(dir);
    } else {
        perror("No se pudo abrir el directorio");
        return EXIT_FAILURE;
    }

    return 0;
}
*/
void construirEstructura(const std::string &dataset1, const std::string &dataset2){
    Build_cBiK KDTreeA, KDTreeB;

    KDTreeA = Build_cBiK(dataset1);
    KDTreeB = Build_cBiK(dataset2);
}

int main(int argc, char **argv) {

    if (argc < 4) {
        printf("%s <dataset1> <dataset2> <EXPERIMENTO> <CONSTRUIR (0 o 1, o 2 para construir y ejecutar)>\n", argv[0]);
        printf("EXPERIMENTO: 0 (HDKD1), 1 (HDKD2), 2 (KAMATA)\n");
        return -1;
    }

    std::string dataset1 = argv[1];
    std::string dataset2 = argv[2];

    int experimentoSeleccionado = std::stoi(argv[3]);
    int construir = std::stoi(argv[4]); // 1 para construir, 0 para no construir, 2 para construir y ejecutar

    if (construir == 1) {
        construirEstructura(dataset1, dataset2);
    } else {
        ejecutaExperimentos(dataset1, dataset2, experimentoSeleccionado);
    }

    return 0;
}