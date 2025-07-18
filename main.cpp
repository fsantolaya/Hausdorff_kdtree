/* 
 * File:   main.cpp
 * Author: Carlos San Juan
 *
 * Created on 01 de octubre de 2020
 *
 * Recibio modificación de : Fernando Santolaya.
 * Obs: se agrego para medir con TimeMesure.h, y utlidades para mejorar la experimentacion..
 */

#include <vector>
#include <list>
#include <iostream>

#include "Build_cBiK.h"
#include "NodeSKQ.h"
#include "TimeMesure.h"

//Para leer archivos
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char **argv) {
    //SET HAUSDORFF 1 - GAUSS
    //Parametros generales
    Cronometer *crono = cCronometer();
   // string datasetNumber = "1G";
    if (argc != 3) {
        std::cerr << "Error: Se requieren exactamente dos argumentos." << std::endl;
        std::cerr << "Uso: " << argv[0] << " <dataSetA> <dataSetB>" << std::endl;
        return 1;  // Indicar error
    }
    std::string dataSetA = argv[1];
    std::string dataSetB = argv[2];

  //  string dataSetA = "set1Conjunto" + datasetNumber + ".csv";
  //  string dataSetB = "set2Conjunto" + datasetNumber + ".csv";
    
    /**************************************************************************/
    /***************************** BUILD cBiK *********************************/
    /**************************************************************************/
    //Build_cBiK KDTreeA = Build_cBiK(dataSetA);
    //Build_cBiK KDTreeB = Build_cBiK(dataSetB);
    
    /**************************************************************************/
    /********************** FIND DISTANCE HAUSDORFF ***************************/
    /**************************************************************************/
    Build_cBiK cBiK = Build_cBiK();
    cBiK.load_cBiK(dataSetA, dataSetB);
    //cBiK.printMapa();
    
    //Ejecuta las pruebas
    //cBiK.bruteForce();
    double time = 0.0;
    start_clock(crono);
    auto directHauss = cBiK.hdKD1();
    time = (double) (stop_clock(crono) * 1000000.0);


    double time2 = 0.0;
    start_clock(crono);
    auto directHauss2 = cBiK.hdKD2();
    time2 = (double) (stop_clock(crono) * 1000000.0);



    //IMPRIMENDO RESULTADOS
    cout << "************************ hdKD1 ********************************" << endl;
    cout << "Distancia H: " << directHauss << endl;
    cout << "Tiempo: " << time << " Microsegundos" << endl;

    cout << "************************ hdKD2 ********************************" << endl;
    cout << "Distancia H: " << directHauss2 << endl;
    cout << "Tiempo: " << time2 << " Microsegundos" << endl;
    
    return 0;
}