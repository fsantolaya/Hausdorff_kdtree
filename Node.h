/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Node.h
 * Author: carlos
 *
 * Created on January 13, 2018, 3:57 AM
 */

#ifndef NODE_H
#define NODE_H

#include <vector>
#include <list>
#include <iostream>
using namespace std;

class Node {
    
//private:
    
    
public:
    int start;
    int end;
    int depth;
    double distance;
    int r1[2];
    int r2[2];
    
    Node();
    Node(int start, int end, int depth, double distance, int r1[], int r2[]);
};

#endif /* NODE_H */

