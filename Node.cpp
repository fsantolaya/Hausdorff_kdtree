/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Node.cpp
 * Author: carlos
 * 
 * Created on January 13, 2018, 3:57 AM
 */

#include "Node.h"
Node::Node() {
}

Node::Node(int start, int end, int depth, double distance, int r1[], int r2[]) {
    this->start = start;
    this->end = end;
    this->depth = depth;
    this->distance = distance;
    this->r1[0] = r1[0];
    this->r1[1] = r1[1];
    this->r2[0] = r2[0];
    this->r2[1] = r2[1];
}

