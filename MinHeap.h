/* 
 * File:   MinHeap.h
 * Author: csanjuan
 *
 * Created on 30 de octubre de 2020, 10:48 AM
 */

#ifndef MINHEAP_H
#define	MINHEAP_H

#include "Node.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
using namespace std;

class MinHeap {
    
private:
    // vector to store heap elements
    vector<Node> A;
    
    int PARENT(int i);	
    int LEFT(int i);	
    int RIGHT(int i);
    void heapify_down(int i);
    void heapify_up(int i);

public:
    MinHeap();
    unsigned int size();
    bool empty();
    void push(Node node);
    void pop();
    Node top();
};

#endif	/* MINHEAP_H */

