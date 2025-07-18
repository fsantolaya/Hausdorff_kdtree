/* 
 * File:   MaxHeap.h
 * Author: csanjuan
 *
 * Created on 30 de octubre de 2020, 10:48 AM
 */

#ifndef MAXHEAP_H
#define	MAXHEAP_H

#include "Node.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>

class MaxHeap {
    
private:
    // vector to store heap elements
    vector<Node> heap;
    
    int PARENT(int i);	
    int LEFT(int i);	
    int RIGHT(int i);
    void heapify_down(int i);
    void heapify_up(int i);

public:
    MaxHeap();
    unsigned int size();
    bool empty();
    void push(Node node);
    void pop();
    Node top();
};

#endif	/* MAXHEAP_H */

