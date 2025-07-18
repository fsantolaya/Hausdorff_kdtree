/* 
 * File:   MaxHeap.cpp
 * Author: csanjuan
 *
 * Created on 16 de octubre de 2020, 08:23 PM
 */

#include "MaxHeap.h"

//Obtenido de https://www.techiedelight.com/min-heap-max-heap-implementation-c/

// Data structure for Min Heap

MaxHeap::MaxHeap() {
    
}

int MaxHeap::PARENT(int i) {
    return (i - 1) / 2;
}

// return left child of A[i]	

int MaxHeap::LEFT(int i) {
    return (2 * i + 1);
}

// return right child of A[i]	

int MaxHeap::RIGHT(int i) {
    return (2 * i + 2);
}

// Recursive Heapify-down algorithm
// the node at index i and its two direct children
// violates the heap property

void MaxHeap::heapify_down(int i) {
    // get left and right child of node at index i
    int left = LEFT(i);
    int right = RIGHT(i);

    int biggest = i;

    // compare A[i] with its left and right child
    // and find biggest value
    if (left < size() && heap[left].distance > heap[i].distance)
        biggest = left;

    if (right < size() && heap[right].distance > heap[biggest].distance)
        biggest = right;

    // swap with child having lesser value and 
    // call heapify-down on the child
    if (biggest != i) {
        swap(heap[i], heap[biggest]);
        heapify_down(biggest);
    }
}

// Recursive Heapify-up algorithm

void MaxHeap::heapify_up(int i) {
    // check if node at index i and its parent violates 
    // the heap property
    if (i && heap[PARENT(i)].distance < heap[i].distance) {
        // swap the two if heap property is violated
        swap(heap[i], heap[PARENT(i)]);

        // call Heapify-up on the parent
        heapify_up(PARENT(i));
    }
}

// return size of the heap

unsigned int MaxHeap::size() {
    return heap.size();
}

// function to check if heap is empty or not

bool MaxHeap::empty() {
    return size() == 0;
}

// insert key into the heap

void MaxHeap::push(Node node) {
    // insert the new element to the end of the vector
    heap.push_back(node);

    // get element index and call heapify-up procedure
    int index = size() - 1;
    heapify_up(index);
}

// function to remove element with lowest priority (present at root)
void MaxHeap::pop() {
    try {
        // if heap has no elements, throw an exception
        if (size() == 0)
            throw out_of_range("Vector<X>::at() : "
                "index is out of range(Heap underflow)");

        // replace the root of the heap with the last element
        // of the vector
        heap[0] = heap.back();
        heap.pop_back();

        // call heapify-down on root node
        heapify_down(0);
    }        // catch and print the exception
    catch (const out_of_range& oor) {
        cout << "\n" << oor.what();
    }
}

// function to return element with lowest priority (present at root)

Node MaxHeap::top() {
    try {
        // if heap has no elements, throw an exception
        if (size() == 0)
            throw out_of_range("Vector<X>::at() : "
                "index is out of range(Heap underflow)");

        // else return the top (first) element
        return heap.at(0); // or return A[0];
    }        // catch and print the exception
    catch (const out_of_range& oor) {
        cout << "\n" << oor.what();
    }
}


