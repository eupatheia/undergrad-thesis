#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quicksort.h"
#include "transforms.h"

// helper function swaps indexes i and j in arr
void swap(struct specimen * arr, int i, int j) {
  struct specimen temp = arr[i];
  arr[i] = arr[j];
  arr[j] = temp;
}

extern void quicksort(struct specimen * arr, int left, int right) {
  int random = rand();
  if (left >= right) {
    return;
  }
  // pivot is a random element in A[left]...A[right]
  int pivot = (rand() % (right + 1 - left)) + left;
  swap(arr, left, pivot);
  int mid = left;
  for (int i = left + 1; i <= right; i++) {
    if (arr[i].quality < arr[left].quality) {
      swap(arr, ++mid, i);
    }
  }
  swap(arr, left, mid);   // arr now partitioned at arr[mid]
  quicksort(arr, left, mid - 1);
  quicksort(arr, mid + 1, right);
}

void insertionSort(struct specimen * arr, int left, int right) {
  int size = right - left + 1;
  for (int i = 1; i < size; i++) {
    struct specimen key = arr[i];  // to be inserted
    // all elements at lower indexes are assumed to be already sorted
    int j = i - 1;
    // look for first element where arr[j] < key
    while (j > -1 && arr[j].quality > key.quality) {
      // shift element up the list
      arr[j + 1] = arr[j];
      j--;
    }
    // insert key in correct spot
    arr[j + 1] = key;
  }
}

void modifiedQuicksort(struct specimen * arr, int left, int right) {
  int random = rand();
  // the following if statement is modified
  // so that the end result is ALMOST sorted
  if (right - left < 25) {
    return;
  }
  // pivot is a random element in A[left]...A[right]
  int pivot = (rand() % (right + 1 - left)) + left;
  swap(arr, left, pivot);
  int mid = left;
  for (int i = left + 1; i <= right; i++) {
    if (arr[i].quality < arr[left].quality) {
      swap(arr, ++mid, i);
    }
  }
  swap(arr, left, mid);   // arr now partitioned at arr[mid]
  modifiedQuicksort(arr, left, mid - 1);
  modifiedQuicksort(arr, mid + 1, right);
}

extern void hybridSort(struct specimen * arr, int left, int right) {
  modifiedQuicksort(arr, left, right);
  insertionSort(arr, left, right);
}