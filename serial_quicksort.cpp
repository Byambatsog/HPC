#include <iostream>     // std::cout
#include <functional>   // std::less
#include <algorithm>    // std::sort
#include <ctype.h>	// isdigit

#include <sys/time.h>	// gettimeofday, struct timeval

#include <stdio.h>	// printf
#include <stdlib.h>	// srand, rand, RAND_MAX

#include <my_timer.h>

template <typename T, class BinaryCompOp>
void selectionSort (const int n, T *arr, const BinaryCompOp &comp)
{
   for (int j = 0; j < (n-1); ++j)
   {
      int idx = j;
      for (int i = j+1; i < n; ++i)
         if ( comp( arr[i], arr[idx] ) )
            idx = i;
      if (j != idx)
         std::swap( arr[j], arr[idx] );
   }
}

template <typename T, class BinaryCompOp>
T* myPartition (T *first, T *last, const T &pivot, const BinaryCompOp &comp)
{
   while (first != last)
   {
      while (comp(*first, pivot))
      {
         ++first;
         if (first == last)
            return first;
      }
      do
      {
         --last;
         if (first == last)
            return first;
      }
      while (!comp(*last, pivot));

      std::swap (*first, *last);
      ++first;
   }

   return first;
}

template <typename T, class BinaryComparisonOp>
T * isSorted (T *first, T *last, const BinaryComparisonOp &comp)
{
   T *next = first+1;
   for (; next != last; ++first, ++next)
      if (comp(*next, *first))
         break;
   return next;
}

template <typename T>
inline T medianOfThree (const T x, const T y, const T z)
{
   if (x < y) {
      if (y < z) return y; // [x, y, z] <-
      if (x < z) return z; // [x, z, y] ->
   }
   if (y < z) {
      if (z < x) return z; // [y, z, x]
      if (x < z) return x; // [y, x, z]
   }
   if (z < x) {
      if (z < y) return y; // [z, y, x]
      if (x < y) return x; // [z, x, y]
   }
   return y;
}

template <typename T>
inline T selectPivot (T *first, T *last)
{
   if (last - first > 3)
   {
      const T lo  = *first;
      const T mid = first[(last - first)/2];
      const T hi  = *(last-1);
      return medianOfThree(lo, mid, hi);
   }
   else
      return first[0];
}

template <typename T>
void printList (T *first, T *last)
{
   const int n = last - first;
   printf("[");
   for(int i = 0; i < (n-1); ++i)
      std::cout << first[i] << ",";
   if (n>0) std::cout << first[n-1];
   printf("]\n");
}

template <typename T, class BinaryComparisonOp>
void quickSortSerial (T *first, T *last, const BinaryComparisonOp &comp, int level)
{
   const size_t cutOff = 10;
   const size_t n = last - first;

   if (n > 1)
      if (n < cutOff)
         // Switch to simple sorting algorithm.
         selectionSort (n, first, comp);
      else
      {
         // Take middle value.
         //const T pivotValue = selectPivot(first,last);

         const T pivotValue = first[(last-first) / 2];

         T *middle = myPartition(first, last, pivotValue, comp);

         // Tricky case when the pivot falls at the head.
         if (middle == first || middle == last)
         {
            while (middle < last and !comp(pivotValue,*middle))
               middle++;

            if (middle == last)
            {
               return;
            }

            const T newPivot = *middle;
            middle = myPartition(first, last, newPivot, comp);
         }

         // Do small side first.
         if (middle - first <= last - middle)
         {
            if (middle != first)
               quickSortSerial (first, middle, comp, level);
            //if (middle != last)
               quickSortSerial (middle, last, comp, level);
         }
         else
         {
            if (middle != last)
               quickSortSerial (middle, last, comp, level);
            //if (middle != first)
               quickSortSerial (first, middle, comp, level);
         }
      }
}


template <typename T, class BinaryComparisonOp>
void quickSort (T *first, T *last, const BinaryComparisonOp &comp)
{
   if (first >= last)
      return;

   quickSortSerial (first, last, comp, 0);
}


template <typename T, class BinaryComparisonOp>
T * BinarySearch (const T value, T *lo, T *hi, const BinaryComparisonOp &comp)
{
   while (lo < hi)
   {
      T *mid = lo + (hi - lo) / 2;
      if (comp(value,*mid))
         hi = mid;
      else
         lo = mid+1;
   }
   return hi;
}

template <typename T>
struct random_value
{
   T operator()(void) const;
   T operator()(const T, const T) const;
};

template <> int random_value<int>::
   operator()(const int lo, const int hi) const
      { return lo + rand() % (hi - lo + 1); }
template <> int random_value<int>::
   operator()(void) const
      { return rand(); }

typedef int ValueType;

int main (int argc, char* argv[])
{
   // Define the size of the array to be sorted. The default is 10.
    int n = 10;
    if (argc > 1)
        if (isdigit(*argv[1]))
            n = atoi( argv[1] );

   // Define the number of tests to run. The default is 100.
    int numTests = 10;
    if (argc > 2)
        if (isdigit(*argv[2]) or *argv[2] == '-')
            numTests = atoi( argv[2] );


    fprintf(stderr,"size = %d, numTests = %d\n", n,  numTests);


    ValueType *a = new ValueType[n];

    typedef std::less<ValueType> comp;

    // Initialize a[] with random numbers (0,1].
    {
      // 1. Seed the pseudo-random generator.
      srand(n);

      random_value<ValueType> random;

      // 2. Loop over a[] and set each value.
      for (int i = 0; i < n; ++i)
      {
         a[i] = random(0,n);
      }

      //Print a[] (if not too large).
      if (1 && n < 50)
        for (int i = 0; i < n; ++i)
           std::cout << i << "," << a[i] << std::endl;
    }

    bool iterate = true;
    if (numTests < 0)
    {
      iterate = false;
      numTests = -numTests;
    }

    // Run the test several times.
    double tSort = 0;
//    for (;;)
//    {
      myTimer_t t_start = getTimeStamp();
      for (int k = 0; k < numTests; ++k)
      {
         // 1. Sort a[].
         quickSort (a, a + n, comp());

         //printf("step=%d, time=%g\n", k, getElapsedTime(t_start, getTimeStamp())*1000.0);

      }
      tSort = getElapsedTime(t_start, getTimeStamp());
//      if (tSort < 0.1 and iterate)
//         numTests *= 2;
//      else
//         break;
    //}
    tSort /= numTests;
    ValueType *lastPtr = isSorted(a, a+n, comp());
    printf("%d, %g, %d, %d, %f\n", n, tSort*1000.0, numTests, lastPtr == (a+n), sizeof(ValueType)*size_t(n)/1024.);

    if (not(iterate))
    {
      printf("isSorted = %d %d\n", lastPtr == (a+n), lastPtr - a);
      ValueType *p = BinarySearch (a[n/2], a, a+n, comp());
      std::cout << "Search: " << a[n/2] << ", " << p-a << ", " << *p << std::endl;
    }

    //Print a[] (if not too large).
    if (1 && n < 50)
      for (int i = 0; i < n; ++i)
         std::cout << i << "," << a[i] << std::endl;

    delete [] a;

    return 0;
}
