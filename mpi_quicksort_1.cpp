#include <iostream>     // std::cout
#include <functional>   // std::less
#include <algorithm>    // std::sort
#include <ctype.h>	// isdigit
#include <math.h>

#include <mpi.h>

#include <sys/time.h>	// gettimeofday, struct timeval

#include <stdio.h>	// printf
#include <stdlib.h>	// srand, rand, RAND_MAX

#include <my_timer.h>

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

template <typename T, class BinaryComparisonOp>
void quickSortSerial (T *first, T *last, const BinaryComparisonOp &comp, int level)
{

    if (first >= last)
        return;

    const size_t cutOff = 10;
    const size_t n = last - first;

    if (n > 1)
        if (n < cutOff)
            // Switch to simple sorting algorithm.
            selectionSort (n, first, comp);
    else {
         // Take middle value.
         const T pivotValue = selectPivot(first,last);
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

int *Partition;
int *nParts;
int numProcs;
int myRank;
int maxParts;

template <typename T>
double sortParallel(T *first, int size, int pr_rank, int max_rank, int rank_index) {

    MPI_Status dtIn;
	int share_pr = pr_rank + pow(2, rank_index); /* Calculate the rank of sharing process*/
	rank_index ++; 				/*Increment the count index*/
	typedef std::less<ValueType> comp;
	if(share_pr > max_rank){    		/*if there is no process to share, run serial quicksort*/

	    //myTimer_t t_start = getTimeStamp();
	    quickSortSerial(first, first + size, comp(), 0);

//		double tSerial = getElapsedTime(t_start, getTimeStamp());
//        printf("myRank=%d, timeSerial=%g\n", pr_rank, tSerial*1000.0);
        //printf("myRank=%d, elementSize=%d\n", pr_rank, size);
        return 0;
	}

    //myTimer_t t_start = getTimeStamp();
    int pivotValue = first[size/2]; 			/* Select the pivot */
    T *partition_pt = myPartition(first, first + size, pivotValue, comp());   	/* partition array */
    int offset = partition_pt - first;

    //double tPartition = getElapsedTime(t_start, getTimeStamp());
    //printf("myRank=%d, time=%g\n", pr_rank, tPartition*1000.0);

    /* Send partition based on size, sort the remaining partitions,
    receive sorted partition */
    double collTime = 0;
	if (offset > size - offset){

        int sizeOfArray = size - offset;

		MPI_Send(&sizeOfArray, 1, MPI::INT, share_pr, 3, MPI_COMM_WORLD);
		MPI_Send((first + offset), sizeOfArray, MPI::INT, share_pr, 1, MPI_COMM_WORLD);

		collTime += sortParallel (first, offset, pr_rank, max_rank, rank_index);

        myTimer_t t_start = getTimeStamp();

		MPI_Recv((first + offset), sizeOfArray, MPI::INT, share_pr, 2, MPI_COMM_WORLD, &dtIn);

		collTime += getElapsedTime(t_start, getTimeStamp());

	} else {

        int sizeOfArray = offset;

        MPI_Send(&sizeOfArray, 1, MPI::INT, share_pr, 3, MPI_COMM_WORLD);
		MPI_Send(first, sizeOfArray, MPI::INT, share_pr, 1, MPI_COMM_WORLD);

		collTime += sortParallel((first + offset), size - offset, pr_rank, max_rank, rank_index);

        myTimer_t t_start = getTimeStamp();

	    MPI_Recv(first, sizeOfArray, MPI::INT, share_pr, 2, MPI_COMM_WORLD, &dtIn);

	    collTime += getElapsedTime(t_start, getTimeStamp());
	}

	return collTime;
}

int main (int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // Define the size of the array to be sorted. The default is 10.
    int n = 10;
    if (argc > 1)
        if (isdigit(*argv[1]))
            n = atoi( argv[1] );

   // Define the number of tests to run. The default is 10.
    int numTests = 10;
    if (argc > 2)
        if (isdigit(*argv[2]) or *argv[2] == '-')
            numTests = atoi( argv[2] );

    int numThreads = 1;

    if (myRank == 0) {

        fprintf(stderr,"Number Objects = %d\n", n);
        fprintf(stderr,"Number Steps   = %d\n", numTests);
        fprintf(stderr,"NumberProcesses= %d\n", numProcs);
    }

    typedef std::less<ValueType> comp;

    if (myRank == 0) {

        double tSort = 0;
        double tColl = 0;

        ValueType *a = new ValueType[n];
        srand(n);
        random_value<ValueType> random;

        for (int i = 0; i < n; ++i) {
            a[i] = random(0,n);
        }

        myTimer_t t_start = getTimeStamp();
        for (int k = 0; k < numTests; ++k) {
            tColl += sortParallel(a, n, 0, numProcs - 1, 0);
        }

        tSort = getElapsedTime(t_start, getTimeStamp());
        tSort /= numTests;
        tColl /= numTests;

        ValueType *lastPtr = isSorted(a, a+n, comp());
        printf("%d, %g, %g, %d, %d, %f\n", n, tSort*1000.0, tColl*1000.0, numTests, lastPtr == (a+n), sizeof(ValueType)*size_t(n)/1024.);

        //Print a[] (if not too large).
        if (1 && n < 51)
            for (int i = 0; i < n; ++i)
                std::cout << i << "," << a[i] << std::endl;

        delete [] a;

    } else {

        ValueType *a;
        MPI_Status msgSt, dtIn;
        int arraySize = 0;
        int indexCount = 0;
        int sourceProc = 0;
        while(pow(2, indexCount) <= myRank) 	/* calculate the indexCount */
            indexCount ++;

        sourceProc = myRank - pow(2, (indexCount - 1));

        for (int k = 0; k < numTests; ++k) {
            MPI_Recv(&arraySize, 1, MPI::INT, sourceProc, 3, MPI_COMM_WORLD, &msgSt);

            a = (int*) malloc(arraySize * sizeof(int));
            MPI_Recv(a, arraySize, MPI::INT, sourceProc, 1, MPI_COMM_WORLD, &dtIn);

            int pivot = a[(arraySize / 2)]; /* Find the pivot */
            sortParallel(a, arraySize, myRank, numProcs, indexCount);                 	/* sort recursively */

            /* send sorted sub array */
            MPI_Send(a, arraySize, MPI::INT, sourceProc, 2, MPI_COMM_WORLD);
            free(a);
        }
    }

    MPI_Finalize();

    return 0;
}






