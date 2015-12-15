//it is not completed

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

unsigned long globalRandState;
#define sGlobalRand(s) (globalRandState = (s))
#define globalRand() (globalRandState = ((1664525L*globalRandState + 1013904223L) & 0xffffffffL))
#define globalRandInt(n) ((globalRand() >> 10) % (n))

#define SMALLER 0
#define GREATER 1



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

static void reallocMPIBuffer(int inSize) {
   int size;
   int *ptr;
   MPI_Buffer_detach(&ptr, &size);

   int newSize = inSize * sizeof(int) + 4 * MPI_BSEND_OVERHEAD;
   if (newSize > size) {
      ptr = (int *)realloc(ptr, newSize);
   }
   MPI_Buffer_attach(ptr, newSize);
}

static void freeMPIBuffer() {
   int size;
   void *ptr;
   MPI_Buffer_detach(&ptr, &size);
   free(ptr);
}

template <typename T>
static int getGlobalPivot(T *first, T *last, MPI_Comm comm, int numProcs) {
    int pivotValue = selectPivot(first, last);
    //int pivotValueProc = globalRandInt(numProcs);   /* from random PE */
    int pivotValueProc = 0;

    /* overwrite pivotValue by that one from random selected processor */
    //fprintf(stderr,"broadcasting source=%d numProcs=%d\n",pivotValueProc, numProcs);
    MPI_Bcast(&pivotValue, 1, MPI_INT, pivotValueProc, comm);
    //fprintf(stderr,"broadcasting is ok\n");
    return pivotValue;
}

/* determine prefix-sum and overall sum over value */
static void countGlobal(int *lengthOfArray, int *prefixSum, int *prefixSumTotal, MPI_Comm comm, int numProcsInPartition) {
   MPI_Scan(lengthOfArray, prefixSum, 2, MPI_INT, MPI_SUM, comm);
   prefixSumTotal[SMALLER] = prefixSum[SMALLER];
   prefixSumTotal[GREATER] = prefixSum[GREATER];
   MPI_Bcast(prefixSumTotal, 2, MPI_INT, numProcsInPartition - 1, comm);
   prefixSum[SMALLER] -= lengthOfArray[SMALLER];
   prefixSum[GREATER] -= lengthOfArray[GREATER];
}

static int calculateNumProcsForSmall(int numProcsInPartition, int *prefixSumTotal) {

    double ratioForSmall = (double)(prefixSumTotal[SMALLER]*numProcsInPartition) / (prefixSumTotal[SMALLER]+prefixSumTotal[GREATER]);

    if (ratioForSmall < 1) {
        return ceil(ratioForSmall);
    } else if (numProcsInPartition - ratioForSmall < 1) {
      return floor(ratioForSmall);
    } else if (ratioForSmall - floor(ratioForSmall) > 0.5) {
      return ceil(ratioForSmall);
    } else {
      return floor(ratioForSmall);
    }
}

template <typename T>
static void distributeItems(T *first, int startProc, int fillLevel, int remaining, int bufferSize, MPI_Comm comm) {
   int sendCount;
   int totalSendCount = 0;

   while (remaining > 0) {
      sendCount = bufferSize - fillLevel;
      if (sendCount > remaining) {
         sendCount = remaining;
      }
      MPI_Bsend(first+totalSendCount, sendCount, MPI_INT, startProc, 88, comm);
      startProc++;
      fillLevel = 0;
      totalSendCount += sendCount;
      remaining -= sendCount;
   }
}

template <typename T, class BinaryComparisonOp>
T * quickSort(T *first, T *last, MPI_Comm comm, int *outSize, const BinaryComparisonOp &comp) {

    int numProcsInPartition, rankInPartition;
    int designatedPartition;
    int designatedPartitionSize;
    int desigantedBufferSize;
    int isLastPEInDesigantedPartition;
    MPI_Status status;

    MPI_Comm_size(comm, &numProcsInPartition);
    MPI_Comm_rank(comm, &rankInPartition);
    //fprintf(stderr,"rankInPartition = %d, numProcsInPartition = %d\n", rankInPartition, numProcsInPartition);

    while (numProcsInPartition > 1) {

        //printf("while starting = %d\n", rankInPartition);


        reallocMPIBuffer(last - first);


        int pivotValue = getGlobalPivot(first, last, comm, numProcsInPartition);

        //fprintf(stderr,"myRank = %d, pivotValue = %d\n", rankInPartition, pivotValue);
        T *middle = myPartition(first, last, pivotValue, comp);

        //fprintf(stderr,"myRank = %d, pivotValue = %d, lengthOfSmaller=%d, lengthOfGreater=%d,\n", rankInPartition, pivotValue, middle - first, last - middle);

        int lengthOfArray[2] = {middle - first, last - middle};
        int prefixSum[2] = {0, 0};
        int prefixSumTotal[2] = {0, 0};
        countGlobal(lengthOfArray, prefixSum, prefixSumTotal, comm, numProcsInPartition);

        //Print a[] (if not too large).
//        if (last-first < 50) {
//            char flname[12];
//            sprintf(flname, "mpiquicksortpartitioned%d.out", rankInPartition);
//            FILE *f = fopen(flname,"w");
//            for (int i = 0; i < (last-first); ++i)
//                fprintf(f,"%d, %d \n", i, first[i]);
//            fclose(f);
//        }

        int numProcsForSmall = calculateNumProcsForSmall(numProcsInPartition, prefixSumTotal);
        int numProcsForGreat = numProcsInPartition - numProcsForSmall;

        /* calculate size of buffers for received small / received large items */
        int averageSizeForSmall = ceil(prefixSumTotal[SMALLER]/(double)numProcsForSmall);
        int averageSizeForGreat = ceil(prefixSumTotal[GREATER]/(double)numProcsForGreat);

        //printf("numProcsForSmall: %i, numProcsForGreat %i\n",numProcsForSmall,numProcsForGreat );
        //printf("averageSizeForSmall: %i, averageSizeForGreat: %i\n",averageSizeForSmall,averageSizeForGreat);



        // Send elements in the left partition.
        int startPE = prefixSum[SMALLER] / averageSizeForSmall;
        int fillLevel = prefixSum[SMALLER] % averageSizeForSmall;
        int numOfElementsToSend = lengthOfArray[SMALLER];
        distributeItems(first, startPE, fillLevel, numOfElementsToSend, averageSizeForSmall, comm);

        // Send elements in the right partition.
        startPE = numProcsForSmall + prefixSum[GREATER] / averageSizeForGreat;
        fillLevel = prefixSum[GREATER] % averageSizeForGreat;
        numOfElementsToSend = lengthOfArray[GREATER];
        distributeItems(first+lengthOfArray[SMALLER], startPE, fillLevel, numOfElementsToSend, averageSizeForGreat, comm);


        // Prepare for receive
        if (rankInPartition < numProcsForSmall) {
            designatedPartitionSize = numProcsForSmall;
            desigantedBufferSize = averageSizeForSmall;
            isLastPEInDesigantedPartition = rankInPartition == numProcsForSmall-1;
	        designatedPartition = SMALLER;
        } else {
            designatedPartitionSize = numProcsForGreat;
            desigantedBufferSize = averageSizeForGreat;
            isLastPEInDesigantedPartition = rankInPartition == numProcsForGreat-1;
	        designatedPartition = GREATER;
        }


        if (desigantedBufferSize > (last - first)) { // Only enlarge, don't shrink.
            first = (T *) realloc(first, desigantedBufferSize*sizeof(int));
        }

        int expectedItems = isLastPEInDesigantedPartition ? prefixSumTotal[designatedPartition]-
            desigantedBufferSize*(designatedPartitionSize-1) : desigantedBufferSize;
        int numOfReceivedItems = 0;
        int totalNumOfReceivedItems = 0;

        // Receive all incoming elements
        while (totalNumOfReceivedItems < expectedItems ) {
            MPI_Recv(first+totalNumOfReceivedItems, desigantedBufferSize, MPI_INT, MPI_ANY_SOURCE, 88, comm, &status);
            MPI_Get_count(&status, MPI_INT, &numOfReceivedItems);
            totalNumOfReceivedItems = totalNumOfReceivedItems + numOfReceivedItems;
            //printf("PE %i received %i items from PE %i\n", rankInPartition, numOfReceivedItems, status.MPI_SOURCE);
        }


        last = first + totalNumOfReceivedItems;
        MPI_Comm_split(comm, rankInPartition < numProcsForSmall , 0, &comm);
        MPI_Comm_size(comm, &numProcsInPartition);
        MPI_Comm_rank(comm, &rankInPartition);

        //printf("group numProcsInPartition=%d, rankInPartition=%d partitionSize=%d\n", numProcsInPartition, rankInPartition, last-first);

    }

//    if (last-first < 50) {
//        char flname[12];
//        sprintf(flname, "mpiquicksortreceived%d.out", rankInPartition);
//        FILE *f = fopen(flname,"w");
//        for (int i = 0; i < (last-first); ++i)
//            fprintf(f,"%d, %d \n", i, first[i]);
//        fclose(f);
//    }

    quickSortSerial(first, last, comp, 0);
    *outSize = last-first;
    return first;
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

   // Define the number of tests to run. The default is 100.
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


    ValueType *Partition = new ValueType[numProcs+1];
    ValueType *nParts = new ValueType[numProcs];

    maxParts = 0;
    Partition[0] = 0;
    for (int i = 0; i < numProcs; ++i) {
        int numParts_i = n / numProcs;
        if (i < n % numProcs)
            numParts_i++;

        Partition[i+1] = Partition[i] + numParts_i;
        if (myRank == 0)
            printf("Partition[%d] = %d\n", i+1, Partition[i+1]);
    }

    for (int i = 0; i < numProcs; ++i) {
        nParts[i] = Partition[i+1] - Partition[i];
        if (myRank == 0)
            printf("nParts[%d] = %d\n", i, nParts[i]);

        maxParts = std::max(maxParts, nParts[i]);
    }

    const int nGlobal = n;
    n = nParts[myRank];


    double timings[numTests];
    int outSize;

    ValueType *a = new ValueType[n];


    if (myRank == 0) {
         /* 1. Seed the pseudo-random generator. */
         //srand(n);
         srand(nGlobal);

         // Do proc 0's initialize first.
        random_value<ValueType> random;

         // 2. Loop over a[] and set each value.
        for (int i = 0; i < n; ++i) {
            a[i] = random(0,nGlobal);
        }

        for (int p = 1; p < numProcs; ++p) {

            int np = nParts[p];

            ValueType *tmp_a = new ValueType[n];

            for (int i = 0; i < np; ++i) {
               /* 2. Loop over a[] and set each value */
               tmp_a[i] = random(0,nGlobal);
            }

            MPI_Send (&tmp_a[0], np, MPI_INT, p, 1, MPI_COMM_WORLD);
        }
    } else {

        ValueType *tmp_a = new ValueType[n];

        MPI_Status status;
        MPI_Recv (&tmp_a[0], n, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        for (int i = 0; i < n; ++i) {
            a[i] = tmp_a[i];
        }
    }

    for (int k = 0; k < numTests; ++k) {

        //MPI_Barrier(MPI_COMM_WORLD);
        ValueType *b = new ValueType[n];

        std::copy(a, a + n, b);

        //int *a = malloc(n*sizeof(int));

        typedef std::less<ValueType> comp;

        // Initialize a[] with random numbers (0,1].
        /* In parallel, this is harder since each process will have a different
         * random sequence. Instead, just let one process do this and send the
         * initial values to everyone else. */

        //printf("Done initializing %d size=%d\n", myRank, n);

        myTimer_t t_start = getTimeStamp();
        b = quickSort(b, b + n, MPI_COMM_WORLD, &outSize, comp());
        double tSort = getElapsedTime(t_start, getTimeStamp());

        MPI_Reduce (&tSort, &timings[k], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (outSize < 50) {
            char flname[32];
            sprintf(flname, "mpiquicksortsorted%d.out", myRank);
            FILE *f = fopen(flname,"w");
            for (int i = 0; i < outSize; ++i)
                //printf("%d, %d \n", i, a[i]);
                fprintf(f,"%d, %d \n", i, b[i]);
            fclose(f);
        }

        ValueType *lastPtr = isSorted(b, b+outSize, comp());
        if(lastPtr != (b+outSize))
            printf("myRank %d isSorted=%d\n", myRank, lastPtr == (b+outSize));

        delete [] b;

        freeMPIBuffer();
        //printf("myRank %d sortTime=%g\n", myRank, tSort*1000.0);

    }

    delete [] a;

    printf("myRank=%d partitionSize=%d\n", myRank, outSize);


    if (myRank == 0) {
        double globalTime = 0;
        for(int i = 0 ; i < numTests ; i++)
            globalTime +=timings[i];

        globalTime /= numTests;
        printf("%d, %g, %d, %f\n", nGlobal, globalTime*1000.0, numTests, sizeof(ValueType)*size_t(nGlobal)/1024.);
    }





    //freeMPIBuffer();
    MPI_Finalize();

    return 0;
}