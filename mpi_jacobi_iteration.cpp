#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <my_timer.h>


template <typename T>
T *Allocate (const int n)
{
   T *ptr = new T[n];
   return ptr;
}

template <typename T>
T *Deallocate (const T *ptr)
{
   if (ptr)
      delete [] ptr;

   return NULL;
}

template <typename T>
T ** AllocateMesh (const int M, const int N)
{
   T **data = new T* [M];
   data[0] = new T [M*N];
   for (int i = 1; i < M; ++i)
      data[i] = data[0] + i*N;

   return data;
}

template <typename T>
T ** DeallocateMesh (T ** ptr)
{
   if (ptr)
   {
      if (ptr[0])
         delete [] ptr[0];

      delete [] ptr;
   }

   return NULL;
}

int numProcs = -1;
int myRank = -1;
int **ProcMap = NULL;
int iProc = -1;
int jProc = -1;
int numProcs_i = -1;
int numProcs_j = -1;
int **iStart = NULL;
int **iEnd   = NULL;
int **jStart = NULL;
int **jEnd   = NULL;
MPI_Comm Comm = MPI_COMM_NULL;

// Send/recv the edges of each grid to the north/south/east/west
// neighbors.
void exchange_boundaries (double **x, int Ni, int Nj)
{
    MPI_Request sendRequest[4];
    MPI_Request recvRequest[4];
    int sendRequests = 0;
    int recvRequests = 0;
    MPI_Status sendStatus[4];
    MPI_Status recvStatus[4];

    double *toSendNorthEdge;
    double *toRecvNorthEdge;
    double *toSendSouthEdge;
    double *toRecvSouthEdge;
    double *toSendWestEdge;
    double *toRecvWestEdge;
    double *toSendEastEdge;
    double *toRecvEastEdge;

    if(jProc!=numProcs_j-1){

        toSendNorthEdge = Allocate<double>(Ni-2);
        for(int i = 0; i < Ni-2; i++)
            toSendNorthEdge[i] = x[i+1][Nj-2];
        toRecvNorthEdge = Allocate<double>(Ni-2);

        MPI_Irecv (toRecvNorthEdge, Ni-2, MPI_DOUBLE, ProcMap[iProc][jProc+1],
                    2, MPI_COMM_WORLD, &recvRequest[recvRequests++]);
        //printf("process %d posted recv from %d\n", ProcMap[iProc][jProc], ProcMap[iProc][jProc+1]);

        MPI_Isend (toSendNorthEdge, Ni-2, MPI_DOUBLE, ProcMap[iProc][jProc+1],
                    1, MPI_COMM_WORLD, &sendRequest[sendRequests++]);
        //printf("process %d posted send %lf to %d\n", ProcMap[iProc][jProc], toSendNorthEdge[0], ProcMap[iProc][jProc+1]);
    }

    if(jProc!=0){

        toSendSouthEdge = Allocate<double>(Ni-2);
        for(int i = 0; i < Ni-2; i++)
            toSendSouthEdge[i] = x[i+1][1];
        toRecvSouthEdge = Allocate<double>(Ni-2);

        MPI_Irecv (toRecvSouthEdge, Ni-2, MPI_DOUBLE, ProcMap[iProc][jProc-1],
                    1, MPI_COMM_WORLD, &recvRequest[recvRequests++]);
        //printf("process %d posted recv from %d\n", ProcMap[iProc][jProc], ProcMap[iProc][jProc-1]);

        MPI_Isend (toSendSouthEdge, Ni-2, MPI_DOUBLE, ProcMap[iProc][jProc-1],
                    2, MPI_COMM_WORLD, &sendRequest[sendRequests++]);
        //printf("process %d posted send %lf to %d\n", ProcMap[iProc][jProc], toSendSouthEdge[0], ProcMap[iProc][jProc-1]);
    }


    if(iProc!=numProcs_i-1){

        toSendWestEdge = Allocate<double>(Nj-2);
        for(int j = 0; j < Nj-2; j++)
            toSendWestEdge[j] = x[Ni-2][j+1];
        toRecvWestEdge = Allocate<double>(Nj-2);

        MPI_Irecv (toRecvWestEdge, Nj-2, MPI_DOUBLE, ProcMap[iProc+1][jProc],
                    2, MPI_COMM_WORLD, &recvRequest[recvRequests++]);
        //printf("process %d posted recv from %d\n", ProcMap[iProc][jProc], ProcMap[iProc+1][jProc]);

        MPI_Isend (toSendWestEdge, Nj-2, MPI_DOUBLE, ProcMap[iProc+1][jProc],
                    1, MPI_COMM_WORLD, &sendRequest[sendRequests++]);
        //printf("process %d posted send %lf to %d\n", ProcMap[iProc][jProc], toSendWestEdge[0], ProcMap[iProc+1][jProc]);
    }


    if(iProc!=0){

        toSendEastEdge = Allocate<double>(Nj-2);
        for(int j = 0; j < Nj-2; j++)
            toSendEastEdge[j] = x[1][j+1];
        toRecvEastEdge = Allocate<double>(Nj-2);

        MPI_Irecv (toRecvEastEdge, Nj-2, MPI_DOUBLE, ProcMap[iProc-1][jProc],
                    1, MPI_COMM_WORLD, &recvRequest[recvRequests++]);
        //printf("process %d posted recv from %d\n", ProcMap[iProc][jProc], ProcMap[iProc-1][jProc]);

        MPI_Isend (toSendEastEdge, Nj-2, MPI_DOUBLE, ProcMap[iProc-1][jProc],
                    2, MPI_COMM_WORLD, &sendRequest[sendRequests++]);
        //printf("process %d posted send %lf to %d\n", ProcMap[iProc][jProc], toSendEastEdge[0], ProcMap[iProc-1][jProc]);
    }

    MPI_Waitall ( sendRequests, sendRequest, sendStatus );

    if(jProc!=numProcs_j-1) Deallocate(toSendNorthEdge);
    if(jProc!=0)            Deallocate(toSendSouthEdge);
    if(iProc!=numProcs_i-1) Deallocate(toSendWestEdge);
    if(iProc!=0)            Deallocate(toSendEastEdge);

    MPI_Waitall ( recvRequests, recvRequest, sendStatus );

    if(jProc!=numProcs_j-1){
        for(int i = 0; i < Ni-2; i++)
            x[i+1][Nj-1] = toRecvNorthEdge[i];
        Deallocate(toRecvNorthEdge);
    }

    if(jProc!=0){
        for(int i = 0; i < Ni-2; i++)
            x[i+1][0] = toRecvSouthEdge[i];
        Deallocate(toRecvSouthEdge);
    }

    if(iProc!=numProcs_i-1){
        for(int j = 0; j < Nj-2; j++)
            x[Ni-1][j+1] = toRecvWestEdge[j];
        Deallocate(toRecvWestEdge);
    }

    if(iProc!=0){
        for(int j = 0; j < Nj-2; j++)
            x[0][j+1] = toRecvEastEdge[j];
        Deallocate(toRecvEastEdge);
    }
}

int main (int argc, char* argv[])
{
   int numThreads = 1;
#ifdef _OPENMP
   numThreads = omp_get_max_threads();
#endif

   int mpi_error;
   mpi_error = MPI_Init (&argc, &argv);

   MPI_Comm_dup (MPI_COMM_WORLD, &Comm);

   //int myRank, numProcs;
   mpi_error = MPI_Comm_rank (Comm, &myRank);
   mpi_error = MPI_Comm_size (Comm, &numProcs);

   //printf("I am %d of %d %d\n", myRank, numProcs, numThreads);

   int N = 10; // 10 x 10 global mesh.
   if (argc > 1)
      if (isdigit(*argv[1]))
         N = atoi(argv[1]);

   int maxIterations = 100; // Maximum # of iterations.
   if (argc > 2)
      if (isdigit(*argv[2]))
         maxIterations = atoi(argv[2]);

   double maxResidual = 1e-4; // Maximum residual before terminating.
   if (argc > 3)
      if (isdigit(*argv[3]))
         maxResidual = atof(argv[3]);

   // Partition the mesh across the processes.

   // Number of partitions in the x and y directions.
   numProcs_i = 1;
   numProcs_j = numProcs;

   // Try to find a nice partition if even or square.
   {
      {
         int imin = 1, jmin = 1, min = 2*numProcs;
         for (int i = 1; i < numProcs; i *= 2)
         {
            if (numProcs % i == 0)
            {
               int ni = i;
               int nj = numProcs / ni;
               if (ni+nj < min)
               {
                  imin = ni;
                  jmin = nj;
                  min = ni+nj;
               }
            }
         }
         numProcs_i = imin;
         numProcs_j = jmin;
      }
   }
//   if (myRank == 0)
//      printf("numProcs i,j = %d, %d, %d\n", numProcs, numProcs_i, numProcs_j);

   // Create a mapping of processes onto the 2d mesh.
   //int **ProcMap = AllocateMesh<int>(numProcs_i, numProcs_j);
   ProcMap = AllocateMesh<int>(numProcs_i, numProcs_j);

   // Where am I in the process grid?
   //int iProc, jProc;

   for (int i = 0; i < numProcs_i; ++i)
      for (int j = 0; j < numProcs_j; ++j)
      {
         int rank = j + i * numProcs_j;
         ProcMap[i][j] = rank;

         if (rank == myRank)
         {
            iProc = i;
            jProc = j;
         }
      }

   // Translate the process coordinates into mesh coordinates.
   iStart = AllocateMesh<int>(numProcs_i, numProcs_j);
   iEnd   = AllocateMesh<int>(numProcs_i, numProcs_j);
   jStart = AllocateMesh<int>(numProcs_i, numProcs_j);
   jEnd   = AllocateMesh<int>(numProcs_i, numProcs_j);

   {
      for (int i = 0; i < numProcs_i; ++i)
         for (int j = 0; j < numProcs_j; ++j)
         {
            int Npts_i = (N-2) / numProcs_i;
            if (i < (N-2) % numProcs_i)
               Npts_i++;

            if (i == 0)
               iStart[i][j] = 1;

            iEnd[i][j] = iStart[i][j] + Npts_i - 1;
            if (i != numProcs_i-1)
               iStart[i+1][j] = iEnd[i][j] + 1;

            int Npts_j = (N-2) / numProcs_j;
            if (j < (N-2) % numProcs_j)
               Npts_j++;

            if (j == 0)
               jStart[i][j] = 1;

            jEnd[i][j] = jStart[i][j] + Npts_j - 1;
            if (j != numProcs_j-1)
               jStart[i][j+1] = jEnd[i][j] + 1;
         }
   }

   int Ni = iEnd[iProc][jProc] - iStart[iProc][jProc] + 3;
   int Nj = jEnd[iProc][jProc] - jStart[iProc][jProc] + 3;

   //printf("rank=%d, iProc,jProc=%d,%d, iStart,iEnd,Ni=%d,%d,%d, jStart,jEnd,Nj=%d,%d,%d\n", myRank, iProc, jProc, iStart[iProc][jProc], iEnd[iProc][jProc], Ni, jStart[iProc][jProc], jEnd[iProc][jProc], Nj);

   double **x = AllocateMesh<double>(Ni, Nj);
   double **xtemp = AllocateMesh<double>(Ni, Nj);

   // x[][] is initially zero everywhere expect ...
   // x[][0] is the lower boundary = 1
   // x[0][] is the left boundary = 1

   for (int i = 0; i < Ni; ++i)
      for (int j = 0; j < Nj; ++j)
         x[i][j] = 0;

   if (iProc == 0)
      for (int j = 0; j < Nj; ++j)
         x[0][j] = 1;

   if (jProc == 0)
      for (int i = 0; i < Ni; ++i)
         x[i][0] = 1;

   // Set xtemp = x so the boundaries are consistent.
   for (int i = 0; i < Ni; ++i)
      for (int j = 0; j < Nj; ++j)
         xtemp[i][j] = x[i][j];

   double mpi_p2p_time = 0;
   double mpi_coll_time = 0;

   myTimer_t total_timer = getTimeStamp();

   // Iterate for some number of steps or until we converge.
   int iteration = 0;
   double residual = 1;

   while (residual > maxResidual and iteration < maxIterations)
   {

      myTimer_t mpi_p2p_timer = getTimeStamp();

      if (numProcs > 1)
      {
         exchange_boundaries (x, Ni, Nj);
      }

      mpi_p2p_time += getElapsedTime( mpi_p2p_timer, getTimeStamp() );

      residual = 0;
      #pragma omp parallel for reduction(+:residual)
      for (int i = 1; i < Ni-1; ++i)
         for (int j = 1; j < Nj-1; ++j)
         {
            xtemp[i][j] = (x[i+1][j] + x[i-1][j] +
                           x[i][j+1] + x[i][j-1]) / 4.0;
            double delta = xtemp[i][j] - x[i][j];
            residual += delta*delta;
         }

      #pragma omp parallel for
      for (int i = 1; i < Ni-1; ++i)
         for (int j = 1; j < Nj-1; ++j)
            x[i][j] = xtemp[i][j];

      myTimer_t mpi_coll_timer = getTimeStamp();

      double localResidual = residual;

      // Add a method to reduce the residual values on each process to a single value
      // and make sure that all processes have the same value.

      double globalResidual;
      MPI_Allreduce(&localResidual, &globalResidual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      mpi_coll_time += getElapsedTime( mpi_coll_timer, getTimeStamp() );

      residual = sqrt(globalResidual);
      iteration++;
      //if (myRank == 0 and iteration % 1 == 0)
         //printf("%d %4d: %e\n", myRank, iteration, residual);
   }

   double total_time = getElapsedTime( total_timer, getTimeStamp() );
   //printf("rank=%d, timers = %f, %f, %f\n", myRank, total_time, mpi_p2p_time, mpi_coll_time);

   double local_times[3] = {1000*total_time, 1000*mpi_p2p_time, 1000*mpi_coll_time};
   double global_times[3];

   MPI_Reduce (local_times, global_times, 3, MPI_DOUBLE, MPI_MAX, 0, Comm);

   if (myRank == 0)
      printf("N = %d Iterations = %d  residual = %e total time = %f p2p_time = %f coll_time = %f Procs = %d %d %d %d\n", N, iteration, residual, global_times[0], global_times[1], global_times[2], numProcs, numProcs_i, numProcs_j, numThreads);

//   if (N < 100)
//   {
//      char flname[12];
//      sprintf(flname, "jacobi%d.out", myRank);
//      FILE *f = fopen(flname,"w");
//      for (int i = 0; i < Ni; ++i)
//         for (int j = 0; j < Nj; ++j)
//            fprintf(f,"%e %e %e\n", (iStart[iProc][jProc]+i-1)/double(N-1), (jStart[iProc][jProc]+j-1)/double(N-1), x[i][j]);
//            //fprintf(f,"%d %d %e\n", (iStart[iProc][jProc]+i-1), (jStart[iProc][jProc]+j-1), x[i][j]);
//      fclose(f);
//   }

   DeallocateMesh(xtemp);
   DeallocateMesh(x);

   MPI_Finalize();

   return 0;
}
