/*
Author: James Root
Class: ECE 6122 R
Last Date Modified: 12/04/23

Description:

Main function to run OpenMPI simulations of integration estimation.

*/

#include <iostream>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <cstring>
#include <cmath>

#define MASTER 0

// this function implements integration of a random point between
// x = 0 and 1 on a defined function (determined by value integral)
double integrate(int integral)
{
    double r;
    // long random(void);
    switch (integral)
    {
    case 1:
        // y = x^2
        r = (((double)rand()) / (RAND_MAX));
        return r * r;
        break;
    case 2:
        // y = e^(-x^2)
        r = (((double)rand()) / (RAND_MAX));
        return exp(-1 * r * r);
        break;

    default:
        return 1000;
        break;
    }
}

// main method for running program
int main(int argc, char *argv[])
{
    // check for input errors (must have one -P and -N)
    if (argc != 5 || !(((strcmp(argv[1], "-P") == 0) && (strcmp(argv[3], "-N") == 0)) || ((strcmp(argv[1], "-N") == 0) && (strcmp(argv[3], "-P") == 0))))
    {
        std::cout << "Input error! Must be of format -P [integral number] -N [number of points]" << std::endl;
        return 1;
    }

    // declare flag values
    int integral;
    int numPoints;
    try
    {
        // std::cout << argv[1] << argv[2] << argv[3] << argv[4] << std::endl;
        // check which order flags are given, assign values appropriately
        if ((strcmp(argv[1], "-P") == 0))
        {
            // P given first
            integral = std::stoi(argv[2]);
            numPoints = std::stoi(argv[4]);
        }
        else
        {
            // N given first
            integral = std::stoi(argv[4]);
            numPoints = std::stoi(argv[2]);
        }
        // verify P value being 1 or 2
        if (integral != 1 && integral != 2) {
            throw std::invalid_argument("received negative value");
        }
    }
    catch (const std::exception &e)
    {
        // we have an issue w stoi, must have positive integers
        std::cout << "Input error! input values must be positive integers. P must be 1 or 2" << std::endl;
        return 1;
    }

    // declare mpi variables
    int taskid, numtasks, rc;
    double homeint, intsum, aveint, integ;
    
    intsum = 0;

    MPI_Status status;
    // initialize MPI
    rc = MPI_Init(&argc, &argv);
    // check for error with initialization
    if (rc != MPI_SUCCESS)
    {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // obtain task numbers and id
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    // set our random seeds
    srandom(taskid);

    // calculate how many points we need to calculate per task
    int p = (int) (numPoints / numtasks);

    // iterate through our points
    for (int i = 0; i < p; i++)
    {
        // cakcukate integral
        homeint = integrate(integral);

        // sum together homeint we had i
        rc = MPI_Reduce(&homeint, &intsum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

        if (taskid == MASTER)
        {
            // std::cout << "process " << i << std::endl;

            integ = intsum / numtasks;
            aveint = ((aveint * i) + integ) / (i + 1);
        }
        // std::cout << "Current Integral Value: " << aveint << std::endl;
    }
    // clear value
    homeint = 0;

    // we need to add on any final operations in the master branch
    // if there isnt even division between numPoints and numTasks
    for (int i = 0; i < numPoints % numtasks; i++) {
        if (taskid == MASTER)
        {
            // compute and sum as we did before
            homeint += integrate(integral);
        }
    }

    // add these final few calculations onto our average
    aveint = aveint * (numPoints - numPoints % numtasks) / (1.0 * numPoints) +  ((1.0 * homeint) / numPoints);

    // print our result (only from master processor)
    if (taskid == MASTER)
    {
        std::cout << "Integral value: " << aveint << std::endl;
    }
    // finalize MPI
    MPI_Finalize();
    return 0;
}
