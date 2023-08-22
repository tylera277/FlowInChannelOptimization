#define CATCH_CONFIG_MAIN // This tells Catch2 to provide its own main function
#include "../catch.hpp"

#include <vector>
#include <iostream>

#include "/Users/starman/Desktop/denver_2023_projects/CFD/NavierStokes/src/simulation/simulation.hpp"

using namespace navier_stokes;

TEST_CASE("pde solver tests", "[unit]")
{


    int imin = 2;
    int imax = 4;
    int jmin = 2;
    int jmax = 4;

    double dxi = 1 / (2.0 / (imax - 1));
    double dyi = 1 / (2.0 /  (jmax - 1));
    
    double viscosity = 1;
    double density = 1;
    double time_step_size = 0.1;

    simulation::Simulation test_sim(
        imin, imax, jmin, jmax,
        dxi, dyi,
        viscosity, density, time_step_size);


    SECTION("createLaplacian")
    {
        /*
         * Testing the laplacian constructor on a specific known example,
         * making sure it is operating correctly before scaling it up to larger problems 
         */
        test_sim.CreateLaplacian();

        double D = (2 * pow(dxi, 2)) + (2 * pow(dyi, 2));
        double Dx = D - pow(dxi ,2);
        double Dy = D - pow(dyi, 2);
        double Dxy = D - pow(dxi, 2) - pow(dyi, 2);

        std::vector<std::vector<double>> truth {
        {1, 0, 0, 0, 0, 0, 0, 0, 0},
        {-pow(dxi,2), Dy, -pow(dxi, 2), 0, -pow(dyi, 2),0 , 0, 0, 0},
        {0, -pow(dxi,2), Dxy, 0, 0, -pow(dyi,2), 0, 0, 0},
        {-pow(dyi,2), 0, 0, Dx, -pow(dxi,2),0, -pow(dyi,2), 0, 0},
        {0, -pow(dyi,2), 0, -pow(dxi,2), D, -pow(dxi,2), 0, -pow(dyi,2), 0},
        {0, 0, -pow(dyi, 2), 0, -pow(dxi,2), Dx, 0, 0, -pow(dyi,2)},
        {0, 0, 0, -pow(dyi,2), 0, 0, Dxy, -pow(dxi, 2), 0},
        {0, 0, 0, 0, -pow(dyi,2), 0, -pow(dxi,2), Dy, -pow(dxi,2)},
        {0, 0, 0, 0, 0, -pow(dyi,2), 0, -pow(dxi,2), Dxy}};


        int trueCounter = 0;
        std::cout<<"D="<<D << std::endl;
        std::cout<<"D="<<D << std::endl;
        for(int i=0; i < (imax - 1)*(jmax - 1); i++)
        {
            for(int j=0; j < (jmax - 1)*(imax - 1); j++)
            {
                if(abs(test_sim.GetLaplacianValue(i, j) - truth[i][j]) < 0.0005)
                {
                    trueCounter++;
                }
            }
        } 
        REQUIRE(trueCounter == pow((imax-1)*(jmax-1),2));
    }

    
}
