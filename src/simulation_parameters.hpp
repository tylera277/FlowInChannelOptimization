// 12Jul2023
// Program housing the parameters for whatever simulation
// I am running at the moment

// #ifndef SIMULATION_PARAMETERS_HPP
// #define SIMULATION_PARAMETERS_HPP

namespace navier_stokes{

    namespace sim_params{

        // The viscosity of the system that I'm simulating
        const double viscosity = 0.01;

        // The density of the system (kg/m3)
        const double density = 1.0;



        // Number of grid points in the x-direction
        const int x_grid_points = 32;

        //         "            "       y-direction
        const int y_grid_points = 32;


        // Length of simulation in the x-direction
        const double length_x = 2;

         // Length of simulation in the y-direction
        const double length_y = 2;


        // Things im putting in this namespace to try out accessing them from here
        const int imin = 2;
        const int jmin = 2;
        const int imax = imin + x_grid_points - 1;
        const int jmax = jmin + y_grid_points - 1;

        // ( this calculation of dx&dy may be incorrect )
        const double dx = length_x / x_grid_points;
        const double dy = length_y / y_grid_points;
        const double dxi = 1 / dx;
        const double dyi = 1 / dy;

        // The total amount of time I want the simulation to run for
        const int total_simulation_time = 3;

        
        const double time_step_size = dx / 1400;


    }// namespace: simulation_parameters
} // namespace: navier_stokes


// #endif