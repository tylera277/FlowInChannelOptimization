
#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "simulation/simulation.hpp"
#include "simulation_parameters.hpp"


using namespace navier_stokes;

int main()
{

    std::string outputFileName = "../outputData/";

    double current_time_step = 0;

    // I need to think of a better way to pass all of these parameters.
    // My current motivation for doing this is to be able to more easily test each unit of a class
    // by being able to create an instance of the class.
    simulation::Simulation sim(
        sim_params::imin, sim_params::imax, sim_params::jmin, sim_params::jmax,
        sim_params::dxi, sim_params::dyi,
        sim_params::viscosity, sim_params::density, sim_params::time_step_size
    );

    std::cout << "Creating laplacian..." << std::endl;
    sim.CreateLaplacian();
    std::cout << "Laplacian created." << std::endl;

    //sim.PrintLaplacian("../outputData/22jul_laplacianv2.csv");
    std::cout << "STEPPER: " << sim_params::time_step_size << std::endl;

    //while(current_time_step < sim_params::total_simulation_time)
    while(current_time_step < (2*sim_params::time_step_size))
    {   
        
        //sim.PrintVelocities(outputFileName + "30jul_u_velocities_one.csv");


        // Applying boundary conditions to the problem domain.
        // Values are entered directly into function simulation.cpp definition,
        // I will change this in the future
        sim.ApplyBoundaryConditions("generic");



        // Compute an intermediate velocity by solving the momentum eqn.
        // but omitting the effect of pressure
        sim.predictor_step();

        sim.PrintStarVelocities(outputFileName + "30jul_us_velocities_one.csv");

        //sim.PrintVelocities(outputFileName + "22jul_star_velocities.csv");

        // Updating the value calculated for velocities where youre now taking into 
        // account the influence of pressure
        sim.corrector_step();



        current_time_step += sim_params::time_step_size;

    }
    sim.PrintStarVelocities(outputFileName + "30jul_us_velocities.csv");
    sim.PrintVelocities(outputFileName + "30jul_u_velocities.csv");
    sim.PrintPressure(outputFileName + "30jul_pressure.csv");


    return 0;
}