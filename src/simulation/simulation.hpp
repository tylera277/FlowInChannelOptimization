// 16Jul2023
// Houses the main simulation body class, which will
// be interfaced via the main.cpp file

#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <exception>
#include <iostream>
#include <cmath>

namespace navier_stokes{

    namespace simulation{


        class Simulation
        {
        
            // X and Y velocity components. 
            //(Values which are based at the faces of each cell.
            // s is for predictor step values, no s is corrector step solution)
            Eigen::MatrixXd us;
            Eigen::MatrixXd vs;

            Eigen::MatrixXd u;
            Eigen::MatrixXd v;
            // Pressure. 
            // (A value which is stored at the center of each cell)
            Eigen::MatrixXd p;

            // Matrix which stores the once computed Laplacian operator
            Eigen::MatrixXd laplacian;

            // 
            Eigen::MatrixXd rhs;


            int imin;
            int imax;
            int jmin;
            int jmax;
            double dxi;
            double dyi;
            double viscosity;
            double density;
            double time_step_size;

            public:
                // Basic constructor which initializes the matrices that are going to be used
                Simulation(int, int, int, int, double, double, double, double, double);


                /**
                 * @brief Sets/Enforces the boundary conditions of the problem that is trying to be solved
                 * 
                 * @param boundary_conditions_name  The name of the particular B.C.'s one wants to utilize
                 *
                 */
                void ApplyBoundaryConditions(std::string boundary_conditions_name);

                void ApplyGeneric_BC();
                void ApplyLidDrivenCavity_BC();

                /**
                 * @brief Calculates the velocity components for the initial predictor step
                 * (u and v star).
                 * REMINDER: the pressure component is ignored for this step
                 *
                 */
                void predictor_step();

                /**
                 * @brief Houses the math for the calculation of u_star, 
                 * needed only for the predictor step
                 *
                 *  TODO: Still need to implement a test for this
                 */
                double u_star(int xPos, int yPos);

                /**
                 * @brief Houses the math for the calculation of v_star, 
                 * needed only for the predictor step
                 *
                 *  TODO: Still need to implement a test for this
                 */
                double v_star(int xPos, int yPos);



                void corrector_step();


                /**
                 * @brief Assigns value into a matrix called the Laplacian,
                 * this can be computed before the program once and then used 
                 * throughout. Depends mainly on the grid mesh.
                 *
                 *  
                 */
                void CreateLaplacian();

                double GetLaplacianValue(int xPos, int yPos);

                void CreateLaplacianV2();
                
                /**
                 * @brief Determines the right hand side (RHS) of the equation which ultimately
                 * calculates the pressure which is then used in the corrector step
                 * 
                 */
                void CalculateRHS();



                void ConvertPressureVectorIntoMatrix(Eigen::MatrixXd pressure_vector);


                void calculateUpdatedVelocities();

                /////////////////////////////////////////
                //////////  UTILITY FUNCTIONS  //////////
                /////////////////////////////////////////

                void PrintPressure(std::string outputFileName);

                void PrintVelocities(std::string outputFileName);

                void PrintStarVelocities(std::string outputFileName);

                void PrintLaplacian(std::string outputFileName);

                void PrintPressureVector(Eigen::MatrixXd pv, std::string outputFileName);

                void PrintMatrixToCSVFile(Eigen::MatrixXd matrix, std::string outputFileName);

        };



    } // namespace: simulation
} //namespace: navier_stokes