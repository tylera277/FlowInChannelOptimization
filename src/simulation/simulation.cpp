

#include "simulation.hpp"

using namespace navier_stokes;


// Basic constructor
// (which looks attrocious at the moment, I know...)
simulation::Simulation::Simulation(int imin, int imax, int jmin, int jmax, 
double dxi, double dyi, double viscosity, double density,
double time_step_size):
imin(imin), imax(imax), jmin(jmin), jmax(jmax),
dxi(dxi), dyi(dyi),
viscosity(viscosity), density(density), time_step_size(time_step_size),

us(imax+2, jmax+2), vs(imax+2, jmax+2),
u(imax+2, jmax+2), v(imax+2, jmax+2),
p(imax+2, jmax+2),
//laplacian( (imax)*(jmax), (imax)*(jmax)),
laplacian((imax-1) * (jmax-1), (imax-1) * (jmax-1)),
rhs( (imax-1)*(jmax-1),1 )
{ }

/////////////////////////////////////////
//////////  BOUNDARY CONDITIONS  /////////
/////////////////////////////////////////

void simulation::Simulation::ApplyBoundaryConditions(std::string bc_name)
{

    if(bc_name == "generic")
    {
        this->ApplyGeneric_BC();
    }
    else if(bc_name == "lid_driven_cavity")
    {
        this->ApplyLidDrivenCavity_BC();
    }       
    else
    {
        throw std::invalid_argument("You did not enter one of the boundary condition options!");
    }
            
}

void simulation::Simulation::ApplyGeneric_BC()
{   

    // Declare the matrices to all be zero, so I dont need to be worried
    // about retrieving undefined values
    p.topLeftCorner(imax+2, jmax+2) = Eigen::MatrixXd::Zero(imax+2, jmax+2);
    u.topLeftCorner(imax+2, jmax+2) = Eigen::MatrixXd::Zero(imax+2, jmax+2);
    v.topLeftCorner(imax+2, jmax+2) = Eigen::MatrixXd::Zero(imax+2, jmax+2);

    // HACK: I will enter these specific boundary conditions in a better way
    // later on
    double u_bot = 0;
    double u_top = 1;
    double v_left = 0;
    double v_right = 0;

    v(int(imax/2), 2) = -50;
    p(int(imax/2), jmax-2) = -1;


    /* for(int it = 1; it < imax; it++)
    {


        u(imin-1, it) = (2 * u_top) - u(imin, it);
        u(imax+1, it) = (2 * u_bot) - u(imax, it);

        v(it, jmin-1) = (2 * v_left) - v(it, jmin);
        v(it, jmax+1) = (2 * v_right) - v(it, jmax); 
        
    } */

    

}

void simulation::Simulation::ApplyLidDrivenCavity_BC()
{


}

/////////////////////////////////////////
//////////  MAIN CALCULATIONS  //////////
/////////////////////////////////////////

void simulation::Simulation::predictor_step()
{
    //std::cout << "BEEPERS \n ";
    for(int i = 1; i <= imax; i++)
    {   
        for(int j = 1; j <= jmax; j++)
        {
            us(i, j) = this->u_star(i, j);
            vs(i, j) = this->v_star(i, j);
        }
    }
}

double simulation::Simulation::u_star(int xPos, int yPos)
{
    int i = xPos;
    int j = yPos;

    double d2u_dx2 = (u(i-1, j) - 2 * u(i, j) + u(i+1, j)) * pow(dxi, 2);
    double d2u_dy2 = (u(i, j-1) - 2 * u(i, j) + u(i, j+1)) * pow(dyi, 2);

    double u_du_dx = u(i,j) * ( (u(i+1, j) - u(i-1, j)) * (0.5 * dxi));
    double v_du_dy = (1.0/4.0) * (v(i-1, j) + v(i, j) + v(i-1, j+1) + v(i, j+1)) * (u(i, j+1) - u(i,j-1))*(0.5*dyi);

    double value =  u(i,j) + time_step_size * ((viscosity)*(d2u_dx2 + d2u_dy2) - (u_du_dx + v_du_dy));

    std::cout << "VAL= " << value << std::endl;
    return value;
}

double simulation::Simulation::v_star(int xPos, int yPos)
{
    int i = xPos;
    int j = yPos;

    double d2v_dx2 = (v(i-1,j) - 2 * v(i,j) + v(i+1,j)) * pow(dxi, 2);
    double d2v_dy2 = (v(i,j-1) - 2*v(i,j) + v(i,j+1)) * pow(dyi, 2);

    double udv_dx = (1.0/4.0) * (u(i, j-1) + u(i, j) + u(i+1, j-1) + u(i+1, j)) * (v(i+1,j)-v(i-1,j))*(0.5*dxi);
    double vdv_dy = v(i,j) * (v(i,j+1) - v(i,j-1)) * (0.5 * dyi);

    double value = v(i,j) + time_step_size * (viscosity*(d2v_dx2 + d2v_dy2) - (udv_dx + vdv_dy));

    return value;
}



void simulation::Simulation::corrector_step()
{
    CalculateRHS();
    
    PrintPressureVector(rhs, "../outputData/30jul_rhs.csv" );
    // PrintLaplacian("../outputData/22jul_laplacianV2.csv");
    std::cout <<"LAP row,cols: " << laplacian.rows() << ", " << laplacian.cols() << std::endl;
    std::cout <<"rhs row,cols: " << rhs.rows() << ", " << rhs.cols() << std::endl;


    std::cout << "Starting solving of linear system of equations..." << std::endl;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(laplacian);
    Eigen::MatrixXd pv = lu.solve(rhs);
    //Eigen::MatrixXd pv = laplacian.fullPivLu().solve(rhs);

    std::cout <<"pv row,cols: " << pv.rows() << ", " << pv.cols() << std::endl;

    std::cout << "Solving completed." << std::endl;

    PrintPressureVector(pv, "../outputData/30jul_pv.csv");
    this->ConvertPressureVectorIntoMatrix(pv);


    calculateUpdatedVelocities();



}

void simulation::Simulation::CreateLaplacian()
{
    // Storage matrix for laplacian values.
    // Initializing all elements to zero so if a place isnt reached during computation,
    // an undefined element isnt being passed on.
    Eigen::MatrixXd temporary_laplacian((imax * jmax), (imax * jmax));
    temporary_laplacian.topLeftCorner((imax * jmax), (imax * jmax)) = Eigen::MatrixXd::Zero((imax * jmax), (imax * jmax));

    // Number of grid points under consideration
    int nx = imax - 1;
    int ny = jmax - 1;

    for(int j = 1; j <= ny; j++)
    {
        for(int i = 1; i <= nx; i++)
        {
            //std::cout << "HMM: " << i << ", " << j << std::endl;
            temporary_laplacian((i+(j-1)*nx), (i+(j-1)*nx)) = 2 * pow(dxi, 2) + 2 * pow(dyi,2);

            for(int ii=(i-1); ii<=(i+1); ii+=2)
            {
                if(ii>0 && ii<=nx)
                {
                    temporary_laplacian((i+(j-1)*nx), (ii+(j-1)*nx)) = -pow(dxi,2);
                }
                 else
                {
                    temporary_laplacian((i+(j-1)*nx), (i+(j-1)*nx)) += -pow(dxi,2);
                }
            }
            for(int jj=(j-1); jj<=(j+1); jj+=2)
            {

                if(jj>0 && jj<=ny)
                {
                    temporary_laplacian((i+(j-1)*nx), (i+(jj-1)*nx)) = -pow(dyi,2);
                }
                 else
                {
                    temporary_laplacian((i+(j-1)*nx), (i+(j-1)*nx)) += -pow(dyi,2);
                }

            }


        }
    }
    temporary_laplacian(1, 1) = 1;
    for(int beep = 2; beep < (imax)*(jmax); beep++)
    {
        temporary_laplacian(1,beep) = 0;
    }


    // Shaving off the excess of the temporary_laplacian matrix
    for(int i = 0; i < ((imax - 1) * (jmax - 1)); i++)
    {
        for(int j = 0; j < ((imax - 1) * (jmax - 1)); j++)
        {
            laplacian(i,j) = temporary_laplacian(i+1, j+1);
        }
    }
}

double simulation::Simulation::GetLaplacianValue(int xPos, int yPos)
{
    return laplacian(xPos, yPos);
}

void simulation::Simulation::CreateLaplacianV2()
{
    
    int nx = imax - 1;
    int ny = jmax - 1;

    double dx = 1 / dxi;
    double dy = 1 / dyi;

    double invDeltaX2 = 1.0 / (dx * dx);
    double invDeltaY2 = 1.0 / (dy * dy);

    for(int i = 0; i < nx; ++i)
    {
        for(int j = 0; j < ny; ++j)
        {
            int currentPoint = i * ny + j;

            laplacian(currentPoint, currentPoint) = -2.0 * (invDeltaX2 + invDeltaY2);
             if(i > 0)
            {
                laplacian(currentPoint, currentPoint - ny) = invDeltaX2;
            }
            if (i < nx - 1)
            {
                laplacian(currentPoint, currentPoint + ny) = invDeltaX2;
            }
            if (j > 0)
            {
                laplacian(currentPoint, currentPoint - 1) = invDeltaY2;
            }
            if(j < ny - 1)
            {
                laplacian(currentPoint, currentPoint + 1) = invDeltaY2;
            } 
        }
    }
    laplacian(0, 0) = 999;

}

void simulation::Simulation::CalculateRHS()
{
    int n = 0;

    for(int j = 1; j < jmax; j++)
    {
        for(int i = 1; i < imax; i++)
        {
            rhs(n,0) = (density/time_step_size) * ((us(i+1, j) - us(i, j)) * dxi +\
            (vs(i, j+1) - vs(i,j))*dyi);
            
            n++;
        }

    }

}

void simulation::Simulation::ConvertPressureVectorIntoMatrix(Eigen::MatrixXd pv)
{
    int n = 0;
    Eigen::MatrixXd temporary_pressure(imax+2, jmax+2);
    for(int j = 1; j < jmax; j++)
    {
        for(int i = 1; i < imax; i++)
        {
            temporary_pressure(i,j) = pv(n,0);

            n++;
        }
    }
    //p = temporary_pressure.transpose();
    p = temporary_pressure;
}

void simulation::Simulation::calculateUpdatedVelocities()
{
    for(int j = jmin; j <= jmax; j++)
    {
        for(int i = imin + 1; i <= imax; i++)
        {
            u(i,j) = us(i,j) - (time_step_size / density) * (p(i,j) - p(i-1,j)) * dxi;
        }
    }
    for(int j = jmin + 1; j <= jmax; j++)
    {
        for(int i = imin; i <= imax; i++)
        {
            v(i,j) = vs(i,j) - (time_step_size / density) * (p(i,j) - p(i,j-1)) * dyi;
        }
    }


}

/////////////////////////////////////////
//////////  UTILITY FUNCTIONS  //////////
/////////////////////////////////////////

void simulation::Simulation::PrintPressure(std::string outputFile)
{
    this->PrintMatrixToCSVFile(p, outputFile);
}

void simulation::Simulation::PrintVelocities(std::string outputFile)
{
    this->PrintMatrixToCSVFile(u, outputFile);
}

void simulation::Simulation::PrintStarVelocities(std::string outputFile)
{
    this->PrintMatrixToCSVFile(us, outputFile);
}

void simulation::Simulation::PrintLaplacian(std::string outputFile)
{
    std::ofstream file;
    // Clear the contents of the file and write afresh for now
    file.open(outputFile);

    for(int i = 0; i < ((imax-1)*(jmax-1)); i++)
    {
        for(int j = 0; j < ((imax-1)*(jmax-1)); j++)
        {
            file << laplacian(i, j) << ",";
        }
        file << "\n";
    }

    file.close();
}

void simulation::Simulation::PrintPressureVector(Eigen::MatrixXd pv, std::string outputFile)
{
    std::ofstream file;
    // Clear the contents of the file and write afresh for now
    file.open(outputFile);

    for(int i = imin; i < ((imax-1)*(jmax-1)); i++)
    {
        file << pv(i, 0) << "\n";
    }

    file.close();
}

void simulation::Simulation::PrintMatrixToCSVFile(Eigen::MatrixXd matrix, std::string outputfile)
{
    std::ofstream file;
    // Clear the contents of the file and write afresh for now
    file.open(outputfile);

    // Used to make it easier to use pandas in python for plotting the results
    for(int filler = imin; filler <= imax; filler++)
    {
        file << 0 << ",";
    }
    file << "\n";

    for(int i = imin; i < imax; i++)
    {
        for(int j = jmin; j < jmax; j++)
        {
            file << matrix(i, j) << ",";
        }
        file << "\n";
    }

    file.close();


}