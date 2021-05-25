#include <robot_model.h>
#include <robot_dynamics.h>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <time.h>

using namespace Eigen;

int main(int argc, char ** argv)
{
    robot_dyn::RobotModel robot;
    robot.set_kinematics_parameters();

    VectorXd param = VectorXd::Zero(robot.Psi_num*robot.dof);
    param << 0.1, 0.3, 0.12, 011, 034, 0.24, 035, 0.223, 0.345, 0.483224, 
             0, 5.83642, 0.690557, 0.1, 3.01686, 16.5103, 1.00592, 0, -0.503815, 2.8116, 
             0, -7.00666, -22.4839, 0, 6.53773, 5.28274, 0.1, -1.3729, -0.538932, -0.397061, 
             -9.09525, 6.03299, -0.318568, -18.3559, -2.18405, -8.17642, -5.54922, -0.341478, 0.886727, 0.1, 
             1.15881, -2.01251, 1.02208, -3.80232, 1.44075, 0.327297, 7.21613, -1.36446, -1.07751, -9.19562, 
             -1.7909, 2.94775, 0.1, 0.365652, 0.427724, -1.47969, 5.30299, -0.20627, 2.60709, 4.31089, 
             -0.864087, 2.22257, -0.239355, 5.57593, -1.51509, 0.1, 0.187228, -0.0657682, -0.337444, -0.111909; 
    robot.set_dynamics_parameters(param);

    MatrixXd Mq = MatrixXd::Zero(robot.dof, robot.dof);
    VectorXd H = VectorXd::Zero(robot.dof);
    VectorXd Gq = VectorXd::Zero(robot.dof);

    clock_t start, finish;
    double totaltime;

    for (unsigned int i=0; i<100; i++)
    {
        VectorXd q = VectorXd::Random(robot.dof);
        VectorXd qDot = VectorXd::Random(robot.dof);

        start = clock();

        Mq = robot_dyn::calcu_InertiaMatrix(&robot, q);
        H = robot_dyn::calcu_CoriolisCentripetal(&robot, q, qDot);
        Gq = robot_dyn::calcu_Gravity(&robot, q);

        finish = clock();
        totaltime = (double)(finish-start)/CLOCKS_PER_SEC;
        std::cout << "totaltime: " << totaltime*1000 << " ms" << std::endl;
    }

    for (unsigned int i=0; i<10; i++)
    {
        VectorXd q = VectorXd::Random(robot.dof);
        VectorXd qDot = VectorXd::Random(robot.dof);

        Mq = robot_dyn::calcu_InertiaMatrix(&robot, q);
        H = robot_dyn::calcu_CoriolisCentripetal(&robot, q, qDot);
        Gq = robot_dyn::calcu_Gravity(&robot, q);

        std::cout << "Mq:" << std::endl << Mq << std::endl;
        std::cout << "H:" << std::endl << H << std::endl;
        std::cout << "Gq:" << std::endl << Gq << std::endl;
    }

    return 0;
}
