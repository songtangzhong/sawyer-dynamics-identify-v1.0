#include <robot_dynamics.h>

namespace robot_dyn
{
MatrixXd calcu_InertiaMatrix(RobotModel *robot, const VectorXd q)
{
    VectorXd qDot = VectorXd::Zero(robot->dof);

    double g = robot->g;
    robot->g = 0;

    MatrixXd Mq = MatrixXd::Zero(robot->dof, robot->dof);

    for (unsigned int i=0; i< robot->dof; i++)
    {
        VectorXd qDDot = VectorXd::Zero(robot->dof);

        qDDot(i) = 1;
        Mq.col(i) = robot->calcu_inv_dyn(q,qDot,qDDot);
    }

    robot->g = g;
    
    return Mq;
}

VectorXd calcu_CoriolisCentripetal(RobotModel *robot, const VectorXd q, const VectorXd qDot)
{
    VectorXd qDDot = VectorXd::Zero(robot->dof);

    double g = robot->g;
    robot->g = 0;

    VectorXd H = VectorXd::Zero(robot->dof);

    H = robot->calcu_inv_dyn(q,qDot,qDDot);

    robot->g = g;

    return H;
}

VectorXd calcu_Gravity(RobotModel *robot, const VectorXd q)
{
    VectorXd qDot = VectorXd::Zero(robot->dof);
    VectorXd qDDot = VectorXd::Zero(robot->dof);

    VectorXd Gq = VectorXd::Zero(robot->dof);

    Gq = robot->calcu_inv_dyn(q,qDot,qDDot);

    return Gq;
}

}
