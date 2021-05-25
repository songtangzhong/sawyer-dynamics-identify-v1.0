#ifndef ROBOT_DYNAMICS_H_
#define ROBOT_DYNAMICS_H_

#include <robot_model.h>
#include <Eigen/Dense>

using namespace Eigen;

namespace robot_dyn
{
MatrixXd calcu_InertiaMatrix(RobotModel *robot, const VectorXd q);

VectorXd calcu_CoriolisCentripetal(RobotModel *robot, const VectorXd q, const VectorXd qDot);

VectorXd calcu_Gravity(RobotModel *robot, const VectorXd q);

}

#endif