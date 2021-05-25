#ifndef DYNAMICS_IDENTIFY_H_
#define DYNAMICS_IDENTIFY_H_

#include <robot_model.h>
#include <Eigen/Dense>

using namespace Eigen;

namespace robot_iden
{
class Fourier
{
public:
    Fourier(robot_dyn::RobotModel robot_);
    ~Fourier();

    // Which can be redefined by user.
    double wf; 
    double Tf; 
    double h;   //sampling period of system
    // End Which

    int N; 
    int num;    //parameters of each joint

    robot_dyn::RobotModel robot; 

    VectorXd q;
    VectorXd qDot;
    VectorXd qDDot;

private:

};

MatrixXd calcu_Ys(robot_dyn::RobotModel *robot, 
    const VectorXd q, const VectorXd qDot, const VectorXd qDDot);

void qr_decompose(robot_dyn::RobotModel *robot);

// ai_1 ... ai_5, bi_1 ... bi_5, qi_0
void generate_fourier_trajectory(const VectorXd x, const double t, Fourier *fourier);

double optimal_object_fun(unsigned n, const double *x, double *grad, void *f_data);

void inequality_constraint(unsigned m, double *result, unsigned n, 
    const double *x, double *grad, void *f_data);

void equality_constraint(unsigned m, double *result, unsigned n, 
    const double *x, double *grad, void *f_data);

}

#endif