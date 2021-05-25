#include <dynamics_identify.h>
#include <robot_model.h>
#include <robot_math.h>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace robot_iden
{
Fourier::Fourier(robot_dyn::RobotModel robot_)
{
    wf = 0.05*2*M_PI; 
    Tf = 20; 
    h = 0.1;  

    N = 5; 
    num = 11; 

    robot = robot_;

    q.resize(robot.dof);
    qDot.resize(robot.dof);
    qDDot.resize(robot.dof);
}

Fourier::~Fourier(){}

MatrixXd calcu_Ys(robot_dyn::RobotModel *robot, 
    const VectorXd q, const VectorXd qDot, const VectorXd qDDot)
{
    MatrixXd Ys = MatrixXd::Zero(robot->dof,robot->Ps_num);

    for (unsigned int i=0; i< robot->Ps_num; i++)
    {
        VectorXd Ps = VectorXd::Zero(robot->Ps_num);

        if ((i%robot->Psi_num)==1 || (i%robot->Psi_num)==2 || (i%robot->Psi_num)==3)
        {
            Ps(0) = 1;
            Ps(i) = 1;
            robot->set_dynamics_parameters(Ps);
            Ys.col(i) = robot->calcu_inv_dyn(q, qDot, qDDot)-Ys.col(0);
        }
        else
        {
            Ps(i) = 1;
            robot->set_dynamics_parameters(Ps);
            Ys.col(i) = robot->calcu_inv_dyn(q, qDot, qDDot);
        }
    }

    return Ys;
}

void qr_decompose(robot_dyn::RobotModel *robot)
{
    unsigned int count = 1000;

    MatrixXd Ys = MatrixXd::Zero(robot->dof,robot->Ps_num);
    MatrixXd W = MatrixXd::Zero(count*robot->dof,robot->Ps_num);

    VectorXd q = VectorXd::Zero(robot->dof);
    VectorXd qDot = VectorXd::Zero(robot->dof);
    VectorXd qDDot = VectorXd::Zero(robot->dof);

    for (unsigned int i=0; i<count; i++)
    {
        for (unsigned int j=0; j< robot->dof; j++)
        {
            q(j) = random(robot->qMin(j), robot->qMax(j));
            qDot(j) = random(robot->qDotMin(j), robot->qDotMax(j));
            qDDot(j) = random(robot->qDDotMin(j), robot->qDDotMax(j));
        }
        
        Ys = calcu_Ys(robot, q, qDot, qDDot);
        W.middleRows(i*robot->dof, robot->dof) = Ys;
    }

    //normalization
    /*double factor;
    for (unsigned int i=0; i< robot->Ps_num; i++)
    {
        factor = W.col(i).norm();
        
        if (factor > 0.0)
        {
            W.col(i) = W.col(i)/factor;
        }
    }*/

    std::cout << W << std::endl;

    HouseholderQR<MatrixXd> qr;
    qr.compute(W);
    MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
    MatrixXd Q = qr.householderQ();

    for (unsigned int i=0; i< robot->Ps_num; i++)
    {
        if (abs(R(i,i))==0) 
        {
            robot->Ps_flag(i) = 0; 
        }
        else
        {
            robot->Ps_flag(i) = 1; 
            robot->Pb_num += 1;
        }
    }

    robot->R1.resize(robot->Pb_num, robot->Pb_num);
    robot->R2.resize(robot->Pb_num, robot->Ps_num-robot->Pb_num);

    int j=-1;
    int k=-1;
    for (unsigned int i=0; i< robot->Ps_num; i++)
    {
        if (robot->Ps_flag(i) == 1)
        {
            j+=1;
            robot->R1.col(j) = R.block(0,i,robot->Pb_num,1);
        }
        else
        {
            k+=1;
            robot->R2.col(k) = R.block(0,i,robot->Pb_num,1);
        }
    }
}

// ai_1 ... ai_5, bi_1 ... bi_5, qi_0
void generate_fourier_trajectory(const VectorXd x, const double t, Fourier *fourier)
{   
    fourier->q = VectorXd::Zero(fourier->robot.dof);
    fourier->qDot = VectorXd::Zero(fourier->robot.dof);
    fourier->qDDot = VectorXd::Zero(fourier->robot.dof);

    double a;
    double b;

    for (unsigned int i=0; i< fourier->robot.dof; i++)
    {
        for (unsigned int j=1; j<= fourier->N; j++)
        {
            a = x(i*(2*fourier->N+1)+j-1);
            b = x(i*(2*fourier->N+1)+fourier->N+j-1);

            fourier->q(i) = fourier->q(i)+a/(fourier->wf*j)*sin(fourier->wf*j*t)-b/(fourier->wf*j)*cos(fourier->wf*j*t);
            fourier->qDot(i) = fourier->qDot(i)+a*cos(fourier->wf*j*t)+b*sin(fourier->wf*j*t);
            fourier->qDDot(i) = fourier->qDDot(i)-fourier->wf*a*j*sin(fourier->wf*j*t)+fourier->wf*b*j*cos(fourier->wf*j*t);
        }

        fourier->q(i) = fourier->q(i)+x(i*(2*fourier->N+1)+2*fourier->N);
    }
}

double optimal_object_fun(unsigned n, const double *x, double *grad, void *f_data)
{
    static unsigned int evaluations = 0;
    ++evaluations;
    std::cout << "evaluations: " << evaluations << " times." << std::endl;

    if (grad){}

    Fourier *fourier = (Fourier *)f_data;

    unsigned int count = (unsigned int)(fourier->Tf/fourier->h);
    double t = -fourier->h;

    MatrixXd W = MatrixXd::Zero((count+1)*fourier->robot.dof,fourier->robot.Ps_num);

    unsigned int x_num = fourier->robot.dof*fourier->num;
    VectorXd xx = VectorXd::Zero(x_num);

    for (unsigned int i=0; i<x_num; i++)
    {
        xx(i)=x[i];
    }
    
    for (unsigned int i=0; i<=count; i++)
    {
        t += fourier->h; 

        generate_fourier_trajectory(xx, t, fourier);
        W.middleRows(i*fourier->robot.dof, fourier->robot.dof) = 
            calcu_Ys(&fourier->robot, fourier->q, fourier->qDot, fourier->qDDot);
    }

    unsigned int row = (count+1)*fourier->robot.dof;
    double sum1;
    double sum2 = 1;

    for (unsigned int i=0; i< fourier->robot.Ps_num; i++)
    {
        if (fourier->robot.Ps_flag(i) == 1)
        {
            sum1 = 0;
            for (unsigned int j=0; j<row; j++)
            {
                sum1 += W(j,i)*W(j,i);
            }

            sum2 = sum2*sum1;
        }
    }  

    return 1/sum2;
}

void inequality_constraint(unsigned m, double *result, unsigned n, 
    const double *x, double *grad, void *f_data)
{
    if (grad){}

    Fourier *fourier = (Fourier *)f_data;

    VectorXd q_ineq = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDot_ineq = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDDot_ineq = VectorXd::Zero(fourier->robot.dof);

    for (unsigned int i=0; i< fourier->robot.dof; i++)
    {
        q_ineq(i) = min(abs(fourier->robot.qMin(i)), abs(fourier->robot.qMax(i)));
        qDot_ineq(i) = min(abs(fourier->robot.qDotMin(i)), abs(fourier->robot.qDotMax(i)));
        qDDot_ineq(i) = min(abs(fourier->robot.qDDotMin(i)), abs(fourier->robot.qDDotMax(i)));
    }

    VectorXd q = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDot = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDDot = VectorXd::Zero(fourier->robot.dof);

    double a;
    double b;

    for (unsigned int i=0; i< fourier->robot.dof; i++)
    {
        for (unsigned int j=1; j<= fourier->N; j++)
        {
            a = x[i*(2*fourier->N+1)+j-1];
            b = x[i*(2*fourier->N+1)+fourier->N+j-1];

            q(i) = q(i)+sqrt(a*a+b*b)/(fourier->wf*j);
            qDot(i) = qDot(i)+sqrt(a*a+b*b);
            qDDot(i) = qDDot(i)+fourier->wf*sqrt(a*a+b*b)*j;
        }

        q(i) = q(i)+abs(x[i*(2*fourier->N+1)+2*fourier->N]);
    }
        
    for (unsigned int i=0; i< fourier->robot.dof; i++)
    {
        result[i] = q(i)-q_ineq(i);
        result[fourier->robot.dof+i] = qDot(i)-qDot_ineq(i);
        result[2*fourier->robot.dof+i] = qDDot(i)-qDDot_ineq(i);
    }
}

void equality_constraint(unsigned m, double *result, unsigned n, 
    const double *x, double *grad, void *f_data)
{
    if (grad){}

    Fourier *fourier=(Fourier *)f_data;

    VectorXd q_T0 = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDot_T0 = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDDot_T0 = VectorXd::Zero(fourier->robot.dof);
    VectorXd q_Tf = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDot_Tf = VectorXd::Zero(fourier->robot.dof);
    VectorXd qDDot_Tf = VectorXd::Zero(fourier->robot.dof);

    double a;
    double b;

    for (unsigned int i=0; i< fourier->robot.dof; i++)
    {
        for (unsigned int j=1; j<= fourier->N; j++)
        {
            a = x[i*(2*fourier->N+1)+j-1];
            b = x[i*(2*fourier->N+1)+fourier->N+j-1];

            q_T0(i) = q_T0(i)-b/(fourier->wf*j);
            qDot_T0(i) = qDot_T0(i)+a;
            qDDot_T0(i) = qDDot_T0(i)+fourier->wf*b*j;

            q_Tf(i) = q_Tf(i)+a/(fourier->wf*j)*sin(fourier->wf*j*fourier->Tf)-b/(fourier->wf*j)*cos(fourier->wf*j*fourier->Tf);
            qDot_Tf(i) = qDot_Tf(i)+a*cos(fourier->wf*j*fourier->Tf)+b*sin(fourier->wf*j*fourier->Tf);
            qDDot_Tf(i) = qDDot_Tf(i)-fourier->wf*a*j*sin(fourier->wf*j*fourier->Tf)+fourier->wf*b*j*cos(fourier->wf*j*fourier->Tf);
        }

        q_T0(i) = q_T0(i)+x[i*(2*fourier->N+1)+2*fourier->N];
        q_Tf(i) = q_Tf(i)+x[i*(2*fourier->N+1)+2*fourier->N];
    }

    for (unsigned int i=0; i< fourier->robot.dof; i++)
    {
        result[i] = q_T0(i);
        result[fourier->robot.dof+i] = qDot_T0(i);
        result[2*fourier->robot.dof+i] = qDDot_T0(i);
        result[3*fourier->robot.dof+i] = q_Tf(i);
        result[4*fourier->robot.dof+i] = qDot_Tf(i);
        result[5*fourier->robot.dof+i] = qDDot_Tf(i);
    }
}

}
