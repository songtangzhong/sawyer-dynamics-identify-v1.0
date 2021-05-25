#ifndef ROBOT_MODEL_H_
#define ROBOT_MODEL_H_

#include <Eigen/Dense>

namespace robot_dyn
{
using namespace Eigen;

class RobotModel
{
public:
    RobotModel();
    ~RobotModel();

    void set_kinematics_parameters();

    // m Pcx Pcy Pcz Ixx Ixy Ixz Iyy Iyz Izz
    void set_dynamics_parameters(const VectorXd param);

    VectorXd calcu_inv_dyn(const VectorXd q, const VectorXd qDot, const VectorXd qDDot);

    unsigned int dof = 7;

    double g = -9.81;

    unsigned int Psi_num = 10;

    unsigned int Ps_num;

    VectorXd qMin; VectorXd qMax;
    VectorXd qDotMin; VectorXd qDotMax;
    VectorXd qDDotMin; VectorXd qDDotMax;

    double qr_threshold = 1e-100;

    VectorXi Ps_flag;

    unsigned int Pb_num = 0;

    MatrixXd R1;
    MatrixXd R2;

private:
    double theta1; double theta2; double theta3; double theta4; double theta5; double theta6; double theta7;
    double d1; double d2; double d3; double d4; double d5; double d6; double d7;
    double alpha1; double alpha2; double alpha3; double alpha4; double alpha5; double alpha6; double alpha7;
    double a1; double a2; double a3; double a4; double a5; double a6; double a7;
    double offset1; double offset2; double offset3; double offset4; double offset5; double offset6; double offset7;

    Vector3d w0; Vector3d w1; Vector3d w2; Vector3d w3; Vector3d w4; Vector3d w5; Vector3d w6; Vector3d w7;

    Vector3d wDot0; Vector3d wDot1; Vector3d wDot2; Vector3d wDot3; Vector3d wDot4; Vector3d wDot5; Vector3d wDot6; Vector3d wDot7;

    Vector3d vDot0; Vector3d vDot1; Vector3d vDot2; Vector3d vDot3; Vector3d vDot4; Vector3d vDot5; Vector3d vDot6; Vector3d vDot7;

    Vector3d vcDot1; Vector3d vcDot2; Vector3d vcDot3; Vector3d vcDot4; Vector3d vcDot5; Vector3d vcDot6; Vector3d vcDot7;

    Vector3d F1; Vector3d F2; Vector3d F3; Vector3d F4; Vector3d F5; Vector3d F6; Vector3d F7;

    Vector3d N1; Vector3d N2; Vector3d N3; Vector3d N4; Vector3d N5; Vector3d N6; Vector3d N7;

    Vector3d f1; Vector3d f2; Vector3d f3; Vector3d f4; Vector3d f5; Vector3d f6; Vector3d f7; Vector3d f8;

    Vector3d n1; Vector3d n2; Vector3d n3; Vector3d n4; Vector3d n5; Vector3d n6; Vector3d n7; Vector3d n8;

    double tau1; double tau2; double tau3; double tau4; double tau5; double tau6; double tau7;

    Vector3d Z;

    Matrix3d R01; Matrix3d R12; Matrix3d R23; Matrix3d R34; Matrix3d R45; Matrix3d R56; Matrix3d R67; Matrix3d R78;
    Matrix3d R10; Matrix3d R21; Matrix3d R32; Matrix3d R43; Matrix3d R54; Matrix3d R65; Matrix3d R76; Matrix3d R87;

    Vector3d P10; Vector3d P21; Vector3d P32; Vector3d P43; Vector3d P54; Vector3d P65; Vector3d P76; Vector3d P87;

    Vector3d Pc1; Vector3d Pc2; Vector3d Pc3; Vector3d Pc4; Vector3d Pc5; Vector3d Pc6; Vector3d Pc7;

    double m1; double m2; double m3; double m4; double m5; double m6; double m7;

    Matrix3d Ic1; Matrix3d Ic2; Matrix3d Ic3; Matrix3d Ic4; Matrix3d Ic5; Matrix3d Ic6; Matrix3d Ic7;

    Matrix3d I1; Matrix3d I2; Matrix3d I3; Matrix3d I4; Matrix3d I5; Matrix3d I6; Matrix3d I7;

    Matrix3d I33 = Matrix3d::Identity();

};

}

#endif