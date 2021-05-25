#include <robot_model.h>
#include <math.h>

namespace robot_dyn
{
RobotModel::RobotModel()
{
    Ps_num = Psi_num*dof;

    qMin.resize(dof); qMax.resize(dof);
    qDotMin.resize(dof); qDotMax.resize(dof);
    qDDotMin.resize(dof); qDDotMax.resize(dof);

    qMin << -3.0503, -2.2736, -3.0426, -3.0439, -2.9761, -2.9761, -3.14;
    qMax << 3.0503, 2.2736, 3.0426, 3.0439, 2.9761, 2.9761, 3.14;
    qDotMin << -1.74, -1.328, -1.957, -1.957, -3.485, -3.485, -4.545;
    qDotMax << 1.74, 1.328, 1.957, 1.957, 3.485, 3.485, 4.545;
    qDDotMin << -10.0, -8.0, -10.0, -10.0, -12.0, -12.0, -12.0;
    qDDotMax << 10.0, 8.0, 10.0, 10.0, 12.0, 12.0, 12.0;

    Ps_flag.resize(Ps_num);
}

RobotModel::~RobotModel(){}

void RobotModel::set_kinematics_parameters()
{
    d1 = 0.317; d2 = 0.1925; d3 = 0.4; d4 = -0.1685; d5 = 0.4; d6 = 0.1363; d7 = 0.13375;
    alpha1 = M_PI_2; alpha2 = M_PI_2; alpha3 = M_PI_2; alpha4 = M_PI_2; alpha5 = M_PI_2; alpha6 = M_PI_2; alpha7 = 0;
    a1 = -0.081; a2 = 0; a3 = 0; a4 = 0; a5 = 0; a6 = 0; a7 = 0;
    offset1 = M_PI; offset2 = -M_PI_2; offset3 = M_PI; offset4 = M_PI; offset5 = M_PI; offset6 = M_PI; offset7 = -M_PI_2;

    P10 << -a1, 0, d1;
    P21 << a2, 0, d2;
    P32 << a3, 0, d3;
    P43 << a4, 0, d4;
    P54 << a5, 0, d5;
    P65 << a6, 0, d6;
    P76 << a7, 0, d7;
    P87 << 0, 0, 0;
}

// m Pcx Pcy Pcz Ixx Ixy Ixz Iyy Iyz Izz
void RobotModel::set_dynamics_parameters(const VectorXd param)
{
    m1 = param(0*dof); m2 = param(1*dof); m3 = param(2*dof); m4 = param(3*dof); m5 = param(4*dof); m6 = param(5*dof); m7 = param(6*dof);
    
    Pc1 << param(0*dof+1), param(0*dof+2), param(0*dof+3);
    Pc2 << param(1*dof+1), param(1*dof+2), param(1*dof+3);
    Pc3 << param(2*dof+1), param(2*dof+2), param(2*dof+3);
    Pc4 << param(3*dof+1), param(3*dof+2), param(3*dof+3);
    Pc5 << param(4*dof+1), param(4*dof+2), param(4*dof+3);
    Pc6 << param(5*dof+1), param(5*dof+2), param(5*dof+3);
    Pc7 << param(6*dof+1), param(6*dof+2), param(6*dof+3);

    I1 <<  param(0*dof+4), -param(0*dof+5), -param(0*dof+6),
          -param(0*dof+5),  param(0*dof+7), -param(0*dof+8),
          -param(0*dof+6), -param(0*dof+8),  param(0*dof+9);
    I2 <<  param(1*dof+4), -param(1*dof+5), -param(1*dof+6),
          -param(1*dof+5),  param(1*dof+7), -param(1*dof+8),
          -param(1*dof+6), -param(1*dof+8),  param(1*dof+9);
    I3 <<  param(2*dof+4), -param(2*dof+5), -param(2*dof+6),
          -param(2*dof+5),  param(2*dof+7), -param(2*dof+8),
          -param(2*dof+6), -param(2*dof+8),  param(2*dof+9);
    I4 <<  param(3*dof+4), -param(3*dof+5), -param(3*dof+6),
          -param(3*dof+5),  param(3*dof+7), -param(3*dof+8),
          -param(3*dof+6), -param(3*dof+8),  param(3*dof+9);
    I5 <<  param(4*dof+4), -param(4*dof+5), -param(4*dof+6),
          -param(4*dof+5),  param(4*dof+7), -param(4*dof+8),
          -param(4*dof+6), -param(4*dof+8),  param(4*dof+9);
    I6 <<  param(5*dof+4), -param(5*dof+5), -param(5*dof+6),
          -param(5*dof+5),  param(5*dof+7), -param(5*dof+8),
          -param(5*dof+6), -param(5*dof+8),  param(5*dof+9);
    I7 <<  param(6*dof+4), -param(6*dof+5), -param(6*dof+6),
          -param(6*dof+5),  param(6*dof+7), -param(6*dof+8),
          -param(6*dof+6), -param(6*dof+8),  param(6*dof+9);

    Ic1 = I1-m1*(Pc1.transpose()*Pc1*I33-Pc1*Pc1.transpose());
    Ic2 = I2-m2*(Pc2.transpose()*Pc2*I33-Pc2*Pc2.transpose());
    Ic3 = I3-m3*(Pc3.transpose()*Pc3*I33-Pc3*Pc3.transpose());
    Ic4 = I4-m4*(Pc4.transpose()*Pc4*I33-Pc4*Pc4.transpose());
    Ic5 = I5-m5*(Pc5.transpose()*Pc5*I33-Pc5*Pc5.transpose());
    Ic6 = I6-m6*(Pc6.transpose()*Pc6*I33-Pc6*Pc6.transpose());
    Ic7 = I7-m7*(Pc7.transpose()*Pc7*I33-Pc7*Pc7.transpose());
}

VectorXd RobotModel::calcu_inv_dyn(const VectorXd q, const VectorXd qDot, const VectorXd qDDot)
{
    theta1 = q(0)+offset1;
    theta2 = q(1)+offset2;
    theta3 = q(2)+offset3;
    theta4 = q(3)+offset4;
    theta5 = q(4)+offset5;
    theta6 = q(5)+offset6;
    theta7 = q(6)+offset7;

    R01 << cos(theta1),  -sin(theta1)*cos(alpha1),  sin(theta1)*sin(alpha1),
           sin(theta1),   cos(theta1)*cos(alpha1), -cos(theta1)*sin(alpha1),
           0,             sin(alpha1),              cos(alpha1);
    R12 << cos(theta2),  -sin(theta2)*cos(alpha2),  sin(theta2)*sin(alpha2),
           sin(theta2),   cos(theta2)*cos(alpha2), -cos(theta2)*sin(alpha2),
           0,             sin(alpha2),              cos(alpha2);
    R23 << cos(theta3),  -sin(theta3)*cos(alpha3),  sin(theta3)*sin(alpha3),
           sin(theta3),   cos(theta3)*cos(alpha3), -cos(theta3)*sin(alpha3),
           0,             sin(alpha3),              cos(alpha3);
    R34 << cos(theta4),  -sin(theta4)*cos(alpha4),  sin(theta4)*sin(alpha4),
           sin(theta4),   cos(theta4)*cos(alpha4), -cos(theta4)*sin(alpha4),
           0,             sin(alpha4),              cos(alpha4);
    R45 << cos(theta5),  -sin(theta5)*cos(alpha5),  sin(theta5)*sin(alpha5),
           sin(theta5),   cos(theta5)*cos(alpha5), -cos(theta5)*sin(alpha5),
           0,             sin(alpha5),              cos(alpha5);
    R56 << cos(theta6),  -sin(theta6)*cos(alpha6),  sin(theta6)*sin(alpha6),
           sin(theta6),   cos(theta6)*cos(alpha6), -cos(theta6)*sin(alpha6),
           0,             sin(alpha6),              cos(alpha6);
    R67 << cos(theta7),  -sin(theta7)*cos(alpha7),  sin(theta7)*sin(alpha7),
           sin(theta7),   cos(theta7)*cos(alpha7), -cos(theta7)*sin(alpha7),
           0,             sin(alpha7),              cos(alpha7);
    R78 << 1, 0, 0,
           0, 1, 0,
           0, 0, 1;

    R10 = R01.transpose();
    R21 = R12.transpose();
    R32 = R23.transpose();
    R43 = R34.transpose();
    R54 = R45.transpose();
    R65 = R56.transpose();
    R76 = R67.transpose();
    R87 = R78.transpose();

    Z << 0, 0, 1;

    w0 << 0, 0, 0; wDot0 << 0, 0, 0; vDot0 << 0, 0, g;
    f8 << 0, 0, 0; n8 << 0, 0, 0;

    w1 = R01*w0+qDot(0)*Z;
    w2 = R12*w1+qDot(1)*Z;
    w3 = R23*w2+qDot(2)*Z;
    w4 = R34*w3+qDot(3)*Z;
    w5 = R45*w4+qDot(4)*Z;
    w6 = R56*w5+qDot(5)*Z;
    w7 = R67*w6+qDot(6)*Z;

    wDot1 = R01*wDot0+R01*w0.cross(qDot(0)*Z)+qDDot(0)*Z;
    wDot2 = R12*wDot1+R12*w1.cross(qDot(1)*Z)+qDDot(1)*Z;
    wDot3 = R23*wDot2+R23*w2.cross(qDot(2)*Z)+qDDot(2)*Z;
    wDot4 = R34*wDot3+R34*w3.cross(qDot(3)*Z)+qDDot(3)*Z;
    wDot5 = R45*wDot4+R45*w4.cross(qDot(4)*Z)+qDDot(4)*Z;
    wDot6 = R56*wDot5+R56*w5.cross(qDot(5)*Z)+qDDot(5)*Z;
    wDot7 = R67*wDot6+R67*w6.cross(qDot(6)*Z)+qDDot(6)*Z;

    vDot1 = R01*(wDot0.cross(P10)+w0.cross(w0.cross(P10))+vDot0);
    vDot2 = R12*(wDot1.cross(P21)+w1.cross(w1.cross(P21))+vDot1);
    vDot3 = R23*(wDot2.cross(P32)+w2.cross(w2.cross(P32))+vDot2);
    vDot4 = R34*(wDot3.cross(P43)+w3.cross(w3.cross(P43))+vDot3);
    vDot5 = R45*(wDot4.cross(P54)+w4.cross(w4.cross(P54))+vDot4);
    vDot6 = R56*(wDot5.cross(P65)+w5.cross(w5.cross(P65))+vDot5);
    vDot7 = R67*(wDot6.cross(P76)+w6.cross(w6.cross(P76))+vDot6);

    vcDot1 = wDot1.cross(Pc1)+w1.cross(w1.cross(Pc1))+vDot1;
    vcDot2 = wDot2.cross(Pc2)+w2.cross(w2.cross(Pc2))+vDot2;
    vcDot3 = wDot3.cross(Pc3)+w3.cross(w3.cross(Pc3))+vDot3;
    vcDot4 = wDot4.cross(Pc4)+w4.cross(w4.cross(Pc4))+vDot4;
    vcDot5 = wDot5.cross(Pc5)+w5.cross(w5.cross(Pc5))+vDot5;
    vcDot6 = wDot6.cross(Pc6)+w6.cross(w6.cross(Pc6))+vDot6;
    vcDot7 = wDot7.cross(Pc7)+w7.cross(w7.cross(Pc7))+vDot7;

    F1 = m1*vcDot1;
    F2 = m2*vcDot2;
    F3 = m3*vcDot3;
    F4 = m4*vcDot4;
    F5 = m5*vcDot5;
    F6 = m6*vcDot6;
    F7 = m7*vcDot7;

    N1 = Ic1*wDot1+w1.cross(Ic1*w1);
    N2 = Ic2*wDot2+w2.cross(Ic2*w2);
    N3 = Ic3*wDot3+w3.cross(Ic3*w3);
    N4 = Ic4*wDot4+w4.cross(Ic4*w4);
    N5 = Ic5*wDot5+w5.cross(Ic5*w5);
    N6 = Ic6*wDot6+w6.cross(Ic6*w6);
    N7 = Ic7*wDot7+w7.cross(Ic7*w7);

    f7 = R87*f8+F7;
    f6 = R76*f7+F6;
    f5 = R65*f6+F5;
    f4 = R54*f5+F4;
    f3 = R43*f4+F3;
    f2 = R32*f3+F2;
    f1 = R21*f2+F1;

    n7 = N7+R87*n8+Pc7.cross(F7)+P87.cross(R87*f8);
    n6 = N6+R76*n7+Pc6.cross(F6)+P76.cross(R76*f7);
    n5 = N5+R65*n6+Pc5.cross(F5)+P65.cross(R65*f6);
    n4 = N4+R54*n5+Pc4.cross(F4)+P54.cross(R54*f5);
    n3 = N3+R43*n4+Pc3.cross(F3)+P43.cross(R43*f4);
    n2 = N2+R32*n3+Pc2.cross(F2)+P32.cross(R32*f3);
    n1 = N1+R21*n2+Pc1.cross(F1)+P21.cross(R21*f2);

    tau7 = n7.transpose()*Z;
    tau6 = n6.transpose()*Z;
    tau5 = n5.transpose()*Z;
    tau4 = n4.transpose()*Z;
    tau3 = n3.transpose()*Z;
    tau2 = n2.transpose()*Z;
    tau1 = n1.transpose()*Z;

    VectorXd tau = VectorXd::Zero(dof);
    tau << tau1, tau2, tau3, tau4, tau5, tau6, tau7;

    return tau;
}

}