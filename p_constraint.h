//
// Created by han on 2020/8/26.
//

#ifndef PENDULUM_P_CONSTRAINT_H
#define PENDULUM_P_CONSTRAINT_H

#include <ifopt/constraint_set.h>
#include <casadi/casadi.hpp>


namespace ifopt{
using namespace casadi;

class p_constraint : public ifopt::ConstraintSet{

public:
    p_constraint(const std::string& name, int num_constr);

    VectorXd GetValues() const override;

    VecBound GetBounds() const override;

    void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override;

private:
    int num_constr_ = 0;
    SX tau = SX::sym("tau");
    SX phi = SX::sym("phi");
    SX theta = SX::sym("theta");
    SX dphi = SX::sym("dphi");
    SX dtheta = SX::sym("dtheta");
    SX ddphi = SX::sym("ddphi");
    SX ddtheta = SX::sym("ddtheta");
    SX g = SX::sym("g");
    Function f1 = Function("f1",{theta, dtheta, ddtheta, phi, dphi, ddphi,tau, g},{-ddtheta + ((tau*2.50325E+7+tau*pow(cos(phi),2.0)*9.17E+6-(dphi*dphi)*sin(phi)*2.50325E+7+tau*pow(cos(theta),2.0)*9.17E+6-(dphi*dphi)*sin(theta)*9.181921E+6+(dphi*dphi)*sin(phi*2.0)*cos(phi)*4.585E+6-(dphi*dphi)*pow(cos(phi),2.0)*sin(phi)*9.17E+6+(dtheta*dtheta)*sin(phi*2.0)*cos(phi)*1.7085E+7-tau*pow(cos(phi),2.0)*pow(cos(theta),2.0)*9.17E+6+(dphi*dphi)*pow(cos(phi),2.0)*sin(theta)*1.5000286E+7-(dphi*dphi)*pow(cos(theta),2.0)*sin(phi)*9.17E+6+(dphi*dphi)*pow(cos(phi),4.0)*sin(theta)*6.727112E+6-(dphi*dphi)*pow(cos(theta),2.0)*sin(theta)*3.363556E+6-(dphi*dphi)*pow(cos(theta),3.0)*sin(theta)*3.363556E+6+g*cos(phi)*sin(phi)*2.5E+7-(dphi*dphi)*cos(theta)*sin(theta)*9.181921E+6+(dtheta*dtheta)*pow(cos(phi),2.0)*cos(theta)*sin(phi)*9.17E+6+(dphi*dphi)*pow(cos(phi),2.0)*cos(theta)*sin(theta)*5.818365E+6+(dphi*dphi)*pow(cos(phi),4.0)*cos(theta)*sin(theta)*3.363556E+6-dphi*dtheta*pow(cos(phi),3.0)*sin(phi)*2.5067112E+7+dphi*dtheta*sin(theta*2.0)*cos(phi)*9.17E+6+(dphi*dphi)*pow(cos(phi),2.0)*pow(cos(theta),2.0)*sin(theta)*1.0090668E+7+(dphi*dphi)*pow(cos(phi),2.0)*pow(cos(theta),3.0)*sin(theta)*6.727112E+6-(dphi*dphi)*pow(cos(phi),4.0)*pow(cos(theta),2.0)*sin(theta)*6.727112E+6-(dphi*dphi)*pow(cos(phi),4.0)*pow(cos(theta),3.0)*sin(theta)*3.363556E+6+g*cos(phi)*pow(sin(phi),2.0)*sin(theta)*9.17E+6-dphi*dtheta*cos(phi)*sin(phi)*6.8428842E+7+(dphi*dphi)*sin(phi*2.0)*cos(phi)*sin(phi)*sin(theta)*1.681778E+6+(dtheta*dtheta)*sin(phi*2.0)*cos(phi)*sin(phi)*sin(theta)*6.266778E+6-dphi*dtheta*cos(phi)*pow(cos(theta),2.0)*sin(phi)*2.5067112E+7-dphi*dtheta*pow(cos(phi),3.0)*cos(theta)*sin(theta)*1.834E+7+(dtheta*dtheta)*pow(cos(phi),2.0)*cos(theta)*pow(sin(phi),2.0)*sin(theta)*3.363556E+6+dphi*dtheta*pow(cos(phi),3.0)*pow(cos(theta),2.0)*sin(phi)*2.5067112E+7-(dphi*dphi)*pow(cos(phi),2.0)*pow(cos(theta),2.0)*pow(sin(phi),2.0)*sin(theta)*3.363556E+6+dphi*dtheta*sin(theta*2.0)*cos(phi)*sin(phi)*sin(theta)*3.363556E+6-dphi*dtheta*pow(cos(phi),3.0)*cos(theta)*sin(phi)*pow(sin(theta),2.0)*6.727112E+6)*-4.0)/(pow(cos(phi),2.0)*1.37744656E+8+pow(cos(phi),4.0)*5.0134224E+7-pow(cos(theta),2.0)*9.9113028E+7+pow(cos(phi),2.0)*pow(cos(theta),2.0)*1.49247252E+8-pow(cos(phi),4.0)*pow(cos(theta),2.0)*5.0134224E+7+pow(cos(phi),2.0)*pow(sin(phi),2.0)*pow(sin(theta),2.0)*1.3454224E+7+pow(cos(phi),2.0)*sin(phi)*sin(theta)*7.336E+7-2.70561273E+8)});
    Function  f2 = Function("f2",{theta, dtheta, ddtheta, phi, dphi, ddphi,tau, g},{-ddphi + (((dphi*dphi)*sin(phi*2.0)*-2.4778257E+7-(dtheta*dtheta)*sin(phi*2.0)*9.2330757E+7-tau*cos(phi)*5.0E+7-g*sin(phi)*1.35105E+8-(dphi*dphi)*pow(cos(phi),3.0)*sin(theta)*3.668E+7-dphi*dtheta*sin(theta*2.0)*4.9556514E+7+(dphi*dphi)*sin(phi*2.0)*pow(cos(phi),2.0)*1.2533556E+7+(dtheta*dtheta)*sin(phi*2.0)*pow(cos(phi),2.0)*4.6703556E+7+(dphi*dphi)*cos(phi)*sin(phi)*5.0E+7+g*pow(cos(phi),2.0)*sin(phi)*6.834E+7+(dphi*dphi)*cos(phi)*sin(theta)*1.834E+7+(dphi*dphi)*cos(phi)*pow(cos(theta),2.0)*sin(phi)*4.9556514E+7+(dtheta*dtheta)*pow(cos(phi),3.0)*cos(theta)*sin(phi)*2.5067112E+7-(dphi*dphi)*pow(cos(phi),3.0)*cos(theta)*sin(theta)*1.834E+7+(dphi*dphi)*cos(phi)*sin(phi)*pow(sin(theta),2.0)*6.727112E+6+(dphi*dphi)*cos(phi)*pow(sin(phi),2.0)*sin(theta)*1.834E+7+dphi*dtheta*pow(cos(phi),2.0)*sin(phi)*1.3668E+8-(dphi*dphi)*pow(cos(phi),3.0)*pow(cos(theta),2.0)*sin(phi)*2.5067112E+7-tau*cos(phi)*sin(phi)*sin(theta)*1.834E+7-(dphi*dphi)*pow(cos(phi),3.0)*sin(phi)*pow(sin(theta),2.0)*1.3454224E+7+dphi*dtheta*sin(theta*2.0)*pow(cos(phi),2.0)*2.5067112E+7-(dtheta*dtheta)*cos(phi)*cos(theta)*sin(phi)*4.9556514E+7+(dphi*dphi)*cos(phi)*cos(theta)*sin(theta)*1.834E+7+dphi*dtheta*pow(cos(phi),2.0)*cos(theta)*sin(theta)*9.9113028E+7-dphi*dtheta*pow(cos(phi),4.0)*cos(theta)*sin(theta)*5.0134224E+7-(dphi*dphi)*pow(cos(phi),3.0)*cos(theta)*sin(phi)*pow(sin(theta),2.0)*6.727112E+6+dphi*dtheta*pow(cos(phi),2.0)*pow(sin(phi),2.0)*sin(theta)*5.0134224E+7+(dphi*dphi)*cos(phi)*cos(theta)*sin(phi)*pow(sin(theta),2.0)*6.727112E+6)*2.0)/(pow(cos(phi),2.0)*1.37744656E+8+pow(cos(phi),4.0)*5.0134224E+7-pow(cos(theta),2.0)*9.9113028E+7+pow(cos(phi),2.0)*pow(cos(theta),2.0)*1.49247252E+8-pow(cos(phi),4.0)*pow(cos(theta),2.0)*5.0134224E+7+pow(cos(phi),2.0)*pow(sin(phi),2.0)*pow(sin(theta),2.0)*1.3454224E+7+pow(cos(phi),2.0)*sin(phi)*sin(theta)*7.336E+7-2.70561273E+8)});
    SX x = SX::sym("x");
    SX xn = SX::sym("xn");
    SX dx = SX::sym("dx");
    SX dxn = SX::sym("dxn");
    Function j_f1 = f1.jacobian();
    Function j_f2 = f2.jacobian();
    std::vector<DM> res_ = {0};

};

}

#endif //PENDULUM_P_CONSTRAINT_H
