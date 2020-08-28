//
// Created by han on 2020/8/26.
//

#include "p_constraint.h"

namespace ifopt {

    p_constraint::p_constraint(const std::string &name, int num_constr) : ConstraintSet(num_constr, name) {
        num_constr_ = num_constr;
    }

    Eigen::VectorXd p_constraint::GetValues() const {
        Eigen::VectorXd g(GetRows());
        Eigen::VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
        // 1. for init and state constraint
        g(0) = x(0);
        g(1) = x(1);
        g(2) = x(3) - 3.14159;
        g(3) = x(4);

        g(4) = x(700);
        g(5) = x(701);
        g(6) = x(703);
        g(7) = x(704);

        // 2. for directCollocation, we use Trapezoidal method
        for (int i = 0; i < 100; ++i) {
            // 1. for theta dtheta
            g(i * 4 + 8 + 0) = (x((i + 1) * 7 + 0) - x(i * 7 + 0) - 0.05 * (x((i + 1) * 7 + 1) + x(i * 7 + 1)));

            // 2. for dtheta ddtheta
            g(i * 4 + 8 + 1) = (x((i + 1) * 7 + 1) - x(i * 7 + 1) - 0.05 * (x((i + 1) * 7 + 2) + x(i * 7 + 2)));

            // 3. for phi dphi
            g(i * 4 + 8 + 2) = (x((i + 1) * 7 + 3) - x(i * 7 + 3) - 0.05 * (x((i + 1) * 7 + 4) + x(i * 7 + 4)));

            // 4. for dphi ddphi
            g(i * 4 + 8 + 3) = (x((i + 1) * 7 + 4) - x(i * 7 + 4) - 0.05 * (x((i + 1) * 7 + 5) + x(i * 7 + 5)));

        }

        // 3. for dynamic equation
        for (int i = 0; i < 101; ++i) {
            /*std::vector<DM> f_arg = {x(i * 7 + 0), x(i * 7 + 1), x(i * 7 + 2), x(i * 7 + 3), x(i * 7 + 4), x(i * 7 + 5),
                                     9.81};
            std::vector<DM> res_1 = f1(f_arg);
            std::vector<DM> res_2 = f2(f_arg);*/

            g(i * 2 + 408 + 0) = 0.0;
            g(i * 2 + 408 + 1) = 0.0;
        }

        return g;
    };

    ifopt::Component::VecBound p_constraint::GetBounds() const {
        VecBound b(GetRows());
        for (int i = 0; i < GetRows(); ++i) {
            b.at(i) = ifopt::Bounds(-0.0001, 0.0001);
        }
        return b;
    }

    void p_constraint::FillJacobianBlock(std::string var_set, Jacobian &jac_block) const {
        if (var_set == "var_set1") {
            VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();

            //jac_block.coeffRef(0, 0) = 2.0*x(0); // derivative of first constraint w.r.t x0
            //jac_block.coeffRef(0, 1) = 1.0;      // derivative of first constraint w.r.t x1
            // 1. for init and state constraint

            jac_block.coeffRef(0, 0) = 1.0;
            jac_block.coeffRef(1, 1) = 1.0;
            jac_block.coeffRef(2, 3) = 1.0;
            jac_block.coeffRef(3, 4) = 1.0;

            jac_block.coeffRef(4, 700) = 1.0;
            jac_block.coeffRef(5, 701) = 1.0;
            jac_block.coeffRef(6, 703) = 1.0;
            jac_block.coeffRef(7, 704) = 1.0;


            // 2. for directCollocation, we use Trapezoidal method
            for (int i = 0; i < 100; ++i) {
                // 1. for theta dtheta
                //g(i*4+8+0) = (x((i+1)*7+0) - x(i*7+0)-0.05*(x((i+1)*7+1) + x(i*7+1)));
                jac_block.coeffRef(i * 4 + 8 + 0, (i + 1) * 7 + 0) = 1;
                jac_block.coeffRef(i * 4 + 8 + 0, i * 7 + 0) = -1;
                jac_block.coeffRef(i * 4 + 8 + 0, (i + 1) * 7 + 1) = -0.05;
                jac_block.coeffRef(i * 4 + 8 + 0, i * 7 + 1) = -0.05;


                // 2. for dtheta ddtheta
                //g(i*4+8+1) = (x((i+1)*7+1) - x(i*7+1)-0.05*(x((i+1)*7+2) + x(i*7+2)));
                jac_block.coeffRef(i * 4 + 8 + 1, (i + 1) * 7 + 1) = 1;
                jac_block.coeffRef(i * 4 + 8 + 1, i * 7 + 1) = -1;
                jac_block.coeffRef(i * 4 + 8 + 1, (i + 1) * 7 + 2) = -0.05;
                jac_block.coeffRef(i * 4 + 8 + 1, i * 7 + 2) = -0.05;

                // 3. for phi dphi
                // g(i*4+8+2) = (x((i+1)*7+3) - x(i*7+3)-0.05*(x((i+1)*7+4) + x(i*7+4)));

                jac_block.coeffRef(i * 4 + 8 + 2, (i + 1) * 7 + 3) = 1;
                jac_block.coeffRef(i * 4 + 8 + 2, i * 7 + 3) = -1;
                jac_block.coeffRef(i * 4 + 8 + 2, (i + 1) * 7 + 4) = -0.05;
                jac_block.coeffRef(i * 4 + 8 + 2, i * 7 + 4) = -0.05;

                // 4. for dphi ddphi
                // g(i*4+8+3) = (x((i+1)*7+4) - x(i*7+4)-0.05*(x((i+1)*7+5) + x(i*7+5)));

                jac_block.coeffRef(i * 4 + 8 + 3, (i + 1) * 7 + 4) = 1;
                jac_block.coeffRef(i * 4 + 8 + 3, i * 7 + 4) = -1;
                jac_block.coeffRef(i * 4 + 8 + 3, (i + 1) * 7 + 5) = -0.05;
                jac_block.coeffRef(i * 4 + 8 + 3, i * 7 + 5) = -0.05;

            }
            /*
            // 3. for dynamic equation
            for (int i = 0; i < 101; ++i) {
                std::vector<DM> f_arg = {x(i * 7 + 0), x(i * 7 + 1), x(i * 7 + 2), x(i * 7 + 3), x(i * 7 + 4),
                                         x(i * 7 + 5), 9.81};
                std::vector<DM> res_1 = f1(f_arg);
                std::vector<DM> res_2 = f2(f_arg);

                //g(i*2+408+0) = res_1[0]->at(0);
                //g(i*2+408+1) = res_2[0]->at(0);

                // for f1
                std::vector<DM> j1_arg = {x(i * 7 + 0), x(i * 7 + 1), x(i * 7 + 2), x(i * 7 + 3), x(i * 7 + 4),
                                          x(i * 7 + 5), 9.81, res_1[0]->at(0)};
                std::vector<DM> j2_arg = {x(i * 7 + 0), x(i * 7 + 1), x(i * 7 + 2), x(i * 7 + 3), x(i * 7 + 4),
                                          x(i * 7 + 5), 9.81, res_2[0]->at(0)};
                std::vector<DM> res_j1 = j_f1(j1_arg);
                std::vector<DM> res_j2 = j_f2(j2_arg);

                for (int j = 0; j < 7; ++j) {
                    jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + j) = res_j1[0].get_elements()[j];
                }

                for (int j = 0; j < 7; ++j) {
                    jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + j) = res_j2[0].get_elements()[j];
                }
            }*/
        }
    }
}
