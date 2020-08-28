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
        double f_data_[8];
        double f_ret_[7];
        for (int i = 0; i < 101; ++i) {

            for (int j = 0; j < 7; ++j) {
                f_data_[j] = x(i * 7 + j);
            }
            f_data_[7] = 9.81;

            f1(f_data_, f_ret_, nullptr, nullptr, 0);
            g(i * 2 + 408 + 0) = f_ret_[0];
            f2(f_data_, f_ret_, nullptr, nullptr, 0);
            g(i * 2 + 408 + 1) = f_ret_[0];
        }

        /*for (int i = 0; i < 101; ++i) {

            g(i * 2 + 408 + 0) = 0.0;
            g(i * 2 + 408 + 1) = 0.0;
        }*/

        return g;
    };

    ifopt::Component::VecBound p_constraint::GetBounds() const {
        VecBound b(GetRows());
        for (int i = 0; i < GetRows(); ++i) {
            b.at(i) = ifopt::Bounds(-0.001, 0.001);
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

            double f_data_[8];
            double f_ret_[7];
            double inter_item_[321];
            // 3. for dynamic equation
            for (int i = 0; i < 101; ++i) {

                for (int j = 0; j < 7; ++j) {
                    f_data_[j] = x(i * 7 + j);
                }
                f_data_[7] = 9.81;

                jac_f1(f_data_, f_ret_, nullptr, inter_item_, 0);


                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 0) = f_ret_[0];
                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 1) = f_ret_[1];
                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 2) = f_ret_[2];
                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 3) = f_ret_[3];
                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 4) = f_ret_[4];
                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 5) = 0.0;
                jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 6) = f_ret_[5];
                //jac_block.coeffRef(i * 2 + 408 + 0, 7 * i + 7) = f_ret_[6];


                jac_f2(f_data_, f_ret_, nullptr, inter_item_, 0);


                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 0) = f_ret_[0];
                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 1) = f_ret_[1];
                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 2) = 0.0;
                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 3) = f_ret_[2];
                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 4) = f_ret_[3];
                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 5) = f_ret_[4];
                jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 6) = f_ret_[5];
                //jac_block.coeffRef(i * 2 + 408 + 1, 7 * i + 7) = f_ret_[6];



            }
        }
    }
}
