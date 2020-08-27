//
// Created by han on 2020/8/27.
//

#ifndef PENDULUM_P_COST_H
#define PENDULUM_P_COST_H

#include <ifopt/cost_term.h>

namespace ifopt {

    class p_cost : public CostTerm {
    public:
        p_cost(const std::string &name) : CostTerm(name) {}

        double GetCost() const override {
            Eigen::VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
            double cost = 0;
            for (int i = 0; i < 101; ++i) {
                cost += std::abs(x((6 + i * 7)));
            }
            return cost;
        };

        void FillJacobianBlock(std::string var_set, Jacobian &jac) const override {
            if (var_set == "var_set1") {
                Eigen::VectorXd x = GetVariables()->GetComponent("var_set1")->GetValues();
                for (int i = 0; i < 101; ++i) {
                    if (x((6 + i * 7)) >= 0) {
                        jac.coeffRef(0, (6 + i * 7)) = 1;
                    } else {
                        jac.coeffRef(0, (6 + i * 7)) = -1;
                    }

                }
            }
        }
    };
}


#endif //PENDULUM_P_COST_H
