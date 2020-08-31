//
// Created by han on 2020/8/26.
//

#ifndef PENDULUM_P_VARIABLES_H
#define PENDULUM_P_VARIABLES_H


#include <ifopt/variable_set.h>


namespace ifopt {

    class p_variables : public VariableSet {
    public:
        p_variables(const std::string &name, int num_var) : ifopt::VariableSet(num_var, name) {

            for (double &i : x_list_) {
                i = 0;
            }
            // for init state
            x_list_[3] = 3.14159;
            num_var_ = num_var;
        }

        void SetVariables(const VectorXd &x) override {
            for (int i = 0; i < num_var_; ++i) {
                x_list_[i] = x(i);
            }
        };

        VectorXd GetValues() const override {
            return Eigen::Map<const Eigen::VectorXd>(x_list_, num_var_);
        };

        VecBound GetBounds() const override {
            VecBound bounds(GetRows());
            for (int i = 0; i < 5001; ++i) {
                bounds.at(0 + 7 * i) = ifopt::Bounds(-6.28, 6.28);
                bounds.at(1 + 7 * i) = ifopt::NoBound;
                bounds.at(2 + 7 * i) = ifopt::NoBound;
                bounds.at(3 + 7 * i) = ifopt::Bounds(-6.28, 6.28);
                bounds.at(4 + 7 * i) = ifopt::NoBound;
                bounds.at(5 + 7 * i) = ifopt::NoBound;
                bounds.at(6 + 7 * i) = ifopt::Bounds(-100, 100);
            }

            return bounds;
        }

    private:
        double x_list_[35007];
        int num_var_;


    };
}


#endif //PENDULUM_P_VARIABLES_H
