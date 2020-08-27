#include <iostream>
#include <memory>

#include "p_variables.h"
#include "p_constraint.h"
#include "p_cost.h"

#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>

int main() {
    ifopt::Problem nlp;
    nlp.AddVariableSet  (std::make_shared<ifopt::p_variables>("var_set1", 707));
    nlp.AddConstraintSet(std::make_shared<ifopt::p_constraint>("constraint1", 610));
    nlp.AddCostSet      (std::make_shared<ifopt::p_cost>("cost1"));

    // Initialize solver and options
    ifopt::IpoptSolver ipopt;
    ipopt.SetOption("hessian_approximation", "limited-memory");
    ipopt.SetOption("bound_frac", 1e-8);
    ipopt.SetOption("bound_push", 1e-8);

    ipopt.SetOption("tol", 1e-3);

    ipopt.SetOption("dual_inf_tol", 1e-3);
    ipopt.SetOption("constr_viol_tol", 1e-5);
    ipopt.SetOption("compl_inf_tol", 1e-3);

    ipopt.SetOption("max_iter", 1000);

    ipopt.SetOption("print_timing_statistics", "yes");

    ipopt.SetOption("nlp_lower_bound_inf", -1e6);
    ipopt.SetOption("nlp_upper_bound_inf", 1e6);


    // Solve
    ipopt.Solve(nlp);

    std::cout << nlp.GetOptVariables()->GetValues().transpose() << std::endl;
}
