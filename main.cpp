#include <iostream>
#include <memory>

#include "p_variables.h"
#include "p_constraint.h"
#include "p_cost.h"

#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>

#include <fstream>

int main() {
    ifopt::Problem nlp;
    nlp.AddVariableSet  (std::make_shared<ifopt::p_variables>("var_set1", 707));
    nlp.AddConstraintSet(std::make_shared<ifopt::p_constraint>("constraint1", 610));
    nlp.AddCostSet      (std::make_shared<ifopt::p_cost>("cost1"));

    // Initialize solver and options
    ifopt::IpoptSolver ipopt;

    //ipopt.SetOption("hessian_approximation", "limited-memory");

    //ipopt.SetOption("jacobian_approximation", "exact");
    ipopt.SetOption("max_iter", 3000);
    ipopt.SetOption("print_level", 5);




    // Solve
    ipopt.Solve(nlp);

    auto re_data = nlp.GetOptVariables()->GetValues();


    std::ofstream myfile;
    myfile.open ("data.txt");
    for (int i = 0; i < 707; ++i) {
        myfile << re_data[i] << std::endl;
    }
    myfile.close();

}
