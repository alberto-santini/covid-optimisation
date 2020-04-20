#include <iostream>
#include "TestLabAssignmentInstance.h"
#include "TestLabAssignmentSolver.h"

int main(int argc, char* argv[]) {
    using namespace covid;

    if(argc != 2) {
        std::cerr << "Usage: " << argv[0] << " [instance]\n";
        std::exit(EXIT_FAILURE);
    }

    const TestLabAssignmentInstance i{argv[1]};
    const TestLabAssignmentSolver s{i};
    const auto sol = s.solve_moo_model();

    if(sol) {
        sol->to_json();
        sol->reg_tests_to_csv();
    }

    return 0;
}
