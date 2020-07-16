#include <iostream>
#include <cstring>
#include "TestLabAssignmentInstance.h"
#include "TestLabAssignmentSolver.h"

int main(int argc, char* argv[]) {
    using namespace covid;

    if(argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [instance] [transshipment]\n";
        std::cerr << "\t[instance]: path to the instance file;\n";
        std::cerr << "\t[transshipment]: if TRUE transshipments of reagent between labs are enabled.\n";
        std::exit(EXIT_FAILURE);
    }

    const bool transshipments = (std::strcmp(argv[2], "TRUE") == 0);

    const TestLabAssignmentInstance i{argv[1]};
    const TestLabAssignmentSolver s{i};
    const auto sol = s.solve_moo_model(transshipments);

    if(sol) {
        sol->to_json();
        sol->reg_tests_to_csv();
    }

    return 0;
}
