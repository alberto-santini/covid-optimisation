//
// Created by alberto on 12/04/2020.
//

#ifndef COVID_TESTLABASSIGNMENTSOLVER_H
#define COVID_TESTLABASSIGNMENTSOLVER_H

#include "TestLabAssignmentInstance.h"
#include "TestLabAssignmentSolution.h"
#include <optional>

namespace covid {
    struct TestLabAssignmentSolver {
        const TestLabAssignmentInstance& i;

        explicit TestLabAssignmentSolver(const TestLabAssignmentInstance& instance) : i{instance} {}

        std::optional<TestLabAssignmentSolution> solve_moo_model() const;
    };
}

#endif //COVID_TESTLABASSIGNMENTSOLVER_H
