//
// Created by alberto on 12/04/2020.
//

#ifndef COVID_TESTLABASSIGNMENTSOLUTION_H
#define COVID_TESTLABASSIGNMENTSOLUTION_H

#include "TestLabAssignmentInstance.h"
#include <iostream>

namespace covid {
    enum class Material {
        SWAB,
        REAGENT
    };

    struct MaterialMovement {
        Material material;
        size_t origin_id;
        size_t destination_id;
        size_t day;
        uint32_t qty;
    };

    struct MaterialStorage {
        Material material;
        size_t facility_id;
        size_t day; // Stored from day to day+1
        uint32_t qty;
    };

    struct TestsAction {
        size_t lab_id;
        size_t day;
        uint32_t qty;
    };

    struct RegionTestsCsvRow {
        std::string reg_name;
        uint32_t reg_id;
        uint32_t day;
        uint32_t tests;
    };

    struct TestLabAssignmentSolution {
        const TestLabAssignmentInstance& instance;
        std::vector<MaterialMovement> reagent_fac_lab_mov;
        std::vector<MaterialMovement> reagent_lab_lab_mov;
        std::vector<MaterialMovement> swab_mov;
        std::vector<MaterialStorage> reagent_fac_store;
        std::vector<MaterialStorage> reagent_lab_store;
        std::vector<MaterialStorage> swab_lab_store;
        std::vector<TestsAction> tests_assigned;
        std::vector<TestsAction> tests_performed;

        std::vector<RegionTestsCsvRow> reg_tests_csv;

        uint32_t untested_swabs_at_end;
        uint32_t total_waiting_swabs;

        explicit TestLabAssignmentSolution(const TestLabAssignmentInstance& instance) :
            instance{instance}, untested_swabs_at_end{0u}, total_waiting_swabs{0u} {}

        void to_json() const;
        void reg_tests_to_csv() const;
    };

    std::ostream& operator<<(std::ostream& out, const RegionTestsCsvRow& r);
    std::ostream& operator<<(std::ostream& out, const TestLabAssignmentSolution& s);
}

#endif //COVID_TESTLABASSIGNMENTSOLUTION_H
