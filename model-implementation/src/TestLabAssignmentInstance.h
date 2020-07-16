//
// Created by alberto on 11/04/2020.
//

#ifndef COVID_TESTLABASSIGNMENTINSTANCE_H
#define COVID_TESTLABASSIGNMENTINSTANCE_H

#include <cstddef>
#include <vector>
#include <cstdint>
#include <string>
#include <iostream>

namespace covid {
    struct TestLabAssignmentInstance {
        std::string inst_file;

        size_t n_labs;
        size_t n_regions;
        size_t n_factories;
        size_t n_days;

        std::vector<std::string> reg_names;

        std::vector<uint32_t> lab_region;
        std::vector<uint32_t> lab_capacity;
        std::vector<uint32_t> lab_start_reagents;
        std::vector<uint32_t> fac_start_reagents;
        std::vector<uint32_t> reg_max_inbound_reagents;
        std::vector<uint32_t> reg_max_inbound_swabs;

        std::vector<std::vector<uint32_t>> fac_day_production;
        std::vector<std::vector<uint32_t>> reg_day_demand;
        std::vector<std::vector<uint32_t>> reg_labs;
        std::vector<std::vector<bool>> lab_lab_compatible;
        std::vector<std::vector<bool>> fac_lab_compatible;
        std::vector<std::vector<float>> lab_lab_distance;
        std::vector<std::vector<float>> fac_lab_distance;

        explicit TestLabAssignmentInstance(const std::string& inst_file);
    };

    std::ostream& operator<<(std::ostream& out, const TestLabAssignmentInstance& i);
}

#endif //COVID_TESTLABASSIGNMENTINSTANCE_H
