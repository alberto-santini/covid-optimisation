//
// Created by alberto on 11/04/2020.
//

#include "TestLabAssignmentInstance.h"
#include <json.hpp>
#include <fstream>
#include <iostream>

namespace covid {
    TestLabAssignmentInstance::TestLabAssignmentInstance(const std::string& inst_file) : inst_file{inst_file} {
        using nlohmann::json;

        std::ifstream ifs{inst_file};
        if(ifs.fail()) {
            std::cerr << "Error opening instance file " << inst_file << "\n";
            std::exit(EXIT_FAILURE);
        }

        json j;

        if(!(ifs >> j)) {
            std::cerr << "Invalid json in instance file " << inst_file << "\n";
            std::exit(EXIT_FAILURE);
        }

        n_labs = j["n_labs"].get<size_t>();
        n_regions = j["n_regions"].get<size_t>();
        n_factories = j["n_factories"].get<size_t>();
        n_days = j["n_days"].get<size_t>();

        if(j.find("reg_names") != j.end()) {
            reg_names = j["reg_names"].get<std::vector<std::string>>();
        } else {
            for(auto reg_id = 0u; reg_id < n_regions; ++reg_id) {
                reg_names.push_back("Region " + std::to_string(reg_id));
            }
        }

        lab_region = j["lab_region"].get<std::vector<std::uint32_t>>();
        lab_capacity = j["lab_capacity"].get<std::vector<std::uint32_t>>();
        lab_start_reagents = j["lab_start_reagents"].get<std::vector<std::uint32_t>>();
        fac_start_reagents = j["fac_start_reagents"].get<std::vector<std::uint32_t>>();
        reg_max_inbound_reagents = j["reg_max_inbound_reagents"].get<std::vector<std::uint32_t>>();
        reg_max_inbound_swabs = j["reg_max_inbound_swabs"].get<std::vector<std::uint32_t>>();

        reg_labs.resize(n_regions);
        reg_day_demand.resize(n_regions);
        for(auto r = 0u; r < n_regions; ++r) {
            for(auto l = 0u; l < n_labs; ++l) {
                if(lab_region[l] == r) {
                    reg_labs[r].push_back(l);
                }
            }
            reg_day_demand[r] = j["reg_day_demand"][r].get<std::vector<uint32_t>>();
        }

        lab_lab_compatible.resize(n_labs);
        lab_lab_distance.resize(n_labs);
        for(auto l1 = 0u; l1 < n_labs; ++l1) {
            lab_lab_compatible[l1] = std::vector<bool>(n_labs);
            lab_lab_distance[l1] = j["lab_lab_distance"][l1].get<std::vector<float>>();

            auto jrow = j["lab_lab_compatible"][l1].get<std::vector<uint32_t>>();
            for(auto l2 = 0u; l2 < n_labs; ++l2) {
                lab_lab_compatible[l1][l2] = (jrow[l2] == 1u);
            }
        }

        fac_lab_compatible.resize(n_factories);
        fac_lab_distance.resize(n_factories);
        fac_day_production.resize(n_factories);
        for(auto r = 0u; r < n_factories; ++r) {
            fac_lab_compatible[r] = std::vector<bool>(n_labs);
            fac_lab_distance[r] = j["fac_lab_distance"][r].get<std::vector<float>>();
            fac_day_production[r] = j["fac_day_production"][r].get<std::vector<uint32_t>>();

            const auto jrow = j["fac_lab_compatible"][r].get<std::vector<uint32_t>>();
            for(auto l = 0u; l < n_labs; ++l) {
                fac_lab_compatible[r][l] = (jrow[l] == 1u);
            }
        }
    }

    std::ostream& operator<<(std::ostream& out, const TestLabAssignmentInstance& i) {
        out << i.n_labs << " labs in " << i.n_regions << " regions\n";
        out << i.n_factories << " factories\n";
        out << i.n_days << " days\n";

        out << "\n--- Regions ---\n";
        for(auto j = 0u; j < i.n_regions; ++j) {
            out << "\tRegion " << j << " (" << i.reg_names[j] << ")\n";
            out << "\t\tMax inbound reagents: " << i.reg_max_inbound_reagents[j] << "\n";
            out << "\t\tMax inbound swabs: " << i.reg_max_inbound_swabs[j] << "\n";
            out << "\t\tLabs: ";
            for(const auto& l : i.reg_labs[j]) {
                out << l << " ";
            }
            out << "\n";
            out << "\t\tDemands: ";
            for(const auto& d : i.reg_day_demand[j]) {
                out << d << " ";
            }
            out << "\n";
        }

        out << "\n--- Labs ---\n";
        for(auto l = 0u; l < i.n_labs; ++l) {
            out << "\tLab " << l << "\n";
            out << "\t\tRegion: " << i.lab_region[l] << "\n";
            out << "\t\tCapacity: " << i.lab_capacity[l] << "\n";
            out << "\t\tStart reagents: " << i.lab_start_reagents[l] << "\n";
        }

        out << "\n--- Factories ---\n";
        for(auto r = 0u; r < i.n_factories; ++r) {
            out << "\tFactory " << r << "\n";
            out << "\t\tDaily production: ";
            for(const auto& p : i.fac_day_production[r]) {
                out << p << " ";
            }
            out << "\n";
            out << "\t\tStart reagents: " << i.fac_start_reagents[r] << "\n";
        }

        return out;
    }
}