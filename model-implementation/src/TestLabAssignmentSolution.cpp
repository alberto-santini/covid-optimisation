//
// Created by alberto on 12/04/2020.
//

#include "TestLabAssignmentSolution.h"
#include <json.hpp>
#include <fstream>
#include <filesystem>

using nlohmann::json;

namespace covid {
    std::ostream& operator<<(std::ostream& out, const TestLabAssignmentSolution& s) {
        out << s.untested_swabs_at_end << "," << s.total_waiting_swabs;

        out << "\n\n--- Labs ---\n";
        for(auto l = 0u; l < s.instance.n_labs; ++l) {
            out << "Lab " << l << "\n";
            for(auto t = 0u; t < s.instance.n_days; ++t) {
                out << "\tDay " << t << "\n";

                if(t == 0u) {
                    out << "\t\tStarted the day with " << s.instance.lab_start_reagents[l] << " units of reagent stored\n";
                } else {
                    for(const auto& st : s.reagent_lab_store) {
                        if(st.facility_id == l && st.day == t - 1u) {
                            out << "\t\tStarted the day with " << st.qty << " units of reagent stored\n";
                        }
                    }
                }

                for(const auto& m : s.reagent_fac_lab_mov) {
                    if(m.destination_id == l && m.day == t) {
                        out << "\t\tReceived " << m.qty << " reagents from Factory " << m.origin_id << "\n";
                    }
                }
                for(const auto& m : s.reagent_lab_lab_mov) {
                    if(m.destination_id == l && m.day == t) {
                        out << "\t\tReceived " << m.qty << " reagents from Lab " << m.origin_id << "\n";
                    }
                    if(m.origin_id == l && m.day == t) {
                        out << "\t\tSent " << m.qty << " reagents to Lab " << m.destination_id << "\n";
                    }
                }
                for(const auto& tst : s.tests_assigned) {
                    if(tst.lab_id == l && tst.day == t) {
                        out << "\t\tAssigned " << tst.qty << " new swabs to test (out of " <<
                               s.instance.reg_day_demand[s.instance.lab_region[l]][t] <<
                               " for its region)\n";
                    }
                }
                for(const auto& m : s.swab_mov) {
                    if(m.origin_id == l && m.day == t) {
                        out << "\t\tSent " << m.qty << " swabs to Lab " << m.destination_id << "\n";
                    }
                    if(m.destination_id == l && m.day == t) {
                        out << "\t\tReceived " << m.qty << " swabs from Lab " << m.origin_id << "\n";
                    }
                }
                for(const auto& tst : s.tests_performed) {
                    if(tst.lab_id == l && tst.day == t) {
                        out << "\t\tPerformed " << tst.qty << " test (capacity = " << s.instance.lab_capacity[l] << ")\n";
                    }
                }
                for(const auto& st : s.reagent_lab_store) {
                    if(st.facility_id == l && st.day == t) {
                        out << "\t\tEnded the day with " << st.qty << " units of reagent stored\n";
                    }
                }
                for(const auto& st : s.swab_lab_store) {
                    if(st.facility_id == l && st.day == t) {
                        out << "\t\tEnded the day with " << st.qty << " swabs stored\n";
                    }
                }
            }
        }

        return out;
    }

    std::ostream& operator<<(std::ostream &out, const RegionTestsCsvRow &r) {
        out << r.reg_name << "," << r.reg_id << "," << r.day << "," << r.tests;
        return out;
    }

    void to_json(json& j, const MaterialMovement& m) {
        j = {{"from", m.origin_id}, {"to", m.destination_id}, {"day", m.day}, {"qty", m.qty}};
    }
    void to_json(json& j, const MaterialStorage& m) {
        j = {{"where", m.facility_id}, {"day", m.day}, {"qty", m.qty}};
    }
    void to_json(json& j, const TestsAction& t) {
        j = {{"where", t.lab_id}, {"day", t.day}, {"qty", t.qty}};
    }

    void TestLabAssignmentSolution::to_json() const {
        namespace fs = std::filesystem;

        uint32_t total_swabs_requested = std::accumulate(instance.reg_day_demand.begin(), instance.reg_day_demand.end(),
                0u, [] (uint32_t tot, const auto& row) {
                    return tot + std::accumulate(row.begin(), row.end(), 0u);
                });

        json j;
        j["untested_swabs_at_end"] = untested_swabs_at_end;
        j["total_waiting_swabs"] = total_waiting_swabs;
        j["avg_wait_time_tested_swabs"] =
                (double)(total_swabs_requested - untested_swabs_at_end) /
                (double)(total_waiting_swabs);
        j["reagent_fac_lab_movements"] = reagent_fac_lab_mov;
        j["reagent_lab_lab_movements"] = reagent_lab_lab_mov;
        j["swab_movements"] = swab_mov;
        j["reagents_stored_at_factories"] = reagent_fac_store;
        j["reagents_stored_at_labs"] = reagent_lab_store;
        j["swabs_stored_at_labs"] = swab_lab_store;
        j["tests_assignemtn_to_labs"] = tests_assigned;
        j["swabs_tested_at_labs"] = tests_performed;

        auto ip = fs::path{instance.inst_file};
        ip.replace_filename("sol-" + ip.filename().string());
        std::ofstream ofs{ip};

        ofs << j.dump(4) << "\n";
    }

    void TestLabAssignmentSolution::reg_tests_to_csv() const {
        namespace fs = std::filesystem;

        auto ip = fs::path{instance.inst_file};
        ip.replace_filename("reg-tests-" + ip.replace_extension("csv").filename().string());
        std::ofstream ofs{ip};

        ofs << "region,region_id,day_from_start,tested_swabs\n";
        for(const auto& row : reg_tests_csv) {
            ofs << row << "\n";
        }
    }
}