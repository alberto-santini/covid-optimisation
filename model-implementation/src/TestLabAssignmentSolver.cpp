//
// Created by alberto on 12/04/2020.
//

#include "TestLabAssignmentSolver.h"
#include <gurobi_c++.h>
#include <cmath>

namespace covid {
    template<typename T> using vec = std::vector<T>;
    template<typename T> using vec2 = vec<vec<T>>;
    template<typename T> using vec3 = vec<vec2<T>>;

    namespace {
        inline vec2<GRBVar> mk_vec2(size_t d1, size_t d2) {
            return vec2<GRBVar>(d1, vec<GRBVar>(d2));
        }
        inline vec3<GRBVar> mk_vec3(size_t d1, size_t d2, size_t d3) {
            return vec3<GRBVar>(d1, mk_vec2(d2, d3));
        }
    }

    std::optional<TestLabAssignmentSolution> TestLabAssignmentSolver::solve_moo_model(bool transshipments) const {
        GRBEnv env{};
        GRBModel model{env};

        // Variables:

        auto fac_lab_day_reagent_mov = mk_vec3(i.n_factories, i.n_labs, i.n_days); // x
        auto lab_lab_day_reagent_mov = mk_vec3(i.n_labs, i.n_labs, i.n_days);      // x
        auto lab_lab_day_swab_mov = mk_vec3(i.n_labs, i.n_labs, i.n_days);         // y
        auto lab_day_swab_store = mk_vec2(i.n_labs, i.n_days);                     // z
        auto lab_day_swab_test = mk_vec2(i.n_labs, i.n_days);                      // w
        auto lab_day_swab_assign = mk_vec2(i.n_labs, i.n_days);                    // u
        auto fac_day_reagent_store = mk_vec2(i.n_factories, i.n_days);             // rho
        auto lab_day_reagent_store = mk_vec2(i.n_labs, i.n_days);                  // rho
        auto lab_day_at_capacity = mk_vec2(i.n_labs, i.n_days);                    // gamma
        auto lab_day_swab_receive = mk_vec2(i.n_labs, i.n_days);                   // gamma^-
        auto lab_day_swab_send = mk_vec2(i.n_labs, i.n_days);                      // gamma^+

        auto fac_lab_day_reagent_mov_bool = mk_vec3(i.n_factories, i.n_labs, i.n_days); // 0 if x==0, 1 if x>0
        auto lab_lab_day_reagent_mov_bool = mk_vec3(i.n_labs, i.n_labs, i.n_days);      // 0 if x==0, 1 if x>0
        auto lab_lab_day_swab_mov_bool = mk_vec3(i.n_labs, i.n_labs, i.n_days);         // 0 if y==0, 1 if y>0

        for(auto r = 0u; r < i.n_factories; ++r) {
            for(auto l = 0u; l < i.n_labs; ++l) {
                for(auto t = 0u; t < i.n_days; ++t) {
                    fac_lab_day_reagent_mov[r][l][t] = model.addVar(0.0,
                            (i.fac_lab_compatible[r][l] ? i.fac_day_production[r][t] : 0u),
                            0.0, GRB_INTEGER);
                    fac_lab_day_reagent_mov_bool[r][l][t] = model.addVar(0.0,
                            (i.fac_lab_compatible[r][l] ? 1u : 0u), 0.0, GRB_BINARY);
                }
            }
            for(auto t = 0u; t < i.n_days; ++t) {
                fac_day_reagent_store[r][t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
            }
        }

        for(auto l = 0u; l < i.n_labs; ++l) {
            for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                for(auto t = 0u; t < i.n_days; ++t) {
                    if(transshipments) {
                        lab_lab_day_reagent_mov[l][l2][t] = model.addVar(0.0,
                                (i.lab_lab_compatible[l][l2] ? i.reg_max_inbound_swabs[i.lab_region[l2]] : 0u),
                                0.0, GRB_INTEGER);

                        lab_lab_day_reagent_mov_bool[l][l2][t] = model.addVar(0.0,
                                (i.lab_lab_compatible[l][l2] ? 1u : 0u), 0.0, GRB_BINARY);
                    }

                    lab_lab_day_swab_mov[l][l2][t] = model.addVar(0.0,
                            (i.lab_lab_compatible[l][l2] ? i.reg_max_inbound_swabs[i.lab_region[l2]] : 0u),
                            0.0, GRB_INTEGER);

                    lab_lab_day_swab_mov_bool[l][l2][t] = model.addVar(0.0,
                            (i.lab_lab_compatible[l][l2] ? 1u : 0u), 0.0, GRB_BINARY);
                }
            }

            for(auto t = 0u; t < i.n_days; ++t) {
                lab_day_swab_store[l][t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
                lab_day_swab_test[l][t] = model.addVar(0.0, i.lab_capacity[l], 0.0, GRB_INTEGER);
                lab_day_swab_assign[l][t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
                lab_day_reagent_store[l][t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER);
                lab_day_at_capacity[l][t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                lab_day_swab_send[l][t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                lab_day_swab_receive[l][t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        // Constraints:

        for(auto j = 0u; j < i.n_regions; ++j) {
            for(auto t = 0u; t < i.n_days; ++t) {
                {
                    GRBLinExpr expr = 0;

                    for(auto r = 0u; r < i.n_factories; ++r) {
                        for(auto l : i.reg_labs[j]) {
                            expr += fac_lab_day_reagent_mov[r][l][t];
                        }
                    }

                    if(transshipments) {
                        for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                            for(auto l : i.reg_labs[j]) {
                                if(i.lab_lab_compatible[l][l2]) {
                                    expr += lab_lab_day_reagent_mov[l2][l][t];
                                }
                            }
                        }
                    }

                    model.addConstr(expr <= i.reg_max_inbound_reagents[j]);
                }

                {
                    GRBLinExpr expr = 0;

                    for(auto l = 0u; l < i.n_labs; ++l) {
                        for(auto l2 : i.reg_labs[j]) {
                            expr += lab_lab_day_swab_mov[l][l2][t];
                        }
                    }

                    model.addConstr(expr <= i.reg_max_inbound_swabs[j]);
                }

                {
                    GRBLinExpr expr = 0;

                    for(auto l : i.reg_labs[j]) {
                        expr += lab_day_swab_assign[l][t];
                    }

                    model.addConstr(expr == i.reg_day_demand[j][t]);
                }
            }
        }

        for(auto l = 0u; l < i.n_labs; ++l) {
            // Case t = 0 is special:
            {
                GRBLinExpr lhs = lab_day_swab_assign[l][0u];
                GRBLinExpr rhs = lab_day_swab_store[l][0u] + lab_day_swab_test[l][0u];

                for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                    if(i.lab_lab_compatible[l2][l]) {
                        lhs += lab_lab_day_swab_mov[l2][l][0u];
                    }
                    if(i.lab_lab_compatible[l][l2]) {
                        rhs += lab_lab_day_swab_mov[l][l2][0u];
                    }
                }

                model.addConstr(lhs == rhs);
            }

            // Case t = 0 is special:
            {
                GRBLinExpr lhs = i.lab_start_reagents[l];
                GRBLinExpr rhs = lab_day_reagent_store[l][0u] + lab_day_swab_test[l][0u];

                for(auto r = 0u; r < i.n_factories; ++r) {
                    if(i.fac_lab_compatible[r][l]) {
                        lhs += fac_lab_day_reagent_mov[r][l][0u];
                    }
                }

                if(transshipments) {
                    for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                        if(i.lab_lab_compatible[l][l2]) {
                            rhs += lab_lab_day_reagent_mov[l][l2][0u];
                        }
                    }
                }

                model.addConstr(lhs == rhs);
                model.addConstr(lab_day_swab_test[l][0u] <= lhs); // TMP
            }

            for(auto t = 1u; t < i.n_days; ++t) {
                {
                    GRBLinExpr lhs = lab_day_swab_store[l][t - 1u] + lab_day_swab_assign[l][t];
                    GRBLinExpr rhs = lab_day_swab_store[l][t] + lab_day_swab_test[l][t];

                    for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                        if(i.lab_lab_compatible[l2][l]) {
                            lhs += lab_lab_day_swab_mov[l2][l][t];
                        }
                        if(i.lab_lab_compatible[l][l2]) {
                            rhs += lab_lab_day_swab_mov[l][l2][t];
                        }
                    }

                    model.addConstr(lhs == rhs);
                }

                {
                    GRBLinExpr lhs = lab_day_reagent_store[l][t - 1u];
                    GRBLinExpr rhs = lab_day_reagent_store[l][t] + lab_day_swab_test[l][t];

                    for(auto r = 0u; r < i.n_factories; ++r) {
                        if(i.fac_lab_compatible[r][l]) {
                            lhs += fac_lab_day_reagent_mov[r][l][t];
                        }
                    }

                    if(transshipments) {
                        for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                            if(i.lab_lab_compatible[l][l2]) {
                                lhs += lab_lab_day_reagent_mov[l2][l][t - 1u];
                                rhs += lab_lab_day_reagent_mov[l][l2][t];
                            }
                        }
                    }

                    model.addConstr(lhs == rhs);
                    model.addConstr(lab_day_swab_test[l][t] <= lhs); // TMP
                }
            }
                         
            for(auto t = 0u; t < i.n_days; ++t) {
                model.addGenConstrIndicator(
                        lab_day_at_capacity[l][t], false, lab_day_swab_test[l][t] <= i.lab_capacity[l] - 1u);

                model.addGenConstrIndicator(
                        lab_day_at_capacity[l][t], false, lab_day_reagent_store[l][t] >= 1);

                {
                    GRBLinExpr expr = 0;

                    for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                        if(i.lab_lab_compatible[l][l2]) {
                            expr += lab_lab_day_swab_mov[l][l2][t];
                        }
                    }

                    model.addGenConstrIndicator(
                            lab_day_at_capacity[l][t], false, expr == 0);

                    model.addGenConstrIndicator(
                            lab_day_swab_send[l][t], true, expr >= 1);

                    model.addGenConstrIndicator(
                            lab_day_swab_send[l][t], false, expr == 0);
                }

                {
                    GRBLinExpr expr = 0;

                    for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                        if(i.lab_lab_compatible[l2][l]) {
                            expr += lab_lab_day_swab_mov[l2][l][t];
                        }
                    }

                    model.addGenConstrIndicator(
                            lab_day_swab_receive[l][t], true, expr >= 1);

                    model.addGenConstrIndicator(
                            lab_day_swab_receive[l][t], false, expr == 0);
                }

                model.addConstr(lab_day_swab_send[l][t] + lab_day_swab_receive[l][t] <= 1);
            }
        }

        for(auto r = 0u; r < i.n_factories; ++r) {
            // Case t = 0 is special:
            {
                GRBLinExpr lhs = i.fac_start_reagents[r] + i.fac_day_production[r][0u];
                GRBLinExpr rhs = fac_day_reagent_store[r][0u];

                for(auto l = 0u; l < i.n_labs; ++l) {
                    if(i.fac_lab_compatible[r][l]) {
                        rhs += fac_lab_day_reagent_mov[r][l][0u];
                    }
                }

                model.addConstr(lhs == rhs);
            }

            for(auto t = 1u; t < i.n_days; ++t) {
                GRBLinExpr lhs = fac_day_reagent_store[r][t - 1u] + i.fac_day_production[r][t];
                GRBLinExpr rhs = fac_day_reagent_store[r][t];

                for(auto l = 0u; l < i.n_labs; ++l) {
                    if(i.fac_lab_compatible[r][l]) {
                        rhs += fac_lab_day_reagent_mov[r][l][t];
                    }
                }

                model.addConstr(lhs == rhs);
            }
        }

        // Setting boolean variables:
        for(auto t = 0u; t < i.n_days; ++t) {
            for(auto l = 0u; l < i.n_labs; ++l) {
                for(auto r = 0u; r < i.n_factories; ++r) {
                    const auto M = std::min(i.fac_day_production[r][t], i.reg_max_inbound_reagents[i.lab_region[l]]);

                    model.addConstr(
                        fac_lab_day_reagent_mov[r][l][t] <=
                        M * fac_lab_day_reagent_mov_bool[r][l][t]);
                }

                for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                    model.addConstr(
                        lab_lab_day_swab_mov[l][l2][t] <=
                        i.reg_max_inbound_swabs[i.lab_region[l2]] * lab_lab_day_swab_mov_bool[l][l2][t]);

                    if(transshipments) {
                        model.addConstr(
                            lab_lab_day_reagent_mov[l][l2][t] <=
                            i.reg_max_inbound_reagents[i.lab_region[l2]] * lab_lab_day_reagent_mov_bool[l][l2][t]);
                    }
                }
            }
        }

        // Constraints with boolean variables:

        for(auto t = 0u; t < i.n_days; ++t) {
            for(auto l1 = 0u; l1 < i.n_labs; ++l1) {
                for(auto l2 = l1 + 1u; l2 < i.n_labs; ++l2) {
                    for(auto r1 = 0u; r1 < i.n_factories; ++r1) {
                        if(!i.fac_lab_compatible[r1][l2]) { continue; }

                        for(auto r2 = r1 + 1u; r2 < i.n_factories; ++r2) {
                            if(!i.fac_lab_compatible[r2][l1]) { continue; }

                            if(i.fac_lab_distance[r1][l1] < i.fac_lab_distance[r1][l2] &&
                               i.fac_lab_distance[r2][l2] < i.fac_lab_distance[r2][l1])
                            {
                                model.addConstr(
                                    fac_lab_day_reagent_mov_bool[r2][l1][t] +
                                    fac_lab_day_reagent_mov_bool[r1][l2][t] <= 1);
                            }
                        }
                    }

                    for(auto ll1 = l2 + 1u; ll1 < i.n_labs; ++ll1) {
                        if(!i.lab_lab_compatible[ll1][l2]) { continue; }

                        for(auto ll2 = ll1 + 1u; ll2 < i.n_labs; ++ll2) {
                            if(!i.lab_lab_compatible[ll2][l1]) { continue; }

                            if(i.lab_lab_distance[l1][ll1] < i.lab_lab_distance[l2][ll1] &&
                               i.lab_lab_distance[l2][ll2] < i.lab_lab_distance[l1][ll2])
                            {
                                if(transshipments) {
                                    model.addConstr(
                                        lab_lab_day_reagent_mov_bool[ll2][l1][t] +
                                        lab_lab_day_reagent_mov_bool[ll1][l2][t] <= 1);

                                    model.addConstr(
                                        lab_lab_day_swab_mov_bool[ll2][l1][t] +
                                        lab_lab_day_swab_mov_bool[ll1][l2][t] <= 1);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Objective functions:

        // All objectives must have the same sense.
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        { // Primary objective
            GRBLinExpr expr = 0;

            for(auto l = 0u; l < i.n_labs; ++l) {
                expr += lab_day_swab_store[l][i.n_days - 1u];
            }

            model.setObjectiveN(expr, 0, 100, 1.0); // Priority: 100
        }

        { // Secondary objective
            GRBLinExpr expr = 0;

            for(auto l = 0u; l < i.n_labs; ++l) {
                for(auto t = 0u; t < i.n_days - 1u; ++t) {
                    expr += lab_day_swab_store[l][t];
                }
            }

            model.setObjectiveN(expr, 1, 50, 1.0); // Priority: 50
        }

        // Solve:

        // Time limit: separate for each objective
        auto primary_obj_env = model.getMultiobjEnv(0);
        primary_obj_env.set(GRB_DoubleParam_TimeLimit, 900.0); // 15 Minutes
        auto secondary_obj_env = model.getMultiobjEnv(1);
        secondary_obj_env.set(GRB_DoubleParam_TimeLimit, 300.0); // 5 Minutes

        // Gap tolerance: 0.1%
        model.set(GRB_DoubleParam_MIPGap, 0.001);
        
        model.optimize();

        const auto status = model.get(GRB_IntAttr_Status);

        if(status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED) {
            std::cerr << "Model infeasible or unbounded!\n";
            return std::nullopt;
        }

        if(status != GRB_OPTIMAL) {
            std::cerr << "Warning: Solution not optimal!\n";
        }

        // Extract solution:

        auto s = TestLabAssignmentSolution{i};

        for(auto r = 0u; r < i.n_factories; ++r) {
            for(auto l = 0u; l < i.n_labs; ++l) {
                if(i.fac_lab_compatible[r][l]) {
                    for(auto t = 0u; t < i.n_days; ++t) {
                        const uint32_t qty = std::round(fac_lab_day_reagent_mov[r][l][t].get(GRB_DoubleAttr_X));
                        if(qty > 0u) {
                            s.reagent_fac_lab_mov.push_back({Material::REAGENT, r, l, t, qty});
                        }
                    }
                }
            }

            for(auto t = 0u; t < i.n_days; ++t) {
                const uint32_t qty = std::round(fac_day_reagent_store[r][t].get(GRB_DoubleAttr_X));
                if(qty > 0u) {
                    s.reagent_fac_store.push_back({Material::REAGENT, r, t, qty});
                }
            }
        }

        for(auto l = 0u; l < i.n_labs; ++l) {
            for(auto l2 = 0u; l2 < i.n_labs; ++l2) {
                if(i.lab_lab_compatible[l][l2]) {
                    for(auto t = 0u; t < i.n_days; ++t) {
                        {
                            const uint32_t qty = std::round(lab_lab_day_swab_mov[l][l2][t].get(GRB_DoubleAttr_X));
                            if(qty > 0u) {
                                s.swab_mov.push_back({Material::SWAB, l, l2, t, qty});
                            }
                        }
                        
                        if(transshipments) {
                            const uint32_t qty = std::round(lab_lab_day_reagent_mov[l][l2][t].get(GRB_DoubleAttr_X));
                            if(qty > 0u) {
                                s.reagent_lab_lab_mov.push_back({Material::REAGENT, l, l2, t, qty});
                            }
                        }
                    }
                }
            }

            for(auto t = 0u; t < i.n_days; ++t) {
                {
                    const uint32_t qty = std::round(lab_day_reagent_store[l][t].get(GRB_DoubleAttr_X));
                    if(qty > 0u) {
                        s.reagent_lab_store.push_back({Material::REAGENT, l, t, qty});
                    }
                }
                {
                    const uint32_t qty = std::round(lab_day_swab_store[l][t].get(GRB_DoubleAttr_X));
                    if(qty > 0u) {
                        s.swab_lab_store.push_back({Material::SWAB, l, t, qty});

                        if(t < i.n_days - 1u) {
                            s.total_waiting_swabs += qty;
                        } else {
                            s.untested_swabs_at_end += qty;
                        }
                    }
                }
                {
                    const uint32_t qty = std::round(lab_day_swab_assign[l][t].get(GRB_DoubleAttr_X));
                    if(qty > 0u) {
                        s.tests_assigned.push_back({l, t, qty});
                    }
                }
                {
                    const uint32_t qty = std::round(lab_day_swab_test[l][t].get(GRB_DoubleAttr_X));
                    if(qty > 0u) {
                        s.tests_performed.push_back({l, t, qty});
                    }
                }
            }
        }

        for(auto j = 0u; j < i.n_regions; ++j) {
            for(auto t = 0u; t < i.n_days; ++t) {
                uint32_t reg_tests = 0u;

                for(auto l : i.reg_labs[j]) {
                    reg_tests += std::round(lab_day_swab_test[l][t].get(GRB_DoubleAttr_X));
                }

                s.reg_tests_csv.push_back({i.reg_names[j], j, t, reg_tests});
            }
        }

        return s;
    }
}
