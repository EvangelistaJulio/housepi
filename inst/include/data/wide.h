/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "data/household.h"
#include "data/reports.h"

//////////////////////////////  data::Wide  //////////////////////////////

namespace housepi::data
{
    struct State;

    /// Wide-format table: one row per member, households and their equivalence classes.
    struct Wide : public KeyHolderWithRequests
    {
        KeyHolder keys_member_val; /// Cache key holder, member values.
        KeyHolder keys_member_vec; /// Cache key holder, member vectors.
        KeyHolder keys_double_vec; /// Cache key holder, double vectors.
        KeyHolder keys_double_mtx; /// Cache key holder, double matrices.

        ReportIndividual report_indiv; /// Member report.
        ReportHousehold report_house;  /// Household report.

        //////  Constructor  //////

        /// Register model-declared keys: data columns `keys`, then cache keys per layout.
        Wide(std::vector<std::string> keys,
             std::vector<std::string> member_keys,
             std::vector<std::string> member_vec_keys,
             std::vector<std::string> double_vec_keys,
             std::vector<std::string> double_mtx_keys)
            : KeyHolderWithRequests(std::move(keys)),
              keys_member_val(std::move(member_keys)),
              keys_member_vec(std::move(member_vec_keys)),
              keys_double_vec(std::move(double_vec_keys)),
              keys_double_mtx(std::move(double_mtx_keys)) {}

        //////  Setter  //////

        /// Build members and households from wide-format table `input`.
        void build(const str_table &input);
        /// Recover every status.
        void recover_all_status() noexcept;
        /// Declare data columns `produced_cols`, produced by observations, missing for every member.
        void add_produced_cols(const std::vector<std::string> &produced_cols);

        //////  Iterators  //////

        auto begin() noexcept { return house_vector_.begin(); }
        auto begin() const noexcept { return house_vector_.begin(); }
        auto end() noexcept { return house_vector_.end(); }
        auto end() const noexcept { return house_vector_.end(); }

        //////  Getter  //////

        /// Household at internal position `H`, or with user id `hh_id`.
        Household &at(size_t H) noexcept { return house_vector_[H]; }
        const Household &at(size_t H) const noexcept { return house_vector_[H]; }
        Household &at(std::string_view hh_id) { return at(index_of_(hh_id)); }

        /// Imported value of column `d` for member `indiv`.
        /// - column produced by an observation: missing.
        double raw(size_t d, const Individual &indiv) const noexcept { return indiv_data_[Index::of(indiv)].vec[d]; }

        /// User id of member `indiv`.
        const std::string &id(const Individual &indiv) const noexcept { return indiv_ids_[Index::of(indiv)]; }
        /// User id of household `house`.
        const std::string &hh_id(const Household &house) const noexcept { return hh_ids_[Index::of(house)]; }
        /// Name of member `indiv`, for error messages.
        std::string str(const Individual &indiv) const { return "individual " + msg::val(id(indiv)); }
        /// Name of household `house`, for error messages.
        std::string str(const Household &house) const { return "household " + msg::val(hh_id(house)); }

        /// Column names known at import, those produced by an observation excluded.
        std::vector<std::string> imported_cols() const;

        //////  Equivalence classes  //////

        /// Declare data columns `data_list` as equivalence-class keys.
        void set_eq_class_data_names(const std::vector<std::string> &data_list);
        /// Group households into equivalence classes, weighted by their frequency.
        /// - equivalence class: households with the same likelihood, computed once;
        /// - `force_unique`: one class per household;
        /// - covariates read through the read-cells of `state`.
        void init_eq_classes(bool force_unique, const State &state);

        /// Number of households in equivalence class `R`.
        size_t eq_class_weight(size_t R) const noexcept { return eq_classes_[R].weight_; }
        /// Equivalence class of household `H`.
        size_t get_eq_class(size_t H) const noexcept { return house_vector_[H].eq_class_; }
        /// Number of equivalence classes.
        size_t eq_class_size() const noexcept { return report_house.distinct; }
        /// Head household of equivalence class `R`.
        Household &eq_class_head(size_t R) noexcept { return at(eq_classes_[R].head_); }
        const Household &eq_class_head(size_t R) const noexcept { return at(eq_classes_[R].head_); }
        /// Member of `house` standing for each member of its class head, by head member index.
        /// - equivalent members pair up in sorted order, so ties pair arbitrarily;
        /// - covariates read through the read-cells of `state`.
        std::vector<MemberIdx> map_head_members(const Household &house, const State &state) const;

        //////  Flags  //////

        /// Throw if every member has missing status.
        void assert_not_all_missing_status() const;
        /// Throw if some member has missing status.
        void assert_no_missing_status() const;

        /// Copy the requests registered on `other`, a store over the same keys.
        void copy_requests(const Wide &other)
        {
            KeyHolderWithRequests::copy_requests(other);
            request_everyone_ = other.request_everyone_;
        }
        /// Request a value for every member, escaped ones included.
        void request_everyone(std::vector<size_t> vec);
        /// Whether key `d` was requested for every member.
        bool is_required_everyone(size_t d) const noexcept;
        /// Throw if a column requested for every member has no value for one.
        /// - values read through the read-cells of `state`;
        /// - `live` columns: skipped, they take their value during the run.
        void assert_required_everyone(const State &state, const std::vector<std::string> &live = {}) const;
        /// Throw if boolean column `d` does not flag exactly one member per household.
        /// - values read through the read-cells of `state`.
        void assert_one_per_house(const State &state, size_t d) const;

        /// Whether no data was provided.
        bool is_empty() const noexcept { return house_vector_.empty(); }

    private:
        //////  Data  //////

        /// Wide-format data of one member.
        struct IndivWide_
        {
            StatusLvl status = status::MISSING; /// Imported status.
            std::vector<double> vec = {};       /// Imported values, by column.
        };

        std::vector<size_t> request_everyone_ = {}; /// Keys requested for every member.
        dict<size_t> hh_id_to_index_ = {};          /// Internal index per household user id.
        std::vector<IndivWide_> indiv_data_;        /// Member data, by member index.
        std::vector<std::string> indiv_ids_ = {};   /// Member user ids, by member index.
        std::vector<std::string> hh_ids_ = {};      /// Household user ids, by household index.
        size_t n_imported_ = 0;                     /// Number of columns known at import.

        std::vector<Household> house_vector_ = {}; /// Every household.

        //////  Equivalence class  //////

        /// Equivalence class: head household and number of households.
        struct EquivalenceClass_
        {
            EquivalenceClass_() = default;
            EquivalenceClass_(const std::vector<size_t> &group)
                : head_(group.front()),
                  weight_(group.size()) {}
            size_t head_ = na::INDEX; /// Head household.
            size_t weight_ = 1;       /// Number of households.
        };
        std::vector<EquivalenceClass_> eq_classes_ = {};      /// Equivalence classes.
        std::unordered_set<size_t> eq_class_covariates_ = {}; /// Covariates hashed into equivalence classes.

        //////  Helpers  //////

        /// Internal index of household with user id `hh_id`.
        size_t index_of_(std::string_view hh_id) const;
        /// Clear the previous data and report.
        void helper_build_clear_data_();
        /// Throw if the headers of `input` are malformed.
        void helper_build_data_headers_(const str_table &input);
        /// Fill the report from the loaded data.
        void build_helper_fill_report_data_();
        /// Initialize equivalence classes, one per household.
        void helper_init_unique_eq_classes_();
    };
}
