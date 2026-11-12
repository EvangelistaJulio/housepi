/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "data/household.h"
#include "data/reports.h"

//////////////////////////////  data::Dataset  //////////////////////////////

namespace housepi::data
{
    struct State;

    /// Dataset: one row per member, and households.
    struct Dataset : public KeyHolderWithRequests
    {
        /// Cache key holder, member values.
        KeyHolder keys_member_val;
        /// Cache key holder, member vectors.
        KeyHolder keys_member_vec;
        /// Cache key holder, double vectors.
        KeyHolder keys_double_vec;
        /// Cache key holder, double matrices.
        KeyHolder keys_double_mtx;

        ReportIndividual report_indiv;
        ReportHousehold report_house;

        //////  Constructor  //////

        /// Register model-declared keys: data columns `keys`, then cache keys per layout.
        Dataset(std::vector<std::string> keys,
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

        /// Build members and households from identifiers `ids` and data columns `cells`.
        void build(const str_table &ids, const dbl_table &cells);
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
        /// * column produced by an observation -> missing.
        double raw(size_t d, const Individual &indiv) const noexcept { return indiv_values_[Index::of(indiv) * size() + d]; }

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

        //////  Flags  //////

        /// Throw if every member has missing status.
        void assert_not_all_missing_status() const;
        /// Throw if some member has missing status.
        void assert_no_missing_status() const;

        /// Copy the requests registered on `other`, a store over the same keys.
        void copy_requests(const Dataset &other)
        {
            KeyHolderWithRequests::copy_requests(other);
            request_everyone_ = other.request_everyone_;
        }
        /// Request a value for every member, escaped ones included.
        void request_everyone(std::vector<size_t> vec);
        /// Whether key `d` was requested for every member.
        bool is_required_everyone(size_t d) const noexcept;
        /// Throw if a column requested for every member has no value for one.
        /// * values read through the read-cells of `state`;
        /// * `drawn` columns -> skipped, they take their value during the run.
        void assert_required_everyone(const State &state, const std::vector<std::string> &drawn = {}) const;

        /// Whether no data was provided.
        bool is_empty() const noexcept { return house_vector_.empty(); }

    private:
        //////  Data  //////

        /// Keys requested for every member.
        std::vector<size_t> request_everyone_ = {};
        /// Internal index per household user id.
        dict<size_t> hh_id_to_index_ = {};
        /// Imported status, by member index.
        std::vector<StatusLvl> indiv_status_ = {};
        /// Imported values, by member index then column.
        std::vector<double> indiv_values_ = {};
        /// Member user ids, by member index.
        std::vector<std::string> indiv_ids_ = {};
        /// Household user ids, by household index.
        std::vector<std::string> hh_ids_ = {};
        /// Number of columns known at import.
        size_t n_imported_ = 0;

        std::vector<Household> house_vector_ = {};

        //////  Helpers  //////

        /// Internal index of household with user id `hh_id`.
        size_t index_of_(std::string_view hh_id) const;
        /// Clear the previous data and report.
        void helper_build_clear_data_();
        /// Throw if the headers of `ids` are malformed, register those of `cells`.
        void helper_build_data_headers_(const str_table &ids, const dbl_table &cells);
        /// Fill the report from the loaded data.
        void build_helper_fill_report_data_();
    };
}
