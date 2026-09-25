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
#include <vector>

#include "tools/data_reports.h"
#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/types.h"

#include "function.h"
#include "household.h"

namespace housepi
{
#define MAKE_HOUSEPI_CACHE_ACCESSOR(TYPE, NAME)                                                                        \
    TYPE &cache_##NAME(const Household &house, size_t e) noexcept { return cached_.current[house.eq_class_].NAME[e]; } \
    const TYPE &cache_##NAME(const Household &house, size_t e) const noexcept { return cached_.current[house.eq_class_].NAME[e]; }

    //////////////////////////////  DataWide  //////////////////////////////

    /// Wide-format data: one row per member, households and their equivalence classes.
    struct DataWide : public KeyModuleWithRequests
    {
        KeyHolder khold_member_val; /// Cache key holder, member values.
        KeyHolder khold_member_vec; /// Cache key holder, member vectors.
        KeyHolder khold_double_vec; /// Cache key holder, double vectors.
        KeyHolder khold_double_mtx; /// Cache key holder, double matrices.

        ReportIndividual report_indiv;   /// Member report.
        ReportHousehold report_house;    /// Household report.
        ReportHouseSizes report_hh_size; /// Household-size report.

        //////  Setter  //////

        /// Build members and households from wide-format table `input`.
        void build(const str_table &input);
        /// Recover every status.
        void recover_all_status() noexcept;
        /// Point every read-cell back at its imported raw value.
        /// - `produced_cols`: data columns produced by observations, declared alongside.
        /// @note Call before rebuilding Variables: wired cells would point at stale ones.
        void reset_readcells(const std::vector<std::string> &produced_cols = {});

        /// Transmission log-likelihood of the head of equivalence class `R`.
        double house_transm_log_lik(size_t R) const { return transm_log_lik(eq_class_head(R)).value(); }

        /// Allocate transmission log-likelihood storage, one per equivalence class.
        /// @note Requires equivalence classes (see init_eq_classes).
        void allocate_transm_log_lik() { house_transm_log_lik_.assign(eq_class_size(), LogProba()); }

        /// Transmission log-likelihood of the equivalence class of `house`.
        LogProba &transm_log_lik(const Household &house) noexcept { return house_transm_log_lik_[house.eq_class_]; }
        const LogProba &transm_log_lik(const Household &house) const noexcept { return house_transm_log_lik_[house.eq_class_]; }

        //////  Iterators  //////

        auto begin() noexcept { return house_vector_.begin(); }
        auto begin() const noexcept { return house_vector_.begin(); }
        auto end() noexcept { return house_vector_.end(); }
        auto end() const noexcept { return house_vector_.end(); }

        //////  Getter  //////

        /// Household at internal position `H`, or with user id `hhid`.
        Household &at(size_t H) noexcept { return house_vector_[H]; }
        const Household &at(size_t H) const noexcept { return house_vector_[H]; }
        Household &at(std::string_view hhid) { return at(index_of_(hhid)); }
        const Household &at(std::string_view hhid) const { return at(index_of_(hhid)); }

        /// Datum `d` of member `indiv`.
        /// - read through the member read-cell: raw value or wired Variable.
        double data(size_t d, const Individual &indiv) const noexcept { return readcell(d, indiv).value(); }
        /// Read-cell `d` of member `indiv`, resolved live by `@D`.
        const Attribute &readcell(size_t d, const Individual &indiv) const noexcept { return indiv_data_[Index::of(indiv)].readcell[d]; }
        /// Point read-cell `d` of member `indiv` at the Variable of a surrogate column.
        void wire_readcell(size_t d, const Individual &indiv, const Attribute *p) noexcept { indiv_data_[Index::of(indiv)].readcell[d].init_pointer(p); }
        /// Finalize every Variable wired into a read-cell.
        /// - evaluation only: simulation re-wires the cells live.
        void finalize();

        /// User id of member `indiv`.
        const std::string &id(const Individual &indiv) const noexcept { return indiv_ids_[Index::of(indiv)]; }
        /// User id of household `house`.
        const std::string &hhid(const Household &house) const noexcept { return hhids_[Index::of(house)]; }
        /// Name of member `indiv`, for error messages.
        std::string str(const Individual &indiv) const { return "individual " + msg::val(id(indiv)); }
        /// Name of household `house`, for error messages.
        std::string str(const Household &house) const { return "household " + msg::val(hhid(house)); }

        /// Column names known at import, those produced by an observation excluded.
        std::vector<std::string> imported_cols() const;

        //////  Equivalence classes  //////

        /// Declare data columns `data_list` as equivalence-class keys.
        void set_eq_class_data_names(const std::vector<std::string> &data_list);
        /// Group households into equivalence classes, weighted by their frequency.
        /// - equivalence class: households with the same likelihood, computed once.
        /// - `force_unique`: one class per household.
        void init_eq_classes(bool force_unique);

        /// Number of households in equivalence class `R`.
        size_t eq_class_weight(size_t R) const noexcept { return eq_classes_[R].weight_; }
        /// Equivalence class of household `H`.
        size_t get_eq_class(size_t H) const noexcept { return house_vector_[H].eq_class_; }
        /// Number of equivalence classes.
        size_t eq_class_size() const noexcept { return report_house.distinct; }
        /// Head household of equivalence class `R`.
        Household &eq_class_head(size_t R) noexcept { return at(eq_classes_[R].head_); }
        const Household &eq_class_head(size_t R) const noexcept { return at(eq_classes_[R].head_); }

        //////  Cache  //////

        /// Initialize cache key holders with model-declared keys.
        void init_cache_key_holders(
            const std::vector<std::string> &member_keys,
            const std::vector<std::string> &member_vec_keys,
            const std::vector<std::string> &double_vec_keys,
            const std::vector<std::string> &double_mtx_keys);

        /// Allocate cached structures according to household size.
        void allocate_all_cache();
        /// Save every household cache.
        void save_all_cached() { cached_.save(); }
        /// Recover every household cache from the last save.
        void recover_old_cached() noexcept { cached_.recover(); }
        /// Save the cache of member `indiv`.
        void save_member_cache(const Household &house, const Individual &indiv);
        /// Recover the cache of member `indiv` from the last save.
        void recover_member_old_cached(const Household &house, const Individual &indiv);

        /// Cached member value.
        MAKE_HOUSEPI_CACHE_ACCESSOR(MemberIdx, member)
        /// Cached member vector.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<MemberIdx>, member_vector)
        /// Cached double vector.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<double>, double_vector)
        /// Cached double matrix.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<std::vector<double>>, double_matrix)

        //////  Flags  //////

        /// Throw if every member has missing status.
        void assert_not_all_missing_status() const;
        /// Throw if some member has missing status.
        void assert_no_missing_status() const;

        /// Copy the requests registered on `other`, a store over the same keys.
        void copy_requests(const DataWide &other)
        {
            KeyModuleWithRequests::copy_requests(other);
            request_everyone_ = other.request_everyone_;
        }
        /// Request a value for every member, escaped ones included.
        void request_everyone(std::vector<size_t> vec);
        /// Whether key `d` was requested for every member.
        bool is_required_everyone(size_t d) const noexcept;
        /// Throw if a column requested for every member has no value for one.
        /// - `live` columns: skipped, they take their value during the run.
        void assert_required_everyone(const std::vector<std::string> &live = {}) const;

        /// Maximum number of infected members observed in one household.
        size_t get_max_hh_cases() const noexcept { return report_hh_size.max_cases.first; }

        /// Whether no data was provided.
        bool is_empty() const noexcept { return house_vector_.empty(); }

    private:
        //////  Data  //////

        /// Wide-format data of one member.
        struct IndivWide_
        {
            StatusLvl status = status::MISSING;   /// Imported status.
            std::vector<double> vec = {};         /// Imported raw values, immutable.
            std::vector<Attribute> readcell = {}; /// Read-cell per column: raw value or wired Variable.
        };

        std::vector<size_t> request_everyone_ = {}; /// Keys requested for every member.
        dict<size_t> hhid_to_index_ = {};           /// Internal index per household user id.
        std::vector<IndivWide_> indiv_data_;        /// Member data, by member index.
        std::vector<std::string> indiv_ids_ = {};   /// Member user ids, by member index.
        std::vector<std::string> hhids_ = {};       /// Household user ids, by household index.
        size_t n_imported_ = 0;                     /// Number of columns known at import.

        std::vector<Household> house_vector_ = {};   /// Every household.
        std::vector<LogProba> house_transm_log_lik_; /// Transmission log-likelihood per equivalence class.

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

        //////  Cached  //////

        /// Model cache per equivalence class, current and last saved.
        struct Cached_
        {
            struct State_
            {
                std::vector<MemberIdx> member = {};
                std::vector<std::vector<MemberIdx>> member_vector = {};
                std::vector<std::vector<double>> double_vector = {};
                std::vector<std::vector<std::vector<double>>> double_matrix = {};
            };
            std::vector<State_> previous, current;

            void save() { previous = current; }
            void recover() noexcept { current.swap(previous); }
        } cached_;

        //////  Helpers  //////

        /// Internal index of household with user id `hhid`.
        size_t index_of_(std::string_view hhid) const;
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

#undef MAKE_HOUSEPI_CACHE_ACCESSOR
