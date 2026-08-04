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
    namespace risks
    {
        constexpr const char *REFERENCE = ".reference"; /// Risk reference choice.
        constexpr double VALUE_REFERENCE = 1.0;         /// Alias for reference value.
    }

#define MAKE_HOUSEPI_CACHE_ACCESSOR(TYPE, NAME)                                                                        \
    TYPE &cache_##NAME(const Household &house, size_t e) noexcept { return cached_.current[house.eq_class_].NAME[e]; } \
    const TYPE &cache_##NAME(const Household &house, size_t e) const noexcept { return cached_.current[house.eq_class_].NAME[e]; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    struct DataWide : public KeyModule
    {
        KeyHolder khold_uint_val;   /// cache key-holder unsigned values
        KeyHolder khold_uint_vec;   /// cache key-holder unsigned vectors
        KeyHolder khold_double_vec; /// cache key-holder double vectors
        KeyHolder khold_double_mtx; /// cache key-holder double matrices

        ReportIndividual report_indiv;   /// individuals report
        ReportHousehold report_house;    /// households report
        ReportHouseSizes report_hh_size; /// house-sizes report

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign data (wide format).
        void build(const str_table &input);
        /// Recover all status.
        void recover_all_status() noexcept;
        /// Re-wire every read-cell back at its imported raw value.
        /// Optionally associate data columns produced by observations.
        /// NOTE: necessary because of stale pointers during Variable re-build.
        void reset_readcells(const std::vector<std::string> &produced_cols = {});

        /// Request data column to be discrete-valued.
        void request_discrete(std::vector<size_t> d) { user_discrete_set_.insert(d.begin(), d.end()); }
        /// Request data column to be strictly greater than zero.
        void request_positive(std::vector<size_t> d) { user_positive_set_.insert(d.begin(), d.end()); }

        /// Whether the modeler requested column 'd' to be discrete-valued.
        bool requires_discrete(size_t d) const { return user_discrete_set_.count(d) > 0; }
        /// Whether the modeler requested column 'd' to be strictly positive.
        bool requires_positive(size_t d) const { return user_positive_set_.count(d) > 0; }

        /// Transmission log-likelihood of head of household equivalence class.
        double house_transm_log_lik(size_t R) const { return transm_log_lik(eq_class_head(R)).value(); }

        /// Allocate per-equivalence-class transmission log-likelihood storage.
        /// NOTE: requires equivalence classes to be initialized (see init_eq_classes).
        void allocate_transm_log_lik() { house_transm_log_lik_.assign(eq_class_size(), LogProba()); }

        /// Access the transmission log-likelihood of a household equivalence class.
        LogProba &transm_log_lik(const Household &house) noexcept { return house_transm_log_lik_[house.eq_class_]; }
        const LogProba &transm_log_lik(const Household &house) const noexcept { return house_transm_log_lik_[house.eq_class_]; }

        /////////////////////////
        //////  ITERATORS  //////
        /////////////////////////

        auto begin() noexcept { return house_vector_.begin(); }
        auto begin() const noexcept { return house_vector_.begin(); }
        auto end() noexcept { return house_vector_.end(); }
        auto end() const noexcept { return house_vector_.end(); }

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Get reference of household by internal absolute indexing position H or by house_id.
        Household &at(size_t H) noexcept { return house_vector_[H]; }
        const Household &at(size_t H) const noexcept { return house_vector_[H]; }
        Household &at(std::string_view house_id) { return at(house_id_to_index_(house_id)); }
        const Household &at(std::string_view house_id) const { return at(house_id_to_index_(house_id)); }

        /// Access wide-format data of given column for an individual.
        /// NOTE: resolves through the per-member read-surface (raw or a wired Variable).
        double data(size_t d, const Individual &indiv) const noexcept { return readcell(d, indiv).value(); }
        /// Read-surface attribute of given column for an individual (for live @D resolution).
        const Attribute &readcell(size_t d, const Individual &indiv) const noexcept { return indiv_data_[Index::of(indiv)].readcell[d]; }
        /// Point a column's read-surface at the Variable of a surrogate column.
        void wire_readcell(size_t d, const Individual &indiv, const Attribute *p) noexcept { indiv_data_[Index::of(indiv)].readcell[d].init_pointer(p); }
        /// Finalize every Variable object in the read-cells.
        /// NOTE: evaluation-only; simulation re-wires the cells live.
        void finalize();

        /// Unique identifier of an individual, as provided in data.
        const std::string &id(const Individual &indiv) const noexcept { return indiv_ids_[Index::of(indiv)]; }
        /// Identifier of a household, as provided in data.
        const std::string &house_id(const Household &house) const noexcept { return house_ids_[Index::of(house)]; }
        /// Names for error messages.
        std::string err_name(const Individual &indiv) const { return "individual " + msg::val(id(indiv)); }
        std::string err_name(const Household &house) const { return "household " + msg::val(house_id(house)); }

        /// Column names known at import, without those produced by an observation.
        std::vector<std::string> imported_cols() const;

        ///////////////////////////////////
        //////  EQUIVALENCE CLASSES  //////
        ///////////////////////////////////

        // Associate data names as relevant for equivalence classes.
        void set_eq_class_data_names(const std::vector<std::string> &data_list);
        /// Initialize frequency weights for equivalence class compression.
        /// NOTE: An equivalence class is a unit of computational redundancy.
        void init_eq_classes(bool force_unique);

        /// Number of households in equivalence class R.
        size_t eq_class_weight(size_t R) const noexcept { return eq_classes_[R].weight_; }
        /// Index to household equivalence class head.
        size_t get_eq_class(size_t H) const noexcept { return house_vector_[H].eq_class_; }
        /// Restricted number of household equivalence classes.
        size_t eq_class_size() const noexcept { return report_house.distinct; }
        /// Head household of equivalence class R.
        Household &eq_class_head(size_t R) noexcept { return at(eq_classes_[R].head_); }
        const Household &eq_class_head(size_t R) const noexcept { return at(eq_classes_[R].head_); }

        /////////////////////
        //////  CACHE  //////
        /////////////////////

        /// Initialize cache key-holder with model declared keys.
        void init_cache_key_holders(
            const std::vector<std::string> &size_keys,
            const std::vector<std::string> &size_vec_keys,
            const std::vector<std::string> &double_vec_keys,
            const std::vector<std::string> &double_mtx_keys);

        /// Initialize cached structeres according to household size.
        void allocate_all_cache();
        /// Store all household cache attributes.
        void save_all_cached() { cached_.save(); }
        /// Recover old values of cached depedencies and attributes.
        void recover_old_cached() noexcept { cached_.recover(); }
        /// Save member cached attributes.
        void save_member_cache(const Household &house, const Individual &indiv);
        /// Recover member old values of cached member depedencies and member attributes.
        void recover_member_old_cached(const Household &house, const Individual &indiv);

        /// Return cached unsigned value.
        MAKE_HOUSEPI_CACHE_ACCESSOR(unsigned, uint)
        /// Return cached unsigned vector.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<unsigned>, uint_vector)
        /// Return cached double vector.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<double>, double_vector)
        /// Return cached double matrix.
        MAKE_HOUSEPI_CACHE_ACCESSOR(std::vector<std::vector<double>>, double_matrix)

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Assert that not all individuals have missing status.
        void assert_not_all_missing_status() const;
        /// Assert no individual has missing status.
        void assert_no_missing_status() const;

        /// Return maximum number of observed infected household members.
        size_t get_max_hh_cases() const noexcept { return report_hh_size.max_cases.first; }

        /// Wheter data has been provided.
        bool is_empty() const noexcept { return house_vector_.empty(); }

    private:
        ////////////////////
        //////  DATA  //////
        ////////////////////

        /// Individual wide-format data.
        struct IndivWide_
        {
            StatusLevel status = status::MISSING;
            std::vector<double> vec = {};         /// immutable imported-raw backing
            std::vector<Attribute> readcell = {}; /// per-column read-surface (raw or wired)
        };

        dict<size_t> hhid_to_index_ = {};         /// Map house-id to internal index
        std::vector<IndivWide_> indiv_data_;      /// External storage for individual data
        std::vector<std::string> indiv_ids_ = {}; /// User individual ids, by individual index
        std::vector<std::string> house_ids_ = {}; /// User house ids, by household index
        size_t n_imported_ = 0;                   /// Number of columns known at import

        std::vector<Household> house_vector_ = {};   /// All households in data
        std::vector<LogProba> house_transm_log_lik_; /// Per-equivalence-class transmission log-likelihood

        std::unordered_set<size_t> user_discrete_set_ = {}; /// modeler requested discrete data
        std::unordered_set<size_t> user_positive_set_ = {}; /// modeler requested positive data

        /////////////////////////////////
        //////  EQUIVALENCE CLASS  //////
        /////////////////////////////////

        struct EquivalenceClass_
        {
            EquivalenceClass_() = default;
            EquivalenceClass_(const std::vector<size_t> &group)
                : head_(group.front()), weight_(group.size()) {}
            size_t head_ = na::INDEX;
            size_t weight_ = 1;
        };
        std::vector<EquivalenceClass_> eq_classes_ = {};      /// Equivalence classes
        std::unordered_set<size_t> eq_class_covariates_ = {}; /// Covariates used for equivalence class hashing

        //////////////////////
        //////  CACHED  //////
        //////////////////////

        struct Cached_
        {
            struct State_
            {
                std::vector<unsigned> uint;
                std::vector<std::vector<unsigned>> uint_vector = {};
                std::vector<std::vector<double>> double_vector = {};
                std::vector<std::vector<std::vector<double>>> double_matrix = {};
            };
            std::vector<State_> previous, current;

            void save() { previous = current; }
            void recover() noexcept { current.swap(previous); }
        } cached_;

        ///////////////////////
        //////  HELPERS  //////
        ///////////////////////

        /// Get housepi household index from house_id.
        size_t house_id_to_index_(std::string_view house_id) const;
        /// Clear previous data report.
        void helper_build_clear_data_();
        /// Assert that data headers are correct.
        void helper_build_data_headers_(const str_table &input);
        /// Use the loaded data to fill-in the respective report.
        void build_helper_fill_report_data_();
        /// Initialize equivalence classes as if all distinct.
        void helper_init_unique_eq_classes_();
    };
}
