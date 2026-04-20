/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
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
#include "tools/types.h"

#include "indiv_storage.h"
#include "household.h"

namespace housepi
{
    namespace rates
    {
        constexpr const char *REFERENCE = ".reference"; /// Rate reference choice.
        constexpr double VALUE_REFERENCE = 1.0;         /// Alias for reference value.
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    struct DataWide : public KeyModule
    {
        KeyHolder khold_size_val;   /// cache key-holder size_t values
        KeyHolder khold_size_vec;   /// cache key-holder size_t vectors
        KeyHolder khold_double_vec; /// cache key-holder double vectors
        KeyHolder khold_double_mtx; /// cache key-holder double matrices

        ReportIndividual report_indiv;   /// individuals report
        ReportHousehold report_house;    /// households report
        ReportHouseSizes report_hh_size; /// house-sizes report

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Initialize cache key-holder with model declared keys.
        void init_cache_key_holders(
            const std::vector<std::string> &size_keys,
            const std::vector<std::string> &size_vec_keys,
            const std::vector<std::string> &double_vec_keys,
            const std::vector<std::string> &double_mtx_keys);

        /// Assign data (wide format).
        void build(const str_table &input);

        /// Set covariates as frequency-relevant by providing data_wide column names.
        void set_frequency_covariates(const std::vector<std::string> &data_list);
        /// Initialize frequency weights for equivalence classes from household data.
        void init_frequencies(bool force_unique);

        /// Initialize cached structeres according to household size.
        void allocate_all_cache();

        /// Store all household cache attributes.
        void save_all_cached() { cached_.save(); }
        /// Recover old values of cached depedencies and attributes.
        void recover_old_cached() noexcept { cached_.recover(); }

        /// Save member cached attributes.
        void save_member_cache(Household &house, size_t i);
        /// Recover member old values of cached member depedencies and member attributes.
        void recover_member_old_cached(Household &house, size_t i);

        /// Raise error log-likelihood is not finite.
        void assert_finite_log_lik(bool is_any_augm) const;

        /// Compute log-likelihood for distinct_head household R (obs + augmented).
        double house_log_lik(size_t R, bool is_any_augm) const;

        /// Recover all status.
        void recover_all_status() noexcept;

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Get reference of household by internal absolute indexing position H or by house_id.
        Household &at(size_t H) { return house_vector_[H]; }
        const Household &at(size_t H) const { return house_vector_[H]; }
        Household &at(std::string_view house_id) { return at(house_id_to_index_(house_id)); }
        const Household &at(std::string_view house_id) const { return at(house_id_to_index_(house_id)); }

        /// Frequency weight of households in the distinc group.
        size_t frequency_weight(size_t R) const { return frequencies_[R].weight_; }
        /// Index to household frequency head.
        size_t get_restricted_index(size_t H) const { return house_vector_[H].ridx_; }
        /// Restricted number of household frequency classes.
        size_t distinct_size() const noexcept { return report_house.distinct; }
        /// Restrict to head of distinct household frequency class.
        Household &distinct_head(size_t R) { return at(frequencies_[R].head_); }
        const Household &distinct_head(size_t R) const { return at(frequencies_[R].head_); }

        /////////////////////
        //////  CACHE  //////
        /////////////////////

        /// Return cached size_t value.
        size_t &cache_size(Household &house, size_t e) { return cached_.current[house.ridx_].size_val[e]; }
        size_t cache_size(const Household &house, size_t e) const { return cached_.current[house.ridx_].size_val[e]; }

        /// Return cached size_t vector.
        std::vector<size_t> &cache_size_vector(Household &house, size_t e) { return cached_.current[house.ridx_].size_vec[e]; }
        const std::vector<size_t> &cache_size_vector(const Household &house, size_t e) const { return cached_.current[house.ridx_].size_vec[e]; }

        /// Return cached double vector.
        std::vector<double> &cache_double_vector(Household &house, size_t e) { return cached_.current[house.ridx_].dbl_vec[e]; }
        const std::vector<double> &cache_double_vector(const Household &house, size_t e) const { return cached_.current[house.ridx_].dbl_vec[e]; }

        /// Return cached double matrix.
        std::vector<std::vector<double>> &cache_double_matrix(Household &house, size_t e) { return cached_.current[house.ridx_].dbl_mtx[e]; }
        const std::vector<std::vector<double>> &cache_double_matrix(const Household &house, size_t e) const { return cached_.current[house.ridx_].dbl_mtx[e]; }

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

        /////////////////////////
        //////  ITERATORS  //////
        /////////////////////////

        auto begin() noexcept { return house_vector_.begin(); }
        auto begin() const noexcept { return house_vector_.begin(); }
        auto end() noexcept { return house_vector_.end(); }
        auto end() const noexcept { return house_vector_.end(); }

    private:
        ////////////////////
        //////  DATA  //////
        ////////////////////

        std::vector<Indiv::Wide> indiv_data_;      /// External storage for individual data
        std::vector<Household> house_vector_ = {}; /// All households in data
        dict<size_t> hhid_to_index_ = {};          /// Map house-id to internal index

        /////////////////////////
        //////  FREQUENCY  //////
        /////////////////////////

        struct Frequency_
        {
            Frequency_() = default;
            Frequency_(const std::vector<size_t> &group) : head_(group.front()),
                                                           weight_(group.size()) {}
            size_t head_ = -1;
            size_t weight_ = 1;
        };
        std::vector<Frequency_> frequencies_ = {};         /// Frequency groups
        std::unordered_set<size_t> frequencies_data_ = {}; /// Frequency-relevant data indexes

        //////////////////////
        //////  CACHED  //////
        //////////////////////

        struct Cached_
        {
            struct State_
            {
                std::vector<size_t> size_val;
                std::vector<std::vector<size_t>> size_vec = {};
                std::vector<std::vector<double>> dbl_vec = {};
                std::vector<std::vector<std::vector<double>>> dbl_mtx = {};
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
        /// Use the fill frequencies for the all distinc case.
        void helper_init_unique_frequencies_();
    };
}
