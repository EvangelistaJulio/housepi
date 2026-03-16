/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "tools/types.h"

#include "tools/data_reports.h"
#include "household.h"

namespace hsp
{
    namespace rates
    {
        constexpr const char *REFERENCE = ".reference";  /// Rate reference choice.
        constexpr const char *MISSING = ".missing_rate"; /// Missing rate.
        constexpr const char *CONTACT = ".contact";      /// Contanct rate.
        constexpr double VALUE_REFERENCE = 1.0;          /// Alias for reference value.
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    namespace key
    {
        struct Registry;
    }

    struct DataWide
    {
        std::unique_ptr<key::Registry> pRegistry;
        std::unique_ptr<key::Registry> p_size_Registry;
        std::unique_ptr<key::Registry> p_size_vec_Registry;
        std::unique_ptr<key::Registry> p_double_vec_Registry;
        std::unique_ptr<key::Registry> p_double_mtx_Registry;

        ReportIndividual report_indiv;
        ReportHousehold report_house;
        ReportHouseSizes report_hh_size;

        std::vector<std::vector<double>> smpl_values = {};

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        DataWide();
        ~DataWide();

        void init_registry(const std::vector<std::string> &keys,
                           const std::vector<std::string> &size_keys,
                           const std::vector<std::string> &size_vec_keys,
                           const std::vector<std::string> &double_vec_keys,
                           const std::vector<std::string> &double_mtx_keys);

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign data (wide format).
        void build(const str_table &input);

        /// Set a new multiplicity-relevant data entry.
        void set_multiplicity_data(size_t d) { multiplicity_data_.emplace(d); }
        /// Build multiplicities from household data.
        void multiplicity_build(bool is_any_augm);

        /// Initialize cached structeres according to household size.
        void allocate_all_cache();

        /// Store all household cache attributes.
        void save_all_cached() { old_ = cached_; }
        /// Recover old values of cached depedencies and attributes.
        void recover_old_cached() noexcept { cached_.swap(old_); }

        /// Save member cached attributes.
        void save_member_cache(Household &house, size_t i);
        /// Recover member old values of cached member depedencies and member attributes.
        void recover_member_old_cached(Household &house, size_t i);

        /// Raise error log-likelihood is not finite.
        void assert_finite_log_lik() const;

        /// Initialize posterior sample of parameters.
        void init_draws(unsigned long sample_size);
        /// Store current sample of log-likelihood for restricted households.
        void store_draws(unsigned long smpl_count);

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        Household &at(size_t H) { return house_vector_[H]; }
        const Household &at(size_t H) const { return house_vector_[H]; }
        Household &at(std::string_view house_id) { return at(house_id_to_index_(house_id)); }
        const Household &at(std::string_view house_id) const { return at(house_id_to_index_(house_id)); }

        /// Number of households in the multiplicity group.
        size_t multiplicity(size_t R) const
        {
            switch (mode_)
            {
            case Mode::Unique:
                return 1;
            case Mode::Multiple:
                return multiplicity_[R].value_;
            }
            __builtin_unreachable();
        }
        /// Index to household multiplicity group.
        size_t get_restricted_index(size_t H) const
        {
            switch (mode_)
            {
            case Mode::Unique:
                return H;
            case Mode::Multiple:
                return house_vector_[H].restricted_;
            }
            __builtin_unreachable();
        }
        /// Number of household multiplicity groups.
        size_t restricted_size() const noexcept
        {
            switch (mode_)
            {
            case Mode::Unique:
                return report_house.total;
            case Mode::Multiple:
                return multiplicity_.size();
            }
            __builtin_unreachable();
        }
        /// Restrict to head of household multiplicity group.
        Household &restricted(size_t R)
        {
            switch (mode_)
            {
            case Mode::Unique:
                return at(R);
            case Mode::Multiple:
                return at(multiplicity_[R].head_);
            }
            __builtin_unreachable();
        }
        const Household &restricted(size_t R) const
        {
            switch (mode_)
            {
            case Mode::Unique:
                return at(R);
            case Mode::Multiple:
                return at(multiplicity_[R].head_);
            }
            __builtin_unreachable();
        }

        /////////////////////
        //////  CACHE  //////
        /////////////////////

        /// Return cached size_t value.
        size_t &cache_size(Household &house, size_t e) { return cached_[house.restricted_].size_[e]; }
        size_t cache_size(const Household &house, size_t e) const { return cached_[house.restricted_].size_[e]; }

        /// Return cached size_t vector.
        std::vector<size_t> &cache_size_vector(Household &house, size_t e) { return cached_[house.restricted_].size_vec_[e]; }
        const std::vector<size_t> &cache_size_vector(const Household &house, size_t e) const { return cached_[house.restricted_].size_vec_[e]; }

        /// Return cached double vector.
        std::vector<double> &cache_double_vector(Household &house, size_t e) { return cached_[house.restricted_].dbl_vec_[e]; }
        const std::vector<double> &cache_double_vector(const Household &house, size_t e) const { return cached_[house.restricted_].dbl_vec_[e]; }

        /// Return cached double matrix.
        std::vector<std::vector<double>> &cache_double_matrix(Household &house, size_t e) { return cached_[house.restricted_].dbl_mtx_[e]; }
        const std::vector<std::vector<double>> &cache_double_matrix(const Household &house, size_t e) const { return cached_[house.restricted_].dbl_mtx_[e]; }

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Assert data wide is not empty.
        void assert_non_empty() const;
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

        std::vector<Household> house_vector_ = {}; /// All households in data
        dict<size_t> hhid_to_index_ = {};          /// Map house-id to internal index

        ////////////////////////////
        //////  MULTIPLICITY  //////
        ////////////////////////////

        enum class Mode : std::uint8_t
        {
            Multiple,
            Unique
        } mode_ = Mode::Unique;

        struct Multiplicity_
        {
            Multiplicity_(const std::vector<size_t> &group) : head_(group.front()),
                                                              value_(group.size()) {}
            size_t head_;
            size_t value_;
        };
        std::unordered_set<size_t> multiplicity_data_ = {}; /// Multiplicity-relevant data
        std::vector<Multiplicity_> multiplicity_ = {};      /// Multiplicity groups

        //////////////////////
        //////  CACHED  //////
        //////////////////////

        struct Cached_
        {
            std::vector<size_t> size_ = {};
            std::vector<std::vector<size_t>> size_vec_ = {};
            std::vector<std::vector<double>> dbl_vec_ = {};
            std::vector<std::vector<std::vector<double>>> dbl_mtx_ = {};
        };
        std::vector<Cached_> cached_ = {};
        std::vector<Cached_> old_ = {};

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
    };
}
