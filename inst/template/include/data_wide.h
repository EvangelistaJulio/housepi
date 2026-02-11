/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "tools/types.h"

#include "tools/data_reports.h"
#include "household.h"

namespace hsp::rates
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

namespace hsp::key
{
    class Registry;
}

struct ReportIndividual;
struct ReportHousehold;
struct ReportHouseSizes;

struct DataWide
{
    std::unique_ptr<hsp::key::Registry> pRegistry;
    std::unique_ptr<hsp::key::Registry> p_size_Registry;
    std::unique_ptr<hsp::key::Registry> p_size_vec_Registry;
    std::unique_ptr<hsp::key::Registry> p_double_vec_Registry;
    std::unique_ptr<hsp::key::Registry> p_double_mtx_Registry;

    std::vector<std::vector<double>> smpl_values = {};
    std::vector<size_t> head_class_house_indexes = {};

    ReportIndividual report_indiv;
    ReportHousehold report_house;
    ReportHouseSizes report_hh_size;

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    DataWide(const std::vector<std::string> &keys,
             const std::vector<std::string> &size_keys,
             const std::vector<std::string> &size_vec_keys,
             const std::vector<std::string> &double_vec_keys,
             const std::vector<std::string> &double_mtx_keys);
    ~DataWide();

    Household &at(size_t H) { return house_vector_[H]; }
    const Household &at(size_t H) const { return house_vector_[H]; }
    Household &at(std::string_view house_id) { return at(house_id_to_index_(house_id)); }
    const Household &at(std::string_view house_id) const { return at(house_id_to_index_(house_id)); }

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Assign data (wide format).
    void build(const hsp::str_table &input);
    /// Build multiplicities.
    /// UPDATE: I should find a way to deal internally with multiplicity.
    void multiplicity_build(const hsp::str_table &input, bool is_any_augm);

    /// Initialize cached structeres according to household size.
    void init_house_cache(Household &house);
    /// Deep copy all household cache attributes in the buffer.
    void buffer_all_cache(Household &house);
    /// Deep copy member cache attributes in the buffer.
    void buffer_member_cache(Household &house, size_t i);
    /// Recover all buffered values of cached depedencies.
    void recover_all_buffer(Household &house);
    /// Recover all buffered values of member cached depedencies.
    void recover_member_buffer(Household &house, size_t i);

    void assert_finite_log_lik() const;

    /// Initialize posterior sample of parameters.
    void init_draws(unsigned long sample_size);
    /// Store current sample of log-likelihood for restricted households.
    void store_draws(unsigned long smpl_count);

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    /// Assert data wide is not empty.
    void assert_non_empty();
    /// Assert that not all individuals have missing status.
    void assert_not_all_missing_status();
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

    Household &restricted(size_t R) { return house_vector_[head_class_house_indexes[R]]; }

private:
    std::vector<Household> house_vector_ = {};
    hsp::dict<size_t> hhid_to_index_ = {};

    /// Get housepi household index from house_id.
    size_t house_id_to_index_(std::string_view house_id) const;
    /// Clear previous data report.
    void helper_build_clear_data_();
    /// Assert that data headers are correct.
    void helper_build_data_headers_(const hsp::str_table &input);
    /// Use the loaded data to fill-in the respective report.
    void build_helper_fill_report_data_();
};
