/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>
#include "tools/keys.h"
#include "tools/parser.h" // UPDATE: this shouldn't be public
#include "tools/random.h"
#include "tools/types.h"
#include "attribute.h"
#include "variable.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Parameter : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Individual;

/// Implementation of Parameters.
struct Parameter
{
    /// NOTE: list_const + list_sample = all parameters in param_vec_
    /// NOTE: list_sample_standard + list_sample_cached = list_sample
    /// NOTE: list_adaptable_sample_standard + list_non_adapt_sample_standard = list_sample_standard
    /// NOTE: list_adaptable_sample_cached + list_non_adapt_sample_cached = list_sample_cached

    std::vector<size_t> list_const = {};  /// list of constant parameters
    std::vector<size_t> list_sample = {}; /// list of parameters to sample

    std::vector<size_t> list_sample_standard = {};
    std::vector<size_t> list_sample_cached = {};

    std::vector<size_t> list_adaptable_sample_standard = {};
    std::vector<size_t> list_non_adapt_sample_standard = {};
    std::vector<size_t> list_adaptable_sample_cached = {};
    std::vector<size_t> list_non_adapt_sample_cached = {};

    std::vector<std::vector<double>> smpl_values = {}; /// values of the sample
    std::vector<std::vector<double>> smpl_priors = {}; /// prior values of the sample

    std::unordered_set<size_t> user_discrete_set = {}; /// user defined discrete parameters

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Parameter : getter  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    /// Get reference to parameter Variable class.
    Variable &at(size_t p) { return param_vec_[p]; }
    /// Get reference to parameter constant Variable class.
    const Variable &at(size_t p) const { return param_vec_[p]; }

    /// Flag to check if there are parameters to sample.
    bool is_any_sample() const { return !list_sample.empty(); }

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Parameter : setter  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    void build(const hsp::str_table &deftable);
    void build_for_expec(const hsp::dbl_table &deftable);
    void make_cache_dependent(size_t p);
    void make_cache_dependent(const hsp::parse::Entry &node);

    /// Randomly initialize values of the parameters to be inferred.
    void user_reset(hsp::rnd_gen &prng);
    /// Load parameters values of the n-th posterrior sample.
    void load_sample_values(size_t n);

    /// Initialize posterior sample of parameters and respective priors.
    void init_draws(unsigned long sample_size);
    /// Store current sample of parameters and respective priors.
    void store_draws(unsigned long smpl_count);

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Parameter : definition helpers  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Auxiliary method to read a definition as a single attribute.
    const Attribute &read_attribute(const hsp::parse::Entry &node, const Individual *member_ptr = nullptr) const;
    /// Auxiliary method to convert a typed string into a numeric attribute.
    const Attribute &read_numeric(const hsp::parse::Entry &node, const Individual *member_ptr = nullptr) const;

    /// Auxiliary method to read a definition as a single attribute (overload for string).
    const Attribute &read_attribute(const std::string &expr, const Individual *member_ptr = nullptr) const
    {
        return read_attribute(hsp::parse::Definition(expr).root, member_ptr);
    }
    /// Auxiliary method to read a definition as a single attribute (overload for string).
    const Attribute &read_numeric(const std::string &expr, const Individual *member_ptr = nullptr) const
    {
        return read_numeric(hsp::parse::Definition(expr).root, member_ptr);
    }

private:
    std::vector<Variable> param_vec_ = {}; /// parameters vector
    mutable Attribute storage_;            /// auxiliary storage attribute used for reading.

    void clear_lists_();
    void helper_build_constant_(const hsp::str_table &deftable, size_t k);
    void helper_build_to_sample_(const hsp::str_table &deftable, size_t k);
};
