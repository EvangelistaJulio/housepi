/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__BINOMIAL__H
#define HOUSEPI__BINOMIAL__H

#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  NAMESPACE : bax  /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace bax
{
    /// Functions for the generations lookup table.
    void init_generations_lookuptable(unsigned max_N);

    unsigned get_num_generations(unsigned n_positive);
    std::vector<unsigned> get_generation_vector(unsigned n_positive, unsigned i);
    unsigned get_generation_multipl(unsigned n_positive, unsigned i);

    /// Recursive likelihood formulation lookuptable.
    void init_escape_patterns_lookuptable(unsigned max_N, bool clean_combinations = true);

    /// Heterogeneous chain binomial log-likelihood for ordered probabilities (seropositives, seronegatives).
    /// @param N_star n. of seropositves
    /// @param commu_esc_proba ordered community escape probabilities
    /// @param p2p_esc_proba ordered person-to-person escape probabilities
    double chain_binomial_heter_log_lik(unsigned N_star,
                                        const dbl_vec &commu_esc_proba,
                                        const matrix<double> &p2p_esc_proba);
}

#endif // HOUSEPI__BINOMIAL__H
