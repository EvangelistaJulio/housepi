/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__PROPOSAL__H
#define HOUSEPI__PROPOSAL__H

#include <memory>
#include <string>
#include <vector>

#include "tools/random.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Proposal : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Attribute;
class Distribution;

/// Base class for Metropolis proposers.
class Proposal
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Proposal : constructor  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Empty proposal constructor.
    Proposal();
    /// Proposal constructor from Distribution pointer.
    Proposal(std::shared_ptr<Distribution> distr,
             double proposal_scale,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);
    /// Proposal constructor from Distribution name.
    Proposal(const std::string &distr_name,
             double proposal_scale,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);
    /// Proposal deconstructor.
    virtual ~Proposal();

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposal : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string name() const { return this->name_; }
    virtual double propose(hsp::rnd_gen &, double) = 0;
    virtual double logratio(double, double) const { return 0.; }

    /// Wheter proposer can be adapted.
    bool is_adaptable() const { return !args_to_adapt_.empty(); }

    /// Return proposer scale.
    double scale() const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposal : setter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    virtual void user_reset() {}
    virtual void associate_location(const Attribute &) {}

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Proposal : adaptation  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double get_accept_rate() const { return accept_rate_; }
    void update_accept_rate(bool is_accepted)
    {
        accept_rate_ = static_cast<double>(accept_rate_ * n_proposed_ + is_accepted) / (1. + n_proposed_);
        ++n_proposed_;
    }
    void reset_accept_rate()
    {
        accept_rate_ = 0.;
        n_proposed_ = 0;
    }
    /// Increase spread if acceptation rate is high and viceversa.
    void adapt(double ratio);

protected:
    std::string name_ = "";                /// name of the proposal
    std::shared_ptr<Distribution> pDistr_; /// distribution associated to proposal

    double user_scale_ = 1.;                         /// scale of the Metropolis proposal
    double accept_rate_ = 0.;                        /// acceptance rate
    unsigned long n_proposed_ = 0;                   /// n. of proposed moves
    std::vector<unsigned short> args_to_adapt_ = {}; /// which arguments should be adapted

    //===== propose : functional methods =====//
    double do_propose_sum_to_previous_(hsp::rnd_gen &prng, double xo);
    double do_propose_times_previous_(hsp::rnd_gen &prng, double xo);
    double do_propose_sample_rand_(hsp::rnd_gen &prng, double);

    //===== logratio : functional methods =====//
    double do_logratio_pdf_(double xold, double xnew) const;
    double do_logratio_symm_(double, double) const { return 0.; }
    double do_logratio_lognorm_symm_(double xold, double xnew) const;
};

#endif // HOUSEPI__PROPOSAL__H
