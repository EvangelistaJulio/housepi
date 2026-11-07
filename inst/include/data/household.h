/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <vector>

#include "data/individual.h"

//////////////////////////////  Household  //////////////////////////////

namespace housepi
{
    namespace data
    {
        struct Dataset;
    }

    /// Members of one household, partitioned by status.
    ///
    /// The user household id lives in data::Dataset, keyed by household index, see Index::of().
    class Household
    {
        /// Reads the opaque storage index.
        friend struct Index;
        /// Owns and mutates its households.
        friend struct data::Dataset;
        /// Changes member status.
        friend class ModelBase;

        /// Range of member indexes sharing a status.
        struct MemberRange_
        {
            const MemberIdx *data_;
            size_t size_;
            const MemberIdx *begin() const noexcept { return data_; }
            const MemberIdx *end() const noexcept { return data_ + size_; }
            const MemberIdx &operator[](size_t i) const noexcept { return data_[i]; }
        };

        /// Range of individuals sharing a status.
        struct ConstIndivRange_
        {
            const Household *house_;
            const MemberIdx *data_;
            size_t size_;
            struct Iterator_
            {
                const Household *house_;
                const MemberIdx *ptr_;
                const Individual &operator*() const noexcept { return house_->indivs_[*ptr_]; }
                Iterator_ &operator++() noexcept
                {
                    ++ptr_;
                    return *this;
                }
                bool operator!=(const Iterator_ &other) const noexcept { return ptr_ != other.ptr_; }
            };
            Iterator_ begin() const noexcept { return {house_, data_}; }
            Iterator_ end() const noexcept { return {house_, data_ + size_}; }
        };

    public:
        //////  Setter  //////

        /// Emplace an individual with `index` and `status`.
        void add_indiv(IndivIdx index, StatusLvl status);

        //////  Getter  //////

        /// Individual at member index `i`.
        Individual &at(MemberIdx i) noexcept { return indivs_[i]; }
        const Individual &at(MemberIdx i) const noexcept { return indivs_[i]; }

        /// Range of member indexes with infected status.
        MemberRange_ infected_members() const noexcept { return {status_vec_.data(), ninf_}; }
        /// Range of member indexes with escaped status.
        MemberRange_ escaped_members() const noexcept { return {status_vec_.data() + ninf_, nsus_}; }

        /// Range of individuals with infected status.
        ConstIndivRange_ infected_indivs() const noexcept { return {this, status_vec_.data(), ninf_}; }
        /// Range of individuals with escaped status.
        ConstIndivRange_ escaped_indivs() const noexcept { return {this, status_vec_.data() + ninf_, nsus_}; }

        /// Number of members with infected status.
        MemberIdx n_infected_status() const noexcept { return ninf_; }
        /// Number of members with escaped status.
        MemberIdx n_escaped_status() const noexcept { return nsus_; }
        /// Number of members with missing status.
        MemberIdx n_missing_status() const noexcept { return hh_size_ - ninf_ - nsus_; }

        //////  Iterators  //////

        auto begin() noexcept { return indivs_.begin(); }
        auto end() noexcept { return indivs_.end(); }
        auto begin() const noexcept { return indivs_.begin(); }
        auto end() const noexcept { return indivs_.end(); }

        //////  Show  //////

        /// Number of household members.
        MemberIdx size() const noexcept { return hh_size_; }

    private:
        /// Internal household index.
        HouseIdx house_index_ = 0;
        MemberIdx hh_size_ = 0;
        /// Number of infected members.
        MemberIdx ninf_ = 0;
        /// Number of escaped members.
        MemberIdx nsus_ = 0;
        std::vector<Individual> indivs_ = {};
        /// Partitioned: [infected | escaped | missing].
        std::vector<MemberIdx> status_vec_ = {};

        /// Move `indiv` to status `new_status`.
        ///
        /// Member statuses are tracked by partition refinement.
        void set_new_status_(Individual &indiv,
                             StatusLvl new_status);
    };

    /// Storage index of `house`, declared in data/individual.h.
    inline HouseIdx Index::of(const Household &house) noexcept { return house.house_index_; }
}
