/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "tools/data_reports.h"
#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/status.h"
#include "tools/string.h"
#include "tools/types.h"

#include "data_long.h"
#include "data_wide.h"
#include "function.h"
#include "household.h"
#include "individual.h"

namespace housepi
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** Household : class ***  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Household::add_member(size_t index,
                               const std::string &id,
                               int status)
    {
        size_t hhindex = size();
        members_.emplace_back();
        ++hhsize_;

        Individual &last = members_.back();
        last.index_ = index;
        last.hhindex_ = hhindex;
        last.id_ = id;

        // new member added as if missing status
        last.status_ = status::MISSING;
        status_vec_.push_back(hhindex);
        if (!status::is_missing(status))
            set_new_status_(last, status);
    }

    void Household::set_new_status_(Individual &member, int snew)
    {
        const size_t i = member.hhindex();
        size_t cold = status::category(member.status());
        size_t cnew = status::category(snew);

        member.status_ = snew; // NOTE: might just change within infected status
        if (cold == cnew)
            return;

        // segment boundaries [infected | escaped | missing]
        size_t bounds[] = {0, ninf_, ninf_ + nsus_, hhsize_};
        // extract old position given the old segment
        size_t pos = std::find(&status_vec_[bounds[cold]], &status_vec_[bounds[cold + 1]], i) -
                     &status_vec_[0];

        // alias for swap between the status segments
        auto segment_swap = [&](size_t a, size_t b)
        { std::swap(status_vec_[a], status_vec_[b]); };

        if (cold == 0 && cnew == 1)
        {
            segment_swap(pos, ninf_ - 1);
            --ninf_;
            ++nsus_;
        }
        else if (cold == 1 && cnew == 0)
        {
            segment_swap(pos, ninf_);
            ++ninf_;
            --nsus_;
        }
        else if (cold == 1 && cnew == 2)
        {
            segment_swap(pos, ninf_ + nsus_ - 1);
            --nsus_;
        }
        else if (cold == 2 && cnew == 1)
        {
            segment_swap(pos, ninf_ + nsus_);
            ++nsus_;
        }
        else if (cold == 0 && cnew == 2)
        {
            segment_swap(pos, ninf_ - 1);
            if (nsus_ > 0)
                segment_swap(ninf_ - 1, ninf_ + nsus_ - 1);
            --ninf_;
        }
        else // if (cold == 2 && cnew == 0)
        {
            segment_swap(pos, ninf_ + nsus_);
            if (nsus_ > 0)
                segment_swap(ninf_ + nsus_, ninf_);
            ++ninf_;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void DataWide::init_cache_key_holders(
        const std::vector<std::string> &size_keys,
        const std::vector<std::string> &size_vec_keys,
        const std::vector<std::string> &double_vec_keys,
        const std::vector<std::string> &double_mtx_keys)
    {
        khold_uint_val.init(size_keys);
        khold_uint_vec.init(size_vec_keys);
        khold_double_vec.init(double_vec_keys);
        khold_double_mtx.init(double_mtx_keys);
    }

    void DataWide::build(const str_table &input)
    {
        khold_.clear_optional();
        khold_uint_val.clear_optional();
        khold_uint_vec.clear_optional();
        khold_double_vec.clear_optional();
        khold_double_mtx.clear_optional();

        helper_build_clear_data_();        // drop previous data and lists
        helper_build_data_headers_(input); // extract headers
        report_indiv.total = read::csv_n_rows(input);

        // pre-allocate external data storage
        indiv_data_.resize(report_indiv.total);

        // store households
        for (size_t i = 0; i < report_indiv.total; ++i)
        {
            std::string id = std::to_string(i + 1);                          // default id
            std::string house_id = std::to_string(house_vector_.size() + 1); // default house_id
            int status = status::MISSING;

            if (input.count("id"))
            {
                id = input.at("id")[i]; // user id
                if (na::is(id))
                    msg::stop(
                        "Individual in row n." + std::to_string(i + 1) + " had missing " + msg::field("id") + ".");
            }
            if (input.count("house_id"))
            {
                house_id = input.at("house_id")[i]; // user house_id
                if (na::is(house_id))
                    msg::stop(
                        "Individual in row n." + std::to_string(i + 1) + " had missing " + msg::field("house_id") + ".");
            }
            if (input.count("status"))
            {
                // note: NA status is a missing status
                if (!na::is(input.at("status")[i]))
                    status = str::toi(input.at("status")[i], "individual status");
            }

            // read individual data into external storage
            IndivWide_ &idata = indiv_data_[i];
            idata.status = status;
            idata.vec.resize(khold_.size());
            for (size_t d = 0; d < idata.vec.size(); ++d)
            {
                const std::string &dname = khold(d);
                double &x = idata.vec[d];

                // placeholder for model-required columns absent at import
                // NOTE: capital check is deferred to the R-side wall (see build_full_)
                x = na::DOUBLE;

                if (input.count(dname))
                {
                    x = str::tod(input.at(dname)[i], "data " + msg::arg(dname));

                    if (user_discrete_set_.count(d) && x != std::floor(x))
                        msg::stop(
                            "Data " + msg::arg(dname) + " for individual in row n." + std::to_string(i + 1) +
                            " was expected to be discrete-valued, as requested by model.");
                    if (user_positive_set_.count(d) && !(x > 0.))
                        msg::stop(
                            "Data " + msg::arg(dname) + " for individual in row n." + std::to_string(i + 1) +
                            " was expected to be greater than zero, as requested by model.");
                }
            }

            // per-column read-surface defaults to the imported raw value; derived
            // columns (household observations) are wired to their Variable later.
            idata.cell.resize(idata.vec.size());
            for (size_t d = 0; d < idata.vec.size(); ++d)
                idata.cell[d].init_numeric(idata.vec[d]);

            // if household index doesn't exist, add new household
            if (!hhid_to_index_.count(house_id))
            {
                // map data index to household number
                size_t H = house_vector_.size();
                hhid_to_index_.emplace(house_id, H);
                // create a new empty household
                house_vector_.emplace_back();
                Household &last = house_vector_.back();
                last.house_id_ = house_id;
                last.index_ = H;
            }
            // extract the household reference
            Household &house = house_vector_.at(hhid_to_index_.at(house_id));

            // assign the individual the respective household
            house.add_member(i, id, idata.status);
        }

        // update report for data
        build_helper_fill_report_data_();
    }

    void DataWide::optimize_cells()
    {
        for (IndivWide_ &idata : indiv_data_)
            for (Attribute &attr : idata.cell)
                attr.optimize();
    }

    size_t DataWide::house_id_to_index_(std::string_view house_id) const
    {
        auto it = hhid_to_index_.find(house_id);
        if (it == hhid_to_index_.end())
            msg::stop(
                "Index " + msg::val(std::string(house_id)) + " did not match any household in the data.");
        return it->second;
    }

    void DataWide::helper_build_clear_data_()
    {
        house_vector_.clear();
        hhid_to_index_.clear();
        indiv_data_.clear();
        report_indiv = ReportIndividual();
        report_house = ReportHousehold();
        report_hh_size = ReportHouseSizes();
    }

    void DataWide::helper_build_data_headers_(const str_table &input)
    {
        std::vector<std::string> headers = read::extract_keys(input);

        // assess the presence of house IDs
        if (!contains(headers, "house_id"))
            msg::warn("Ignored household structure: missing data column " + msg::field("house_id") + ".");

        // assess individual IDs are unique
        if (contains(headers, "id") &&
            input.at("id").size() != str::unique(input.at("id")).size())
            msg::stop("Provided individual IDs are not unique.");

        // remove standard headers
        headers = headers - std::vector<std::string>({"id", "house_id", "status"});

        // add optional keys to registry
        for (const std::string &datastring : headers)
            khold_.add(datastring);
    }

    void DataWide::build_helper_fill_report_data_()
    {
        report_house.total = house_vector_.size();
        report_house.distinct = report_house.total;

        // set report counts
        for (const Household &house : house_vector_)
        {
            // observed quantities in-house
            size_t house_size = house.size();
            size_t n_infec = house.n_infected_status();
            size_t n_miss = house.n_missing_status();
            size_t n_esc = house.n_escaped_status();

            report_indiv.infectious += n_infec;
            report_indiv.missing += n_miss;
            report_indiv.escaped += n_esc;

            report_house.no_cases += (n_infec == 0);
            report_house.at_least_1_case += (n_infec > 0);
            report_house.at_least_2_cases += (n_infec > 1);

#define UPDATE_REPORT_SIZE(Report, Value, Symbol) \
    if (Report.first Symbol Value)                \
    {                                             \
        Report.first = Value;                     \
        Report.second = 1;                        \
    }                                             \
    else if (Report.first == Value)               \
        ++Report.second;

            // max number of house cases
            UPDATE_REPORT_SIZE(report_hh_size.max_cases, n_infec, <);
            // max number of missing status
            UPDATE_REPORT_SIZE(report_hh_size.max_miss, n_miss, <);
            // minimum house number of members
            UPDATE_REPORT_SIZE(report_hh_size.min_size, house_size, >);
            // maximum house number of members
            UPDATE_REPORT_SIZE(report_hh_size.max_size, house_size, <);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** DataWide : frequency weights ***  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        struct FrequencyMapper_
        {
            FrequencyMapper_(const std::unordered_set<size_t> data_list, const DataWide &data_wide)
                : data_list_(data_list), data_wide_(data_wide) {}

            /// Compute frequency-relevant hash for each household
            std::map<size_t, std::vector<size_t>> define_hash_groups(
                const std::vector<Household> &house_vector)
            {
                std::map<size_t, std::vector<size_t>> output;

                for (size_t H = 0; H < house_vector.size(); ++H)
                {
                    const Household &house = house_vector[H];

                    // sorted member indices
                    std::vector<size_t> order(house.size());
                    std::iota(order.begin(), order.end(), 0);
                    std::sort(order.begin(), order.end(),
                              [this, &house](size_t i, size_t j)
                              { return member_less_(house.member(i), house.member(j)); });

                    // equal house sizes share seeding
                    size_t seed = std::hash<size_t>{}(house.size());
                    // combine hashes over sorted members
                    for (size_t i : order)
                    {
                        const Individual &member = house.member(i);
                        hash_combine_(seed, std::hash<int>{}(member.status()));
                        for (size_t d : data_list_)
                            hash_combine_(seed, std::hash<double>{}(data_wide_.data(d, member)));
                    }
                    output[seed].emplace_back(H);
                }

                return output;
            }

        private:
            /// List of data relevant for frequency groups
            const std::unordered_set<size_t> data_list_;
            /// Wide-format data storage
            const DataWide &data_wide_;

            /// Comparison functor: sort members by (status, data[0], data[1], ...)
            bool member_less_(const Individual &a, const Individual &b) const
            {
                if (a.status() != b.status())
                    return a.status() < b.status();
                for (size_t d : data_list_)
                    if (data_wide_.data(d, a) != data_wide_.data(d, b))
                        return data_wide_.data(d, a) < data_wide_.data(d, b);
                return false;
            }

            /// Using golden ratio fractional part to ensure bit distribution.
            /// Based on @cite james_boostcontainerhash_2022
            void hash_combine_(size_t &seed, size_t h)
            {
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
                seed ^= h + 0x9e3779b97f4a7c15ULL + (seed << 12) + (seed >> 4);
#else
                seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
#endif
            }
        };
    }

    void DataWide::set_eq_class_data_names(const std::vector<std::string> &data_list)
    {
        for (size_t d : khold(data_list))
            eq_class_covariates_.emplace(d);
    }

    void DataWide::init_eq_classes(bool force_unique)
    {
        if (force_unique)
        {
            helper_init_unique_eq_classes_();
            return;
        }

        // compute frequency-relevant hash for each household
        std::map<size_t, std::vector<size_t>> hash_groups =
            FrequencyMapper_(eq_class_covariates_, *this).define_hash_groups(house_vector_);

        // update number of distinct frequency groups
        report_house.distinct = hash_groups.size();

        // init all unique households and skip if it's the case
        if (report_house.total == report_house.distinct)
        {
            helper_init_unique_eq_classes_();
            return;
        }

        // compute frequency heads and weights
        eq_classes_.clear();
        eq_classes_.reserve(report_house.distinct);
        for (const auto &[hash, group] : hash_groups)
        {
            const size_t refer_restrict = eq_classes_.size();
            for (size_t H : group)
                at(H).eq_class_ = refer_restrict;

            eq_classes_.emplace_back(group);
        }
    }

    void DataWide::helper_init_unique_eq_classes_()
    {
        // restore distinct count
        report_house.distinct = report_house.total;

        // skip if already unique equivalence classes
        if (eq_classes_.size() == report_house.total)
            return;

        eq_classes_.clear();
        eq_classes_.resize(report_house.total);
        for (size_t H = 0; H < report_house.total; ++H)
        {
            eq_classes_[H].head_ = H;
            at(H).eq_class_ = H;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : cache attributes  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void DataWide::allocate_all_cache()
    {
        cached_.previous.resize(eq_class_size());
        cached_.current.resize(eq_class_size());

        for (size_t R = 0; R < eq_class_size(); ++R)
        {
            Household &house = eq_class_head(R);
            size_t hh_size = house.size();
            Cached_::State_ &prev = cached_.previous[R];
            Cached_::State_ &curr = cached_.current[R];

            // size_t
            size_t s_size = khold_uint_val.size();
            curr.uint_val.assign(s_size, 0);
            prev.uint_val.assign(s_size, 0);

            // size_t vector
            size_t s_size_vec = khold_uint_vec.size();
            curr.uint_vec.resize(s_size_vec);
            prev.uint_vec.resize(s_size_vec);
            for (size_t e = 0; e < s_size_vec; ++e)
            {
                curr.uint_vec[e].assign(hh_size, 0);
                prev.uint_vec[e].assign(hh_size, 0);
            }

            // double vector
            size_t s_double_vec = khold_double_vec.size();
            curr.dbl_vec.resize(s_double_vec);
            prev.dbl_vec.resize(s_double_vec);
            for (size_t e = 0; e < s_double_vec; ++e)
            {
                curr.dbl_vec[e].assign(hh_size, 0);
                prev.dbl_vec[e].assign(hh_size, 0);
            }

            // double matrix
            size_t s_double_mtx = khold_double_mtx.size();
            curr.dbl_mtx.resize(s_double_mtx);
            prev.dbl_mtx.resize(s_double_mtx);
            for (size_t e = 0; e < s_double_mtx; ++e)
            {
                curr.dbl_mtx[e].assign(hh_size, std::vector<double>(hh_size, 0.));
                prev.dbl_mtx[e].assign(hh_size, std::vector<double>(hh_size, 0.));
            }
        }
    }

    void DataWide::save_member_cache(Household &house, size_t i)
    {
        Cached_::State_ &prev = cached_.previous[house.eq_class_];
        Cached_::State_ &curr = cached_.current[house.eq_class_];
        prev.uint_val = curr.uint_val;

        for (size_t e = 0; e < curr.uint_vec.size(); ++e)
            prev.uint_vec[e][i] = curr.uint_vec[e][i];

        for (size_t e = 0; e < curr.dbl_vec.size(); ++e)
            prev.dbl_vec[e][i] = curr.dbl_vec[e][i];

        for (size_t e = 0; e < curr.dbl_mtx.size(); ++e)
            for (size_t j = 0; j < house.size(); ++j)
            {
                prev.dbl_mtx[e][i][j] = curr.dbl_mtx[e][i][j];
                if (j != i)
                    prev.dbl_mtx[e][j][i] = curr.dbl_mtx[e][j][i];
            }
    }

    void DataWide::recover_member_old_cached(Household &house, size_t i)
    {
        Cached_::State_ &prev = cached_.previous[house.eq_class_];
        Cached_::State_ &curr = cached_.current[house.eq_class_];
        curr.uint_val.swap(prev.uint_val);

        for (size_t e = 0; e < prev.uint_vec.size(); ++e)
            curr.uint_vec[e][i] = prev.uint_vec[e][i];

        for (size_t e = 0; e < prev.dbl_vec.size(); ++e)
            curr.dbl_vec[e][i] = prev.dbl_vec[e][i];

        for (size_t e = 0; e < prev.dbl_mtx.size(); ++e)
        {
            curr.dbl_mtx[e][i].swap(prev.dbl_mtx[e][i]);
            for (size_t j = 0; j < house.size(); ++j)
                if (j != i)
                    curr.dbl_mtx[e][j][i] = prev.dbl_mtx[e][j][i];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  DataWide : posterior draws  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void DataWide::assert_not_all_missing_status() const
    {
        if (report_indiv.total == report_indiv.missing)
            msg::stop(
                "All individuals in the data had missing status.");
    }

    void DataWide::assert_no_missing_status() const
    {
        if (report_indiv.missing > 0)
            msg::stop(
                "Individual with missing status are incompatible with current model.");
    }

    void DataWide::recover_all_status() noexcept
    {
        for (Household &house : house_vector_)
            for (Individual &member: house)
                house.set_new_status_(member, indiv_data_[Index::of(member)].status);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** DataLong : struct ***  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void DataLong::build(DataWide &data_wide, const str_table &input)
    {
        khold_.clear_optional();

        if (input.empty())
            return;

        helper_build_clear_data_();
        helper_build_data_headers_(data_wide, input);

        // data-long allocation
        size_t n_indiv = data_wide.report_indiv.total;
        indiv_datalist.resize(n_indiv);

        // individual map for id -> index
        dict<size_t> id_to_index = {};
        for (const Household &house : data_wide)
            for (const Individual &member : house)
                id_to_index.emplace(member.id(), Index::of(member));

        // load individual data list
        size_t n_rows = read::csv_n_rows(input);
        for (size_t r = 0; r < n_rows; ++r)
        {
            const std::string &id = input.at("id")[r];
            if (na::is(id))
                msg::stop(
                    "Invalid id found in row n." + std::to_string(r + 1) + ".");
            // UPDATE: these errors should be more self-explicative
            if (!id_to_index.count(id))
                msg::stop(
                    "Id " + msg::val(id) + " from data long-format had no match in loaded data wide-format.");

            // append (since vector is initialized to empty)
            auto &idat = indiv_datalist[id_to_index.at(id)];
            idat.vec.resize(khold_.size());
            for (size_t d = 0; d < khold_.size(); ++d)
            {
                const std::string &dstr = khold(d);
                idat.vec[d].emplace_back(
                    str::tod(input.at(dstr)[r], "data long format " + msg::arg(dstr)));
            }
        }
    }

    void DataLong::helper_build_clear_data_()
    {
        // data sequence headers
        headers.clear();
        indiv_datalist.clear();
    }

    void DataLong::helper_build_data_headers_([[maybe_unused]] const DataWide &data_wide, const str_table &input)
    {
#ifndef NDEBUG
        if (data_wide.is_empty())
            msg::stop(
                "Tried to define data long format before wide format.");
#endif

        headers = read::extract_keys(input);
        // check for individual IDs
        if (!contains(headers, "id"))
            msg::stop(
                "Individual identifier column " + msg::arg("id") + " was missing.");
        // remove 'id'
        headers = headers - std::vector<std::string>({"id"});

        // add optional keys to registry
        for (const std::string &datastring : headers)
            khold_.add(datastring);
    }
}
