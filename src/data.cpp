/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
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
#include "tools/read.h"
#include "tools/status.h"
#include "tools/string.h"
#include "tools/types.h"
#include "tools/vector.h"

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

    void Household::add_indiv(IndivIdx index, StatusLvl status)
    {
        const MemberIdx m = hhsize_++;
        indivs_.emplace_back();

        Individual &last = indivs_.back();
        last.indiv_index_ = index;
        last.member_index_ = m;

        // new member added as if missing status
        last.status_ = status::MISSING;
        status_vec_.push_back(m);
        if (!status::is_missing(status))
            set_new_status_(last, status);
    }

    void Household::set_new_status_(Individual &indiv, StatusLvl snew)
    {
        const MemberIdx i = indiv.member_index_;
        size_t cold = status::category(indiv.status());
        size_t cnew = status::category(snew);

        indiv.status_ = snew; // NOTE: might just change within infected status
        if (cold == cnew)
            return;

        // segment starts [infected | escaped | missing]
        const size_t starts[] = {0, ninf_, static_cast<size_t>(ninf_ + nsus_)};
        // find old position given the old segment
        // NOTE: segments are small, raw scan profiled better
        size_t pos = starts[cold];
        while (status_vec_[pos] != i)
            ++pos;

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

    void DataWide::init_cache_key_holders(const std::vector<std::string> &member_keys,
                                          const std::vector<std::string> &member_vec_keys,
                                          const std::vector<std::string> &double_vec_keys,
                                          const std::vector<std::string> &double_mtx_keys)
    {
        khold_member_val.init(member_keys);
        khold_member_vec.init(member_vec_keys);
        khold_double_vec.init(double_vec_keys);
        khold_double_mtx.init(double_mtx_keys);
    }

    void DataWide::build(const str_table &input)
    {
        khold_.clear_optional();
        khold_member_val.clear_optional();
        khold_member_vec.clear_optional();
        khold_double_vec.clear_optional();
        khold_double_mtx.clear_optional();

        helper_build_clear_data_();        // drop previous data and lists
        helper_build_data_headers_(input); // extract headers
        report_indiv.total = read::csv_n_rows(input);
        if (report_indiv.total > std::numeric_limits<IndivIdx>::max())
            msg::stop(
                "Data exceeded the maximum of " +
                std::to_string(std::numeric_limits<IndivIdx>::max()) + " individuals.");

        // pre-allocate external data storage
        indiv_data_.resize(report_indiv.total);
        indiv_ids_.resize(report_indiv.total);

        // store households
        for (size_t i = 0; i < report_indiv.total; ++i)
        {
            std::string id = std::to_string(i + 1);                      // default id
            std::string hhid = std::to_string(house_vector_.size() + 1); // default hhid
            StatusLvl status = status::MISSING;

            if (input.count("id"))
            {
                id = input.at("id")[i]; // user id
                if (na::is(id))
                    msg::stop(
                        "Individual in row n." + std::to_string(i + 1) +
                        " had missing " + msg::field("id") + ".");
            }
            if (input.count("hhid"))
            {
                hhid = input.at("hhid")[i]; // user hhid
                if (na::is(hhid))
                    msg::stop(
                        "Individual in row n." + std::to_string(i + 1) +
                        " had missing " + msg::field("hhid") + ".");
            }
            if (input.count("status"))
            {
                // note: NA status is a missing status
                if (!na::is(input.at("status")[i]))
                {
                    const int s = str::toi(input.at("status")[i], "individual status");
                    if (s < 0)
                        msg::stop(
                            "Individual in row n." + std::to_string(i + 1) +
                            " had a negative status.\nUse " + msg::val("NA") +
                            " to declare a missing status.");
                    if (s > status::MAX_INFECTED)
                        msg::stop(
                            "Data exceeded the maximum of " +
                            std::to_string(status::MAX_INFECTED) + " infected statuses.");
                    status = static_cast<StatusLvl>(s);
                }
            }

            // read individual data into external storage
            IndivWide_ &idata = indiv_data_[i];
            idata.status = status;
            idata.vec.resize(khold_.size());
            // per-column read-cell defaults to the imported raw data value
            idata.readcell.resize(idata.vec.size());
            for (size_t d = 0; d < idata.vec.size(); ++d)
            {
                const std::string &dname = khold(d);
                double &x = idata.vec[d];

                // placeholder for model-required columns absent at import
                // NOTE: the check is deferred to the R-side wall (see build_full_)
                x = na::DOUBLE;

                if (input.count(dname))
                    x = str::tod(input.at(dname)[i], "data " + msg::arg(dname));

                idata.readcell[d].name(dname, ErrType::Data);
                idata.readcell[d].init_numeric(x);

                // a missing cell carries no verdict: absence is reported at
                // the run wall, not as a violated requirement
                if (na::is(x))
                    continue;
                try
                {
                    assert_requests(d, idata.readcell[d]);
                }
                catch (const std::exception &e)
                {
                    msg::stop(
                        std::string(e.what()) + "\nCheck " + idata.readcell[d].str() +
                        " for individual in row n." + std::to_string(i + 1) + ".");
                }
            }

            // if household index doesn't exist, add new household
            if (!hhid_to_index_.count(hhid))
            {
                // map data index to household number
                size_t H = house_vector_.size();
                if (H > std::numeric_limits<HouseIdx>::max())
                    msg::stop(
                        "Data exceeded the maximum of " +
                        std::to_string(std::numeric_limits<HouseIdx>::max()) + " households.");
                hhid_to_index_.emplace(hhid, H);
                // create a new empty household
                house_vector_.emplace_back();
                house_vector_.back().house_index_ = static_cast<HouseIdx>(H);
                hhids_.emplace_back(hhid);
            }
            // extract the household reference
            Household &house = house_vector_.at(hhid_to_index_.at(hhid));
            if (house.size() == std::numeric_limits<MemberIdx>::max())
                msg::stop(
                    "Household " + msg::val(hhid) + " exceeded the maximum of " +
                    std::to_string(std::numeric_limits<MemberIdx>::max()) + " members.");

            // assign the individual the respective household
            indiv_ids_[i] = id;
            house.add_indiv(static_cast<IndivIdx>(i), idata.status);
        }

        // columns known at import: model declared keys plus the user ones. Columns
        // produced by an observation are appended later (see reset_readcells).
        n_imported_ = khold_.size();

        // update report for data
        build_helper_fill_report_data_();
    }

    std::vector<std::string> DataWide::imported_cols() const
    {
        std::vector<std::string> all = khold_.all();
        all.resize(n_imported_);
        return all;
    }

    void DataWide::reset_readcells(const std::vector<std::string> &produced_cols)
    {
        for (const std::string &name : produced_cols)
            khold_.add(name);
        const size_t total = khold_.size();

        for (IndivWide_ &idata : indiv_data_)
        {
            idata.vec.resize(total, na::DOUBLE);
            idata.readcell.resize(total);
            for (size_t d = 0; d < total; ++d)
            {
                idata.readcell[d].name(khold(d), ErrType::Data);
                idata.readcell[d].init_numeric(idata.vec[d]);
            }
        }
    }

    void DataWide::finalize()
    {
        for (IndivWide_ &idata : indiv_data_)
            for (Attribute &attr : idata.readcell)
                attr.finalize();
    }

    size_t DataWide::index_of_(std::string_view hhid) const
    {
        auto it = hhid_to_index_.find(hhid);
        if (it == hhid_to_index_.end())
            msg::stop(
                "Index " + msg::val(std::string(hhid)) + " did not match any household in the data.");
        return it->second;
    }

    void DataWide::helper_build_clear_data_()
    {
        house_vector_.clear();
        hhid_to_index_.clear();
        indiv_data_.clear();
        indiv_ids_.clear();
        hhids_.clear();
        report_indiv = ReportIndividual();
        report_house = ReportHousehold();
        report_hh_size = ReportHouseSizes();
    }

    void DataWide::helper_build_data_headers_(const str_table &input)
    {
        std::vector<std::string> headers = read::extract_keys(input);

        // assess the presence of house IDs
        if (!vect::contains(headers, "hhid"))
            msg::warn("Ignored household structure.\n" +
                      msg::alert("Missing data column " + msg::field("hhid") + "."));

        // assess individual IDs are unique
        if (vect::contains(headers, "id") &&
            input.at("id").size() != str::unique(input.at("id")).size())
            msg::stop("Provided individual IDs are not unique.");

        // remove standard headers
        vect::remove(headers, {"id", "hhid", "status"});

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
                              { return indiv_less_(house.at(i), house.at(j)); });

                    // equal house sizes share seeding
                    size_t seed = std::hash<size_t>{}(house.size());
                    // combine hashes over sorted members
                    for (size_t i : order)
                    {
                        const Individual &indiv = house.at(i);
                        hash_combine_(seed, std::hash<StatusLvl>{}(indiv.status()));
                        for (size_t d : data_list_)
                            hash_combine_(seed, std::hash<double>{}(data_wide_.data(d, indiv)));
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
            /// NOTE: a missing cell is a class of its own, ordered last, so the
            /// equivalence it induces stays transitive (see CovariateLess_).
            bool indiv_less_(const Individual &a, const Individual &b) const
            {
                if (a.status() != b.status())
                    return a.status() < b.status();
                for (size_t d : data_list_)
                {
                    const double x = data_wide_.data(d, a), y = data_wide_.data(d, b);
                    const bool na_x = na::is(x), na_y = na::is(y);
                    if (na_x != na_y)
                        return na_y;
                    if (!na_x && x != y)
                        return x < y;
                }
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
                at(H).eq_class_ = static_cast<HouseIdx>(refer_restrict);

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
            at(H).eq_class_ = static_cast<HouseIdx>(H);
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

            // member value
            size_t s_member = khold_member_val.size();
            curr.member.assign(s_member, 0);
            prev.member.assign(s_member, 0);

            // member vector
            size_t s_member_vec = khold_member_vec.size();
            curr.member_vector.resize(s_member_vec);
            prev.member_vector.resize(s_member_vec);
            for (size_t e = 0; e < s_member_vec; ++e)
            {
                curr.member_vector[e].assign(hh_size, 0);
                prev.member_vector[e].assign(hh_size, 0);
            }

            // double vector
            size_t s_double_vec = khold_double_vec.size();
            curr.double_vector.resize(s_double_vec);
            prev.double_vector.resize(s_double_vec);
            for (size_t e = 0; e < s_double_vec; ++e)
            {
                curr.double_vector[e].assign(hh_size, 0);
                prev.double_vector[e].assign(hh_size, 0);
            }

            // double matrix
            size_t s_double_mtx = khold_double_mtx.size();
            curr.double_matrix.resize(s_double_mtx);
            prev.double_matrix.resize(s_double_mtx);
            for (size_t e = 0; e < s_double_mtx; ++e)
            {
                curr.double_matrix[e].assign(hh_size, std::vector<double>(hh_size, 0.));
                prev.double_matrix[e].assign(hh_size, std::vector<double>(hh_size, 0.));
            }
        }
    }

    void DataWide::save_member_cache(const Household &house,
                                     const Individual &indiv)
    {
        size_t i = Index::member(indiv);
        Cached_::State_ &prev = cached_.previous[house.eq_class_];
        Cached_::State_ &curr = cached_.current[house.eq_class_];
        prev.member = curr.member;

        for (size_t e = 0; e < curr.member_vector.size(); ++e)
            prev.member_vector[e][i] = curr.member_vector[e][i];

        for (size_t e = 0; e < curr.double_vector.size(); ++e)
            prev.double_vector[e][i] = curr.double_vector[e][i];

        for (size_t e = 0; e < curr.double_matrix.size(); ++e)
            for (size_t j = 0; j < house.size(); ++j)
            {
                prev.double_matrix[e][i][j] = curr.double_matrix[e][i][j];
                if (j != i)
                    prev.double_matrix[e][j][i] = curr.double_matrix[e][j][i];
            }
    }

    void DataWide::recover_member_old_cached(const Household &house,
                                             const Individual &indiv)
    {
        size_t i = Index::member(indiv);
        Cached_::State_ &prev = cached_.previous[house.eq_class_];
        Cached_::State_ &curr = cached_.current[house.eq_class_];
        curr.member.swap(prev.member);

        for (size_t e = 0; e < prev.member_vector.size(); ++e)
            curr.member_vector[e][i] = prev.member_vector[e][i];

        for (size_t e = 0; e < prev.double_vector.size(); ++e)
            curr.double_vector[e][i] = prev.double_vector[e][i];

        for (size_t e = 0; e < prev.double_matrix.size(); ++e)
        {
            curr.double_matrix[e][i].swap(prev.double_matrix[e][i]);
            for (size_t j = 0; j < house.size(); ++j)
                if (j != i)
                    curr.double_matrix[e][j][i] = prev.double_matrix[e][j][i];
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
            for (Individual &indiv : house)
                house.set_new_status_(indiv, indiv_data_[Index::of(indiv)].status);
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
            for (const Individual &indiv : house)
                id_to_index.emplace(data_wide.id(indiv), Index::of(indiv));

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

    void DataLong::helper_build_data_headers_(const DataWide &data_wide,
                                              const str_table &input)
    {
        if (data_wide.is_empty())
            msg::stop(
                "Tried to define data long format before wide format.");

        headers = read::extract_keys(input);
        // check for individual IDs
        if (!vect::contains(headers, "id"))
            msg::stop(
                "Individual identifier column " + msg::arg("id") + " was missing.");
        // remove 'id'
        vect::remove(headers, {"id"});

        // add optional keys to registry
        for (const std::string &datastring : headers)
            khold_.add(datastring);
    }
}
