/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/loglike.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "data_wide.h"
#include "function.h"
#include "tools/data_reports.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Household : class ***  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace hsp
{
    void Household::init_member(size_t index,
                                const std::string &id,
                                const Indiv::Wide *data_ptr)
    {
        size_t hhindex = size();
        members_.emplace_back();

        Individual &last = members_.back();
        last.index_ = index;
        last.hhindex_ = hhindex;
        last.id_ = id;
        last.status_ = data_ptr->status;
        last.pwide_ = data_ptr;

        // update status lists within the household
        if (last.is_infected_status())
            infected_list_.emplace_back(hhindex);
        else if (last.is_susceptible_status())
            healthy_list_.emplace_back(hhindex);
        else
            missing_list_.emplace_back(hhindex);
    }

    void Household::helper_set_status_(Individual &member, int new_status)
    {
        int old_status = member.status();
        if (new_status == old_status)
            return;

        member.status_ = new_status; // change member status
        size_t i = member.hhindex(); // get within house index

        if (status::is_infected(old_status))
            infected_list_.erase(std::find(infected_list_.begin(), infected_list_.end(), i));
        else if (status::is_susceptible(old_status))
            healthy_list_.erase(std::find(healthy_list_.begin(), healthy_list_.end(), i));
        else
            missing_list_.erase(std::find(missing_list_.begin(), missing_list_.end(), i));

        if (status::is_infected(new_status))
            infected_list_.emplace_back(i);
        else if (status::is_susceptible(new_status))
            healthy_list_.emplace_back(i);
        else
            missing_list_.emplace_back(i);
    }

    void Household::raise_set_non_infec_status_(int new_status) const
    {
        throw std::invalid_argument(
            "Tried to set as infected status " + std::to_string(new_status) +
            " to " + err_name() + ".\nInfected status must be an integer greater than 0.");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    DataWide::DataWide() = default;
    DataWide::~DataWide() = default;

    void DataWide::init_registry(
        const std::vector<std::string> &keys,
        const std::vector<std::string> &size_keys,
        const std::vector<std::string> &size_vec_keys,
        const std::vector<std::string> &double_vec_keys,
        const std::vector<std::string> &double_mtx_keys)
    {
        pKeys = std::make_unique<KeyHolder>(keys);
        pKeys_size_val = std::make_unique<KeyHolder>(size_keys);
        pKeys_size_vec = std::make_unique<KeyHolder>(size_vec_keys);
        pKeys_double_vec = std::make_unique<KeyHolder>(double_vec_keys);
        pKeys_double_mtx = std::make_unique<KeyHolder>(double_mtx_keys);
    }

    void DataWide::build(const str_table &input)
    {
        pKeys->clear();
        pKeys_size_val->clear();
        pKeys_size_vec->clear();
        pKeys_double_vec->clear();
        pKeys_double_mtx->clear();

        if (input.empty())
            return;

        // initialize multiplicity-relevant data entries to only model dependent
        set_multiplicity_data(pKeys->get_all());

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
                    throw std::invalid_argument(
                        "Individual in row n." + std::to_string(i + 2) + " had missing 'id'.");
            }
            if (input.count("house_id"))
            {
                house_id = input.at("house_id")[i]; // user house_id
                if (na::is(house_id))
                    throw std::invalid_argument(
                        "Individual in row n." + std::to_string(i + 2) + " had missing 'house_id'.");
            }
            if (input.count("status"))
            {
                // note: NA status is a missing status
                if (!na::is(input.at("status")[i]))
                    status = str::toi(input.at("status")[i], "individual status");
            }

            // read individual data into external storage
            Indiv::Wide &idata = indiv_data_[i];
            idata.status = status;
            idata.vec.resize(pKeys->size());
            for (size_t d = 0; d < idata.vec.size(); ++d)
            {
                const std::string datastring = pKeys->as_string(d);
                idata.vec[d] = str::tod(input.at(datastring)[i], "data '" + datastring + "'");
            }

            // if household index doesn't exist, add new household
            if (!hhid_to_index_.count(house_id))
            {
                // map data index to household number
                hhid_to_index_.emplace(house_id, house_vector_.size());
                // create a new empty household
                house_vector_.emplace_back();
                Household &last = house_vector_.back();
                last.house_id_ = house_id;
            }
            // extract the household reference
            Household &house = house_vector_.at(hhid_to_index_.at(house_id));

            // assign the individual the respective household
            house.init_member(i, id, &idata);
        }

        // update report for data
        build_helper_fill_report_data_();
    }

    size_t DataWide::house_id_to_index_(std::string_view house_id) const
    {
        auto it = hhid_to_index_.find(house_id);
        if (it != hhid_to_index_.end())
            return it->second;
        throw std::invalid_argument(
            "Index '" + std::string(house_id) + "' did not match any household in the data.");
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
        if (std::find(headers.begin(), headers.end(), "house_id") == headers.end())
            msg::warning("Ignored household structure. Missing data column 'house_id'.");

        // assess individual IDs are unique
        if (std::find(headers.begin(), headers.end(), "id") != headers.end() &&
            input.at("id").size() != str::unique(input.at("id")).size())
            msg::error("Provided individual IDs are not unique.");

        // remove standard headers
        headers = headers - std::vector<std::string>({"id", "house_id", "status"});

        // pedantic check for declared data keys
        str::assert_no_missing_entries(headers, pKeys->get_all());

        // add optional keys to registry
        for (const std::string &datastring : headers)
            pKeys->add(datastring);
    }

    void DataWide::build_helper_fill_report_data_()
    {
        report_house.total = house_vector_.size();

        // set report counts
        for (const Household &house : house_vector_)
        {
            // observed quantities in-house
            size_t house_size = house.size();
            size_t n_infec = house.n_infected_status();
            size_t n_miss = house.n_missing_status();
            size_t n_susc = house.n_susceptible_status();

            report_indiv.infectious += n_infec;
            report_indiv.missing += n_miss;
            report_indiv.susceptible += n_susc;

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
            UPDATE_REPORT_SIZE(report_hh_size.max_cases, house_size, <);
        }
    }

    namespace
    {
        /// Using golden ratio fractional part to ensure bit distribution.
        /// Ref:
        /// https://www.boost.org/doc/libs/latest/libs/container_hash/doc/html/hash.html#notes_hash_combine
        inline void hash_combine(size_t &seed, size_t h)
        {
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
            seed ^= h + 0x9e3779b97f4a7c15ULL + (seed << 12) + (seed >> 4);
#else
            seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
#endif
        }
    }

    void DataWide::set_multiplicity_data(const std::vector<std::string> &data_list)
    {
        for (size_t d : pKeys->from_string(data_list))
            multiplicity_data_.emplace(d);
    }

    void DataWide::init_multiplicity(bool is_any_augm)
    {
        mode_ = Mode::Unique;
        multiplicity_.clear();

        // initialize household indexes
        for (Household &house : house_vector_)
            house.ridx_ = hhid_to_index_[house.house_id()];

        // NOTE: multiplicity incompatible with sampling augmented data
        if (house_vector_.empty() || is_any_augm)
            return;

        // List of data relevant for multiplicity groups
        // FIXME: not all data contribute to this!
        // probably just group rates and model dependencies
        const std::unordered_set<size_t> data_list = multiplicity_data_;

        // comparison functor: sort members by (status, data[0], data[1], ...)
        auto member_less = [data_list](const Individual &a, const Individual &b)
        {
            if (a.status() != b.status())
                return a.status() < b.status();
            for (size_t d : data_list)
                if (a.data(d) != b.data(d))
                    return a.data(d) < b.data(d);
            return false;
        };

        // compute multiplicity-relevant hash for each household
        std::map<size_t, std::vector<size_t>> hash_groups;

        for (size_t H = 0; H < report_house.total; ++H)
        {
            const Household &house = at(H);

            // sorted member indices
            std::vector<size_t> order(house.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&house, &member_less](size_t i, size_t j)
                      { return member_less(house.member(i), house.member(j)); });

            // combine hashes over sorted members
            size_t seed = std::hash<size_t>{}(house.size());
            for (size_t i : order)
            {
                const Individual &member = house.member(i);
                hash_combine(seed, std::hash<int>{}(member.status()));
                for (size_t d : data_list)
                    hash_combine(seed, std::hash<double>{}(member.data(d)));
            }
            hash_groups[seed].emplace_back(H);
        }

        // skip if all unique households
        if (house_vector_.size() == hash_groups.size())
            return;

        // build multiplicity structures
        mode_ = Mode::Multiple;
        multiplicity_.reserve(hash_groups.size());

        for (const auto &[hash, group] : hash_groups)
        {
            const size_t refer_restrict = multiplicity_.size();
            for (size_t H : group)
                at(H).ridx_ = refer_restrict;

            multiplicity_.emplace_back(group);
        }
        multiplicity_.shrink_to_fit();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : cache attributes  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void DataWide::allocate_all_cache()
    {
        cached_.previous.resize(restricted_size());
        cached_.current.resize(restricted_size());

        for (size_t R = 0; R < restricted_size(); ++R)
        {
            Household &house = restricted(R);
            size_t hh_size = house.size();
            Cached_::State_ &prev = cached_.previous[R];
            Cached_::State_ &curr = cached_.current[R];

            // size_t
            size_t s_size = pKeys_size_val->size();
            curr.size_val.assign(s_size, 0);
            prev.size_val.assign(s_size, 0);

            // size_t vector
            size_t s_size_vec = pKeys_size_vec->size();
            curr.size_vec.resize(s_size_vec);
            prev.size_vec.resize(s_size_vec);
            for (size_t e = 0; e < s_size_vec; ++e)
            {
                curr.size_vec[e].assign(hh_size, 0);
                prev.size_vec[e].assign(hh_size, 0);
            }

            // double vector
            size_t s_double_vec = pKeys_double_vec->size();
            curr.dbl_vec.resize(s_double_vec);
            prev.dbl_vec.resize(s_double_vec);
            for (size_t e = 0; e < s_double_vec; ++e)
            {
                curr.dbl_vec[e].assign(hh_size, 0);
                prev.dbl_vec[e].assign(hh_size, 0);
            }

            // double matrix
            size_t s_double_mtx = pKeys_double_mtx->size();
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
        Cached_::State_ &prev = cached_.previous[house.ridx_];
        Cached_::State_ &curr = cached_.current[house.ridx_];
        prev.size_val = curr.size_val;

        for (size_t e = 0; e < curr.size_vec.size(); ++e)
            prev.size_vec[e][i] = curr.size_vec[e][i];

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
        Cached_::State_ &prev = cached_.previous[house.ridx_];
        Cached_::State_ &curr = cached_.current[house.ridx_];
        curr.size_val.swap(prev.size_val);

        for (size_t e = 0; e < prev.size_vec.size(); ++e)
            curr.size_vec[e][i] = prev.size_vec[e][i];

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

    void DataWide::assert_finite_log_lik() const
    {
        for (size_t R = 0; R < restricted_size(); ++R)
        {
            const Household &house = restricted(R);

            for (size_t i : house.infected_hhindexes())
                if (house.member(i).paugm_->log_lik.is_inf())
                    throw std::runtime_error(
                        "Initial augmented log-likelihood computation returned infinity for " +
                        house.member(i).err_name() + ".");

            if (house.obs_log_lik_.is_inf())
                throw std::runtime_error(
                    "Initial log-likelihood computation returned infinity for " + house.err_name() + ".");
        }
    }

    double DataWide::house_log_lik(size_t R) const
    {
        const Household &house = restricted(R);
        double x = house.obs_log_lik();
        for (size_t i : house.infected_hhindexes())
            x += house.member(i).paugm_->log_lik.value();
        return x;
    }

    void DataWide::assert_non_empty() const
    {
        if (report_house.total == 0)
            throw std::length_error(
                "Unable to proceed: no data was found.");
    }

    void DataWide::assert_not_all_missing_status() const
    {
        if (report_indiv.total == report_indiv.missing)
            throw std::length_error(
                "All individuals in the data had missing status.");
    }

    void DataWide::assert_no_missing_status() const
    {
        if (report_indiv.missing > 0)
            throw std::length_error(
                "Individual with missing status are incompatible with current model.");
    }

    void DataWide::recover_all_status() noexcept
    {
        for(Household &house: house_vector_)
            for(Individual &member : house)
                house.helper_set_status_(member, member.pwide_->status);
    }
}
