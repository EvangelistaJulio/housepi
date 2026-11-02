/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "data.h"
#include "data/dataset.h"
#include "data/household.h"
#include "data/individual.h"
#include "data/long_dataset.h"
#include "data/reports.h"
#include "data/state.h"
#include "data/status.h"
#include "function.h"
#include "model.h"

namespace housepi
{
    //////////////////////////////  Data  //////////////////////////////

    Data::Data(const ModelBase &model)
        : store_(std::make_shared<Store_>(model))
    {
        // the transmission model reads its data keys per member
        set_eq_class_keys(Source::Model, model.data_keys());
    }

    Data::Store_::Store_(const ModelBase &model)
        : dataset(model.data_keys(),
                  model.house_member_keys(),
                  model.house_member_vec_keys(),
                  model.house_double_vec_keys(),
                  model.house_double_mtx_keys()),
          long_dataset(model.long_keys()) {}

    void Data::build_dataset(const str_table &input)
    {
        dataset().build(input);
        state.reset_readcells(dataset());
        // definitions push their keys again over the new data
        for (Source source : {Source::GroupRate, Source::ContactRate, Source::Latent})
            set_eq_class_keys(source, {});
    }

    //////  Data : equivalence classes  //////

    namespace
    {
        /// Bits of the cell of column `d` for `indiv`, as compared in partition keys.
        std::uint64_t cell_bits_(const data::State &state, size_t d, const Individual &indiv)
        {
            const double x = state.readcell(d, indiv).value();
            std::uint64_t bits;
            std::memcpy(&bits, &x, sizeof bits);
            return bits;
        }

        /// Whether member `a` sorts before `b`, by status then cell bits of columns `data_list`.
        bool indiv_less_(const std::vector<size_t> &data_list,
                         const data::State &state,
                         const Individual &a,
                         const Individual &b)
        {
            if (a.status() != b.status())
                return a.status() < b.status();
            for (size_t d : data_list)
            {
                const std::uint64_t x = cell_bits_(state, d, a);
                const std::uint64_t y = cell_bits_(state, d, b);
                if (x != y)
                    return x < y;
            }
            return false;
        }

        /// Member indices of `house` sorted by indiv_less_().
        std::vector<MemberIdx> sorted_members_(const std::vector<size_t> &data_list,
                                               const data::State &state,
                                               const Household &house)
        {
            std::vector<MemberIdx> order(house.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&](MemberIdx i, MemberIdx j)
                      { return indiv_less_(data_list, state, house.at(i), house.at(j)); });
            return order;
        }

        /// Households of each equivalence class, keyed on their sorted members.
        /// * key: status, then cell bits of columns `data_list`, per sorted member;
        /// * NA cell -> equivalent to NA only.
        std::map<std::vector<std::uint64_t>, std::vector<size_t>> partition_houses_(
            const std::vector<size_t> &data_list,
            const data::State &state,
            const data::Dataset &dataset)
        {
            std::map<std::vector<std::uint64_t>, std::vector<size_t>> output;
            for (size_t H = 0; H < dataset.report_house.total; ++H)
            {
                const Household &house = dataset.at(H);
                std::vector<std::uint64_t> key;
                for (MemberIdx i : sorted_members_(data_list, state, house))
                {
                    const Individual &indiv = house.at(i);
                    key.push_back(indiv.status());
                    for (size_t d : data_list)
                        key.push_back(cell_bits_(state, d, indiv));
                }
                output[std::move(key)].emplace_back(H);
            }
            return output;
        }
    }

    void Data::init_eq_classes(bool force_singleton)
    {
        const size_t total = dataset().report_house.total;
        std::vector<size_t> &keys = state.eq_class_keys_;
        keys.clear();
        for (const std::vector<std::string> &names : eq_class_keys_)
            for (const std::string &name : names)
            {
                // a column absent from the data tells no household apart
                const size_t d = dataset().find(name);
                if (!na::is(d) && !vect::contains(keys, d))
                    keys.push_back(d);
            }

        std::map<std::vector<std::uint64_t>, std::vector<size_t>> partition;
        if (!force_singleton)
            partition = partition_houses_(keys, state, dataset());

        std::vector<data::State::EquivalenceClass_> &classes = state.eq_classes_;
        classes.clear();
        state.class_of_.resize(total);

        // no two households equivalent -> one class per household
        if (force_singleton || partition.size() == total)
        {
            classes.resize(total);
            for (size_t H = 0; H < total; ++H)
            {
                classes[H].rep_ = H;
                state.class_of_[H] = static_cast<EqClassIdx>(H);
            }
            return;
        }

        // first household of each class as its representative
        classes.reserve(partition.size());
        for (const auto &[key, houses] : partition)
        {
            for (size_t H : houses)
                state.class_of_[H] = static_cast<EqClassIdx>(classes.size());
            classes.emplace_back(houses);
        }
    }

    std::vector<MemberIdx> Data::map_rep_members(const Household &house) const
    {
        const std::vector<size_t> &keys = state.eq_class_keys_;
        const Household &rep = eq_class_rep(state.eq_class(house));
        const std::vector<MemberIdx> rep_order = sorted_members_(keys, state, rep);
        const std::vector<MemberIdx> house_order = sorted_members_(keys, state, house);

        std::vector<MemberIdx> output(house.size());
        for (size_t k = 0; k < house.size(); ++k)
            output[rep_order[k]] = house_order[k];
        return output;
    }

    //////////////////////////////  Household  //////////////////////////////

    void Household::add_indiv(IndivIdx index, StatusLvl status)
    {
        const MemberIdx m = hh_size_++;
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

        indiv.status_ = snew; // may stay within the infected category
        if (cold == cnew)
            return;

        // segment starts [infected | escaped | missing]
        const size_t starts[] = {0, ninf_, static_cast<size_t>(ninf_ + nsus_)};
        // find old position given the old segment
        // linear scan: segments are small
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

}

namespace housepi::data
{
    //////////////////////////////  data::Dataset  //////////////////////////////

    void Dataset::build(const str_table &input)
    {
        clear_optional_();

        helper_build_clear_data_();        // drop previous data and lists
        helper_build_data_headers_(input); // extract headers
        report_indiv.total = count_rows(input);
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
            std::string id = std::to_string(i + 1);                       // default id
            std::string hh_id = std::to_string(house_vector_.size() + 1); // default hh_id
            StatusLvl status = status::MISSING;                           // default status
            // individual context for reading errors
            std::string err_ctx = " for individual in row n." + std::to_string(i + 1);

            if (input.count("id"))
            {
                id = input.at("id")[i]; // user id
                if (na::is(id))
                    msg::stop("Missing identifier " + msg::field("id") + err_ctx + ".");
                err_ctx = " for individual " + msg::val(id);
            }
            if (input.count("hh_id"))
            {
                hh_id = input.at("hh_id")[i]; // user hh_id
                if (na::is(hh_id))
                    msg::stop("Missing household identifier " + msg::field("hh_id") + err_ctx + ".");
            }
            if (input.count("status"))
            {
                // note: NA status is a missing status
                if (!na::is(input.at("status")[i]))
                {
                    const int s = str::toi(input.at("status")[i], "individual status");
                    if (s < 0)
                        msg::stop(
                            "Negative status found " + err_ctx +
                            "\nUse " + msg::val("NA") + " to declare a missing status.");
                    if (s > status::MAX_INFECTED)
                        msg::stop(
                            "Data exceeded the maximum of " +
                            std::to_string(status::MAX_INFECTED) + " infected statuses.");
                    status = static_cast<StatusLvl>(s);
                }
            }

            // read individual data into external storage
            IndivRow_ &idata = indiv_data_[i];
            idata.status = status;
            idata.vec.resize(size());
            for (size_t d = 0; d < idata.vec.size(); ++d)
            {
                const std::string &dname = consult(d);
                double &x = idata.vec[d];

                // model-required column absent at import: missing, reported at the run wall
                x = na::DOUBLE;

                if (input.count(dname))
                    x = str::tod(input.at(dname)[i], "data " + msg::arg(dname));

                // a missing cell carries no verdict, its absence is reported at the run wall
                if (na::is(x))
                    continue;
                Attribute cell;
                cell.name(dname, ErrType::Data);
                cell.init_numeric(x);
                try
                {
                    assert_requests(d, cell);
                }
                catch (const std::exception &e)
                {
                    msg::stop(
                        std::string(e.what()) + "\nCheck " + cell.str() + err_ctx + ".");
                }
            }

            // if household index doesn't exist, add new household
            if (!hh_id_to_index_.count(hh_id))
            {
                // map data index to household number
                size_t H = house_vector_.size();
                if (H > std::numeric_limits<HouseIdx>::max())
                    msg::stop(
                        "Data exceeded the maximum of " +
                        std::to_string(std::numeric_limits<HouseIdx>::max()) + " households.");
                hh_id_to_index_.emplace(hh_id, H);
                // create a new empty household
                house_vector_.emplace_back();
                house_vector_.back().house_index_ = static_cast<HouseIdx>(H);
                hh_ids_.emplace_back(hh_id);
            }
            // extract the household reference
            Household &house = house_vector_.at(hh_id_to_index_.at(hh_id));
            if (house.size() == std::numeric_limits<MemberIdx>::max())
                msg::stop(
                    "Household " + msg::val(hh_id) + " exceeded the maximum of " +
                    std::to_string(std::numeric_limits<MemberIdx>::max()) + " members.");

            // assign the individual the respective household
            indiv_ids_[i] = id;
            house.add_indiv(static_cast<IndivIdx>(i), idata.status);
        }

        // columns known at import: model-declared keys plus the user ones
        // columns produced by an observation are appended later, see add_produced_cols
        n_imported_ = size();

        // update report for data
        build_helper_fill_report_data_();
    }

    std::vector<std::string> Dataset::imported_cols() const
    {
        return std::vector<std::string>(all().begin(), all().begin() + n_imported_);
    }

    void Dataset::add_produced_cols(const std::vector<std::string> &produced_cols)
    {
        for (const std::string &name : produced_cols)
            add_(name);
        for (IndivRow_ &idata : indiv_data_)
            idata.vec.resize(size(), na::DOUBLE);
    }

    size_t Dataset::index_of_(std::string_view hh_id) const
    {
        auto it = hh_id_to_index_.find(hh_id);
        if (it == hh_id_to_index_.end())
            msg::stop(
                "Index " + msg::val(std::string(hh_id)) + " did not match any household in the data.");
        return it->second;
    }

    void Dataset::helper_build_clear_data_()
    {
        house_vector_.clear();
        hh_id_to_index_.clear();
        indiv_data_.clear();
        indiv_ids_.clear();
        hh_ids_.clear();
        report_indiv = ReportIndividual();
        report_house = ReportHousehold();
    }

    void Dataset::helper_build_data_headers_(const str_table &input)
    {
        std::vector<std::string> headers = extract_keys(input);

        // assess the presence of house IDs
        if (!vect::contains(headers, "hh_id"))
            msg::warn("Ignored household structure.\n" +
                      msg::alert("Missing data column " + msg::field("hh_id") + "."));

        // assess individual IDs are unique
        if (vect::contains(headers, "id") &&
            input.at("id").size() != str::unique(input.at("id")).size())
            msg::stop("Provided individual IDs are not unique.");

        // remove standard headers
        vect::remove(headers, {"id", "hh_id", "status"});

        // add optional keys to registry
        for (const std::string &datastring : headers)
            add_(datastring);
    }

    void Dataset::build_helper_fill_report_data_()
    {
        report_house.total = house_vector_.size();

        // set report counts
        for (const Household &house : house_vector_)
        {
            // observed quantities in-house
            size_t n_infec = house.n_infected_status();
            size_t n_miss = house.n_missing_status();
            size_t n_esc = house.n_escaped_status();

            report_indiv.infectious += n_infec;
            report_indiv.missing += n_miss;
            report_indiv.escaped += n_esc;

            report_house.no_cases += (n_infec == 0);
            report_house.at_least_1_case += (n_infec > 0);
            report_house.at_least_2_cases += (n_infec > 1);
        }
    }

    //////  data::Dataset : flags  //////

    void Dataset::assert_not_all_missing_status() const
    {
        if (report_indiv.total == report_indiv.missing)
            msg::stop(
                "All individuals in the data had missing status.");
    }

    void Dataset::assert_no_missing_status() const
    {
        if (report_indiv.missing > 0)
            msg::stop(
                "Individual with missing status are incompatible with current model.");
    }

    void Dataset::request_everyone(std::vector<size_t> vec)
    {
        vect::append_unique(request_everyone_, vec);
    }

    bool Dataset::is_required_everyone(size_t d) const noexcept
    {
        return vect::contains(request_everyone_, d);
    }

    void Dataset::assert_required_everyone(const State &state,
                                           const std::vector<std::string> &live) const
    {
        for (size_t d : request_everyone_)
            if (!vect::contains(live, consult(d)))
                for (const Household &house : house_vector_)
                    for (const Individual &indiv : house)
                        if (na::is(state.readcell(d, indiv).value()))
                            msg::stop(
                                "Data " + msg::val(consult(d)) + " has no value for " + str(indiv) + " for " +
                                str(house) + ".\nThe model request it for every member, escaped ones included.");
    }

    void Dataset::recover_all_status() noexcept
    {
        for (Household &house : house_vector_)
            for (Individual &indiv : house)
                house.set_new_status_(indiv, indiv_data_[Index::of(indiv)].status);
    }

    //////////////////////////////  data::State  //////////////////////////////

    void State::reset_readcells(const Dataset &dataset)
    {
        const size_t total = dataset.size();
        readcell_.resize(dataset.report_indiv.total);
        for (const Household &house : dataset)
            for (const Individual &indiv : house)
            {
                std::vector<Attribute> &cells = readcell_[Index::of(indiv)];
                cells.resize(total);
                for (size_t d = 0; d < total; ++d)
                {
                    cells[d].name(dataset.consult(d), ErrType::Data);
                    cells[d].init_numeric(dataset.raw(d, indiv));
                }
            }
    }

    void State::finalize_readcells()
    {
        for (std::vector<Attribute> &cells : readcell_)
            for (Attribute &cell : cells)
                cell.finalize();
    }

    void State::allocate_transm_log_lik()
    {
        transm_log_lik_.assign(eq_classes_.size(), LogProba());
        evaluated_.assign(eq_classes_.size(), false);
    }

    void State::allocate_cache(const Dataset &dataset)
    {
        evaluated_.assign(eq_classes_.size(), false);
        previous_.resize(eq_classes_.size());
        current_.resize(eq_classes_.size());

        for (size_t R = 0; R < eq_classes_.size(); ++R)
        {
            size_t hh_size = dataset.at(eq_classes_[R].rep_).size();
            Cache_ &prev = previous_[R];
            Cache_ &curr = current_[R];

            // member value
            size_t s_member = dataset.keys_member_val.size();
            curr.member.assign(s_member, 0);
            prev.member.assign(s_member, 0);

            // member vector
            size_t s_member_vec = dataset.keys_member_vec.size();
            curr.member_vector.resize(s_member_vec);
            prev.member_vector.resize(s_member_vec);
            for (size_t e = 0; e < s_member_vec; ++e)
            {
                curr.member_vector[e].assign(hh_size, 0);
                prev.member_vector[e].assign(hh_size, 0);
            }

            // double vector
            size_t s_double_vec = dataset.keys_double_vec.size();
            curr.double_vector.resize(s_double_vec);
            prev.double_vector.resize(s_double_vec);
            for (size_t e = 0; e < s_double_vec; ++e)
            {
                curr.double_vector[e].assign(hh_size, 0);
                prev.double_vector[e].assign(hh_size, 0);
            }

            // double matrix
            size_t s_double_mtx = dataset.keys_double_mtx.size();
            curr.double_matrix.resize(s_double_mtx);
            prev.double_matrix.resize(s_double_mtx);
            for (size_t e = 0; e < s_double_mtx; ++e)
            {
                curr.double_matrix[e].assign(hh_size, std::vector<double>(hh_size, 0.));
                prev.double_matrix[e].assign(hh_size, std::vector<double>(hh_size, 0.));
            }
        }
    }

    void State::save_member_cache(const Household &house,
                                  const Individual &indiv)
    {
        size_t i = Index::member(indiv);
        Cache_ &prev = previous_[eq_class(house)];
        Cache_ &curr = current_[eq_class(house)];
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

    void State::recover_member_cache(const Household &house,
                                     const Individual &indiv)
    {
        size_t i = Index::member(indiv);
        Cache_ &prev = previous_[eq_class(house)];
        Cache_ &curr = current_[eq_class(house)];
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

    //////////////////////////////  data::LongDataset  //////////////////////////////

    void LongDataset::build(Dataset &dataset,
                            const str_table &input)
    {
        clear_optional_();

        if (input.empty())
            return;

        helper_build_clear_data_();
        helper_build_data_headers_(dataset, input);

        // data-long allocation
        size_t n_indiv = dataset.report_indiv.total;
        indiv_datalist.resize(n_indiv);

        // individual map for id -> index
        dict<size_t> id_to_index = {};
        for (const Household &house : dataset)
            for (const Individual &indiv : house)
                id_to_index.emplace(dataset.id(indiv), Index::of(indiv));

        // load individual data list
        size_t n_rows = count_rows(input);
        for (size_t r = 0; r < n_rows; ++r)
        {
            const std::string &id = input.at("id")[r];
            if (na::is(id))
                msg::stop(
                    "Invalid id found in row n." + std::to_string(r + 1) + ".");
            if (!id_to_index.count(id))
                msg::stop(
                    "Id " + msg::val(id) + " from the long dataset has no match in the dataset.");

            // append (since vector is initialized to empty)
            auto &idat = indiv_datalist[id_to_index.at(id)];
            idat.resize(size());
            for (size_t d = 0; d < size(); ++d)
            {
                const std::string &dstr = consult(d);
                idat[d].emplace_back(
                    str::tod(input.at(dstr)[r], "long dataset column " + msg::arg(dstr)));
            }
        }
    }

    void LongDataset::helper_build_clear_data_()
    {
        // data sequence headers
        headers.clear();
        indiv_datalist.clear();
    }

    void LongDataset::helper_build_data_headers_(const Dataset &dataset,
                                                 const str_table &input)
    {
        if (dataset.is_empty())
            msg::stop(
                "Long dataset attached before the dataset.");

        headers = extract_keys(input);
        // check for individual IDs
        if (!vect::contains(headers, "id"))
            msg::stop(
                "Individual identifier column " + msg::arg("id") + " was missing.");
        // remove 'id'
        vect::remove(headers, {"id"});

        // add optional keys to registry
        for (const std::string &datastring : headers)
            add_(datastring);
    }
}
