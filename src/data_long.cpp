/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/nan.h"
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "data_long.h"
#include "data_wide.h"
#include "household.h"
#include "individual.h"
#include "tools/data_reports.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  *** DataLong : struct ***  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

DataLong::DataLong(const std::vector<std::string> &keys,
                   DataWide *ptr)
    : pRegistry(std::make_unique<hsp::key::Registry>(keys)),
      pDataWide_(ptr) {}
DataLong::~DataLong() = default;

void DataLong::build(const hsp::str_table &input)
{
    pRegistry->clear();

    if (input.empty())
        return;

    helper_build_clear_data_();
    helper_build_data_headers_(input);

    // data-long allocation
    size_t n_indiv = pDataWide_->report_indiv.total;
    indiv_datalist.resize(n_indiv);

    // individual map for id -> index
    hsp::dict<size_t> id_to_index = {};
    size_t c = 0;
    for (const Household &house : *pDataWide_)
        for (const Individual &member : house)
            id_to_index.emplace(member.id(), c++);

    // load individual data list
    size_t n_rows = hsp::read::csv_n_rows(input);
    for (size_t r = 0; r < n_rows; ++r)
    {
        const std::string &id = input.at("id")[r];
        if (hsp::na::is(id))
            throw std::invalid_argument(
                "Invalid id found in row n." + std::to_string(r + 2) + ".");
        // UPDATE: these errors should be more self-explicative
        if (!id_to_index.count(id))
            throw std::invalid_argument(
                "Id '" + id + "' from data long format had no match in loaded data wide format.");

        // append (since vector is initialized to empty)
        auto &idat = indiv_datalist[id_to_index.at(id)];
        idat.resize(pRegistry->size());
        for (size_t d = 0; d < pRegistry->size(); ++d)
        {
            const std::string &datastring = pRegistry->as_string(d);
            idat[d].emplace_back(
                hsp::str::tod(input.at(datastring)[r], "data long format '" + datastring + "'"));
        }
    }

    // assign pointers to relative data-long to all individuals
    c = 0;
    for (Household &house : *pDataWide_)
        for (Individual &member : house)
            member.pdatalong_ = &indiv_datalist.at(c++);
}

void DataLong::helper_build_clear_data_()
{
    // data sequence headers
    headers.clear();
    indiv_datalist.clear();
}

void DataLong::helper_build_data_headers_(const hsp::str_table &input)
{
    if (pDataWide_->is_empty())
        throw std::invalid_argument(
            "Tried to define data long format before wide format.");

    headers = hsp::read::extract_keys(input);
    // check for individual IDs
    if (std::find(headers.begin(), headers.end(), "id") == headers.end())
        throw std::invalid_argument(
            "Individual identifier column 'id' was missing.");
    // remove 'id'
    headers = headers - std::vector<std::string>({"id"});

    // pedantic check for declared data keys
    hsp::str::assert_no_missing_entries(headers, pRegistry->get_all());

    // add optional keys to registry
    for (const std::string &datastring : headers)
        if (!pRegistry->count(datastring))
            pRegistry->add(datastring);
}
