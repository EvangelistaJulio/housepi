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

namespace hsp
{
    DataLong::DataLong() = default;
    DataLong::~DataLong() = default;

    void DataLong::init_registry(const std::vector<std::string> &keys)
    {
        pKeys = std::make_unique<KeyHolder>(keys);
    }

    void DataLong::build(DataWide &data_wide, const str_table &input)
    {
        pKeys->clear_optional();

        if (input.empty())
            return;

        helper_build_clear_data_();
        helper_build_data_headers_(data_wide, input);

        // data-long allocation
        size_t n_indiv = data_wide.report_indiv.total;
        indiv_datalist.resize(n_indiv);

        // individual map for id -> index
        dict<size_t> id_to_index = {};
        size_t c = 0;
        for (const Household &house : data_wide)
            for (const Individual &member : house)
                id_to_index.emplace(member.id(), c++);

        // load individual data list
        size_t n_rows = read::csv_n_rows(input);
        for (size_t r = 0; r < n_rows; ++r)
        {
            const std::string &id = input.at("id")[r];
            if (na::is(id))
                throw std::invalid_argument(
                    "Invalid id found in row n." + std::to_string(r + 2) + ".");
            // UPDATE: these errors should be more self-explicative
            if (!id_to_index.count(id))
                throw std::invalid_argument(
                    "Id '" + id + "' from data long-format had no match in loaded data wide-format.");

            // append (since vector is initialized to empty)
            auto &idat = indiv_datalist[id_to_index.at(id)];
            idat.vec.resize(pKeys->size());
            for (size_t d = 0; d < pKeys->size(); ++d)
            {
                const std::string &datastring = pKeys->as_string(d);
                idat.vec[d].emplace_back(
                    str::tod(input.at(datastring)[r], "data long format '" + datastring + "'"));
            }
        }

        // assign pointers to relative data-long to all individuals
        // NOTE: adopted paradigm return empty on dangling pointer
        // UPDATE: just store indexes and access indiv_datalist by index...
        c = 0;
        for (Household &house : data_wide)
            for (size_t i = 0; i < house.size(); ++i, ++c)
                if (indiv_datalist.at(c).vec.size() > 0)
                    house.member(i).plong_ = &indiv_datalist.at(c);
    }

    void DataLong::helper_build_clear_data_()
    {
        // data sequence headers
        headers.clear();
        indiv_datalist.clear();
    }

    void DataLong::helper_build_data_headers_(const DataWide &data_wide, const str_table &input)
    {
        if (data_wide.is_empty())
            throw std::invalid_argument(
                "Tried to define data long format before wide format.");

        headers = read::extract_keys(input);
        // check for individual IDs
        if (std::find(headers.begin(), headers.end(), "id") == headers.end())
            throw std::invalid_argument(
                "Individual identifier column 'id' was missing.");
        // remove 'id'
        headers = headers - std::vector<std::string>({"id"});

        // pedantic check for declared data keys
        str::assert_no_missing_entries(headers, pKeys->get_all());

        // add optional keys to registry
        for (const std::string &datastring : headers)
            pKeys->add(datastring);
    }
}
