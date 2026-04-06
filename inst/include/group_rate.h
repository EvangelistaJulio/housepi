/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"

#include "indiv_storage.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    struct KeyHolder;
    struct DataWide;
    struct Parameter;
    class Individual;

    struct GroupRate
    {
    public:
        std::unique_ptr<KeyHolder> pKeys; /// pointer to registry

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        GroupRate();
        ~GroupRate();

        void init_registry(const std::vector<std::string> &keys);

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Update model report and set group rates associate to new individual.
        /// NOTE: this has to run even if input is empty to initialize existing rates to reference.
        void build(DataWide &data_wide, const Parameter &param, const str_table &deftable);

    private:
        std::vector<Indiv::Group> indiv_group_rates_; /// External storage for individual group rates

        /// Verify that group rates definitions are correct.
        std::vector<std::string> helper_assert_definitions_(DataWide &data_wide, const str_table &input);
        /// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
        void helper_init_all_rates_(const Parameter &param, Individual &member, const str_table &these_rates);
    };
}
