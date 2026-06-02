/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/types.h"

#include "indiv_storage.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;
    struct Parameter;
    class Individual;

    struct GroupRate : public KeyModule
    {
    public:
        /// Update model report and set group rates associate to new individual.
        /// NOTE: this has to run even if input is empty to initialize existing rates to reference.
        void build(DataWide &data_wide, const Parameter &param, const str_table &deftable);
        /// Optimize all variables (see Variable::optimize).
        void optimize();

    private:
        std::vector<Indiv::Group> unique_groups_; /// Unique group rate definitions

        /// Verify that group rates definitions are correct and pre-process theme.
        std::vector<std::string> helper_process_definitions_(DataWide &data_wide, const str_table &input);
        /// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
        void helper_init_unique_rate_(const Parameter &param, const Individual &member, const str_table &these_rates);
    };
}
