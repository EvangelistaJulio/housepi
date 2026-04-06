/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include "data_wide.h"
#include "data_long.h"

namespace hsp
{
    struct DataStore
    {
        DataWide data_wide;
        DataLong data_long;

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        DataStore() = default;
        ~DataStore() = default;

        /// NOTE: I never copy DataStore.
        DataStore(const DataStore &) = delete;
        DataStore &operator=(const DataStore &) = delete;
    };
}
