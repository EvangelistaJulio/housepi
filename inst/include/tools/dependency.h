/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "types.h"

namespace housepi::dependency
{
    /// Helper to build definition order and assert for acyclic dependencies.
    void walk_tree(const std::string &entry,
                   const str_table &tree,
                   std::vector<std::string> &visited,
                   std::vector<std::string> &recursion_stack,
                   std::vector<std::string> &ordered);

    /// Topological order of all definitions dependencies.
    std::vector<std::string> order_dependencies(const str_table &graph);

    /// Topological order of definitions conditioned on 'start' (start included).
    std::vector<std::string> order_conditioned_on(const std::string &start, const str_table &graph);
}
