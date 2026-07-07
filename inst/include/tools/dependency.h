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

namespace housepi::dpdcy
{
    /// Helper to build definition order and assert for acyclic dependencies.
    void walk_graph(const std::string &entry,
                    const str_table &graph,
                    std::vector<std::string> &visited,
                    std::vector<std::string> &recursion_stack,
                    std::vector<std::string> &ordered);

    /// Topological order of all definitions dependencies.
    std::vector<std::string> order_dependencies(const str_table &graph);

    /// Return 'start' and its transitive consumers walking the inverse 'graph'.
    /// Consumers depend on 'start' directly or through intermediaries.
    /// Ordered respects producer before consumer ('start' is always the first element).
    std::vector<std::string> transitive_consumers(const std::string &start, const str_table &graph);

    /// Return 'start' and its transitive producers walking the forward 'graph'.
    /// (the names it reads, directly or through intermediaries)
    /// Ordered respects producer before consumer.
    std::vector<std::string> transitive_producers(const std::string &start, const str_table &graph);
}
