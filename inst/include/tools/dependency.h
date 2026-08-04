/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <set>
#include <string>
#include <vector>

#include "messages.h"
#include "types.h"

/// NOTE:
/// Producer of 'A' : what 'A' reads (upstream, its inputs)
/// Consumer of 'A' : what reads 'A' (downstream, its dependents)
/// Dependency graph : directed acyclic graph, consumer node -> its producers
namespace housepi::dpdcy
{
    /// Producer graph: every node maps to the set of nodes it reads.
    /// Producers are unique and never read by position, only traversed.
    using prod_graph = dict<std::set<std::string, std::less<>>>;

    /// @brief Walk the edges of 'graph' from 'entry'.
    /// Asserts every node is known and that no cycle is closed.
    /// @param entry the starting node for the walk.
    /// @param graph the dependency graph to traverse.
    /// @param visited nodes already visited during the walk.
    /// @param recursion_stack nodes in the current recursion path.
    /// @param ordered the post-order traversal result.
    /// @param err error-type to raise errors.
    void walk_graph(const std::string &entry,
                    const prod_graph &graph,
                    std::vector<std::string> &visited,
                    std::vector<std::string> &recursion_stack,
                    std::vector<std::string> &ordered,
                    ErrType err);

    /// @brief Topological order of every node of the graph, producers first.
    /// @param graph the dependency graph to traverse.
    /// @param err error-type to raise errors.
    /// @return every node in the graph, ordered by producer before consumer.
    std::vector<std::string> order_producers_first(const prod_graph &graph,
                                                   ErrType err);

    /// @brief Return 'start' and its transitive producers.
    /// Producers are read by 'start' directly or through intermediaries.
    /// Ordered producer before consumer, so 'start' comes last.
    /// @param start the starting node for the walk.
    /// @param graph the dependency graph to traverse.
    /// @param err error-type to raise errors.
    std::vector<std::string> transitive_producers(const std::string &start,
                                                  const prod_graph &graph,
                                                  ErrType err);

    /// @brief Return 'start' and its transitive consumers.
    /// Consumers read 'start' directly or through intermediaries.
    /// Ordered producer before consumer, so 'start' comes first.
    /// @param start the starting node for the walk.
    /// @param graph the dependency graph to traverse.
    /// @param err error-type to raise errors.
    std::vector<std::string> transitive_consumers(const std::string &start,
                                                  const prod_graph &graph,
                                                  ErrType err);
}
