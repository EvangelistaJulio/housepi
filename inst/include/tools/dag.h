/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <set>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"

/// Traverse directed acyclic graphs of definitions.
///
/// * producer of `A`: what `A` reads, upstream;
/// * consumer of `A`: what reads `A`, downstream;
/// * graph: each consumer node maps to its producers.
namespace housepi::dag
{
    /// Map every node to the set of nodes it reads.
    /// Producers are unique, traversed but never read by position.
    using prod_graph = dict<std::set<std::string, std::less<>>>;

    /// Walk the edges of `graph` from `entry`, appending nodes to `ordered` in post-order.
    ///
    /// * unknown node -> error;
    /// * closed cycle -> error.
    /// @param visited Nodes visited during the walk.
    /// @param recursion_stack Nodes on the current recursion path.
    /// @param err Entity type named in errors.
    void walk_graph(const std::string &entry,
                    const prod_graph &graph,
                    std::vector<std::string> &visited,
                    std::vector<std::string> &recursion_stack,
                    std::vector<std::string> &ordered,
                    ErrType err);

    /// Return every node of `graph` in topological order, producers first.
    /// @param err Entity type named in errors.
    std::vector<std::string> order_producers_first(const prod_graph &graph,
                                                   ErrType err);

    /// Return `start` and its transitive producers, producers first.
    ///
    /// * producers: read by `start` directly or through intermediaries;
    /// * `start` comes last.
    /// @param err Entity type named in errors.
    std::vector<std::string> transitive_producers(const std::string &start,
                                                  const prod_graph &graph,
                                                  ErrType err);

    /// Return `start` and its transitive consumers, producers first.
    ///
    /// * consumers: read `start` directly or through intermediaries;
    /// * `start` comes first.
    /// @param err Entity type named in errors.
    std::vector<std::string> transitive_consumers(const std::string &start,
                                                  const prod_graph &graph,
                                                  ErrType err);
}
