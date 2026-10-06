/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/dag.h"
#include "tools/dictionary.h"
#include "augmented/definition.h"

//////////////////////////////  augm::graph  //////////////////////////////

/// Dependency analysis over definition rows, by name.
namespace housepi::augm::graph
{
    /// Build the latent producer graph of `rows`, for a simulation or inference run.
    /// - node: definition name;
    /// - edges: latent names it reads, across the source columns of the run, see DefRow::source_cols().
    dag::prod_graph producers(const DefTable &def,
                              const std::vector<size_t> &rows,
                              bool simulation);

    /// Add an edge from each row of `rows` to the definitions `targets` whose data column it reads.
    ///
    /// - node: every row of `rows`, with or without edges;
    /// - row reading its own column: no edge.
    /// @note Edges follow the source columns of a simulation run, see DefRow::source_cols().
    void add_data_edges(dag::prod_graph &pgraph,
                        const DefTable &def,
                        const std::vector<size_t> &rows,
                        const std::vector<std::string> &targets);

    /// Map each model name to its transitive consumers among other model nodes of `order`.
    str_table model_consumers(const DefTable &def,
                              const std::vector<size_t> &order);

    /// List latent names of `order` that no sink reaches.
    ///
    /// Sinks:
    /// - transmission names `transm`;
    /// - observations of `order`, except `absent` ones;
    /// - latents a household observation reduces over.
    ///
    /// Such a latent is an unused parameter.
    /// Marginalizing it integrates its model to one: its draw and density are dropped.
    std::vector<std::string> dead_names(const DefTable &def,
                                        const std::vector<size_t> &order,
                                        const std::vector<std::string> &transm,
                                        const std::vector<std::string> &absent);

    /// List names of `order` that transmission names `transm` consume, `transm` included.
    ///
    /// Consumption is direct or through intermediaries.
    std::vector<std::string> transm_producers(const DefTable &def,
                                              const std::vector<size_t> &order,
                                              const std::vector<std::string> &transm);
}
