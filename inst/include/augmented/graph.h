/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/dependency.h"
#include "tools/types.h"

#include "augmented/definition.h"

//////////////////////////////  augm::graph  //////////////////////////////

/// Dependency analysis over definition rows, by name.
namespace housepi::augm::graph
{
    /// Build the latent producer graph of `rows`.
    /// - node: definition name;
    /// - edges: @A names it reads, across source columns.
    dpdcy::prod_graph producers(const DefTable &def,
                                const std::vector<size_t> &rows);

    /// Add edges from a definition to the member observations of `rows` it reads as data.
    ///
    /// An observation is drawn before its reader.
    void add_member_obs_edges(dpdcy::prod_graph &pgraph,
                              const DefTable &def,
                              const std::vector<size_t> &rows);

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
