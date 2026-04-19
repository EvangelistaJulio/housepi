/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "tools/nan.h"
#include "tools/string.h"
#include "tools/types.h"

#include "registry.h"
#include "function.h"
#include "distribution.h"
#include "model.h"

HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Custom)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Arithmetic)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(Transcendental)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(ContinuosDistrib)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(DiscreteDistrib)
HOUSEPI_DEFINE_REGISTRY_INSTANCE_(InfecProfile)

namespace housepi
{
    ModelRegistry &ModelRegistry::instance()
    {
        static ModelRegistry _inst;
        return _inst;
    }
}

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// *** FunctionRegistry : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<std::string> FunctionRegistry::get_options() const
    {
        std::vector<std::string> entries;
        entries.reserve(factory_map_.size());
        for (const auto &[key, val] : factory_map_)
            entries.push_back(key);
        return entries;
    }

    std::unique_ptr<Function> FunctionRegistry::create(std::string_view opt, bool is_prim) const
    {
        auto it = factory_map_.find(opt);
        if (it == factory_map_.end())
            return nullptr;

        auto output = it->second();
        if (is_prim)
            output->set_eval_as_prim();
        // evaluate as function by default

        return output;
    }

    std::unique_ptr<Function> compose_func_registries(std::string_view opt,
                                                      const std::string &err_type,
                                                      const std::vector<std::string> &registries,
                                                      bool is_prim)
    {
        if (registries.empty())
            throw std::invalid_argument(
                "No registries were provided for '" + err_type + "'.");

        // FIXME: this should be an error
        if (na::is(opt))
            return nullptr;

        std::vector<const FunctionRegistry *> registry_ptrs = {};
        for (const std::string &r : registries)
            if (r == "custom")
                registry_ptrs.push_back(&CustomRegistry::instance());

            else if (r == "transcendental")
                registry_ptrs.push_back(&TranscendentalRegistry::instance());

            else if (r == "arithmetic")
                registry_ptrs.push_back(&ArithmeticRegistry::instance());

            else if (r == "continuos")
                registry_ptrs.push_back(&ContinuosDistribRegistry::instance());

            else if (r == "discrete")
                registry_ptrs.push_back(&DiscreteDistribRegistry::instance());

            else if (r == "infecprofile")
                registry_ptrs.push_back(&InfecProfileRegistry::instance());

            else
                throw std::invalid_argument(
                    "Unknown " + err_type + " registry '" + r + "'.");

        std::vector<std::string> options = {};
        for (const FunctionRegistry *ptr : registry_ptrs)
        {
            if (ptr->count(opt))
                return ptr->create(opt, is_prim);
            const std::vector<std::string> &y = ptr->get_options();
            options.reserve(options.size() + y.size());
            options.insert(options.end(), y.begin(), y.end());
        }

        throw std::invalid_argument(
            "Unknown " + err_type + " '" + std::string(opt) + "'." +
            str::find_closest(opt, options));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// *** DistributionRegistry : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<Function> DistributionRegistry::create(std::string_view opt, bool is_prim) const
    {
        auto it = factory_map_.find(opt);
        if (it == factory_map_.end())
            return nullptr;

        auto output = it->second();
        if (is_prim)
            output->set_eval_as_cdf();
        else
            output->set_eval_as_pdf();

        return output;
    }

    std::unique_ptr<Distribution> compose_distr_registries(std::string_view opt,
                                                           const std::string &err_type,
                                                           const std::vector<std::string> &registries,
                                                           bool is_cumul)
    {
        if (registries.empty())
            throw std::invalid_argument(
                "No registries were provided for '" + err_type + "'.");

        // FIXME: this should be an error
        if (na::is(opt))
            return nullptr;

        std::vector<const DistributionRegistry *> registry_ptrs = {};
        for (const std::string &r : registries)
            if (r == "continuos")
                registry_ptrs.push_back(&ContinuosDistribRegistry::instance());
            else if (r == "discrete")
                registry_ptrs.push_back(&DiscreteDistribRegistry::instance());
            else if (r == "infecprofile")
                registry_ptrs.push_back(&InfecProfileRegistry::instance());
            else
                throw std::invalid_argument(
                    "Unknown " + err_type + " registry '" + r + "'.");

        std::vector<std::string> options = {};
        for (const DistributionRegistry *r : registry_ptrs)
        {
            if (r->count(opt))
                return std::unique_ptr<Distribution>(
                    static_cast<Distribution *>(r->create(opt, is_cumul).release()));
            const std::vector<std::string> &y = r->get_options();
            options.reserve(options.size() + y.size());
            options.insert(options.end(), y.begin(), y.end());
        }

        throw std::invalid_argument(
            "Unknown " + err_type + " '" + std::string(opt) + "'." +
            str::find_closest(opt, options));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// *** ModelRegistry : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<ModelBase> ModelRegistry::create(std::string_view name) const
    {
        auto it = factory_map_.find(name);
        if (it == factory_map_.end())
            throw std::invalid_argument(
                "Model registration failed: unknown model '" + std::string(name) + "'.\n");
        return it->second();
    }
}
