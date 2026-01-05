/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__FUNCTION__REGISTRY__H
#define HOUSEPI__FUNCTION__REGISTRY__H

#include <memory>
#include <functional>

#include "tools/types.h"
#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** FunctionRegistry : class *** //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class FunctionRegistry
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  FunctionRegistry : constructor  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    FunctionRegistry(const std::string &name)
        : name_(name) {}
    virtual ~FunctionRegistry() = default;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  FunctionRegistry : getter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    const std::string &name() const { return name_; }
    virtual str_vec get_options() const;
    virtual std::unique_ptr<Function> get(const std::string &opt, bool is_primitive = false) const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  FunctionRegistry : setter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    void register_func()
    {
        static_assert(std::is_base_of<Function, T>::value,
                      "T must be derived from Function");
        factory_map_[T().name()] = []()
        { return std::make_unique<T>(); };
    }

private:
    std::string name_;
    u_dict<std::function<std::unique_ptr<Function>()>> factory_map_ = {};
};

/// Separate functions registries.
#define DECLARE_HOUSEPI_FUNCTION_REGISTRY(ClassName, Description) \
    class ClassName : public FunctionRegistry                     \
    {                                                             \
    public:                                                       \
        static ClassName &inst();                                 \
                                                                  \
    private:                                                      \
        ClassName() : FunctionRegistry(Description) {}            \
    };
DECLARE_HOUSEPI_FUNCTION_REGISTRY(ArithmeticRegistry, "arithmetic function")
DECLARE_HOUSEPI_FUNCTION_REGISTRY(TranscendentalRegistry, "transcendental function")
DECLARE_HOUSEPI_FUNCTION_REGISTRY(UserFuncRegistry, "user function")

/// Register a Function derived class T to FungRegistry derived class R.
template <typename R, typename T>
struct StaticFunctionRegistatror
{
    StaticFunctionRegistatror()
    {
        static_assert(std::is_base_of<FunctionRegistry, R>::value,
                      "R must be derived from FunctionRegistry.");
        R::inst().template register_func<T>();
    }
};
#define REGISTER_HOUSEPI_FUNCTION(R, T)        \
    namespace                                  \
    {                                          \
        static StaticFunctionRegistatror<R, T> \
            housepi_register_##T##_to_##R;     \
    }

#define REGISTER_USER_FUNCTION(T) REGISTER_HOUSEPI_FUNCTION(UserFuncRegistry, T)

/// Compose registries to get a function from them.
std::unique_ptr<Function> compose_func_registries(const std::string &opt,
                                                  const std::string &name,
                                                  const std::vector<const FunctionRegistry *> &registries,
                                                  bool is_prim = false);

#endif // HOUSEPI__FUNCTION__REGISTRY__H
