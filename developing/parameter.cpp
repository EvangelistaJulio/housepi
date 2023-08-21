// parameter.cpp //

#include "parameter.hpp"

using namespace std;

///////////////////////////////////////////////////////////
// Parameter_Dictionary : Class
///////////////////////////////////////////////////////////

/* add_param */
void Parameter_Dictionary::add_param(const std::string k, Parameter new_param)
{
    params.emplace(k, new_param); // FIXME: check here for double parameters
    if (new_param.get_infer() == true)
    {
        to_be_inferred.push_back(k);
    }
}

/* print */
void Parameter_Dictionary::print()
{
    std::cout << "_NAME_" << TAB << "_VALUE_" << TAB << "_INFER_" << endl;
    for (const auto &entry : params)
    {
        Parameter param = entry.second;
        std::cout << entry.first << TAB << param.get_value();
        std::cout << TAB << (param.get_infer() > 0 ? "true" : "false") << endl;
    }
}

///////////////////////////////////////////////////////////
//////// Augmented : class definition /////////////////////
///////////////////////////////////////////////////////////

// FIXME : this  isn't efficient
kwargs& Augmented::address_args()
{
	for (auto const& entry : model_args_names)
	{
		set_arg(entry.first, glob_params.get_value(entry.second));
	}
	return args;
}