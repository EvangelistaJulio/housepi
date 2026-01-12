/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__INDIVIDUAL__H
#define HOUSEPI__INDIVIDUAL__H

#include "tools/types.h"
#include "tools/operators.h"
#include "tools/global.h"
#include "attribute.h"
#include "augmented_map.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Individual : class *** /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Forward declaration.
class Household;

/// Class for individual household study member.
class Individual : public AugmentedMap
{
public:
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  Individual : constructor  /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	Individual() {}
	Individual(size_t index,
			   size_t hhindex,
			   const std::string &id,
			   int status,
			   const u_dict<double> &data)
	{
		index_ = index;
		hhindex_ = hhindex;
		id_ = id;
		status_ = status;
		data_map_ = data;
	}

	/// Sugar to tell if the corresponding members are the same.
	bool operator==(const Individual &other) const { return this->index_ == other.index_; }
	/// Sugar to tell if the corresponding members are different
	bool operator!=(const Individual &other) const { return this->index_ != other.index_; }
	/// Sugar to call an Individual object as house index.
	operator size_t() const { return hhindex_; }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : getter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	const std::string &id() const { return id_; }
	std::string err_name() const override final { return "individual '" + id_ + "'"; }
	int status() const { return status_; }
	bool is_infected_status() const { return hsp::status::is_infected(status_); }
	bool is_healthy_status() const { return hsp::status::is_susceptible(status_); }
	bool is_missing_status() const { return hsp::status::is_missing(status_); }
	double data(const std::string &d) const { return safe_at(data_map_, d, "data"); }
	const dbl_vec &data_long(const std::string &d) const
	{
		// return a reference to a stable empty vector when no datalong is attached
		static const dbl_vec EMPTY_DBL_VEC = dbl_vec();
		return (pdatalong_ == nullptr) ? EMPTY_DBL_VEC : pdatalong_->at(d);
	}
	double group_rate(const std::string &g) const { return safe_at(group_rate_map_, g, "group rate").get(); }
	const Attribute &get_group_rate_attribute(const std::string &g) const { return safe_at(group_rate_map_, g, "group rate"); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : setter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Initialize pointer to data long format.
	void init_pointer_to_data_long(u_dict<dbl_vec> *y) { pdatalong_ = y; }
	/// Initialize group rate according to given attribute.
	void init_group_rate(const std::string &g, const Attribute &attr)
	{
		group_rate_map_[g].transfer(attr);
		group_rate_map_[g].set_name(g, "group rate");
	}
	/// Initialize group rate to reference.
	void init_group_rate_to_reference(const std::string &g)
	{
		group_rate_map_[g].transfer(Attribute(hsp::rates::VALUE_REFERENCE));
		group_rate_map_[g].set_name(g, "group rate");
	}

private:
	std::string id_ = hsp::na::STRING;	// user id
	size_t index_ = 0;					// housepi internal index
	size_t hhindex_ = 0;				// household relative index
	int status_ = hsp::status::MISSING; // status of the individual

	// NOTE: allowd data variables are numerical only
	u_dict<double> data_map_ = {};			// data
	u_dict<Attribute> group_rate_map_ = {}; // group rates

	u_dict<dbl_vec> *pdatalong_ = nullptr; // data long format associated

	friend class Household;
};

#endif // HOUSEPI__INDIVIDUAL__H
