/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__INDIVIDUAL__H
#define HOUSEPI__INDIVIDUAL__H

#include <string>
#include <vector>

#include "tools/operators.h"
#include "tools/types.h"

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

	Individual();
	Individual(size_t index, size_t hhindex, const std::string &id,
			   int status, const hsp::dict<double> &data)
	{
		index_ = index;
		hhindex_ = hhindex;
		id_ = id;
		status_ = status;
		data_map_ = data;
	}
	~Individual() = default;

	/// Sugar to tell if the corresponding members are the same.
	bool operator==(const Individual &other) const { return this->index_ == other.index_; }
	/// Sugar to tell if the corresponding members are different
	bool operator!=(const Individual &other) const { return this->index_ != other.index_; }
	/// Sugar to call an Individual object as house index.
	operator size_t() const { return hhindex_; }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : getter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Individidaul unique identifier as in data.
	const std::string &id() const { return id_; }
	std::string err_name() const override final { return "individual '" + id_ + "'"; }
	int status() const { return status_; }
	double data(const std::string &d) const { return hsp::safe_at(data_map_, d, "data"); }

	/// Return associated data long format.
	/// NOTE: return a reference to a static empty vector when no datalong is attached
	const std::vector<double> &data_long(const std::string &d) const
	{
		if (pdatalong_ != nullptr)
			return pdatalong_->at(d);

		static const std::vector<double> EMPTY_DBL_VEC = std::vector<double>();
		return EMPTY_DBL_VEC;
	}
	double group_rate(const std::string &g) const { return hsp::safe_at(group_rate_map_, g, "group rate").get(); }
	const Attribute &get_group_rate_attribute(const std::string &g) const { return hsp::safe_at(group_rate_map_, g, "group rate"); }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : flags  /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	bool is_infected_status() const;
	bool is_healthy_status() const;
	bool is_missing_status() const;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : setter  ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	/// Initialize pointer to data long format.
	void init_pointer_to_data_long(hsp::dict<std::vector<double>> *y) { pdatalong_ = y; }
	/// Initialize group rate according to given attribute.
	void init_group_rate(const std::string &g, const Attribute &attr);
	/// Initialize group rate to reference.
	void init_group_rate_to_reference(const std::string &g);

private:
	std::string id_; // user id
	size_t index_;	 // housepi internal index
	size_t hhindex_; // household relative index
	int status_;	 // status of the individual

	// NOTE: allowd data variables are numerical only
	hsp::dict<double> data_map_ = {};			 // data
	hsp::dict<Attribute> group_rate_map_ = {}; // group rates

	hsp::dict<std::vector<double>> *pdatalong_ = nullptr; // data long format associated

	friend class Household;
};

#endif // HOUSEPI__INDIVIDUAL__H
