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

#include "tools/keys.h"
#include "tools/operators.h"
#include "tools/types.h"

#include "attribute.h"
#include "augmented_map.h"

namespace hsp::status
{
	constexpr int INFECTED = 1;	   /// Infected status (>=)
	constexpr int SUSCEPTIBLE = 0; /// Susceptible status (==)
	constexpr int MISSING = -1;	   /// Missing status (<=)

	/// Wheter individual status is infected.
	inline bool is_infected(int s) { return s >= INFECTED; }
	/// Wheter individual status is susceptible.
	inline bool is_susceptible(int s) { return s == SUSCEPTIBLE; }
	/// Wheter individual status is missing.
	inline bool is_missing(int s) { return s <= MISSING; }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Individual : class *** /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Class for individual household study member.
class Individual : public AugmentedMap
{
	friend class Household;
	friend struct GroupRate;

public:
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  Individual : constructor  /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	Individual()
	{
		id_ = "";
		index_ = 0;
		hhindex_ = 0;
		status_ = hsp::status::MISSING;
		data_map_ = {};
	}
	
	// NOTE: allowed data variables are numerical only
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
	double group_rate(GroupRateKey g) const { return group_rate_vec_[static_cast<size_t>(g)].get(); }
	const Attribute &get_group_rate_attribute(size_t g) const { return group_rate_vec_[g]; }

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  Individual : flags  /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	bool is_infected_status() const { return hsp::status::is_infected(status_); }
	bool is_healthy_status() const { return hsp::status::is_susceptible(status_); }
	bool is_missing_status() const { return hsp::status::is_missing(status_); }

	/// Initialize pointer to data long format.
	void init_pointer_to_data_long(hsp::dict<std::vector<double>> *y) { pdatalong_ = y; }

private:
	std::string id_;									  /// user id
	size_t index_;										  /// housepi internal index
	size_t hhindex_;									  /// household relative index
	int status_;										  /// status of the individual
	hsp::dict<double> data_map_ = {};					  /// data
	std::vector<Attribute> group_rate_vec_ = {};		  /// group rates
	hsp::dict<std::vector<double>> *pdatalong_ = nullptr; /// data long format associated
};

#endif // HOUSEPI__INDIVIDUAL__H
