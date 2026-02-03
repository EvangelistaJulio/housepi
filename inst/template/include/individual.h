/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

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
class Individual final : public AugmentedContainer
{
	friend class Household;
	friend struct DataLong;
	friend struct GroupRate;

public:
	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	Individual() = default;

	/// NOTE: allowed data variables are numerical only.
	Individual(size_t index, size_t hhindex, const std::string &id,
			   int status, const std::vector<double> &data)
	{
		index_ = index;
		hhindex_ = hhindex;
		id_ = id;
		status_ = status;
		data_vec_ = data;
	}
	~Individual() = default;

	/////////////////////////
	//////  OPERATORS  //////
	/////////////////////////

	/// Sugar to tell if the corresponding members are the same.
	bool operator==(const Individual &other) const noexcept { return index_ == other.index_; }
	/// Sugar to tell if the corresponding members are different
	bool operator!=(const Individual &other) const noexcept { return index_ != other.index_; }

	//////////////////////
	//////  GETTER  //////
	//////////////////////

	int status() const noexcept { return status_; }

	double data(size_t d) const { return data_vec_[d]; }
	const std::vector<double> &get_data_wide_vec() const noexcept { return data_vec_; }

	/// Return associated data long format.
	/// NOTE: return a reference to a static empty vector when no datalong is attached.
	const std::vector<double> &data_long(size_t d) const
	{
		if (pdatalong_ != nullptr)
			return (*pdatalong_)[d];

		static const std::vector<double> EMPTY_DATA_LONG_VEC = std::vector<double>();
		return EMPTY_DATA_LONG_VEC;
	}

	/// Return associated data long format.
	/// NOTE: return a reference to a static empty vector when no datalong is attached
	const std::vector<std::vector<double>> *get_data_long_ptr() const noexcept { return pdatalong_; }

	double group_rate(size_t g) const { return group_rate_vec_[g].get(); }
	const Attribute &get_group_rate_attribute(size_t g) const { return group_rate_vec_[g]; }

	/////////////////////
	//////  FLAGS  //////
	/////////////////////

	bool is_infected_status() const noexcept { return hsp::status::is_infected(status_); }
	bool is_healthy_status() const noexcept { return hsp::status::is_susceptible(status_); }
	bool is_missing_status() const noexcept { return hsp::status::is_missing(status_); }

	////////////////////
	//////  SHOW  //////
	////////////////////

	/// Individidaul unique identifier as in data.
	const std::string &id() const noexcept { return id_; }
	size_t hhindex() const noexcept { return hhindex_; }
	std::string err_name() const { return "individual '" + id_ + "'"; }

private:
	std::string id_ = "";									/// user id
	size_t index_ = 0;										/// housepi internal index
	size_t hhindex_ = 0;									/// household relative index
	int status_ = hsp::status::MISSING;						/// status of the individual
	std::vector<double> data_vec_ = {};						/// data
	std::vector<Attribute> group_rate_vec_ = {};			/// group rates
	std::vector<std::vector<double>> *pdatalong_ = nullptr; /// data long format associated
};
