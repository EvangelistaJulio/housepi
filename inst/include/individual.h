/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/status.h"

#include "indiv_storage.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Individual : class *** /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
	/// Class for individual household study member.
	/// Acts as a lightweight viewer over externally-stored sub-objects
	/// (data, augmented, group rates) owned by their respective domain classes.
	class Individual final
	{
		friend class Household;
		friend struct DataLong;
		friend struct DataWide;
		friend struct GroupRate;
		friend struct Augmented;
		friend class MCMCBase;
		friend class MCMC;

	public:
		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		Individual() = default;

		/////////////////////////
		//////  OPERATORS  //////
		/////////////////////////

		/// Equal to operator.
		bool operator==(const Individual &other) const noexcept { return index_ == other.index_; }
		/// Not equal to operator.
		bool operator!=(const Individual &other) const noexcept { return !(*this == other); }

		//////////////////////
		//////  SETTER  //////
		//////////////////////

		/// Unique identifier as provided in data.
		const std::string &id() const noexcept { return id_; }
		/// Associated status.
		int status() const noexcept { return status_; }

		/// Set numeric value of augmented variable. If not allowed, raise exception.
		void set_augm_value(size_t a, double x) { paugm_->vec[a].set(x); }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		double data(size_t d) const { return pwide_->vec[d]; }

		/// Return associated data long format.
		/// NOTE: return a reference to a static empty vector when no datalong is attached.
		const std::vector<double> &data_long(size_t d) const
		{
			if (plong_)
				return plong_->vec[d];

			static const std::vector<double> EMPTY_DATA_LONG_VEC = std::vector<double>();
			return EMPTY_DATA_LONG_VEC;
		}

		/// Return associated data long format.
		/// NOTE: return a reference to a static empty vector when no datalong is attached
		const std::vector<std::vector<double>> *get_data_long_ptr() const noexcept { return &plong_->vec; }

		double aug(size_t a) const { return paugm_->vec[a].value(); }
		Variable &get_augm_variable(size_t a) { return paugm_->vec[a]; }
		const Variable &get_augm_variable(size_t a) const { return paugm_->vec[a]; }

		double group_rate(size_t g) const { return pgroup_rate_->vec[g].value(); }
		const Attribute &get_group_rate_attribute(size_t g) const { return pgroup_rate_->vec[g]; }

		/////////////////////
		//////  FLAGS  //////
		/////////////////////

		bool is_infected_status() const noexcept { return status::is_infected(status_); }
		bool is_susceptible_status() const noexcept { return status::is_susceptible(status_); }
		bool is_missing_status() const noexcept { return status::is_missing(status_); }

		////////////////////
		//////  SHOW  //////
		////////////////////

		size_t hhindex() const noexcept { return hhindex_; }
		std::string err_name() const { return "individual '" + id() + "'"; }

	private:
		/// Identity (owned)
		size_t index_ = 0;	  /// housepi internal index
		size_t hhindex_ = 0;  /// household relative index
		std::string id_ = ""; /// user id

		/// status of the individual
		int status_ = status::MISSING;

		/// External sub-object pointers
		const Indiv::Wide *pwide_ = nullptr;		/// wide-format data (see DataWide)
		const Indiv::Long *plong_ = nullptr;		/// long-format data (see DataLong)
		const Indiv::Group *pgroup_rate_ = nullptr; /// group rates (see GroupRate)
		Indiv::Augm *paugm_ = nullptr;				/// augmented variables (see Augmented)
	};
}
