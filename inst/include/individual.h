/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "function.h"
#include "augm_container.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// *** Individual : class *** /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
	/// Class for individual household study member.
	class Individual final : public AugmContainer
	{
		friend class Household;
		friend struct DataLong;
		friend struct GroupRate;

	public:
		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		Individual() = default;
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

		double data(size_t d) const { return data_vec_[d]; }

		/// Return associated data long format.
		/// NOTE: return a reference to a static empty vector when no datalong is attached.
		const std::vector<double> &data_long(size_t d) const
		{
			if (pdatalong_)
				return (*pdatalong_)[d];

			static const std::vector<double> EMPTY_DATA_LONG_VEC = std::vector<double>();
			return EMPTY_DATA_LONG_VEC;
		}

		/// Return associated data long format.
		/// NOTE: return a reference to a static empty vector when no datalong is attached
		const std::vector<std::vector<double>> *get_data_long_ptr() const noexcept { return pdatalong_; }

		double group_rate(size_t g) const { return group_rate_vec_[g].value(); }
		const Attribute &get_group_rate_attribute(size_t g) const { return group_rate_vec_[g]; }

		////////////////////
		//////  SHOW  //////
		////////////////////

		size_t hhindex() const noexcept { return hhindex_; }
		std::string err_name() const { return "individual '" + id_ + "'"; }

	private:
		size_t index_ = 0;										/// housepi internal index
		size_t hhindex_ = 0;									/// household relative index
		std::vector<double> data_vec_ = {};						/// data
		std::vector<Attribute> group_rate_vec_ = {};			/// group rates
		std::vector<std::vector<double>> *pdatalong_ = nullptr; /// data long format associated
	};
}
