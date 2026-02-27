/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <Rcpp.h>
#include "model.h"

inline void invoke_destructor_(ModelBase *p) { delete p; }
using ModelXPtr = Rcpp::XPtr<ModelBase, Rcpp::PreserveStorage,
							 invoke_destructor_>;

