/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 * 
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 * See the COPYRIGHT file for details.
 * 
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the COPYING file for details.
 *
 * tools/messages.h
 *
 */

#ifndef HOUSEPI__TOOLS__MESSAGES__H
#define HOUSEPI__TOOLS__MESSAGES__H

#include <string>

namespace hsp
{
    namespace msg
    {
        /// Throw an error message.
        [[noreturn]] void error(const std::string &txt);
        /// Show a message to stdout and flush.
        void info(const std::string &txt);
        /// Display a warning message if not already displayed before.
        void warning(const std::string &txt);
    }
}

#endif // HOUSEPI__TOOLS__MESSAGES__H
