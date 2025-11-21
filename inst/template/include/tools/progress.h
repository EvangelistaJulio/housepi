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
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * tools/progress.h
 *
 */

#ifndef HOUSEPI__TOOLS__PROGRESS__H
#define HOUSEPI__TOOLS__PROGRESS__H

#include <string>
#include <memory>

namespace hsp
{
    namespace progress
    {
        /// Class for custom progress bar.
        class Bar
        {
        public:
            Bar();
            ~Bar();
            void set_text(const std::string &text);
            void set_thintime(unsigned long thintime);
            void open(unsigned long total);
            void update();
            void crash() const;
            void close();

        private:
            class Impl; // Pointer to Implementation
            std::unique_ptr<Impl> pImpl_;
        };
    }
}

#endif // HOUSEPI__TOOLS__PROGRESS__H
