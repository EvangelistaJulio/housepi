/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <memory>

namespace hsp::progress
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
