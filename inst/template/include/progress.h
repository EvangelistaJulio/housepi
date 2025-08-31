/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * progress.h
 *
 */

#ifndef HOUSEPI__PROGRESS__H
#define HOUSEPI__PROGRESS__H

#include <string>
#include <memory>

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Progress Bar Namespace  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Class for custom progress bar. */
namespace progress
{
    class Bar
    {
    public:
        //===== CONSTRUCTOR =====//
        Bar();
        ~Bar();

        //===== SETTER =====//
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

#endif // HOUSEPI__PROGRESS__H
