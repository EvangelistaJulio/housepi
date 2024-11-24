#include <Rcpp.h>
using namespace Rcpp;

// cli C API
https://cli.r-lib.org/reference/progress-c.html

// chatgpt solution
void print_message_with_theme(std::string message) {
    using Rcpp;
    // Check if rstudioapi is available
    Function isAvailable("rstudioapi::isAvailable");
    bool rstudio_available = as<bool>(isAvailable());

    std::string color_code;
    if (rstudio_available) {
        // Call rstudioapi::getThemeInfo()
        Function getThemeInfo("rstudioapi::getThemeInfo");
        List theme_info = getThemeInfo();

        // Check the "dark" field
        bool is_dark = false;
        if (theme_info.containsElementNamed("dark")) {
            is_dark = theme_info["dark"];
        }

        // Assign color based on theme
        color_code = is_dark ? "\033[38;5;15m" : "\033[38;5;0m";
    } else {
        // Fallback to a default color
        color_code = "\033[38;5;7m"; // Default gray
    }

    // Print the styled message
    std::string styled_message = color_code + message + "\033[0m"; // Reset after message
    Rcout << styled_message << std::endl;
}


#include <R_ext/Print.h>
#include <RcppArmadillo.h>

#include <progress_bar.hpp>

// Adapted from: Othman El Hammouchi (2024)
// https://othmanelhammouchi.com/posts/r-parallel-progress/

class CliProgressBar : public ProgressBar
{
public:
  CliProgressBar() { reset(); }
  ~CliProgressBar() {}

  void display()
  {
    Rcpp::Rcout << "\033[37mRunning simulations ";
  }

  void update(float progress)
  {
    _update_ticks_display(progress);
    if (_ticks_displayed >= _max_ticks)
    {
      _finalize_display();
    }
  }

  void end_display()
  {
    update(1);
    reset();
  }

  void reset()
  {
    Rcpp::Environment cli = Rcpp::Environment::namespace_env("cli");
    Rcpp::Function console_width = cli["console_width"];

    // From fiddling around with it, it seems that dividing the console
    // width by 2 produces the best display.
    _max_ticks = std::floor(Rcpp::as<int>(console_width()) / 2);
    _ticks_displayed = 0;
    _finalized = false;
  }

private:
  int _max_ticks;       // the total number of ticks to print
  int _ticks_displayed; // the nb of ticks already displayed
  bool _finalized;

  void _update_ticks_display(float progress)
  {
    int nb_ticks = _compute_nb_ticks(progress);
    if (nb_ticks > _ticks_displayed)
    {
      Rcpp::Rcout << "\r" << std::flush;
      _ticks_displayed = nb_ticks;
      _display_ticks(progress);
    }
  }

  void _finalize_display()
  {
    if (_finalized)
      return;

    Rcpp::Rcout << std::endl;
    _finalized = true;
  }

  int _compute_nb_ticks(float progress)
  {
    return int(progress * _max_ticks);
  }

  void _display_ticks(double progress)
  {
    Rcpp::Rcout << "\033[37mRunning simulations ";
    for (int i = 0; i < _ticks_displayed; ++i)
    {
      Rcpp::Rcout << "\033[32m\u25A0" << std::flush;
    }
    for (int i = 0; i < (_max_ticks - _ticks_displayed); i++)
    {
      Rcpp::Rcout << " " << std::flush;
    }
    Rcpp::Rcout << "\033[37m | " << int(progress * 100) << "%";
  }
};