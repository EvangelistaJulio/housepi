#include <Rcpp.h>
#include "juke_box.h"

void JukeBox::print_stdout(const std::string &txt)
{
	Rcpp::Rcout << txt;
	Rcpp::Rcout.flush();
}
