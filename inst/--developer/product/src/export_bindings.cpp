#include <Rcpp.h>
#include "juke_box.h"

RCPP_MODULE(jukebox_rcpp)
{
  using namespace Rcpp;

  class_<JukeBox>("JukeBox")
      .constructor()
      .method("pezzali", &JukeBox::pezzali)
      .method("collins", &JukeBox::collins);
}
