#include "juke_box.h"

#define MOD

#ifdef MOD
void JukeBox::pezzali()
{
	print_stdout("(1) Io, di risposte non ne ho");	
}
#else
void JukeBox::pezzali()
{
  print_stdout("(2) C'è un tempo per i baci sperati, desiderati");	
}
#endif