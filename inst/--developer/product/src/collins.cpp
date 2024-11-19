#include "juke_box.h"

//#define MOD

#ifdef MOD
void JukeBox::collins()
{
	print_stdout("(A) Come stop your crying");
}
#else
void JukeBox::collins()
{
  print_stdout("(B) I can feel it coming in the air tonight");
}
#endif