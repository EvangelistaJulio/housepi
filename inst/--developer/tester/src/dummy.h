
#ifndef DUMMY__H
#define DUMMY__H

#include <string.h>
#include <Rcpp.h>

using std::string;

class JukeBox
{

private:
	void print_stdout(const string &txt)
	{
		Rcpp::Rcout << txt;
		Rcpp::Rcout.flush();
	}

public:
	void collins()
	{
		string out = "Come stop your crying\nIt will be alright\n";
		out += "Just take my hand\n Hold it tight\nI will protect you from all around you\n";
		out += "I will be here\nDon't you cry\n";
		print_stdout(out);
	}

	void pezzali()
	{
		string out = "Io, di risposte non ne ho\nMai avute, mai ne avrò\n";
		out += "Di domande ne ho quante ne vuoi\nE tu, neanche tu mi fermerai\n";
		out += "Neanche tu ci riuscirai\nIo non sono quel tipo di uomo e non lo sarò mai\n";
		print_stdout(out);	
	}
}

#endif // DUMMY__H