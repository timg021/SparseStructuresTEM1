//@@@@@@ start TEG code
#include <string>
#include <vector>
#include "XA_file.h"

#define TEG_MULTITHREADED 1 //if defined, multithreaded execution is used

// TEG class for keeping a count across multiple threads (e.g. for counting the number of active threads)
class Counter_Obj
{
	static int counter;
	static bool isUpdated;
public:
	Counter_Obj() { counter++; }
	~Counter_Obj() { counter--; }
	int GetCount() { return counter;  }
	void SetUpdated(bool status) { isUpdated = status; }
	bool GetUpdated() { return isUpdated; }
	void SetTerminate() { counter = -100000; }
};

int autosliccmd(std::vector<std::string> params, std::vector<xar::Pair> defocus, std::vector<std::string> fileout);
int autosliccmd1(std::vector<std::string> params);

//@@@@@ end TEG code

