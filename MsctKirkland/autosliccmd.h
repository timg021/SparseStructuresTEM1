#include <string>
#include <vector>

#define TEG_MULTITHREADED 1 //if defined, multithreaded execution is used

class Counter_Obj
{
	static unsigned int counter;
	static bool isUpdated;
public:
	Counter_Obj() { counter++; }
	~Counter_Obj() { counter--; }
	unsigned int GetCount() { return counter;  }
	void SetUpdated(bool status) { isUpdated = status; }
	bool GetUpdated() { return isUpdated; }
};

int autosliccmd(std::vector<std::string> params);

