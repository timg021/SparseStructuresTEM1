#include <string>

class Counter_Obj
{
private: unsigned int* pthread_counter_in;
public:
	Counter_Obj(unsigned int* pthread_counter) { (*pthread_counter)++; pthread_counter_in = pthread_counter; }
	~Counter_Obj() { (*pthread_counter_in)--; }
};

const size_t numaslicpars(27); // number of input parameter lines of autosliccmd() function

extern int autosliccmd(std::string params[numaslicpars], unsigned int* pthread_counter);

