
#ifndef POWELL_ALGORITHM_FILE_HEADER
#define POWELL_ALGORITHM_FILE_HEADER

class engine//优化目标继承此类，定义其中的虚函数
{
public:
	engine();
	virtual ~engine();
	
	void register_engine();
	void release_engine();
	
	virtual double compute_energy(double* current_value) = 0;
};

class optimization
{
public:
	optimization();
	virtual ~optimization();

	virtual double optimize(double* value, int value_count, double ftol = 0.00000001) = 0;
};

class optimization_powell
{
public:
	optimization_powell();
	virtual ~optimization_powell();

	virtual double optimize(double* value, int value_count, double ftol = 0.00000001);//value: input and output
};



#endif