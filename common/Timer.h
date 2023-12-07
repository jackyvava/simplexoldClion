//////////////////////////////////////////////////////////////////////////
// Timer
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Timer_h__
#define __Timer_h__
#include <chrono>
#include <ctime>
#include <string>
#include <iostream>
#include <iomanip>

template<class T> class Timer
{
public:
	void Reset()
	{
		start=std::chrono::system_clock::now();
	}

	T Elapse()
	{
		std::chrono::time_point<std::chrono::system_clock> end=std::chrono::system_clock::now();
		std::chrono::duration<T,std::ratio<1,1000> > elapse=end-start;
		return elapse.count();
	}

	T Elapse_And_Reset()
	{
		T elapse=Elapse();
		Reset();
		return elapse;
	}

	void Elapse_And_Output(const std::string& message)
	{
		T elapse=Elapse();
		std::cout<<"Timer for "<<std::setw(16)<<message<<": \t"<<elapse<<" ms"<<std::endl;
	}

	void Elapse_And_Output_And_Reset(const std::string& message)
	{
		Elapse_And_Output(message);
		Reset();
	}

protected:
	std::chrono::time_point<std::chrono::system_clock> start;
};
#endif