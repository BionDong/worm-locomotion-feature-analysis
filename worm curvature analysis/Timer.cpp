#include "stdafx.h"
#include "Timer.h"
#include <iostream>
static std::ofstream fout;

myTimer::myTimer(void)
{
	QueryPerformanceFrequency(&freq);//获取主机CPU时钟频率
	startCount.QuadPart = 0;
	endCount.QuadPart = 0;
	totalTime = 0;
	fout.open("TimeRecord.txt");
}

myTimer::~myTimer(void)
{
	fout.close();
	fout.clear();
}

void myTimer::StartTimer(const char* name)
{
	QueryPerformanceCounter(&startCount);//开始计时
	fout << name << "(milli second): ";
}

void myTimer::StopTimer()
{
	QueryPerformanceCounter(&endCount);//停止计时

	dbTime = (endCount.QuadPart - startCount.QuadPart) / (double)freq.QuadPart * 1000;//获取时间差 毫秒

	totalTime += dbTime;
	fout << dbTime << "		Total Time: " << totalTime << std::endl;
}

void myTimer::GetTime()
{
	QueryPerformanceCounter(&currentCount);

	cTime = (double)currentCount.QuadPart/(double)freq.QuadPart * 1000; // 此刻时间 毫秒

	CurrentTime = currentCount.QuadPart * (1000.0 / freq.QuadPart);
}