#include "stdafx.h"
#include "Timer.h"
#include <iostream>
static std::ofstream fout;

myTimer::myTimer(void)
{
	QueryPerformanceFrequency(&freq);//��ȡ����CPUʱ��Ƶ��
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
	QueryPerformanceCounter(&startCount);//��ʼ��ʱ
	fout << name << "(milli second): ";
}

void myTimer::StopTimer()
{
	QueryPerformanceCounter(&endCount);//ֹͣ��ʱ

	dbTime = (endCount.QuadPart - startCount.QuadPart) / (double)freq.QuadPart * 1000;//��ȡʱ��� ����

	totalTime += dbTime;
	fout << dbTime << "		Total Time: " << totalTime << std::endl;
}

void myTimer::GetTime()
{
	QueryPerformanceCounter(&currentCount);

	cTime = (double)currentCount.QuadPart/(double)freq.QuadPart * 1000; // �˿�ʱ�� ����

	CurrentTime = currentCount.QuadPart * (1000.0 / freq.QuadPart);
}