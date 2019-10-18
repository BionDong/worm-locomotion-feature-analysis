/*
* Copyright 2019 Xianke Dong et al <xianke.dong@mail.mcgill.ca>
* This file is the image processing algorithm of the RoboWorm system.
* Specifically, it is developed to analyse the length/width/volume/
* curvature of a free moving C. elegans, and target an interested
* part on the worm body.
*
* This algorithm is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY for a perticular use. For the most up 
* to date version of this software, see:
* https://github.com/BionDong/worm-locomotion-feature-analysis
*
* NOTE: If you use any portion of this code in your research, kindly cite:
* Dong, Xianke, Pengfei Song, and Xinyu Liu. "An Automated 
* Microfluidic System for Morphological Measurement and Size-Based 
* Sorting of C. Elegans." IEEE transactions on nanobioscience (2019).
* and our following works on C. elegans microrobot.
*/

/*
*  Created on: Oct. 10, 2019
*      Author: Xianke (Bion) Dong
*/

/*
* This is the main file for the image processing algorithm to analyse.
* The core algorithms for the morphorlogic feature analysis of free
* moving C. elegans are enclosed in the class named "WormTrack". In
* this code, the usesage of the WormTrack is illustrated by analysing
* 49 images for continuous worm moving.
*/



#include "stdafx.h"
#include "Timer.h" //to analysis the performance of the algorithm
#include "WormTrack.h" //the class I developed for worm shape analysis

using namespace std;
using namespace cv;

int RefineCurvature(seg* segment, vector<Point>& cur);
int DrawImage(seg& segment, Mat& I, int partnum, vector<Point> cur);

/*functions for curvature feature analysis*/
int CurvatureAnalysis(vector<double>* curvature, vector<Point>& cur);
int RefineCurvature(seg* segment, vector<Point>& cur);
int FindMax(vector<double>* input, int start, int end);
int FindMin(vector<double>* input, int start, int end);
int FindPerpPoint(vector<Point>* input, Point x, Point targent, Point &result, int startindex, int endindex);
void Smooth1DSequence(const vector<double>* input, vector<double>& output, double sigma);

/*initialized parameters*/
Point m_A = Point(0, 0);
Point m_B = Point(2040, 1088);
Point head = Point(0, 0);
Point tail = Point(0, 0);

int main()
{
	myTimer Timer; // to record time cost
	seg segment; // to save analysis results
	
	// initiate parameters
	segment.head = head;
	segment.tail = tail;

	for (int j = 1; j <= 49; j++) // 49 continous images for instance
	{
		/*load images*/
		String n, b;
		n = to_string(j);
		b = n + ".jpg";
		String name(b);

	    Mat I = imread(name, 0); // for worm shape analysis
		Mat I_color = imread(name, 1); // to record the analysis result

		/*set parameters*/
		int guassian = 1;
		int threshold = 70;
		int morpology = 1;
		bool switchht = 0;
		int partnum = 120;

		/*create instant for worm shape analysis*/
		WormTrack Track(I, 1, guassian, threshold, morpology, partnum, m_A, m_B, &segment, switchht);
		
		Timer.StartTimer("time cost for one image"); // record time cost

		/*start image processing*/
		Track.PreProcess();
		Track.Contour();
		Track.Analysis();
		Track.Curvaturecalculation();
		Track.Lengthcalculation();
		Track.Widthcalculation();
		Track.Volumecalculation();
		Track.Segment();

		Timer.StopTimer();
		
		segment = Track.segment; // copy results

		/*feature analysis*/
		vector<Point> cur;
		RefineCurvature(&segment, cur); //2ms

		/*draw analysis results*/
		DrawImage(segment, I_color, partnum, cur);

		/*send notifications*/
		cout << "Pic " << j << " : length " << segment.length << " pixels; width " 
			<< segment.width << " pixels; volume: " << segment.volume << " cube pixels." << endl;

		//save results
		String a;
		a = n + "_.jpg";
		String c(a);
		imwrite(c, I_color);

		//update ROI
		head = segment.head;
		tail = segment.tail;
		m_A = Track.segment.topleft - Point(50, 50);
		m_B = Track.segment.botright + Point(50, 50);

		Track.~WormTrack();
	}

	waitKey();
	return 0;
}


int RefineCurvature(seg* segment, vector<Point>& result)
{
	if (segment->curvatureC.size() == 0)
	{
		return 0;
	}

	vector<Point> curC;
	vector<double> curvatureA, curvatureB, scurvA, scurvB;
	CurvatureAnalysis(&(segment->curvatureC), curC);

	if (curC.size() < 5)
	{
		return 0;
	}

	curvatureA.insert(curvatureA.begin(), (segment->curvatureA).begin() + 1, (segment->curvatureA).end() - 1);
	curvatureB.insert(curvatureB.begin(), (segment->curvatureB).begin() + 1, (segment->curvatureB).end() - 1);

	Smooth1DSequence(&(curvatureA), scurvA, 3); // 1 point difference to the curvature
	Smooth1DSequence(&(curvatureB), scurvB, 3); // 1 point difference to the curvature


	int a, b;  // search from a to b
	a = 0;

	int length = ((segment->Partnum / curC.size()) > 1 ? (segment->Partnum / curC.size()) : 1) * 1.5;

	Point backward, forward, tangent;
	int step = 10; // to find perpendicular vector on boundaries
	for (int i = 0; i < curC.size(); i++)
	{
		b = (a + length) < (segment->center.size() - 1) ? (a + length) : (segment->center.size() - 1);

		if (curC[i].x == 0)
		{
			result.push_back(curC[i]);
		}
		else if (curC[i].x == 1)
		{
			int start = (curC[i].y - length) > 0 ? (curC[i].y - length) : 0;
			int end = (curC[i].y + length) < scurvA.size() ? (curC[i].y + length) : scurvA.size();

			Point perp;
			int index = segment->tabA[1 + FindMax(&scurvA, start, end)].y;

			backward = segment->contourA[(index - step) > 0 ? (index - step) : 0];
			forward = segment->contourA[(index + step) < (segment->contourA.size() - 1) ? (index + step) : (segment->contourA.size() - 1)];
			tangent = forward - backward;

			result.push_back(Point(1, FindPerpPoint(&(segment->center), segment->contourA[index], tangent, perp, a, b)));
		}
		else
		{
			int start = (curC[i].y - length) > 0 ? (curC[i].y - length) : 0;
			int end = (curC[i].y + length) < scurvB.size() ? (curC[i].y + length) : scurvB.size();

			Point perp;
			int index = segment->tabB[1 + FindMin(&scurvB, start, end)].y;

			backward = segment->contourB[(index - step) > 0 ? (index - step) : 0];
			forward = segment->contourB[(index + step) < (segment->contourB.size() - 1) ? (index + step) : (segment->contourB.size() - 1)];
			tangent = forward - backward;

			result.push_back(Point(-1, FindPerpPoint(&(segment->center), segment->contourB[index], tangent, perp, a, b)));
		}


		a = result[i].y + length * 0.3;
	}

	return 1;
}

int FindPerpPoint(vector<Point>* input, Point x, Point targent, Point &result, int startindex, int endindex)
{
	if (abs((*input)[0].x - x.x) + abs((*input)[0].y - x.y) <1)
	{
		result = x;
		return 0;
	}
	else if (abs((*input)[input->size() - 1].x - x.x) + abs((*input)[input->size() - 1].y - x.y) <1)
	{
		result = x;
		return input->size() - 1;
	}

	int mindot = INT_MAX;
	int temp, index = startindex;
	Point pt;

	startindex = startindex >= 0 ? startindex : 0;
	startindex = startindex <= input->size() - 1 ? startindex : input->size() - 1;
	endindex = endindex <= input->size() - 1 ? endindex : input->size() - 1;

	// not accurate, since temp is a dot product of two vectors.

	// a1/b1 * a2/b2 = -1 -> temp = a1a2 + b1b2
	for (int i = startindex; i <= endindex; i++)
	{
		pt = (*input)[i] - x;
		temp = pt.x * targent.x + pt.y * targent.y;
		temp = abs(temp);
		if (mindot > temp)
		{
			mindot = temp;
			index = i;
		}
	}

	result = (*input)[index];
	return index;
}


int FindMax(vector<double>* input, int start, int end)
{
	double Max = -100;
	int indexMax = start;

	for (int i = start; i < end; i++)
	{
		if ((*input)[i] > Max)
		{
			indexMax = i;
			Max = (*input)[i];
		}
	}
	return indexMax;
}

int FindMin(vector<double>* input, int start, int end)
{
	double Min = 100;
	int indexMin = start;

	for (int i = start; i < end; i++)
	{
		if ((*input)[i] < Min)
		{
			indexMin = i;
			Min = (*input)[i];
		}
	}
	return indexMin;
}

// for curvature analysis
int CurvatureAnalysis(vector<double>* curvature, vector<Point>& cur)
{
	if (curvature->size() == 0)
	{
		return 0;
	}

	// smooth sequence
	vector<double> curv, scurv, super_scurv;
	curv.insert(curv.begin(), curvature->begin() + 1, curvature->end() - 1);

	Smooth1DSequence(&curv, scurv, 1);   // 1 point difference to the curvature
	Smooth1DSequence(&curv, super_scurv, 3); // 1 point difference to the curvature

	float omit =0.08;
	for (int i = 0; i < curv.size()*omit; i++)
	{
		scurv[i] = 0;
		super_scurv[i] = 0;
	}
	for (int i = curv.size()*(1-omit); i < curv.size(); i++)
	{
		scurv[i] = 0;
		super_scurv[i] = 0;
	}


	vector<int> pre_Zeros, Zeros;
	vector<int>	pre_Max, pre_Min;
	Point Headpeak = Point(0, 0);
	Point Tailpeak = Point(0, 0);


	// find the zeros
	for (int i = 0; i < scurv.size() - 1; i++)
	{
		if (scurv[i] == 0 || scurv[i + 1] * scurv[i] <0)  // find all points on zeros
		{
			pre_Zeros.push_back(i);

		}

	}

	// find peaks in the center body
	for (int i = 0; i < pre_Zeros.size() - 1; i++)
	{
		if (pre_Zeros[i + 1] - pre_Zeros[i] > curv.size() / 6 && scurv[pre_Zeros[i] + 1] > 0) // maximum
		{
			double Max = 0;
			int indexMax = pre_Zeros[i] + 1;
			for (int j = pre_Zeros[i] + 1; j < pre_Zeros[i + 1] - 1; j++)
			{
				if (Max < super_scurv[j])  // super_scure is more smooth
				{
					Max = super_scurv[j];
					indexMax = j;
				}
			}
			pre_Max.push_back(indexMax);

			// store zeros
			if (Zeros.size() == 0)
			{
				Zeros.push_back(pre_Zeros[i]);
				Zeros.push_back(pre_Zeros[i + 1]);
			}
			else
			{
				if (Zeros[Zeros.size() - 1] != pre_Zeros[i])
				{
					Zeros.push_back(pre_Zeros[i]);
					Zeros.push_back(pre_Zeros[i + 1]);
				}
				else
				{
					Zeros.push_back(pre_Zeros[i + 1]);
				}
			}

		}

		if (pre_Zeros[i + 1] - pre_Zeros[i] > curv.size() / 6 && scurv[pre_Zeros[i] + 1] < 0) // 
		{
			double Min = 100;
			int indexMin = pre_Zeros[i] + 1;
			for (int j = pre_Zeros[i] + 1; j < pre_Zeros[i + 1] - 1; j++)
			{
				if (Min > super_scurv[j])  // super_scure is more smooth
				{
					Min = super_scurv[j];
					indexMin = j;
				}
			}
			pre_Min.push_back(indexMin);
			// store zeros
			if (Zeros.size() == 0)
			{
				Zeros.push_back(pre_Zeros[i]);
				Zeros.push_back(pre_Zeros[i + 1]);
			}
			else
			{
				if (Zeros[Zeros.size() - 1] != pre_Zeros[i])
				{
					Zeros.push_back(pre_Zeros[i]);
					Zeros.push_back(pre_Zeros[i + 1]);
				}
				else
				{
					Zeros.push_back(pre_Zeros[i + 1]);
				}
			}

		}
	}


	// find minimum and maximum before the first zero point

	if (Zeros.size() != 0)
	{
		if (Zeros[0]>0)
		{
			if (scurv[Zeros[0] - 1] >0)
			{
				double Max = 0;
				int indexMax = 1;
				for (int j = 1; j < Zeros[0] - 1; j++)
				{
					if (Max < super_scurv[j])  // super_scure is more smooth
					{
						Max = super_scurv[j];
						indexMax = j;
					}
				}
				Headpeak = Point(1, indexMax);
			}
			else
			{
				double Min = 100;
				int indexMin = 0;
				for (int j = 1; j < Zeros[0] - 1; j++)
				{
					if (Min > super_scurv[j])  // super_scure is more smooth
					{
						Min = super_scurv[j];
						indexMin = j;
					}
				}
				Headpeak = Point(-1, indexMin);
			}
		}
	}


	// find minimum and maximum after the last zero point
	if (Zeros.size() != 0)
	{
		if (Zeros[Zeros.size() - 1] < scurv.size() - 1)
		{
			if (scurv[Zeros[Zeros.size() - 1] + 1] >0)
			{
				double Max = 0;
				int indexMax = Zeros[Zeros.size() - 1] + 1;
				for (int j = Zeros[Zeros.size() - 1] + 1; j < scurv.size() - 2; j++)
				{
					if (Max < super_scurv[j])  // super_scure is more smooth
					{
						Max = super_scurv[j];
						indexMax = j;
					}
				}
				Tailpeak = Point(1, indexMax);
			}
			else
			{
				double Min = 100;
				int indexMin = Zeros[Zeros.size() - 1] + 1;
				for (int j = Zeros[Zeros.size() - 1] + 1; j < scurv.size() - 2; j++)
				{
					if (Min > super_scurv[j])  // super_scure is more smooth
					{
						Min = super_scurv[j];
						indexMin = j;
					}
				}
				Tailpeak = Point(-1, indexMin);
			}

		}
	}


	//store headpeak and tailpeak
	if (Headpeak != Point(0, 0) && Zeros.size() != 0)
	{
		if (Zeros[0] - Headpeak.y > scurv.size() / 8)
		{
			if (Headpeak.x > 0)
			{
				pre_Max.insert(pre_Max.begin(), Headpeak.y);
			}
			else
			{
				pre_Min.insert(pre_Min.begin(), Headpeak.y);
			}
		}
	}

	if (Tailpeak != Point(0, 0) && Zeros.size() != 0)
	{
		if (Tailpeak.y - Zeros[Zeros.size() - 1] > scurv.size() / 8)
		{
			if (Tailpeak.x > 0)
			{
				pre_Max.push_back(Tailpeak.y);
			}
			else
			{
				pre_Min.push_back(Tailpeak.y);
			}
		}
	}

	for (int i = 0; i < pre_Max.size(); i++)
	{
		cur.push_back(Point(1, pre_Max[i] + 1));
	}
	for (int i = 0; i < pre_Min.size(); i++)
	{
		cur.push_back(Point(-1, pre_Min[i] + 1));
	}
	for (int i = 0; i < Zeros.size(); i++)
	{
		cur.push_back(Point(0, Zeros[i] + 1));
	}

	if (cur.size() != 0)
	{
		for (int i = 0; i < cur.size() - 1; i++)
		{ // times
			for (int j = 0; j < cur.size() - i - 1; j++)
			{ // position
				if (cur[j].y > cur[j + 1].y)
				{
					Point temp = cur[j];
					cur[j] = cur[j + 1];
					cur[j + 1] = temp;
				}
			}
		}
	}

	return 1;
}

void Smooth1DSequence(const vector<double>* input, vector<double>& output, double sigma)
{
	int *kernel, klength, normfactor;

	int ll, ul, x;
	double n;
	ll = (int)(-3 * sigma) - 1;
	ul = (int)(3 * sigma) + 1;
	klength = ul - ll + 1;
	kernel = (int*)malloc(klength * sizeof(int));

	normfactor = 0;
	n = exp(-1.0*ll*ll / (2 * sigma*sigma));
	for (x = 0; x < klength; x++)
	{
		(kernel)[x] = (int)(exp(-1.0*(x + ll)*(x + ll) / (2 * sigma*sigma)) / n + 0.5);
		normfactor += (kernel)[x];
	}

	/// filtering
	int j, k, ind, anchor;
	double sum;

	double *pt;
	int length = input->size();
	pt = (double*)malloc(length * sizeof(Point));
	anchor = klength / 2;
	for (j = 0; j < length; j++)
	{
		sum = 0;
		for (k = 0; k < klength; k++)
		{
			ind = j + k - anchor;
			ind = ind > 0 ? ind : 0;
			ind = ind < length ? ind : (length - 1);
			sum = sum + (*input)[ind] * kernel[k];
		}

		pt[j] = (double)(1.0*sum / normfactor);

	}
	output = vector<double>(pt, pt + length);

	free(pt);

}

int DrawImage(seg& segment, Mat& I, int partnum, vector<Point> cur)
{
	//parameters
	float offsite = 0.2;
	float extension = 1;

	if (segment.center.size() == 0)
	{
		return 0;
	}

	vector<Point> newCenter, newContourA, newContourB;
	for (int i = 0; i < segment.center.size(); i++)
	{
		newCenter.push_back(segment.center[i] - m_A);
	}

	for (int i = 0; i < segment.contourA.size(); i++)
	{
		newContourA.push_back(segment.contourA[i] - m_A);
	}

	for (int i = 0; i < segment.contourB.size(); i++)
	{
		newContourB.push_back(segment.contourB[i] - m_A);
	}

	// create small ROI     0.5 ms
	Mat roi = I(cv::Rect(m_A, m_B));
	Mat overlay;
	roi.copyTo(overlay);

	// generate pattern
	vector<Point> controlA, controlB;
	for (int i = 0; i < partnum; i++)
	{
		if (i % 3 == 1)
		{
			controlA.push_back(Point(i, 0));
			controlB.push_back(Point(i, 0));
		}
		else
		{
			controlA.push_back(Point(i, 255));
			controlB.push_back(Point(i, 255));
		}
	}

	// draw pattern on small ROI      5 ms
	vector<Point> SegContours;

	Point Int, Ext;
	vector<Point> SegCountours, BoundA, BoundB;
	for (int i = 0; i < segment.Partnum; i++)      // combine the conencted countours to a single one
	{
		if (controlA[i].y != 0)
		{
			do{
				Int = offsite*(newContourA[segment.segAB[i].x] - newCenter[i]) + newCenter[i];
				SegCountours.push_back(Int);

				Ext = extension*(newContourA[segment.segAB[i].x] - newCenter[i]) + newCenter[i];
				BoundA.push_back(Ext);

				i++;
				if (i == segment.Partnum)
					break;
			} while (controlA[i].y == controlA[i - 1].y);

			Int = offsite*(newContourA[segment.segAB[i].x] - newCenter[i]) + newCenter[i];
			SegCountours.push_back(Int);    // store one more point to make the pattern accurate

			Ext = extension*(newContourA[segment.segAB[i].x] - newCenter[i]) + newCenter[i];
			BoundA.push_back(Ext);
			SegCountours.insert(SegCountours.end(), BoundA.rbegin(), BoundA.rend());

			fillConvexPoly(overlay, SegCountours, Scalar(0, 255, 0), 8, 0);

			SegCountours.clear();
			BoundA.clear();
		}
	}

	for (int i = 0; i < segment.Partnum; i++)
	{
		if (controlB[i].y != 0)
		{
			do{
				Int = offsite*(newContourB[segment.segAB[i].y] - newCenter[i]) + newCenter[i];
				SegCountours.push_back(Int);

				Ext = extension*(newContourB[segment.segAB[i].y] - newCenter[i]) + newCenter[i];

				BoundB.push_back(Ext);
				i++;
				if (i == segment.Partnum)
					break;
			} while (controlB[i].y == controlB[i - 1].y);

			Int = offsite*(newContourB[segment.segAB[i].y] - newCenter[i]) + newCenter[i];
			SegCountours.push_back(Int);    // store one more point to make the pattern accurate

			Ext = extension*(newContourB[segment.segAB[i].y] - newCenter[i]) + newCenter[i];

			BoundB.push_back(Ext);
			SegCountours.insert(SegCountours.end(), BoundB.rbegin(), BoundB.rend());

			fillConvexPoly(overlay, SegCountours, Scalar(0, 255, 0), 8, 0);

			SegCountours.clear();
			BoundB.clear();
		}
	}

	// draw left and right sides      3.5 ms
	for (int i = 0; i < segment.contourA.size() - 1; i++)
	{
		line(overlay, newContourA[i], newContourA[i + 1], Scalar(0, 255, 255), 2, 8);
	}
	for (int i = 0; i < segment.contourB.size() - 1; i++)
	{
		line(overlay, newContourB[i], newContourB[i + 1], Scalar(0, 255, 255), 2, 8);
	}

	// draw head and tail
	circle(overlay, segment.head - m_A, 5, Scalar(0, 0, 255), 3, 8, 0);
	circle(overlay, segment.tail - m_A, 3, Scalar(255, 0, 0), 3, 8, 0);

	// plot result
	for (int i = 0; i < cur.size(); i++)
	{
		if (cur[i].x == 0)
		{
			circle(overlay, segment.center[cur[i].y] - m_A, 5, Scalar(0, 0, 255), 2, 8, 0);
		}
		else if (cur[i].x == -1)
		{
			circle(overlay, segment.center[cur[i].y] - m_A, 5, Scalar(0, 255, 255), 2, 8, 0);
		}
		else
		{
			circle(overlay, segment.center[cur[i].y] - m_A, 5, Scalar(0, 255, 255), 2, 8, 0);
		}
	}

	// blend image    2 ms
	float alpha = 0.5;
	addWeighted(overlay, alpha, roi, 1 - alpha, 0, roi);

	return 1;
}

