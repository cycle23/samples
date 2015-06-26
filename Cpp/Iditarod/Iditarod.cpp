// Iditarod.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <sstream>

// day 1: 8:00 AM
// hh:mm xM, DAY n

// hh = hour
// mm = minutes
// x = 'A' or 'P'
// n = int < 100 (no leading zeros)
// string len = 15 (if n < 10) or 16 (if n >= 10)

// round to the nearest minute, 0.5 rounds up

using namespace std;

// 12:00 AM, DAY d" refers to midnight between DAY d-1 and DAY d
// ie: midnight day 2 = exactly 16 hours

// 12:00 PM, DAY d" refers to noon on DAY d

// times = 1 to 50 elements
int avgMinutes(vector<string> times) {
	int totalMinutes = 0;
	int i = 0;
	for (; i < times.size(); i++) {
		int hour;
		char colon;
		int min;
		char x;
		string mComma;
		string Day;
		int day;
		//if (4 == sscanf_s(times[i].c_str(), "%d:%d %cM, DAY %d", &hour, &min, &x, 1, &day)) {
		istringstream iss(times[i]);
		if ((iss >> hour >> colon >> min >> x >> mComma >> Day >> day) 
			&& (colon == ':' && !mComma.compare("M,") && !Day.compare("DAY"))) {
			if (x == 'P' && hour < 12) {
				hour += 12;
			}
			if (x == 'A' && hour == 12) {
				hour = 0;
			}
			if (day == 1) {
				// remove the first hours if we are in a 1 day scenario
				hour -= 8;
			}
			else {
				// account for rest of day 1
				totalMinutes += 16 * 60;
			}
			// add total days for every day past day 2
			if (day > 2) {
				totalMinutes += (day - 2) * 24 * 60;
			}
			// add the remaining time in the final day
			totalMinutes += ((hour * 60) + min);
		}
		else {
			cout << "could not parse: " << times[i] << endl;
			cout << "hour " << hour << ", min " << min << ", colon " << colon << ", Day " << Day << endl;
		}
	}
	// round up..
	return round((double)totalMinutes / (double)i);
}


int _tmain(int argc, _TCHAR* argv[])
{
	int avg = avgMinutes({"02:00 PM, DAY 19", "02:00 PM, DAY 20", "01:58 PM, DAY 20"});
	cout << avg << endl;
	return 0;
}

