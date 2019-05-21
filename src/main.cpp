/* 
 * Double Squareful - Numbertheoretic program to find sequences related to 
 * prime factors / squarefree numbers.
 * Copyright (C) 2018-19 Bernd Zemann - zb4ng@arcor.de
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include <NTL/ZZ.h>

#include "global.h"
#include "sqf.h"
#include "sqfmarmet.h"
#include "primetester.h"
#include "factor.h"
#include "checkresult.h"

using namespace std;
using namespace NTL;

void callSQF(){
    SQF oSQF;
    oSQF.exec();
}
void callSQFMarmet(){
    sqfMarmet oSQFMarmet;
    oSQFMarmet.exec();
}
void checkNTLPrimeFunc(){
	PrimeTester oPTest;
	// oPTest.testProbPrime(100000000);
//	oPTest.testProbPrime(10000000);
// 	oPTest.testNextPrime(1000000);
	
    BNUM exp, limit;    

    exp = 9;
    limit = power_long(10,exp);
    cout << "10^" << exp << "=" << limit << endl;
	
	oPTest.testCheckNextPrimeForComposite	(limit);
}
void testOverFlow(){
	UBNUM n;
	for(n = 2;;n *= 2){
		cout << n << endl;
	}
}
void testFactor(){
	Factor oFactor;	
	UBNUM n;
	cout << "Factor test" << endl;
	UBNUM nStart = 717688804748, cnt = 8;
	
	double dStartTime = GetTime();
	for (n = nStart; n < nStart + cnt; n++){
		cout << n << " = " << oFactor.prnFactor(n) << endl;
		UBNUM rad = oFactor.radical();
		cout << "radical: " << rad << endl;
		cout << "omega = " <<  oFactor.LittleOmega(n / rad) << endl;
		// oFactor.doFactor(n);
	}
	cout << GetTime() - dStartTime << " secs" << endl;
	return;
	dStartTime = GetTime();
	for (n = nStart; n < nStart + cnt; n++){
		oFactor.doFactorNP(n);
	}
	cout << GetTime() - dStartTime << " secs" << endl;
	return;

	string s,sNP;
	for (n = nStart; n < nStart + cnt; n++){
		s = oFactor.prnFactor(n);
		sNP = oFactor.prnFactorNP(n);
		if (s != sNP){
			cout << "Error for n = " << n << endl;
			cout << s << endl;
			cout << sNP << endl;
			return;
		}
	}
}
void callCheckResult(){
	CheckResult cr("doublesqf_list.txt", ChDoubleSquareFree, 5);
	cr.run();
}
int main(int argc, char **argv) {
	string sVersion = "Version: 1.0.1, Date: 21.05.2019 ";
	cout << sVersion << endl;
	// callSQF();	
 	callSQFMarmet();
// 	checkNTLPrimeFunc();
//	testOverFlow();
// 	testFactor();
	// callCheckResult();
    return 0;
}
