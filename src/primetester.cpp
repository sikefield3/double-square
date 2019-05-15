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

#include "primetester.h"

bool PrimeTester::isPrimeTrial(BNUM n) const{
	if (n == 2 or n == 3 or n == 5){
		return true;
	}
	if (n % 2 == 0 or n % 3 == 0 or n % 5 == 0){
		return false;
	}
	int nStep = 2;
	for (BNUM j=7;j*j <= n;j+=nStep){
		if (n % j == 0){
			return false;
		}
		nStep = nStep * 2 % 6;
	}
	return true;
}
// we assume p == 2 or P % 6 == ???
// we don't check for limits!
BNUM PrimeTester::NextPrime_PT(BNUM p) const{
	if (p < 2){
		return 2;
	}
	p += p % 2 + 1;
	BNUM j;
	for (j = p;not isPrimeTrial(j);j += 2);
	return j;	
}

void PrimeTester::testProbPrime(BNUM N) const{
	bool bProbPrime, bPrimeTrial;	
	for (BNUM j = 3; j <= N; j++){
		if ((bProbPrime = ProbPrime(j,3)) != (bPrimeTrial = isPrimeTrial(j))){
			cout << j << ": Probprime = " << bProbPrime << ", Prime Trial = " << bPrimeTrial << endl;
		}
	}
	cout << "Finished: " << N << endl;
	// Time tests
	double dStartTime = GetTime();
	BNUM cnt = 1;
	for (BNUM j = 3; j <= N; j++){
		if (isPrimeTrial(j)){
			cnt ++;
		}
	}
	cout << "Time for trial division: " << GetTime() - dStartTime << " secs" << endl;
	cout << "Number of primes: " << cnt << endl;

	dStartTime = GetTime();
	cnt = 1;
	for (BNUM j = 3; j <= N; j++){
		if (ProbPrime(j)){
			cnt ++;
		}
	}
	cout << "Time for Miller witness: " << GetTime() - dStartTime << " secs" << endl;
	cout << "Number of primes: " << cnt << endl;
	
}
void PrimeTester::testNextPrime(BNUM N) const{
	BNUM j;
	BNUM p = 1;
	double dStartTime = GetTime();
	for (j=0;j<= N;j++)	{
		p = NextPrime_PT(p);
		// cout << j << ": " << p << endl;
	}
	cout << j << ": " << p << endl;
	cout << "Time: " << GetTime() - dStartTime << " secs" << endl;
}
// we assume that the NTL function prime never returns a composite
// we find: test with limit = 10^9, NumTrials default = 10 -> no false primes (3247 sec)
BNUM PrimeTester::testCheckNextPrimeForComposite(BNUM NLimit, int NumTrials ) const{
	std::function<BNUM (BNUM, int)> funcNextPrime;
	if (NumTrials > 0){
		      funcNextPrime = [this] (BNUM parP, int parNumTrials) {return NextPrime (parP, parNumTrials);};
	} else {
		      funcNextPrime = [this] (BNUM parP, int parNumTrials) {return NextPrime (parP);};
	}
	double dStartTime = GetTime();
	BNUM p = 1;
	cout << "testCheckNextPrimeForComposite" << endl;
	while (p < NLimit){
		p = funcNextPrime(p+1, NumTrials);
		if (not isPrimeTrial(p)){
			cout << p << " is a composite !" << endl << "We must get a random seed!" << endl;
			cout << "Problem: NTL uses the current time as seed and if we wanted a fixed seed," << endl;
			cout << "Shoup had to implement it" << endl;
			cout << "Time: " << GetTime() - dStartTime << " secs" << endl;
			return p;
		}
	}
	cout << "testCheckNextPrimeForComposite: allrighty!" << endl;
	cout << "Time: " << GetTime() - dStartTime << " secs" << endl;
	return NLimit;
}
