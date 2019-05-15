/*
 * Double Squareful - Numbertheoretic program to find sequences related to 
 * prime factors / squarefree numbers.
 * Copyright (C) 2018-19 Bernd Zemann - zb4ng@arcor.de
 * 
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2018  Bernd Zemann <email>
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

#include "sqf.h"

SQF::SQF()
{

}
std::__cxx11::string SQF::str(const NTL::ZZ& n) const
{   
    ostringstream oSTR;
    oSTR << n;
    return oSTR.str();
}
void SQF::makePrimorial(BNUM n){
	BNUM p = 2;
	prm = 2;
	double dStartTime = GetTime();
	while (p < n){
		p = NextPrime(p+1);
		prm *= p;		
	}
	cout << GetTime() - dStartTime << endl;
	// cout << prm << endl << "#Bytes: " << NumBytes(prm) << endl;
	prmLimit = n;
}
ZZ SQF::cbrt(const ZZ& n) const{ // -> m, s.t. m^3 <= n < (m+1)^3
	// cout << NumBits(n);
	ZZ rt = power_ZZ (2,(NumBits(n) - 1) / 3);
	// cout << rt << endl;
	ZZ smd = rt / 2;
	ZZ testrt;
	while (smd != 0){
		testrt = rt + smd;
		if (power (testrt,3) <= n){
			rt = testrt;
		}
		smd /= 2;
	}
	return rt;
}
BNUM SQF::cbrt(BNUM n) const{
	return conv <int>( cbrt(ZZ(n)));	
}
void SQF::testcbrt() {
	cout << 1 << " cbrt " << cbrt(1) << endl;
	cout << 5 << " cbrt " << cbrt(5) << endl;
	cout << 8 << " cbrt " << cbrt(8) << endl;
	cout << 9 << " cbrt " << cbrt(9) << endl;
	cout << 27 << " cbrt " << cbrt(27) << endl;
}
void SQF::testsquarefree()
{
	const BNUM testary [11] = {3,4,7,8, 34, 36, 18, 125, 54,210, 242};
	for (int j = 0; j<11;j++){
		cout << testary[j] << " squarefree " << isSquareFree(testary[j]) << endl;
	}
}
void SQF::testdoublesquarefree(){
	const BNUM testary [15] = {3,4,8, 34, 36, 18, 125, 54,210,900,72,144,27000};
	for (int j = 0; j<15;j++){
		cout << testary[j] << " doublesquarefree " << isDoubleSquareFree(ZZ(testary[j])) << endl;
	}	
	for (ZZ n (2595) ; n < 2610;n++){
        cout << n << " doublesquarefree " << isDoubleSquareFree(n) << endl;
    }
}
bool SQF::isSquareFree(const ZZ& n) const{
	if (n % 4 == 0){
		return false;
	}
	ZZ rmd ( (n % 2) ? n : n/2);
	BNUM d = 3;
    ZZ cbrt_rmd = cbrt(rmd);
	while (d <= cbrt_rmd){
		if (rmd % d == 0){
			if (rmd % (d*d) == 0){
				return false;
			}
			rmd /= d;
		}
		d += 2;
	}
	ZZ rt = SqrRoot(rmd);
	return (rt*rt != rmd);
}
bool SQF::isSquareFree(BNUM n) const{	
	return isSquareFree(ZZ(n));
}
// CAVEAT: we assume n < prmLimit^3 (!!!)
// SAGE: [(1, 36), (2, 675), (3, 625974)]
// [0 36 675 625974 364498649]
bool SQF::isDoubleSquareFree_ver1(const ZZ& n) const{
    if (n == 1150){
        // it gets in in there twice !!! ???
        // you can also use "print str(n)" in GDB !
        cout << "n: " << str(n) << endl;
    }
        
	ZZ rad = GCD(prm, ZZ(n));
	ZZ sqfac = n/rad;
	if (sqfac == 1 or n == sqfac){ // n == sqfac means n has no prime divisor < prmlimit 
		return true;
	}
	ZZ gcdrad = GCD(sqfac, rad);
	if (gcdrad > prmLimit){
		return false;
	}
	if (ProbPrime(gcdrad)){
		while(divide(sqfac, sqfac, gcdrad));
		// all div of sqfac are > prmlimit
		ZZ rt = SqrRoot(sqfac);
		return sqfac == 1 or rt*rt != sqfac;
	} else {
		return false;
	}
}
// CAVEAT: we assume n < prmLimit^2 (!!!)
bool SQF::isDoubleSquareFree(const ZZ& n) const{
    ZZ n1 = GCD (prm, n);
    ZZ n2 = GCD (prm, n/n1); // n2 is the product of the "square" prime div. (< prlimit)
    if (n2 == 1){
        return true;        
    }
    return ProbPrime(n2); // return true if n2 is prime    
}
void SQF::loop(const ZZ& nMax, std::function<bool(const ZZ&)> funcFree){
	int L=0;
	ZZ n(3), nSeqMin, nk;
	bool bIsSequence, bBreak;
	int k, j, i, nCurSeqLen;
	BNUM nLong;
	vZZ.kill();
	vZZ.append(ZZ(0));
	
	while(n<nMax){
		nLong = conv<int> (n);
		if (nLong >= 6){
			bBreak = true;
		}
		// look behind
		for (k = 0, bIsSequence = true, nSeqMin = n;k <= L and bIsSequence;k++){
			nk = n-k;
			bIsSequence = not funcFree(nk);
			nSeqMin = (bIsSequence) ? nk : nSeqMin;
		}
		if (bIsSequence or k>1){ // found sth
			// look ahead
			for (j = 1, bIsSequence = true;bIsSequence;j++){
				bIsSequence = not funcFree(n+j);
			}
			nCurSeqLen = k+j-2 - (L>0);
			if (nCurSeqLen > L){
				for (i=0;i<nCurSeqLen-L;i++){
					vZZ.append(nSeqMin);
				}
				L = nCurSeqLen;			
			}
			n += j-1;			
		} else {
			n += (L > 0) ? L : 1;
		}
	}
}
void SQF::exec(){
	makePrimorial(50000);
    
	// testsquarefree();
    // testdoublesquarefree();
    
    
    double dStartTime = GetTime();
	ZZ cnt = power_ZZ (10,8);
	loop (cnt, [this] (const ZZ& n) {return isSquareFree(n);}); // find a way to use only the function name as parameter
    // loop (cnt, [this] (const ZZ& n) {return isDoubleSquareFree(n);}); // find a way to use only the function name as parameter
	cout << vZZ << endl;
    cout << GetTime() - dStartTime << " secs" << endl;
	
}
