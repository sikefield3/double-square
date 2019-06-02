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

/*
 * An adapted version of the algorithm described in the following paper 
 * was used: 
 * Louis Marmet.  First occurrences of square-free gaps and an algorithm for
 * their computation. arXiv preprint arXiv:1210.3829, 2012.
*/

#ifndef SQFMARMET_H
#define SQFMARMET_H

#define SQFMARMET_USEVANILLAARRAY 0
#define SQFMARMET_USESTDARRAY 1
#define SQFMARMET_USESTDVECTOR 2
#define SQFMARMET_CONTAINERTYPE SQFMARMET_USESTDVECTOR

// #define SQFMARMET_POINTERS

#include <functional>
#include<array>
#include<vector>
#include<algorithm>
#include <cassert>
#include <string>
#include <sstream>

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/tools.h>

#include "global.h"
#include "sqfsorter.h"
#include "factor.h"

using namespace std;
using namespace NTL;

class sqfMarmet
{	
	sqfSorter m_oSqfSorter;
    Factor m_oFactor;
	
	static const size_t m_NP2minLen = 21;
	static const int m_improvement = 5;
	size_t m_cntAryFilledLen; // to store the # of proper elements in p2_ary
	int m_cntPrimeFactors = 1; // atm only 1 or 2 is supported, in general we'd need a tuple of exponents
	int m_k2, m_k3; // see [1]
	int m_gaplen = 1;
	int m_exp = 2; // exponent e: we are looking for sequences of numbers n = k x p^e, where p is prime
	// e=2 for the original squarefree gaps
	int m_head;
	UBNUM m_N = 0; // should be used as N in [1]
	UBNUM m_NLimit = 0;
	static constexpr int mc_bitveclen = 8 * sizeof(BITVECTOR);
	static constexpr int mc_bitvecZeroPos = mc_bitveclen / 2;
		
	// [1] = Marmet
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDVECTOR)	
	static const size_t m_cntAryLen = 10000000;
	vector<UBNUM> p2_ary, nqsf_ary, NP2min_ary, pow2_ary, next_ary;
	vector<BNUM> spRemainders;
	
#endif // (SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVECTOR)
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDARRAY)
	static const size_t m_cntAryLen = 10000;
	array<UBNUM, m_cntAryLen> p2_ary, nqsf_ary, next_ary ;
	array<UBNUM, m_NP2minLen> NP2min_ary;
	array<BNUM, m_NP2minLen> spRemainders ;
	array<UBNUM, mc_bitveclen> pow2_ary;
#endif // (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDARRAY)
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY)
	static const size_t m_cntAryLen = 1000000;
	UBNUM* p2_ary = new UBNUM  [m_cntAryLen];
	UBNUM* nqsf_ary = new UBNUM  [m_cntAryLen];
	UBNUM* next_ary = new UBNUM  [m_cntAryLen];	
	UBNUM NP2min_ary[m_NP2minLen];
	BNUM spRemainders[m_NP2minLen];
	UBNUM pow2_ary[mc_bitveclen];
#endif // #if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY)
	
private:
	enum enBatchMode {BMTimeLimit, BMNumberLimit };
	
	function<UBNUM(int, UBNUM, UBNUM)> func_getFirstGap();
	function<void(int, UBNUM)> func_fillNumArrays();
	function<bool(UBNUM)> func_ValueChecker();
	
	void initNumArrays();
	void fillpow2();
	void fillNumArrays();
	void fillNumArrays(int gaplen, UBNUM startVal);
	void fillNumArrays_improV(int gaplen, BNUM startVal);
	void fillNumArrays_improVII(int gaplen, BNUM startVal);
	UBNUM primeLimit(int e) const;
	UBNUM power (UBNUM n, int e) const;
	void fillNP2Min();
	void fillNPexpMin(int e);
	void fillNPexpMinDouble();
	UBNUM getLimitDoubleSquare(BNUM n) const;
	void fillNumArraysDouble(int gaplen, BNUM startVal);
	UBNUM getNextSqFul(UBNUM N, int offset = 0);
	UBNUM sortVII(UBNUM N, int offset = 0);
	UBNUM getNextSqFulGap_improV(UBNUM N, UBNUM endVal);
	UBNUM getNextSqFulGap_improVII(UBNUM N);	
	UBNUM prepareSmallPrimes(UBNUM startVal);
	UBNUM prepareSmallPrimesDouble(UBNUM startVal);
	BNUM getNextSqFulGap(BNUM N, BNUM LMin, bool bJump);
	BNUM getFirstGap_ver01(int l);
	BNUM getFirstGap(int l, BNUM startVal, UBNUM endVal);
	UBNUM getFirstGap_improV(int l, UBNUM startVal, UBNUM endVal);
	UBNUM getFirstGap_improVII(int l, UBNUM startVal, UBNUM endVal);
	void swapNeighbours(int j);
	void initNRemainders(UBNUM N);
	void add2Remainders(BNUM diff);
	void setNRemainders(UBNUM N);
	void initNSQFAry(BNUM startVal);
	void sortAry();
	void sortAlgo();
    bool checkValDouble (UBNUM n);
    bool checkValPower (UBNUM n);
	bool checkSequence (int gaplen, UBNUM N);
	void batchCalc(UBNUM startVal, UBNUM endVal, int mins, enBatchMode batchMode, int startNr);
	void prnTimeStamp(string sText) const;
	// debug functions
	void prv  (size_t cnt = 10)  const;
	std::__cxx11::string vec(size_t cnt = 10) const;
public:
    sqfMarmet();
	void exec();
	void batchCalc(UBNUM startVal, UBNUM endVal, int startNr = 0);
	void batchCalcTime(UBNUM startVal, int mins, int startNr = 0);
    
};

#endif // SQFMARMET_H
