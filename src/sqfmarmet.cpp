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

#include "sqfmarmet.h"


sqfMarmet::sqfMarmet() {
    initNumArrays();
}
void sqfMarmet::initNumArrays() {
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDVECTOR)	
    p2_ary.reserve(m_cntAryLen);
    nqsf_ary.reserve(m_cntAryLen);
	if (m_improvement > 5){
		next_ary.reserve(m_cntAryLen);
	}	
    NP2min_ary.reserve(m_NP2minLen);
    spRemainders.reserve(m_NP2minLen);
    pow2_ary.reserve(mc_bitveclen);
#endif // #if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDVECTOR)
}

std::__cxx11::string sqfMarmet::vec(size_t cnt) const {
    ostringstream oSTR;
    for (size_t j = 0; j < cnt; j++) {
        oSTR << p2_ary[j]  << " ";
    }
    oSTR << endl;
    for (size_t j = 0; j < cnt; j++) {
        oSTR << nqsf_ary[j] << " ";
    }
    oSTR << endl;
    return oSTR.str();
}
// this method returns the fillNumArrays method according to m_improvement
function<void(int, UBNUM)> sqfMarmet::func_fillNumArrays(){
	if (m_cntPrimeFactors > 1){
		static_assert (m_improvement == 5);
		return [this] (int l, UBNUM startVal) {return fillNumArraysDouble(l, startVal);};
	}
	switch (m_improvement){
		case 1:
		case 2:
			return [this] (int l, UBNUM startVal) {return fillNumArrays(l, startVal);};
		case 3:
		case 4:
		case 5:
			return [this] (int l, UBNUM startVal) {return fillNumArrays_improV(l, startVal);};
		case 6:
		case 7:
			return [this] (int l, UBNUM startVal) {return fillNumArrays_improVII(l, startVal);};
	}	
}
// this method returns the getFirstGap method according to m_improvement
function<UBNUM(int, UBNUM, UBNUM)> sqfMarmet::func_getFirstGap(){
	switch (m_improvement){
		case 1:
		case 2:
			return [this] (int l, UBNUM startVal, UBNUM  endVal) {return getFirstGap(l, startVal, endVal);};
		case 3:
		case 4:
		case 5:
			return [this] (int l, UBNUM startVal, UBNUM  endVal) {return getFirstGap_improV(l, startVal, endVal);};
		case 6:
		case 7:
			return [this] (int l, UBNUM startVal, UBNUM  endVal) {return getFirstGap_improVII(l, startVal, endVal);};
	}
}
function<bool(UBNUM)> sqfMarmet::func_ValueChecker(){
	if (m_cntPrimeFactors > 1){
		static_assert (m_improvement == 5);
		return [this] (UBNUM n) {return checkValDouble(n);};
	}
	// plain powers
	return [this] (UBNUM n) {return checkValPower(n);};	
}
void sqfMarmet::fillNumArrays() {
    fillNP2Min();
    UBNUM p = 0;
    for (int j = 0; j < m_cntAryLen; j++) {
        p = NextPrime(p+1);
        nqsf_ary[j] = p2_ary[j] = p*p;
    }
    m_cntAryFilledLen = m_cntAryLen;
}
// fill arrays/vectors for single squareful algo. starting at "startVal"
// -> values in nqsf_ary >= startVal
void sqfMarmet::fillNumArrays(int gaplen, UBNUM startVal) 
{    
	m_gaplen = gaplen;
	m_cntAryFilledLen = m_cntAryLen;
	if (startVal == 1){
		fillNumArrays();
		return;
	}
	fillNP2Min();
    UBNUM p, p2;
    int j, k;
    p = 0;
    for (j = 0; j < m_cntAryLen; j++) {
        p = NextPrime(p+1);
        p2 = p*p;
        p2_ary[j] = p2;
        nqsf_ary[j] = (startVal / p2 + 1) * p2;
        // ensure that arrays are ascending
        k = j;
        while (k > 0 and nqsf_ary[k] < nqsf_ary[k-1]) {
            swapNeighbours(k-1);
            k--;
        }
    }
}
void sqfMarmet::fillNumArrays_improV(int gaplen, BNUM startVal)
{
    UBNUM p,p2,pLimit;
    int j,k;
    m_gaplen = gaplen;    
	fillNPexpMin(m_exp);
    fillpow2();
    m_k2 = NP2min_ary[gaplen]-2;
    p = prepareSmallPrimes(startVal);
	pLimit = primeLimit(m_exp);
	m_NLimit = power(pLimit, m_exp);
    for (j = m_k2; j < m_cntAryLen and p < pLimit; j++) {
        p = NextPrime(p+1);
		p2 = power(p, m_exp);
        p2_ary[j] = p2;
		nqsf_ary[j] = (startVal / p2 + 1) * p2;
        // ensure that arrays are ascending
        k = j;
        while (k > 0 and nqsf_ary[k] < nqsf_ary[k-1]) {
            swapNeighbours(k-1);
            k--;
        }
    }
    cout << "last prime: " << p << ", limit: " << m_NLimit << endl;
	m_cntAryFilledLen = m_cntAryLen;
    return;
}

void sqfMarmet::fillNumArrays_improVII(int gaplen, BNUM startVal){
    UBNUM p,p2,pLimit;
    int j,k;
#ifndef SQF_STARTVAL
	if (startVal > 1){
		throw invalid_argument ("For now only starting value = 1 is supported!");
	}
#endif //  SQF_STARTVAL
    m_gaplen = gaplen;
	m_cntAryFilledLen = m_cntAryLen;
	fillNPexpMin(m_exp);
    fillpow2();
    m_k2 = NP2min_ary[gaplen]-2;
	m_k3 = NP2min_ary[gaplen]-3; // not yet used
    p = prepareSmallPrimes(startVal);
	pLimit = primeLimit(m_exp);
	m_NLimit = power(pLimit, m_exp);
	m_head = m_k2;
	next_ary[m_head] = m_k2;
    for (j = m_k2; j < m_cntAryLen and p < pLimit; j++) {
        p = NextPrime(p+1);
		p2 = power(p, m_exp);
		nqsf_ary[j] = p2_ary[j] = p2;        
#ifdef SQF_STARTVAL
		nqsf_ary[j] = (startVal / p2 + 1) * p2;
		UBNUM tempnsqf = nqsf_ary[j];
		int i;
		// we have to check for the end of the chained list!
		for (i = m_head; tempnsqf>=nqsf_ary[next_ary[i]];){
			i=next_ary[i];
		}
		next_ary[i] = j;
		for (;i < j;i++){ 
			next_ary[i]++;
		}
#else // // SQF_STARTVAL
		next_ary[j] = j + 1;
#endif // SQF_STARTVAL
	}    
    return;
}
UBNUM sqfMarmet::power(UBNUM n, int e) const{
	UBNUM res;
	int j;
	for (res=n, j=0;j < e-1;j++,res *= n);
	return res;
}
// get a (not too sharp) limit for the prime s.t. p^e doesn't give an overflow
// 
UBNUM sqfMarmet::primeLimit(int e) const{
	int explim = (NextPowerOfTwo(MAXUBNUM / 4)+1) / e;
	return to_long(power2_ZZ(explim));
}

void sqfMarmet::fillpow2() {
    int j;
    BITVECTOR pw2;
    for (j = 0, pw2 = 1; j < mc_bitveclen; j++, pw2 *= 2) {
        pow2_ary[j] = pw2;
    }
}
// improvement IV: make an array for the moduli of the smallest primes
UBNUM sqfMarmet::prepareSmallPrimes(UBNUM startVal) {
    BNUM p, p2;
    p = 0;
    for (int j = 0; j < m_k2; j++) {
        p = NextPrime(p+1);
        p2 = power(p, m_exp);
        p2_ary[j] = p2;
        nqsf_ary[j] = 0;
        spRemainders[j] = startVal % (p2);
    }
    return p;
}
UBNUM sqfMarmet::prepareSmallPrimesDouble(UBNUM startVal) {
	for (int j = 0; j < m_k2; j++) {
        nqsf_ary[j] = 0;
        spRemainders[j] = startVal % (p2_ary[j]);		
	}
}
// fill NP2min for improvement II (only single squarefree !)
// we assume NP2minLen < 25
// TODO adjust this to powers other than 2
void sqfMarmet::fillNP2Min() {
	UBNUM val = 1;
	NP2min_ary[0] = 0;
	for (int j = 1; j < m_NP2minLen; j++) {
		NP2min_ary[j] = val;
        // cout << j << " "  << NP2min_ary[j] << ", ";
		val += ((j % 4)*(j % 9)) ? 1 : 0; // is val squarefree
	}
}
// supposed to work for 2 <= gaplen <= 25 and 2 <= e <= 39
void sqfMarmet::fillNPexpMin(int e){
	UBNUM modpow2 = power_long(2,e);
	UBNUM modpow3 = power_long(3,e);
	UBNUM val = 1;
	NP2min_ary[0] = 0;
	for (int j = 1; j < m_NP2minLen; j++) {
		NP2min_ary[j] = val;
        // cout << j << " "  << NP2min_ary[j] << ", ";
		val += ((j % modpow2)*(j % modpow3)) ? 1 : 0; // is val squarefree
	}
}

void sqfMarmet::fillNPexpMinDouble(){	
	NP2min_ary[0] = 0;
	for (int j = 1; j < m_NP2minLen; j++) {
		NP2min_ary[j] = j;
	}
}

// get the limit for the numbers p*q that have to be computed for the array (double squareful)
UBNUM sqfMarmet::getLimitDoubleSquare(BNUM n) const {
	if (n >= 1000000 - 100 and n <= 1000000 + 100){
		return n * 5;
	}
   	if (n >= 10000000 - 100 and n <= 10000000 + 100){		
		return llround(n * 5.6) - 1;
	}
    return n; // simple version limit -> to improve
}
// startVal not used for now
void sqfMarmet::fillNumArraysDouble(int gaplen, BNUM startVal = 1) {
	// startVal = 1;
	m_gaplen = gaplen;
    BNUM maxpq = getLimitDoubleSquare(m_cntAryLen);
    BNUM maxp = SqrRoot(maxpq);
    BNUM p2;
    size_t j;
    j = m_cntAryFilledLen = 0;
    for (BNUM p=2; p <= maxp and j < m_cntAryLen; p = NextPrime(p+1)) {
        p2 = p*p;
        for (BNUM q = NextPrime(p+1); p*q < maxpq and j < m_cntAryLen; q = NextPrime(q+1)) {
            p2_ary[j] = p2*q*q; // p^2*q^2
            j++;
        }
    }
    if (j >= m_cntAryLen){
		cerr <<  "Array overflow! Please adjust your limit!" << endl;
		throw "Array overflow! Please adjust your limit!";
	}
    m_cntAryFilledLen = j;
	// sortAry();
	initNSQFAry(startVal);
	sortAlgo();
        
	fillNPexpMinDouble();
    fillpow2();
    m_k2 = NP2min_ary[gaplen]-2;
    int p = prepareSmallPrimesDouble(startVal);	
	m_NLimit = p2_ary[m_cntAryFilledLen-1];
	cout << m_cntAryFilledLen << ", limit: " << m_NLimit << endl;	
}
// TODO sort both arrays according to nqsf_ary
void sqfMarmet::sortAry(){
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY)
	
#else
	sort(p2_ary.begin(), p2_ary.begin() + m_cntAryFilledLen);
#endif	
}
// sort p2_ary and nqsf_ary simultanounsly using bubble sort
// TODO bubble sort is too slow!
void sqfMarmet::sortAlgo(){
	m_oSqfSorter.sortSynced(&p2_ary, &nqsf_ary, m_cntAryFilledLen);
	return;
}
void sqfMarmet::initNSQFAry(BNUM startVal) {
	for (size_t i = 0;i < m_cntAryFilledLen;i++){
		UBNUM p2 = p2_ary[i];
		nqsf_ary[i] = (startVal / p2 + 1) * p2;
	}
}

// swap the contents of the array elements j and j+1 in both arrays in parallel
void sqfMarmet::swapNeighbours(int j) {
#ifdef SQFMARMET_POINTERS
#if SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY
	UBNUM* ptrP2Ary = &(p2_ary[j]);
	UBNUM* ptrNqsfAry = &(nqsf_ary[j]);
#else // SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY
	UBNUM* ptrP2Ary = &(p2_ary.data()[j]);
	UBNUM* ptrNqsfAry = &(nqsf_ary.data()[j]);
#endif	// SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY
	UBNUM nTmp = *ptrP2Ary;
	*ptrP2Ary = *(ptrP2Ary+1);
	*(ptrP2Ary+1) = nTmp;
	nTmp = *ptrNqsfAry;
	*ptrNqsfAry = *(ptrNqsfAry+1);
	*(ptrNqsfAry+1) = nTmp;
#else	// SQFMARMET_POINTERS	
	UBNUM nTmp = p2_ary[j];
	p2_ary[j] = p2_ary[j+1];
	p2_ary[j+1] = nTmp;
	nTmp = nqsf_ary[j];
	nqsf_ary[j] = nqsf_ary[j+1];
	nqsf_ary[j+1] = nTmp;
#endif	// SQFMARMET_POINTERS
}

// find next squareful number after N
// this is, where the arrays are sorted
UBNUM sqfMarmet::getNextSqFul(UBNUM N, int offset)
{    
	int lastj = m_cntAryFilledLen-1;
    while (nqsf_ary[offset] <= N) {
        nqsf_ary[offset] += p2_ary[offset];
        // move first element to the right, until array is sorted again
        int j = offset;
        while (j < lastj and nqsf_ary[j] > nqsf_ary[j+1]) {
            swapNeighbours(j);
            j++;
        }
    }
    return nqsf_ary[offset];
}
// replacement for "getNextSqFul" using chained list
UBNUM sqfMarmet::sortVII(UBNUM N, int offset){
	// offset = m_head ?;
	int i, lasti, nexthead;
	int lastj = m_cntAryLen - 1;
	while (nqsf_ary[m_head] <= N) {
		nqsf_ary[m_head] += p2_ary[m_head];
		UBNUM tempnsqf = nqsf_ary[m_head];		
		nexthead = next_ary[m_head];		
		// assert (tempnsqf >= nqsf_ary[nexthead]);
		lasti = m_head;
		for (i=lasti=m_head; tempnsqf>=nqsf_ary[next_ary[i]];){
			 lasti=i;
			 i=next_ary[i];
		}
		next_ary[m_head] = next_ary[i];
		// if (lasti != m_head){
		if (i != m_head){
			next_ary[i] = m_head;
			m_head = nexthead;
		}		
	}
	return nqsf_ary[m_head];
}
// improvement II of single squareful
// this version skips some squareful numbers depending on L
BNUM sqfMarmet::getNextSqFulGap(BNUM N, BNUM LMin, bool bJump) {
    // cout << N;
    bJump = bJump and LMin >= 2;
    if (bJump and nqsf_ary[NP2min_ary[LMin-1]] > nqsf_ary[0]+LMin) {
        N = nqsf_ary[NP2min_ary[LMin-1]] - LMin + 1;
        // cout << " -> " << N << endl;
        N--;
    }
    // cout << endl;
    return getNextSqFul(N);
}
// gap length is now a member variable
UBNUM sqfMarmet::getNextSqFulGap_improV(UBNUM N, UBNUM endVal) {
	initNRemainders(N);
	UBNUM nLoopLimit = (endVal <= N) ? m_NLimit : min(m_NLimit, endVal);
    for (;m_N < nLoopLimit;) {
        // improvement V
		N = m_N;
        while (nqsf_ary[m_k2 + 1] - N > m_gaplen - 1) {
            N = getNextSqFul(N, m_k2); // maybe we must use nqsf_ary[m_k2 + 1] - 1 as argument ?
        }
        setNRemainders(N);        

        // check if the numbers around N are squareful
        // 1st check the values below k1
        int posN = m_gaplen; // has to be bigger , but how small can we make it ?
        UBNUM Nlow = m_N - posN;
        UBNUM Nhigh = m_N + posN;
        UBNUM NIt;
		BNUM rem, p2;
		BITVECTOR bitv = 0;
        for (int j = 0; j < m_k2; j++) {
            p2 = p2_ary[j];
            rem = (spRemainders[j] - posN) % p2;
			rem += p2 * (rem < 0); // ensure it's >= 0
            for (NIt = Nlow + p2 - rem; NIt <= Nhigh; NIt += p2) {
                bitv |= pow2_ary[(mc_bitvecZeroPos + NIt) - N]; // do this with pointer!!
            }
        }
        // 2. to get the start of the sequence of squareful numbers, we check the numbers N-1, N-2... until we find a squarefree number
        int bitpos;
        for (bitpos = mc_bitvecZeroPos - 1; bitpos >= 0; bitpos--) {
            if ((bitv & pow2_ary[bitpos]) == 0) { // TODO  put this into for (...)
                break;
            }
        }
        UBNUM N0 = m_N - mc_bitvecZeroPos + bitpos + 1;
        // N0 is the first number of the gap
        // 3. check if the subsequent numbers are squareful while sorting the array
		

        UBNUM NNew = m_N, N1 = NNew - 1;
		bitpos = mc_bitvecZeroPos + N1 - m_N;
		for (bool bSquareful = true; bSquareful;){
			N1++;
			bitpos++;
			if (N1 == NNew){
				NNew = getNextSqFul(N1, m_k2);
				continue;
			}			
			bSquareful = (bitv & pow2_ary[bitpos]);
		}
		N1--;
		
        if (N1 - N0 >= m_gaplen - 1) {            
			return N0;
        }
		setNRemainders(NNew);
    }
    // we have found no gap in the given interval    
    return 0;
}
void sqfMarmet::initNRemainders(UBNUM N){
	m_N = N;
    for(int j = 0; j < m_k2; j++) {
        spRemainders[j] = m_N % p2_ary[j];
    }
}
void sqfMarmet::add2Remainders(BNUM diff) {
    for(int j = 0; j < m_k2; j++) {
        spRemainders[j] += diff;
		spRemainders[j] %= p2_ary[j];
    }
}
void sqfMarmet::setNRemainders(UBNUM N){
	BNUM diff = N - m_N;
	m_N = N;
    for(int j = 0; j < m_k2; j++) {
        spRemainders[j] += diff;
		spRemainders[j] %= p2_ary[j];
    }	
}
// + chained list
UBNUM sqfMarmet::getNextSqFulGap_improVII(UBNUM N){
	initNRemainders(N);    
    for (;m_N < m_NLimit;) {
		N = m_N;
        // while (nqsf_ary[m_k2 + 1] - N > m_gaplen - 1) { use next_ary
		while (nqsf_ary[next_ary[m_head]] - N > m_gaplen - 1){
            N = sortVII(N, m_k2); // maybe we must use nqsf_ary[m_k2 + 1] - 1 as argument ?
        }
        setNRemainders(N);        

        // check if the numbers around N are squareful
        // 1st check the values below k1
        int posN = m_gaplen; // has to be bigger , but how small can we make it ?
        UBNUM Nlow = m_N - posN;
        UBNUM Nhigh = m_N + posN;
        UBNUM NIt;
		BNUM rem, p2;
		BITVECTOR bitv = 0;
        for (int j = 0; j < m_k2; j++) {
            p2 = p2_ary[j];
            rem = (spRemainders[j] - posN) % p2;
			rem += p2 * (rem < 0); // ensure it's >= 0
            for (NIt = Nlow + p2 - rem; NIt <= Nhigh; NIt += p2) {
                bitv |= pow2_ary[(mc_bitvecZeroPos + NIt) - N]; // do this with pointer!!
            }
        }
        // 2. to get the start of the sequence of squareful numbers, we check the numbers N-1, N-2... until we find a squarefree number
        int bitpos;
        for (bitpos = mc_bitvecZeroPos - 1; bitpos >= 0; bitpos--) {
            if ((bitv & pow2_ary[bitpos]) == 0) { // TODO  put this into for (...)
                break;
            }
        }
        UBNUM N0 = m_N - mc_bitvecZeroPos + bitpos + 1;
        // N0 is the first number of the gap
        // 3. check if the subsequent numbers are squareful while sorting the array
		

        UBNUM NNew = m_N, N1 = NNew - 1;
		bitpos = mc_bitvecZeroPos + N1 - m_N;
		for (bool bSquareful = true; bSquareful;){
			N1++;
			bitpos++;
			if (N1 == NNew){
				NNew = sortVII(N1, m_k2);
				continue;
			}			
			bSquareful = (bitv & pow2_ary[bitpos]);
		}
		N1--;
		
        if (N1 - N0 >= m_gaplen - 1) {            
			return N0;
        }
		setNRemainders(NNew);
    }
    return 0;
}

BNUM sqfMarmet::getFirstGap_ver01(int l) {
    BNUM N,N2, exp, limit;

    exp = 10;
    limit = power_long(10,exp);

    N = N2 = 2;
    int curgap = 1;
    for (BNUM j=0; j < limit; j++) {
        N2 = getNextSqFul(N);
        if (N2 - N == 1) {
            curgap++;
            if(curgap == l) {
                return N2-l+1;
            }
        } else {
            curgap = 1;
        }
        N = N2;
    }
    cout << "last N = " << N << endl;
    return -1;
}
BNUM sqfMarmet::getFirstGap(int l, BNUM startVal, UBNUM endVal) {
    BNUM N,N2, exp, limit;
    BNUM NLimit;

    exp = 10;
    limit = power_long(10,exp);
    NLimit = p2_ary[m_cntAryLen-1];
    cout << NLimit << endl;

    N = N2 = startVal;
    int curgap = 1;
    for (BNUM j=0; j < limit and N < NLimit; j++) {
        N2 = getNextSqFulGap(N, l - curgap, true);
        if (N2 - N == 1) {
            curgap++;
            if(curgap == l) {
                return N2-l+1;
            }
        } else {
            curgap = 1;
        }
        N = N2;
    }
    return -1;
}
UBNUM sqfMarmet::getFirstGap_improV(int l, UBNUM startVal, UBNUM  endVal) {
	return getNextSqFulGap_improV(startVal, endVal);
}
UBNUM sqfMarmet::getFirstGap_improVII(int l, UBNUM startVal, UBNUM  endVal){
	return getNextSqFulGap_improVII(startVal);
}
bool sqfMarmet::checkSequence (int gaplen, UBNUM N){
    auto checkFunc = func_ValueChecker();
	UBNUM n = N;	
	for(bool bOk = true ; bOk;n++){
		cout << n << " = " << m_oFactor.prnFactor(n) << endl;		
        bOk = checkFunc(n);
	}
	int gl = n - N - 1;
	if (gl < gaplen){
		cout << "Failure: gap length only " << gl << endl;
		return false;
	}
	if (gl > gaplen){
		cout << "Found bigger gap: " << gl << endl;
	}
	return true;
}
bool sqfMarmet::checkValPower (UBNUM n){
    return m_oFactor.hasPower(n, m_exp);
}
bool sqfMarmet::checkValDouble (UBNUM n){
    return m_oFactor.IsDoubleSquareFul(n);
}

void sqfMarmet::batchCalc(UBNUM startVal, UBNUM endVal, int mins, enBatchMode batchMode, int startNr){
	UBNUM gapNum;
	auto gpFunc = func_getFirstGap();
	double dStartTime = GetTime();
	prnTimeStamp("Starting at");
	endVal = (batchMode == BMTimeLimit) ? startVal : endVal;
	do {
		gapNum = gpFunc(m_gaplen, startVal, endVal);
		if (gapNum != 0){
			startNr++;
			cout << startNr << ". " << gapNum  << ": ~< 10^" << to_string(gapNum).length() << endl;
			if (not checkSequence(m_gaplen, gapNum)){
				cout << "Error in sequence !";
				break;
			}
		}
	} while (gapNum != 0 and (not (batchMode == BMTimeLimit and (GetTime() - dStartTime) > 60 * mins)));
	prnTimeStamp("Finished at");
	
}
// get all the sequences in the given interval
// arrays are supposed to be filled beforehand and parameters defined
void sqfMarmet::batchCalc(UBNUM startVal, UBNUM endVal, int startNr){
	cout << endl;
	cout << "Batch Mode: " << "starting value = " << startVal << ", end value = " << endVal << endl;
	batchCalc(startVal, endVal, 0, BMNumberLimit, startNr);
	return;
}
// get all the sequences above startVal, stop after "mins" minutes
// arrays are supposed to be filled beforehand and parameters defined
void sqfMarmet::batchCalcTime(UBNUM startVal, int mins, int startNr){
	cout << endl;
	cout << "Batch Mode: " << "starting value = " << startVal << ", running at least for " << mins << " minutes" << endl;	
	batchCalc(startVal, startVal, mins, BMTimeLimit, startNr);
	return;
}
void sqfMarmet::prnTimeStamp(string sText) const{
    time_t result = time(nullptr);
    cout << sText << " " << asctime(localtime(&result)) << endl;	
}
void sqfMarmet::exec() {
	double dStartTime;

#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDVECTOR)
    cout << "Using vectors!" << endl;
#endif // #if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDVECTOR)
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDARRAY)	
    cout << "Using arrays!" << endl;
#endif // #if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USESTDARRAY)
#if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY)
	cout << "Using plain vanilla arrays!" << endl;	
#endif // #if (SQFMARMET_CONTAINERTYPE == SQFMARMET_USEVANILLAARRAY)
	cout << "Filling arrays!" << endl;
// ***********************************************************************************************//
// ************************************INPUT PARAMETERS FOR CALCULATION **************************//	
// ***********************************************************************************************//
	m_cntPrimeFactors = 1; // 1: for numbers with divisors p^e, 2: for numbers with divisors p^e*q^e
	m_exp = 3; // exponent of the prime power p^e 
	int gaplength = 7; // length of the sequence having the desired prime factors
	UBNUM startVal =    2287355532961372; // Find sequences starting above startVal
	UBNUM endVal   =    1; // Terminating if no sequence is found below endVal
	int timeInMins = 120; // In batch mode: Stop if a sequence is found and the program has been running for more than timeInMins minutes.
	int startNr = 544; // Counter used in output.
	bool doBatch = true; // batch mode ?
// ***********************************************************************************************//
// ************************************INPUT PARAMETERS FOR CALCULATION **************************//	
// ***********************************************************************************************//

    cout << "gaplength: " << gaplength << ", exponent: " << m_exp << ", number of prime factors: " << m_cntPrimeFactors << endl;
	
	auto fillFunc = func_fillNumArrays();

	dStartTime = GetTime();
	fillFunc(gaplength, startVal);
	cout << "Time to fill array/vector: " << GetTime() - dStartTime << " secs" << endl;
	
	if (doBatch){
		batchCalcTime(startVal, timeInMins, startNr);
		cout << endl;
		return;
	}
	
    dStartTime = GetTime();
    cout << "algorithm/improvement: " << m_improvement << endl;
	cout << "gap length: " << gaplength << endl;
	cout << "exponent: " << m_exp << endl;
	
	UBNUM gapNum;	
	auto gpFunc = func_getFirstGap();
	gapNum = gpFunc(gaplength, startVal, endVal);
	
	int gapNumLen = to_string(gapNum).length();
	if (gapNum > 0) {
		cout << gapNum << ": ~< 10^" << gapNumLen << endl;
	} else {
		cout << "No gap of length " << gaplength << " found in the interval [" << startVal << "," << endVal << "]." << endl;
 	}
    cout << GetTime() - dStartTime << " secs" << endl;
}
