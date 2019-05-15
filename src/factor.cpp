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

#include "factor.h"
Factor::Factor(){
	init();
}

void Factor::init(bool bFirstTime){
	if (bFirstTime){
		prime_ary.reserve(m_vectorSize);
		exp_ary.reserve(m_vectorSize);
	}
	prime_ary.clear();
	exp_ary.clear();
}
void Factor::add(UBNUM p, int e){
	prime_ary.push_back(p);
	exp_ary.push_back(e);
	
}
void Factor::doFactor(UBNUM n){
	init();
	ZZ zz_n = ZZ(n);
	m_n = n;
	int e = NumTwos(zz_n);
	if (e>0){
		add(2,e);
		MakeOdd(zz_n);
		n = to_long(zz_n);
	}
	UBNUM m,d,dLim;
	for (m=n, d=3, dLim = SqrRoot(m);d<=dLim; d += 2){
		if (m % d == 0){
			for (e=0; m % d == 0; m /= d, e++);
			add(d,e);
			dLim = SqrRoot(m);
		}
	}
	if (m>1){
		add(m,1);
	}
}
// this function is much too slow (~factor 50)
void Factor::doFactorNP(UBNUM n){
	init();
	UBNUM p;
	UBNUM m,dLim,e;
	for (p=2, m=n, dLim=SqrRoot(m); p<=dLim;p = NextPrime(p+1)){
		if (m % p == 0){
			for (e=0; m % p == 0; m /= p, e++);
			add(p,e);
			dLim = SqrRoot(m);
		}
	}
	if (m>1){
		add(m,1);
	}
}
string Factor::str() const{
	ostringstream os;
	for (size_t j = 0;j < prime_ary.size(); j++){
		if (j > 0){
			os << " * ";
		}
		os << prime_ary[j];
		int e = exp_ary[j];
		if (e>1){
			os << "^" << exp_ary[j];
		}
	}
	return os.str();	
}
string Factor::prnFactor(UBNUM n){
	doFactor(n);
	return str();
}
string Factor::prnFactorNP(UBNUM n){
	doFactorNP(n);
	return str();
}
int Factor::LittleOmega(UBNUM n){
	doFactor(n);
	return prime_ary.size();
}
int Factor::LittleOmega() const{
	return prime_ary.size();
}
UBNUM Factor::radical(UBNUM n){
	doFactor(n);
	return radical();
}
UBNUM Factor::radical() const{	
	// return accumulate (prime_ary.begin(), prime_ary.end(), 1, multiplies<UBNUM>());
	UBNUM res = 1;
	for (auto pr : prime_ary){
		res *= pr;
	}
	return res;
}
bool Factor::IsDoubleSquareFul(UBNUM n){
	doFactor(n);
	return IsDoubleSquareFul();
}
bool Factor::IsDoubleSquareFul() const{
	Factor oFc;
	return oFc.LittleOmega(m_n / radical()) >= 2;
}
