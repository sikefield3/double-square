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

#ifndef FACTOR_H
#define FACTOR_H

#include <vector>
#include <string>
#include <sstream>
#include <numeric>
#include <functional>

#include <NTL/ZZ.h>

#include "global.h"

using namespace std;
using namespace NTL;

class Factor
{
	UBNUM m_n;
	static const size_t m_vectorSize = 100;
	vector<UBNUM> prime_ary, exp_ary;
private:
	void init(bool bFirstTime = false);
	void add(UBNUM p, int e);
	string str() const;
public:
	Factor();
	void doFactor(UBNUM n);
	void doFactorNP(UBNUM n);
	string prnFactor(UBNUM n);
	string prnFactorNP(UBNUM n);
	int LittleOmega(UBNUM n);
	int LittleOmega() const;
	UBNUM radical(UBNUM n);
	UBNUM radical() const;
	bool IsDoubleSquareFul(UBNUM n);
	bool IsDoubleSquareFul() const;
	
};

#endif // FACTOR_H
