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

#ifndef SQF_H
#define SQF_H

#include <functional>
#include <string>
#include <sstream>

#include <NTL/ZZ.h>
#include <NTL/vector.h>

#include "global.h"

using namespace std;
using namespace NTL;


#define BNUM long long int

class SQF
{
private:
	ZZ prm;
	BNUM prmLimit;
	Vec<ZZ> vZZ; // the resulting sequence is stored here
private:
	void makePrimorial(BNUM n);
	bool isSquareFree(const ZZ& n) const;
	bool isSquareFree(BNUM n) const;
    bool isDoubleSquareFree_ver1(const ZZ& n) const; // -> false <=> there are primes p,q with p != q and p^2q^2 | n
	bool isDoubleSquareFree(const ZZ& n) const; // -> false <=> there are primes p,q with p != q and p^2q^2 | n
	ZZ	 cbrt(const ZZ& n) const; // -> m, s.t. m^3 <= n < (m+1)^3
	BNUM cbrt(BNUM n) const;
	void testcbrt();
	void testsquarefree();
	void testdoublesquarefree();
	void loop(const ZZ& nMax, std::function<bool(const ZZ&)> funcFree);
    string str(const NTL::ZZ& n) const;
public:
	SQF();
	void exec();
};

#endif // SQF_H
