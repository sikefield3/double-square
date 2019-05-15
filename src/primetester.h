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

#ifndef PRIMETESTER_H
#define PRIMETESTER_H

#include <functional>

#include <NTL/ZZ.h>

#include "global.h"

using namespace NTL;
using namespace std;

// test NTL functions ProbPrime, NextPrime against trial division
class PrimeTester
{
private:
	bool isPrimeTrial(BNUM n) const;
	BNUM NextPrime_PT(BNUM p) const;
public:
	void testProbPrime(BNUM N) const;
	void testNextPrime(BNUM N) const;
	BNUM testCheckNextPrimeForComposite(BNUM NLimit, int NumTrials = -1) const;
};

#endif // PRIMETESTER_H
