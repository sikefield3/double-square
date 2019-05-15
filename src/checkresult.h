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

#ifndef CHECKRESULT_H
#define CHECKRESULT_H

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <locale>

#include "global.h"
#include "factor.h"

using namespace std;
// using namespace NTL;

enum enCheckType {ChDoubleSquareFree};

class CheckResult
{	
	string m_fname;
	ifstream m_ifs;
	int m_gaplen;
	enCheckType m_cht;
private:
	string getNextLine();
	string getNextNumberLine();
	bool isDoubleSquareFree(UBNUM n);
public:	
	CheckResult(string parFName, enCheckType cht, int gaplen);
	void run();
	void checkDoubleSquareFree();
};

#endif // CHECKRESULT_H
