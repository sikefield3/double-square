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

#include "checkresult.h"

CheckResult::CheckResult(string parFName, enCheckType cht, int gaplen){
	m_fname = parFName;
	m_cht = cht;
	m_gaplen = gaplen;
}
void CheckResult::run(){
	m_ifs.open (m_fname, ifstream::in);
	if (m_ifs.fail()){
		cout << m_fname << " not available!" << endl;
		exit (0);
	}
	if (m_cht == ChDoubleSquareFree){
		cout << "CheckResult : Start!" << endl;
		checkDoubleSquareFree();
	}
}
string CheckResult::getNextLine(){
	string s;
	while (getline(m_ifs, s)){
		if (s.empty()){
			continue;
		}
		else {
			return s;
		}
	}
	return "";
}
string CheckResult::getNextNumberLine(){
	string s, sRead;	
	bool fndalpha = true;
	while(not (sRead = getNextLine()).empty() and fndalpha){
		s = sRead;
		fndalpha = false;
		for (char& ch : s){			
			fndalpha = isalpha(ch);
			if (fndalpha){
				break;
			}
		}
		if (not fndalpha) {
			return s;
		}
	}
	return "";
}
void CheckResult::checkDoubleSquareFree(){
	int nr = 1, curgaplen = m_gaplen, nr_found;
	size_t pos_fst;
	string s;
	UBNUM val, oldval;
	Factor oFactor;
	
	while (not (s = getNextNumberLine()).empty()){
		// check for "."		
		pos_fst = s.find_first_of(".");
		if (pos_fst != string::npos and curgaplen > 0){
			if (curgaplen > 0 and curgaplen < m_gaplen){
				cout << "Sequence with too small length: " << endl << "nr = " << nr << ", val = " << val << endl;
				exit(0);
			}
			if (nr != (nr_found = stoi(s))){
				cout << "Incorrect number: expected: " << nr << ", got: " << nr_found << endl;
				nr = nr_found;				
			}
			s = s.substr(pos_fst+1);
			val = stoll(s);
			if (oFactor.IsDoubleSquareFul(val-1)){
				cout << "Sequence starting earlier !" << endl;
				exit(0);
			}
			nr++;
			curgaplen = 0;			
		} else if (curgaplen > 1) {
			oldval = val;
			val = stoll (s);
			if (val != oldval + 1){
				cout << "No sequence !" << endl;
				exit (0);
			}
		}
		if (curgaplen <= m_gaplen and not oFactor.IsDoubleSquareFul(val)){
				cout << val << " has not the required property !" << endl;
				exit(0);			
		}
		curgaplen++;
	}
	cout << "Done: " << nr-1 << endl;
}

