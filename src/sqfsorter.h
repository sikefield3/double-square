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

#ifndef SQFSORTER_H
#define SQFSORTER_H

#include <functional>
#include <iostream>
#include <vector>

#include "global.h"

using namespace std;

class sqfSorter
{
public:
    sqfSorter();
	void sortSynced(vector<UBNUM>* p2_ary, vector<UBNUM>* nqsf_ary, size_t cnt, bool lessThan = true);
private:
	function<bool(size_t, size_t)> fc_cmp;
	size_t m_cnt;
	vector<UBNUM>* m_ptr_p2_ary, *m_ptr_nqsf_ary;
	
	void heapsort();
	void heapify();
	void siftUp(size_t ixStart, size_t ixEnd);
	void siftDown(size_t ixStart, size_t ixEnd);
	int iParent(int n);
	int iLeftChild(int n);
	void swap_ary(size_t i, size_t j);
	
	void bubbleSort();
};

#endif // SQFSORTER_H
