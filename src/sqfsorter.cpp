/*
 * Double Squareful - Numbertheoretic program to find sequences related to 
 * prime factors / squarefree numbers.
 * Copyright (C) 2018-19 Bernd Zemann - zb4ng@arcor.de
 * 
 * <one line to give the program's name and a brief idea of what it does.>
 * Copyright (C) 2019  berzem <email>
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

#include "sqfsorter.h"

sqfSorter::sqfSorter()
{

}
// the indices 0,...,cnt-1 define the part of the arrays to be sorted
// less_than=true -> resulting array is ascending
void sqfSorter::sortSynced(vector<UBNUM>* p2_ary, vector<UBNUM>* nqsf_ary, size_t cnt, bool less_than){
	if (less_than) {
		fc_cmp = [this] (size_t i,size_t j) {return (*m_ptr_nqsf_ary)[i] < (*m_ptr_nqsf_ary)[j];};
	} else {
		fc_cmp = [this] (size_t i,size_t j) {return (*m_ptr_nqsf_ary)[i] > (*m_ptr_nqsf_ary)[j];};
	}
	m_cnt = cnt;
	// TODO I assume that the following assignments are shallow, but this has to be checked!
	m_ptr_p2_ary = p2_ary;
	m_ptr_nqsf_ary = nqsf_ary;
	heapsort();
}

void sqfSorter::heapsort(){
	heapify();
	size_t ixEnd = m_cnt - 1;
	while (ixEnd > 0){
		swap_ary(ixEnd, 0);
		ixEnd--;
		siftDown(0, ixEnd); // TODO not clear if we shoud use siftUp or siftDown !
	}
}
void sqfSorter::heapify(){
	size_t ixEnd = 1;
	while (ixEnd < m_cnt){
		siftUp(0, ixEnd);
		ixEnd++;
	}
}
void sqfSorter::siftUp(size_t ixStart, size_t ixEnd){
	size_t ixChild = ixEnd;	
	while (ixChild > ixStart){
		size_t ixParent = iParent(ixChild);
		// if((*m_ptr_nqsf_ary)[ixParent] < (*m_ptr_nqsf_ary)[ixChild]){
		if(fc_cmp(ixParent, ixChild)){
			swap_ary(ixParent, ixChild);
			ixChild = ixParent;
		} else {
			return;
		}
	}
}
void sqfSorter::siftDown(size_t ixStart, size_t ixEnd){
	size_t ixRoot = ixStart;
	while (iLeftChild(ixRoot) <= ixEnd){
		size_t ixChild = iLeftChild(ixRoot);
		size_t ixSwap = ixRoot;
		if(fc_cmp(ixSwap, ixChild)){
			ixSwap = ixChild;
		}
		if (ixChild + 1 <= ixEnd and fc_cmp(ixSwap, ixChild+1)){
			ixSwap = ixChild + 1;
		}
		if (ixSwap == ixRoot){
			return;
		} else{
			swap_ary(ixRoot, ixSwap);
			ixRoot = ixSwap;
		}
	}
}
int sqfSorter::iParent(int n){
	return (n-1)/2;
}
int sqfSorter::iLeftChild(int n){
	return 2*n+1;
}

void sqfSorter::swap_ary(size_t i, size_t j){
	UBNUM nTmp = (*m_ptr_p2_ary)[i];
    (*m_ptr_p2_ary)[i] = (*m_ptr_p2_ary)[j];
    (*m_ptr_p2_ary)[j] = nTmp;
    nTmp = (*m_ptr_nqsf_ary)[i];
    (*m_ptr_nqsf_ary)[i] = (*m_ptr_nqsf_ary)[j];
    (*m_ptr_nqsf_ary)[j] = nTmp;

}
void sqfSorter::bubbleSort(){		
	bool bNotFinished = true;
	int j;
	while (bNotFinished){
		for (j = 0, bNotFinished = false;j < m_cnt-1; j++){
			if ((*m_ptr_nqsf_ary)[j] > (*m_ptr_nqsf_ary)[j+1]){
				bNotFinished = true;
				swap_ary(j,j+1);
			}
		}
	}
}
