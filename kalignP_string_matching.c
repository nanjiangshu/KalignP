/*
	kalignP_string_matching.c
	
    Copyright (C) 2010 Nanjiang Shu <nanjiang@sbc.su.se>
    Department of Biochemistry and Biophysics 
    Stockholm University

    Description: multiple sequence alignment supporting external supplied
    position specific gap penalties

    This software is derived from 
    Kalign version 2.03, Copyright (C) 2006 Timo Lassmann
    http://msa.cgb.ki.se/
    timolassmann@gmail.com

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

	A copy of this license is in the COPYING file.
*/

#include <string.h>

int byg_detect(int* text,int n)/*{{{*/
{/*return 1 (true) if text is DNA 
   return 0 (false) if text is Amino acid seq
   */
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	int mb = 1;
	//char *unique_aa = "EFILPQXZ";//permissiv
	//ABCDEFGHIJKLMNOPQRSTUVWXYZ
	char *unique_aa = "BDEFHIJKLMNOPQRSVWYZ";//restrictive, from the 26 char alphabet, A, C, G, T, U and X which are probable nucleic acids are removed, annotated by Nanjiang, 2010-09-28 
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	for (i= 0;i < 20;i++){
        /*bug found 2010-09-28, when unique_aa[i] == 'J', aacode['J'] == -1,
         * the operation T[aacode['J']] |= 1 will result in invalid write.
         * Therefore, the condition when unique_aa[i] == 'J' should be treated
         * specially*/
        if (aacode[unique_aa[i]-65] >= 0) {   /*this condition added on 2010-09-28*/
            T[(int)aacode[unique_aa[i]-65]] |= 1;
        }
	}
	for (i = 0;i < n;i++){
	//	fprintf(stderr,"%d\n",text[i]);
		if(text[i] != -1){
			s <<= 1; /*shift left by 1 position, equivalent to x2 */
			s |= 1;
			Tc = T[text[i]];
			s &= Tc;
			if(s & mb){
				return 0;
			}
		}
	}
	return 1;
}/*}}}*/
int check_identity(char* n,char*m)/*{{{*/
{ /*check whether string n == m */
	int len_n;
	int len_m;
	int i;
	
	len_n = strlen(n);
	len_m = strlen(m);
	if(len_m != len_n){
		return -1;
	}
	for (i = 0; i < len_n;i++){
		if(n[i] != m[i]){
			return -1;
		}
	}
	return 1;
	
}/*}}}*/
int byg_count(char* pattern,char*text)/*{{{*/
{ /*count the number of pattern matched in text*/
	int Tc;
	int count = 0;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	
	int m = strlen(pattern);
	int n = strlen (text);
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}

	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			count++;
		}
	}
	return count;
}/*}}}*/
int byg_end(char* pattern,char*text)/*{{{*/
{ /*match the pattern in text and return the end index of the substring*/
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	
	int m = strlen(pattern);
	int n = strlen (text);
	int mb = (1 << (m-1));

	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}

	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		if(!text[i]){
			return -1;
		}
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			return i+1;
		}
	}
	return -1;
}/*}}}*/
int byg_start(char* pattern,char*text)/*{{{*/
{ /*similar to strstr(), search "pattern" in "text", return the starting index of the matched substring, return -1 if not matched*/
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	
	int m = strlen(pattern);
	int n = strlen(text);
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}

	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			return i-m+1;
		}
	}
	return -1;
}/*}}}*/

