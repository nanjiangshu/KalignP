/*
	kalignP_hirschberg.h
	
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

/* ==== consensus - consensus alignemnt */
int* hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
int* hirsch_pp_dyn_new(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
int* hirsch_pp_dyn_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path);
int* hirsch_pp_dyn_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path);
int* hirsch_pp_dyn_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path);
int* hirsch_pp_dyn_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path);

struct states* foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* foward_hirsch_pp_dyn_new(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* foward_hirsch_pp_dyn_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm);
struct states* foward_hirsch_pp_dyn_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm);
struct states* foward_hirsch_pp_dyn_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm);
struct states* foward_hirsch_pp_dyn_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm);

struct states* backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* backward_hirsch_pp_dyn_new(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* backward_hirsch_pp_dyn_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm);
struct states* backward_hirsch_pp_dyn_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm);
struct states* backward_hirsch_pp_dyn_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm);
struct states* backward_hirsch_pp_dyn_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm);

int* hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int* hirsch_align_two_pp_vector_new(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int* hirsch_align_two_pp_vector_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int* hirsch_align_two_pp_vector_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int* hirsch_align_two_pp_vector_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int* hirsch_align_two_pp_vector_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);

struct states* backward_hirsch_pp_dynm(const float* prof1,const float* prof2,struct hirsch_mem* hm);

struct states* backward_hirsch_pp_dyn0(const float* prof1,const float* prof2,struct hirsch_mem* hm);

/* ==== consensus - sequence alignemnt */
int* hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip);
int* hirsch_ps_dyn_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm, int* hirsch_path,int sip);
int* hirsch_ps_dyn_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm, int* hirsch_path,int sip);
int* hirsch_ps_dyn_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, int sip);
int* hirsch_ps_dyn_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, float* weightArray1, int sip);
int* hirsch_ps_dyn_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, float* weightArray1, int sip);
int* hirsch_ps_dyn_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, float* weightArray1, int sip);

struct states* foward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip);
struct states* foward_hirsch_ps_dyn_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int sip);
struct states* foward_hirsch_ps_dyn_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int sip);
struct states* foward_hirsch_ps_dyn_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, int sip);
struct states* foward_hirsch_ps_dyn_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip);
struct states* foward_hirsch_ps_dyn_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip);
struct states* foward_hirsch_ps_dyn_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip);

struct states* backward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip);
struct states* backward_hirsch_ps_dyn_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int sip);
struct states* backward_hirsch_ps_dyn_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int sip);
struct states* backward_hirsch_ps_dyn_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, int sip);
struct states* backward_hirsch_ps_dyn_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip);
struct states* backward_hirsch_ps_dyn_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip);
struct states* backward_hirsch_ps_dyn_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip);

int* hirsch_align_two_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);
int* hirsch_align_two_ps_vector_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);
int* hirsch_align_two_ps_vector_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);
int* hirsch_align_two_ps_vector_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, int sip);
int* hirsch_align_two_ps_vector_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, float* weightArray1, int sip);
int* hirsch_align_two_ps_vector_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, float* weightArray1, int sip);
int* hirsch_align_two_ps_vector_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, float* weightArray1, int sip);


/* ==== sequence - sequence alignemnt */
int* hirsch_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path);
int* hirsch_ss_dyn_new(float**subm, const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path);

struct states* foward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm);
struct states* foward_hirsch_ss_dyn_new(float**subm,const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm);

struct states* backward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm);
struct states* backward_hirsch_ss_dyn_new(float**subm,const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm);

int* hirsch_align_two_ss_vector(float**subm,const int* seq1,const int* seq2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
int* hirsch_align_two_ss_vector_new(float**subm,const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);

