/*
    kalignP_hirschberg.c
    
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

#include "kalignP.h"
#include "kalignP_hirschberg.h"
#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)
//#include <emmintrin.h>

int** hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength)/*{{{*/
    /*strength: param->gap_inc, default = 0.0*/
{
    struct hirsch_mem* hm = 0;
    int i,j,g,a,b,c;
    int len_a;
    int len_b;
    float** profile = 0;

    profile = malloc(sizeof(float*)*numprofiles);
    for ( i = 0;i< numprofiles;i++){
        profile[i] = 0;
    }

    map = malloc(sizeof(int*)*numprofiles);
    for ( i = 0;i < numprofiles;i++){
        map[i] = 0;
    }
    
    for (i = 0; i < numseq; i ++){  /*set the number of residues at each position of the sequence  2010-10-19 */
        aln->nRes[i] = malloc(sizeof(unsigned int) * (aln->sl[i]+2));
        for (j = 0; j < aln->sl[i]+2; j ++){
            aln->nRes[i][j] = 1;
        }
    }


    hm = hirsch_mem_alloc(hm,1024);/*initially allocated size is 1024, if the max(lena,lenb) is > 1024, it will be re-allocated*/

    fprintf(stderr,"\nAlignment:\n");
#ifdef DEBUG/*{{{*/
    fprintf(stderr,"Tree:\n");
    for (i = 0; i < numseq; i ++)
    {
        fprintf(stderr,"Tree:%d: %3d %3d %3d\n", i, tree[i*3], tree[i*3+1], tree[i*3+2]);
    }
#endif /*}}}*/
    float* weightArray_a = 0;
    float* weightArray_b = 0;

    for (i = 0; i < (numseq-1);i++){
        weightArray_a = 0;
        weightArray_b = 0;
        a = tree[i*3];
        b = tree[i*3+1];
        c = tree[i*3+2];
        fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)numseq * 100);
#ifdef DEBUG/*{{{*/
        fprintf(stderr,"\nAligning:%d %d->%d	done:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
        fprintf(stderr,"\ngpo=%f\ngpe=%f\ntgpe=%f\n",gpo,gpe,tgpe);
#endif /*}}}*/
        len_a = aln->sl[a];
        len_b = aln->sl[b];

        g = (len_a > len_b)? len_a:len_b; /*max(len_a, len_b)*/
        map[c] = malloc(sizeof(int) * (g+2));
        if(g > hm->size){
            hm = hirsch_mem_realloc(hm,g);
        }

        for (j = 0; j < (g+2);j++){
            map[c][j] = -1;
        }

        if (a < numseq){
            /*profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);*/
            profile[a] = make_profile_new(profile[a],aln->s[a],aln->gpo[a], aln->gpe[a], aln->tgpe[a], len_a,submatrix);
        }else{ /*set_gap_penalties for consensus sequences, e.g. for A, B, C, E four sequences, when compared (AB) to C, the gap penalties for (AB) should be reset*/
            /*set_gap_penalties(profile[a],len_a,aln->nsip[b],strength,aln->nsip[a]);[>this can be modified<]*/
            weightArray_a = malloc(sizeof(float)*(len_a+2));
            set_gap_penalties_new(profile[a],len_a,aln->nsip[b],strength,aln->nsip[a], weightArray_a);/*this can be modified*/
            //smooth_gaps(profile[a],len_a,window,strength);
            //increase_gaps(profile[a],len_a,window,strength);
        }
        if (b < numseq){
            /*profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);*/
            profile[b] = make_profile_new(profile[b],aln->s[b],aln->gpo[b], aln->gpe[b], aln->tgpe[b], len_b,submatrix);
        }else{		
            weightArray_b = malloc(sizeof(float)*(len_b+2));
            /*set_gap_penalties(profile[b],len_b,aln->nsip[a],strength,aln->nsip[b]);*/
            set_gap_penalties_new(profile[b],len_b,aln->nsip[a],strength,aln->nsip[b], weightArray_b);
            //smooth_gaps(profile[b],len_b,window,strength);
            //increase_gaps(profile[b],len_b,window,strength);
        }

        hm->starta = 0;
        hm->startb = 0;
        hm->enda = len_a;
        hm->endb = len_b;
        hm->len_a = len_a;
        hm->len_b = len_b;

        hm->f[0].a = 0.0;
        hm->f[0].ga =  -FLOATINFTY;
        hm->f[0].gb = -FLOATINFTY;
        hm->b[0].a = 0.0;
        hm->b[0].ga =  -FLOATINFTY;
        hm->b[0].gb =  -FLOATINFTY;
#ifdef DEBUG
        fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
#endif 
        if(a < numseq){
            if(b < numseq){
                /*map[c] = hirsch_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);[>submatrix: default=gon250mt, sequences-sequences<]*/
                map[c] = hirsch_ss_dyn_new(submatrix,aln->s[a],aln->s[b],aln->gpo[a]+1,aln->gpe[a]+1, aln->tgpe[a]+1, aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm, map[c]);/*submatrix: default=gon250mt, sequences-sequences*/
            }else{
                hm->enda = len_b;
                hm->endb = len_a;
                hm->len_a = len_b;
                hm->len_b = len_a;
                /*map[c] = hirsch_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);[>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new(profile[b],aln->s[a],aln->gpo[a]+1, aln->gpe[a]+1, aln->tgpe[a]+1, hm,map[c],aln->nsip[b]);[>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new1(profile[b],aln->s[a],aln->gpo[a]+1, aln->gpe[a]+1, aln->tgpe[a]+1, hm,map[c],aln->nRes[b], aln->nsip[b]);[>profile-sequence<]*/
                map[c] = hirsch_ps_dyn_new2(profile[b],aln->s[a],aln->gpo[a]+1, aln->gpe[a]+1, aln->tgpe[a]+1, hm,map[c],aln->nRes[b], weightArray_b, aln->nsip[b]);/*profile-sequence*/
                /*map[c] = hirsch_ps_dyn_new2_test1(profile[b],aln->s[a],aln->gpo[a]+1, aln->gpe[a]+1, aln->tgpe[a]+1, hm,map[c],aln->nRes[b], weightArray_b, aln->nsip[b]);[>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new0(profile[b],aln->s[a],aln->gpo[a]+1, aln->gpe[a]+1, aln->tgpe[a]+1, hm,map[c],aln->nsip[b]);[>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new2_test2(profile[b],aln->s[a],aln->gpo[a]+1, aln->gpe[a]+1, aln->tgpe[a]+1, hm,map[c],aln->nRes[b], weightArray_b, aln->nsip[b]);[>profile-sequence<]*/
                map[c] = mirror_hirsch_path(map[c],len_a,len_b); /*??*/
            }
        }else{
            if(b < numseq){
                /*map[c] = hirsch_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]); [>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new(profile[a],aln->s[b],aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm,map[c],aln->nsip[a]); [>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new1(profile[a],aln->s[b],aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm,map[c],aln->nRes[a], aln->nsip[a]); [>profile-sequence<]*/
                map[c] = hirsch_ps_dyn_new2(profile[a],aln->s[b],aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm,map[c],aln->nRes[a], weightArray_a, aln->nsip[a]); /*profile-sequence*/
                /*map[c] = hirsch_ps_dyn_new0(profile[a],aln->s[b],aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm,map[c],aln->nsip[a]); [>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new2_test1(profile[a],aln->s[b],aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm,map[c],aln->nRes[a], weightArray_a, aln->nsip[a]); [>profile-sequence<]*/
                /*map[c] = hirsch_ps_dyn_new2_test2(profile[a],aln->s[b],aln->gpo[b]+1, aln->gpe[b]+1, aln->tgpe[b]+1, hm,map[c],aln->nRes[a], weightArray_a, aln->nsip[a]); [>profile-sequence<]*/
            }else{
                if(len_a < len_b){
                    /*map[c] = hirsch_pp_dyn(profile[a],profile[b],hm,map[c]); [>profile-profile<]*/
                    /*map[c] = hirsch_pp_dyn_new(profile[a],profile[b],hm,map[c]); [>profile-profile<]*/
                    /*map[c] = hirsch_pp_dyn_new0(profile[a],profile[b],aln->nRes[a],aln->nRes[b], weightArray_a, weightArray_b, hm,map[c]); [>profile-profile<]*/
                    map[c] = hirsch_pp_dyn_new2(profile[a],profile[b],aln->nRes[a],aln->nRes[b], weightArray_a, weightArray_b, hm,map[c]); /*profile-profile*/
                    /*map[c] = hirsch_pp_dyn_new2_test1(profile[a],profile[b],aln->nRes[a],aln->nRes[b], weightArray_a, weightArray_b, hm,map[c]); [>profile-profile<]*/
                    /*map[c] = hirsch_pp_dyn_new2_test2(profile[a],profile[b],aln->nRes[a],aln->nRes[b], aln->nsip[a], aln->nsip[b], weightArray_a, weightArray_b, hm,map[c]); [>profile-profile<]*/
                }else{
                    hm->enda = len_b;
                    hm->endb = len_a;
                    hm->len_a = len_b;
                    hm->len_b = len_a;
                    /*map[c] = hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);*/
                    /*map[c] = hirsch_pp_dyn_new(profile[b],profile[a],hm,map[c]);*/
                    /*map[c] = hirsch_pp_dyn_new0(profile[b],profile[a],aln->nRes[b], aln->nRes[a],weightArray_b, weightArray_a, hm, map[c]);*/
                    map[c] = hirsch_pp_dyn_new2(profile[b],profile[a],aln->nRes[b], aln->nRes[a],weightArray_b, weightArray_a, hm, map[c]);
                    /*map[c] = hirsch_pp_dyn_new2_test1(profile[b],profile[a],aln->nRes[b], aln->nRes[a],weightArray_b, weightArray_a, hm, map[c]);*/
                    /*map[c] = hirsch_pp_dyn_new2_test2(profile[b],profile[a],aln->nRes[b], aln->nRes[a],aln->nsip[b], aln->nsip[a], weightArray_b, weightArray_a, hm, map[c]);*/
                    map[c] = mirror_hirsch_path(map[c],len_a,len_b);
                }
            }
        }

        map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);/*map[c] will be re-assigned, and the size of map[c] changed from "max(len_a,len_b) +2" to (len_a+len_b+2)*/

        if(i != numseq-2){
            profile[c] = malloc(sizeof(float)*64*(map[c][0]+2)); /*map[c][0] = length of consensus seq*/
            /*profile[c] = update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);*/
            /*profile[c] = update_new(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);*/
            profile[c] = update_new2(profile, map[c],aln, a, b, c );
#ifdef DEBUG_PROFILE1/*{{{*/
            int ii = 0;
            int cnt = 0;
            fprintf(stderr,"\nprofile[%d](a):gpo, size = %d\n", a, len_a);
            for (ii = 0; ii < map[c][0]+2; ii++ ) {
                if (cnt > 0 && cnt < len_a+1){
                    if (!map[c][ii] || (map[c][ii] & 2)){ /*aligned or gap_B*/
                        fprintf(stderr,"%6.2f ", profile[a][64*cnt+55]);
                        cnt ++;
                    }else{
                        fprintf(stderr,"%6s ", "-");
                    }
                }else{
                    fprintf(stderr,"%6.2f ", profile[a][64*cnt+55]);
                    cnt ++;
                }
            }
            fprintf(stderr,"\n");

            fprintf(stderr,"profile[%d](b):gpo, size = %d\n", b, len_b);
            cnt = 0;
            for (ii = 0; ii < map[c][0]+2; ii++ ) {
                if (cnt > 0 && cnt < len_b+1){
                    if (!map[c][ii] || (map[c][ii] & 1)){ /*aligned or gap_A*/
                        fprintf(stderr,"%6.2f ", profile[b][64*cnt+55]);
                        cnt ++;
                    }else{
                        fprintf(stderr,"%6s ", "-");
                    } 
                }else{
                    fprintf(stderr,"%6.2f ", profile[b][64*cnt+55]);
                    cnt ++;
                }
            }
            fprintf(stderr,"\n");

            fprintf(stderr,"profile[%d](c):gpo, size = %d\n", c, map[c][0]);
            for (ii = 0; ii < map[c][0]+2; ii++ ) {
                fprintf(stderr,"%6.2f ", profile[c][64*ii+55]);
            }
            fprintf(stderr,"\n");
#endif /*}}}*/
#ifdef DEBUG_PROFILE /*{{{*/
            fprintf(stderr,"\nnRes[%d](a):nRes, size = %d\n", a, len_a);
            int ii = 0;
            int cnt = 0;
            for (ii = 1; ii < map[c][0]+1; ii++ ) {
                if (!map[c][ii] || (map[c][ii] & 2)){ /*aligned or gap_B*/
                    fprintf(stderr,"%6d ", aln->nRes[a][cnt+1]);
                    cnt++;
                }else{
                    fprintf(stderr,"%6s ", "-");
                } 
            }
            fprintf(stderr,"\n");

            fprintf(stderr,"nRes[%d](b):nRes, size = %d\n", b, len_b);
            cnt = 0;
            for (ii = 1; ii < map[c][0]+1; ii++ ) {
                if (!map[c][ii] || (map[c][ii] & 1)){ /*aligned or gap_A*/
                    fprintf(stderr,"%6d ", aln->nRes[b][cnt+1]);
                    cnt ++;
                }else{
                    fprintf(stderr,"%6s ", "-");
                } 
            }
            fprintf(stderr,"\n");

            fprintf(stderr,"nRes[%d](c):nRes, size = %d\n", c, map[c][0]);
            for (ii = 1; ii < map[c][0]+1; ii++ ) {
                fprintf(stderr,"%6d ", aln->nRes[c][ii]);
            }
            fprintf(stderr,"\n");

#endif /*}}}*/
        }

        aln->sl[c] = map[c][0];

        aln->nsip[c] = aln->nsip[a] + aln->nsip[b]; /*nsip[c] updated here*/
        aln->sip[c] = malloc(sizeof(int)*(aln->nsip[a] + aln->nsip[b]));
        g =0;
        for (j = aln->nsip[a];j--;){
            aln->sip[c][g] = aln->sip[a][j];
            g++;
        }
        for (j = aln->nsip[b];j--;){
            aln->sip[c][g] = aln->sip[b][j];
            g++;
        }

        free(profile[a]);
        free(profile[b]);
        if (weightArray_a){
            free(weightArray_a);
        }
        if (weightArray_b){
            free(weightArray_b);
        }
    }
    fprintf(stderr,"\r%8.0f percent done\n",100.0);
    free(profile);
    hirsch_mem_free(hm);
    for (i = 32;i--;){
        free(submatrix[i]);
    }
    free(submatrix);
    return map;
}/*}}}*/

int** hirschberg_alignment_against_a(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength)/*{{{*/
{
    struct hirsch_mem* hm = 0;
    int i,j,g,a,b,c;
    int len_a;
    int len_b;
    float** profile = 0;

    profile = malloc(sizeof(float*)*numprofiles);
    for ( i = 0;i< numprofiles;i++){
        profile[i] = 0;
    }

    map = malloc(sizeof(int*)*numprofiles);
    for ( i = 0;i < numprofiles;i++){
        map[i] = 0;
    }
    
    hm = hirsch_mem_alloc(hm,1024);

    fprintf(stderr,"\nAlignment:\n");

    for (i = 0; i < (numseq-1);i++){
        a = tree[i*3];
        b = tree[i*3+1];
        c = tree[i*3+2];
        fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)numseq * 100);
#ifdef DEBUG
        fprintf(stderr,"\nAligning:%d %d->%d	done:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
#endif 
        len_a = aln->sl[a];
        len_b = aln->sl[b];

        
        g = (len_a > len_b)? len_a:len_b;
        map[c] = malloc(sizeof(int) * (g+2));
        if(g > hm->size){
            hm = hirsch_mem_realloc(hm,g);
        }

        for (j = 0; j < (g+2);j++){
            map[c][j] = -1;
        }

        if (a < numseq){
            profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
        }else{
            set_gap_penalties(profile[a],len_a,aln->nsip[b],0,aln->nsip[a]);
            //smooth_gaps(profile[a],len_a,window,strength);
            
            //increase_gaps(profile[a],len_a,window,strength);
        }
        if (b < numseq){
            profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
        }else{		
            set_gap_penalties(profile[b],len_b,aln->nsip[a],0,aln->nsip[b]);
            //smooth_gaps(profile[b],len_b,window,strength);
            //increase_gaps(profile[b],len_b,window,strength);
        }
        
        hm->starta = 0;
        hm->startb = 0;
        hm->enda = len_a;
        hm->endb = len_b;
        hm->len_a = len_a;
        hm->len_b = len_b;
        
        hm->f[0].a = 0.0;
        hm->f[0].ga =  -FLOATINFTY;
        hm->f[0].gb = -FLOATINFTY;
        hm->b[0].a = 0.0;
        hm->b[0].ga =  -FLOATINFTY;
        hm->b[0].gb =  -FLOATINFTY;
#ifdef DEBUG
        fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
#endif 
        if(a < numseq){
            if(b < numseq){
                map[c] = hirsch_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
            }else{
                hm->enda = len_b;
                hm->endb = len_a;
                hm->len_a = len_b;
                hm->len_b = len_a;
                map[c] = hirsch_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
                map[c] = mirror_hirsch_path(map[c],len_a,len_b);
            }
        }else{
            if(b < numseq){
                map[c] = hirsch_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
            }else{
                if(len_a < len_b){
                    map[c] = hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
                }else{
                    hm->enda = len_b;
                    hm->endb = len_a;
                    hm->len_a = len_b;
                    hm->len_b = len_a;
                    map[c] = hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
                    map[c] = mirror_hirsch_path(map[c],len_a,len_b);
                }
            }
        }
        
        map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

        if(i != numseq-2){
            profile[c] = malloc(sizeof(float)*64*(map[c][0]+2));
            profile[c] = update_only_a(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
        }
            
        aln->sl[c] = map[c][0];
    
        aln->nsip[c] = aln->nsip[a] + aln->nsip[b];
        aln->sip[c] = malloc(sizeof(int)*(aln->nsip[a] + aln->nsip[b]));
        g =0;
        for (j = aln->nsip[a];j--;){
            aln->sip[c][g] = aln->sip[a][j];
            g++;
        }
        for (j = aln->nsip[b];j--;){
            aln->sip[c][g] = aln->sip[b][j];
            g++;
        }

        free(profile[a]);
        free(profile[b]);
    }
    fprintf(stderr,"\r%8.0f percent done\n",100.0);
    free(profile);
    hirsch_mem_free(hm);
    for (i = 32;i--;){
        free(submatrix[i]);
    }
    free(submatrix);
    return map;
}/*}}}*/


/*==== sequence - sequence alignment ==== */
int* hirsch_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta; /*the middle position of the effective alignment region for seq1, 2010-09-28 by Nanjiang*/
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){
        return hirsch_path;
    }
    if(hm->startb  >= hm->endb){
        return hirsch_path;
    }


    hm->enda = mid;

#ifdef DEBUG
    fprintf(stderr,"Forward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
#endif
    hm->f = foward_hirsch_ss_dyn(subm,seq1,seq2,hm);

    hm->starta = mid;
    hm->enda = old_cor[1];

#ifdef DEBUG
    fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
#endif 
    hm->b = backward_hirsch_ss_dyn(subm,seq1,seq2,hm);


    hirsch_path = hirsch_align_two_ss_vector(subm,seq1,seq2,hm,hirsch_path,input_states,old_cor);
    return hirsch_path;
}/*}}}*/
int* hirsch_ss_dyn_new(float**subm, const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
/*sequence to sequence alignment with position specific gap penalties*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta; /*the middle position of the effective alignment region for seq1, 2010-09-28 by Nanjiang*/
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};/*store the original positions*/

    if(hm->starta  >= hm->enda){
        return hirsch_path;
    }
    if(hm->startb  >= hm->endb){
        return hirsch_path;
    }

    hm->enda = mid;

    if ((!gpoArray1) && (!gpeArray1) && (!tgpeArray1) && (!gpoArray2) && (!gpeArray2) && (!tgpeArray2) ){
        hm->f = foward_hirsch_ss_dyn(subm,seq1,seq2,hm);
    }else{
        hm->f = foward_hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1,gpeArray1,tgpeArray1,gpoArray2,gpeArray2,tgpeArray2, hm);
    }
#ifdef DEBUG
    fprintf(stderr,"\n%s: Forward. a:%d-%d	b: %d-%d\n",__FUNCTION__, hm->starta,hm->enda,hm->startb,hm->endb);
    int i;
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif

    hm->starta = mid;
    hm->enda = old_cor[1]; /*recover the original 'enda'*/

    if ((!gpoArray1) && (!gpeArray1) && (!tgpeArray1) && (!gpoArray2) && (!gpeArray2) && (!tgpeArray2) ){
        hm->b = backward_hirsch_ss_dyn(subm,seq1,seq2,hm);
    }else{
        hm->b = backward_hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1,gpeArray1,tgpeArray1,gpoArray2,gpeArray2,tgpeArray2, hm);
    }
#ifdef DEBUG
    fprintf(stderr,"\n%s: Backward. a: %d-%d	b: %d-%d\n",__FUNCTION__, hm->starta,hm->enda,hm->startb,hm->endb);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif 

    if ((!gpoArray1) && (!gpeArray1) && (!tgpeArray1) && (!gpoArray2) && (!gpeArray2) && (!tgpeArray2) ){
        hirsch_path = hirsch_align_two_ss_vector(subm,seq1,seq2,hm,hirsch_path,input_states,old_cor);
    }else{ 
        hirsch_path = hirsch_align_two_ss_vector_new(subm,seq1,seq2,gpoArray1,gpeArray1,tgpeArray1,gpoArray2,gpeArray2,tgpeArray2, hm,hirsch_path,input_states,old_cor);
    }
    return hirsch_path;
}/*}}}*/

int* hirsch_align_two_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;
    
    
    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

#ifdef DEBUG
    fprintf(stderr,"function hirsch_align_two_ss_vector:\n");
#endif
    
    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2]; /*the middle position of sequence b*/
    float sub = 0.0;
    
    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++){
    for(i = old_cor[2]; i < old_cor[3];i++){
    
        sub = abs(middle -i);
        sub /= 1000; 
#ifdef DEBUG
        fprintf(stderr,"%d-%d	%f\n",hm->startb,hm->endb,sub);
#endif
        if(f[i].a+b[i].a-sub > max){
            max = f[i].a+b[i].a-sub;
#ifdef DEBUG
            fprintf(stderr,"aligned->aligned:%f + %f = %f\n",f[i].a,b[i].a,f[i].a+b[i].a);
#endif
            transition = 1;
            c = i;
        }
        if(f[i].a+b[i].ga-gpo-sub > max){
            max = f[i].a+b[i].ga-gpo-sub;
#ifdef DEBUG
            fprintf(stderr,"aligned->gap_a:%f + %f +%f = %f\n",f[i].a,b[i].ga,-gpo-sub,f[i].a+b[i].ga-gpo-sub);
#endif
            transition = 2;
            c = i;
        }
        if(f[i].a+b[i].gb -gpo-sub > max){
            max = f[i].a+b[i].gb - gpo-sub;
#ifdef DEBUG
            fprintf(stderr,"aligned->gap_b:%f + %f +%f = %f\n",f[i].a,b[i].gb,-gpo-sub,f[i].a+b[i].gb-gpo-sub);
#endif
            transition = 3;
            c = i;
        }
        if(f[i].ga+b[i].a - gpo-sub > max){
            max = f[i].ga+b[i].a - gpo-sub;
#ifdef DEBUG
            fprintf(stderr,"gap_a->aligned:%f + %f + %f(gpo) = %f\n",f[i].ga,b[i].a,-gpo-sub,f[i].ga+b[i].a-gpo-sub);
#endif
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            if(f[i].gb+b[i].gb - tgpe-sub > max){
                max = f[i].gb+b[i].gb -tgpe-sub;
#ifdef DEBUG
                fprintf(stderr,"gap_b->gap_b:%f + %f +%f(gpe) =%f \n",f[i].gb, b[i].gb, -tgpe-sub,f[i].gb+b[i].gb-tgpe-sub);
#endif
                transition = 6;
                c = i;
            }
        }else{
            if(f[i].gb+b[i].gb - gpe -sub> max){
                max = f[i].gb+b[i].gb - gpe-sub;
#ifdef DEBUG
                fprintf(stderr,"gap_b->gap_b:%f + %f +%f(gpe) =%f \n",f[i].gb, b[i].gb, -gpe-sub,f[i].gb+b[i].gb-gpe-sub);
#endif
                transition = 6;
                c = i;
            }
        }
        if(f[i].gb+b[i].a - gpo-sub > max){
            max = f[i].gb+b[i].a - gpo-sub;
#ifdef DEBUG
            fprintf(stderr,"gap_b->aligned:%f + %f + %f(gpo) = %f\n",f[i].gb,b[i].a,-gpo-sub,f[i].gb+b[i].a-gpo-sub);
#endif
            transition = 7;
            c = i;
        }
    }
    //i = hm->endb;
    i = old_cor[3];
    sub = abs(middle -i);
    sub /= 1000; 
    
    if(f[i].a+b[i].gb-gpo-sub > max){
        max = f[i].a+b[i].gb - gpo-sub;
#ifdef DEBUG
        fprintf(stderr,"aligned->gap_b:%f + %f +%f = %f\n",f[i].a,b[i].gb,-gpo-sub,f[i].a+b[i].gb-gpo-sub);
#endif
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        if(f[i].gb+b[i].gb -tgpe-sub > max){
            max = f[i].gb+b[i].gb - tgpe-sub;
#ifdef DEBUG
            fprintf(stderr,"gap_b->gap_b:%f + %f +%f(gpe) =%f \n",f[i].gb, b[i].gb, -tgpe-sub,f[i].gb+b[i].gb-tgpe-sub);
#endif
            transition = 6;
            c = i;
        }	
    }else{
        if(f[i].gb+b[i].gb - gpe-sub > max){
            max = f[i].gb+b[i].gb - gpe-sub;
#ifdef DEBUG
            fprintf(stderr,"gap_b->gap_b:%f + %f +%f(gpe) =%f \n",f[i].gb, b[i].gb, -gpe-sub,f[i].gb+b[i].gb-gpe-sub);
#endif
            transition = 6;
            c = i;
        }
    }
    
    
#ifdef DEBUG
    fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
#endif
    
    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1
            
            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;
            
#ifdef DEBUG
            /*fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);*/
            /*fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);*/
            fprintf(stderr,"Forward:\n");
#endif
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
#ifdef DEBUG
            fprintf(stderr,"Using this for start:%f	%f	%f\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
#endif
            
            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            
            hm->startb = old_cor[2];
            hm->endb = c-1;
#ifdef DEBUG
            fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

#ifdef DEBUG
            fprintf(stderr,"Backward:\n");
#endif
            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];
    
#ifdef DEBUG
            fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            break;
        case 2:// a -> ga = 2
            
            hirsch_path[old_cor[4]] = c;
#ifdef DEBUG
            fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
#endif
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            
            
            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            
            hm->startb = old_cor[2];
            hm->endb = c-1;
#ifdef DEBUG
            fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];
#ifdef DEBUG
            fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            break;
        case 3:// a -> gb = 3
            
            hirsch_path[old_cor[4]] = c;
#ifdef DEBUG
            fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
#endif
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            
            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            
            hm->startb = old_cor[2];
            hm->endb = c-1;
#ifdef DEBUG
            fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];
    
#ifdef DEBUG
            fprintf(stderr,"Following last: %d\n",c+1);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
#ifdef DEBUG
            fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
#endif

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;
            
            hm->starta = old_cor[0];
            hm->enda = old_cor[4];
            
            hm->startb = old_cor[2];
            hm->endb = c-1;
#ifdef DEBUG
            fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];
    
#ifdef DEBUG
            fprintf(stderr,"Following last: %d\n",c+1);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            break;
        case 6://gb->gb = 6;
            
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;
            
            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
#ifdef DEBUG
            fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];
    
#ifdef DEBUG
            fprintf(stderr,"Following last: %d\n",c+1);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            break;
        case 7://gb->a = 7;
            
            hirsch_path[old_cor[4]+1] = c+1;
#ifdef DEBUG
            fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
#endif
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;
            
            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
#ifdef DEBUG
            fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];
    
#ifdef DEBUG
            fprintf(stderr,"Following last: %d\n",c+1);
#endif
            hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
            break;
    }
        
    return hirsch_path;
}/*}}}*/
int* hirsch_align_two_ss_vector_new(float**subm,const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[] )/*{{{*/
{
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;
    //transition code:
    // a -> a  = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga-> ga = 4
    // ga-> a  = 5
    // gb-> gb = 6;
    // gb-> a  = 7;

    int mid_a = old_cor[4]; 
    mid_a ++;

#ifdef DEBUG
    fprintf(stderr,"\nFunction:%s\n", __FUNCTION__);
#endif
    
    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2]; /*the middle position of sequence b*/
    float sub = 0.0;
    
    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++){
    for(i = old_cor[2]; i < old_cor[3];i++){
        sub = abs(middle -i); /*position related value ?? what is it for*/
        sub /= 1000; 
#ifdef DEBUG
        fprintf(stderr,"startb:endb %d-%d	sub=%f\n",hm->startb,hm->endb,sub);
#endif
        if(f[i].a+b[i].a-sub > max){
            max = f[i].a+b[i].a-sub;
#ifdef DEBUG
            fprintf(stderr,"i=%d\t aligned->aligned:%f + %f = %f\n",i,f[i].a,b[i].a,f[i].a+b[i].a);
#endif
            transition = 1;
            c = i;
        }
        /*if(f[i].a+b[i].ga-gpo-sub > max){*/
        if(f[i].a+b[i].ga-gpoArray1[mid_a]-sub > max){
            max = f[i].a+b[i].ga-gpoArray1[mid_a]-sub;
#ifdef DEBUG
            fprintf(stderr,"i=%d\t aligned->gap_a:%f + %f +%f = %f\n",i,f[i].a,b[i].ga,-gpoArray1[mid_a]-sub,f[i].a+b[i].ga-gpoArray1[mid_a]-sub);
#endif
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb -gpo-sub > max){*/
        if(f[i].a+b[i].gb -gpoArray2[i]-sub > max){
            max = f[i].a+b[i].gb - gpoArray2[i]-sub;
#ifdef DEBUG
            fprintf(stderr,"i=%d\t aligned->gap_b:%f + %f +%f = %f\n",i,f[i].a,b[i].gb,-gpoArray2[i]-sub,f[i].a+b[i].gb-gpoArray2[i]-sub);
#endif
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a - gpo-sub > max){*/
    if(f[i].ga+b[i].a - gpoArray1[mid_a-1]-sub > max){
        max = f[i].ga+b[i].a - gpoArray1[mid_a-1]-sub;
#ifdef DEBUG
        fprintf(stderr,"i=%d\t gap_a->aligned:%f + %f + %f(gpo) = %f\n",i, f[i].ga,b[i].a,-gpoArray1[mid_a-1]-sub,f[i].ga+b[i].a-gpoArray1[mid_a-1]-sub);
#endif
        transition = 5;
        c = i;
    }

    if(hm->startb == 0){
        /*if(f[i].gb+b[i].gb - tgpe-sub > max){*/
        if(f[i].gb+b[i].gb - tgpeArray2[i]-sub > max){
            max = f[i].gb+b[i].gb -tgpeArray2[i]-sub;
#ifdef DEBUG
            fprintf(stderr,"i=%d\t gap_b->gap_b:%f + %f +%f(tgpe) =%f \n",i, f[i].gb, b[i].gb, -tgpeArray2[i]-sub,f[i].gb+b[i].gb-tgpeArray2[i]-sub);
#endif
            transition = 6;
            c = i;
        }
    }else{
        /*if(f[i].gb+b[i].gb - gpe -sub> max){*/
        if(f[i].gb+b[i].gb - gpeArray2[i] -sub> max){
            max = f[i].gb+b[i].gb - gpeArray2[i]-sub;
#ifdef DEBUG
            fprintf(stderr,"i=%d\t gap_b->gap_b:%f + %f +%f(gpe) =%f \n",i, f[i].gb, b[i].gb, -gpeArray2[i]-sub,f[i].gb+b[i].gb-gpeArray2[i]-sub);
#endif
            transition = 6;
            c = i;
        }
    }
    /*if(f[i].gb+b[i].a - gpo-sub > max){*/
    if(f[i].gb+b[i].a - gpoArray2[i-1]-sub > max){
        max = f[i].gb+b[i].a - gpoArray2[i-1]-sub;
#ifdef DEBUG
        fprintf(stderr,"i=%d\t gap_b->aligned:%f + %f + %f(gpo) = %f\n",i,f[i].gb,b[i].a,-gpoArray2[i-1]-sub,f[i].gb+b[i].a-gpoArray2[i-1]-sub);
#endif
        transition = 7;
        c = i;
    }
}
//i = hm->endb;
i = old_cor[3];
sub = abs(middle -i);
sub /= 1000; 

/*if(f[i].a+b[i].gb-gpo-sub > max){*/
if(f[i].a+b[i].gb-gpoArray2[i]-sub > max){
    max = f[i].a+b[i].gb - gpoArray2[i]-sub;
#ifdef DEBUG
    fprintf(stderr,"i=%d\t aligned->gap_b:%f + %f +%f = %f\n",i, f[i].a,b[i].gb,-gpoArray2[i]-sub,f[i].a+b[i].gb-gpoArray2[i]-sub);
#endif
    transition = 3;
    c = i;
}
if(hm->endb == hm->len_b){
    /*if(f[i].gb+b[i].gb -tgpe-sub > max){*/
    if(f[i].gb+b[i].gb -tgpeArray2[i]-sub > max){
        max = f[i].gb+b[i].gb - tgpeArray2[i]-sub;
#ifdef DEBUG
        fprintf(stderr,"i=%d\t gap_b->gap_b:%f + %f +%f(tgpe) =%f \n",i, f[i].gb, b[i].gb, -tgpeArray2[i]-sub,f[i].gb+b[i].gb-tgpeArray2[i]-sub);
#endif
        transition = 6;
        c = i;
    }	
}else{
    /*if(f[i].gb+b[i].gb - gpe-sub > max){*/
    if(f[i].gb+b[i].gb - gpeArray2[i]-sub > max){
        max = f[i].gb+b[i].gb - gpeArray2[i]-sub;
#ifdef DEBUG
        fprintf(stderr,"i=%d\t gap_b->gap_b:%f + %f +%f(gpe) =%f \n",i, f[i].gb, b[i].gb, -gpeArray2[i]-sub,f[i].gb+b[i].gb-gpeArray2[i]-sub);
#endif
        transition = 6;
        c = i;
    }
}

#ifdef DEBUG
fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
#endif

j = hirsch_path[0];
switch(transition){
    case 1: //a -> a = 1
        
        hirsch_path[old_cor[4]] = c;
        hirsch_path[old_cor[4]+1] = c+1;
        
#ifdef DEBUG
        /*fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);*/
        /*fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);*/
        fprintf(stderr,"Forward:\n");
#endif
        //foward:
        hm->f[0].a = input_states[0];
        hm->f[0].ga = input_states[1];
        hm->f[0].gb = input_states[2];
        hm->b[0].a = 0.0;
        hm->b[0].ga = -FLOATINFTY;
        hm->b[0].gb = -FLOATINFTY;
#ifdef DEBUG
        fprintf(stderr,"Using this for start:%f	%f	%f\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
#endif
        
        hm->starta = old_cor[0];
        hm->enda = old_cor[4]-1;
        
        hm->startb = old_cor[2];
        hm->endb = c-1;
#ifdef DEBUG
        fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path);

#ifdef DEBUG
        fprintf(stderr,"Backward:\n");
#endif
        //backward:
        hm->starta = old_cor[4]+1;
        hm->enda = old_cor[1];
        hm->startb = c+1;
        hm->endb = old_cor[3];
        hm->f[0].a = 0.0;
        hm->f[0].ga = -FLOATINFTY;
        hm->f[0].gb = -FLOATINFTY;
        hm->b[0].a = input_states[3];
        hm->b[0].ga = input_states[4];
        hm->b[0].gb = input_states[5];

#ifdef DEBUG
        fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        break;
    case 2:// a -> ga = 2
        
        hirsch_path[old_cor[4]] = c;
#ifdef DEBUG
        fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
#endif
        //foward:
        hm->f[0].a = input_states[0];
        hm->f[0].ga = input_states[1];
        hm->f[0].gb = input_states[2];
        hm->b[0].a = 0.0;
        hm->b[0].ga = -FLOATINFTY;
        hm->b[0].gb = -FLOATINFTY;
        
        
        hm->starta = old_cor[0];
        hm->enda = old_cor[4]-1;
        
        hm->startb = old_cor[2];
        hm->endb = c-1;
#ifdef DEBUG
        fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        

        //backward:
        hm->starta = old_cor[4];
        hm->enda = old_cor[1];
        hm->startb = c+1;
        hm->endb = old_cor[3];
        hm->f[0].a = -FLOATINFTY;
        hm->f[0].ga = 0.0;
        hm->f[0].gb = -FLOATINFTY;
        hm->b[0].a = input_states[3];
        hm->b[0].ga = input_states[4];
        hm->b[0].gb = input_states[5];
#ifdef DEBUG
        fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        break;
    case 3:// a -> gb = 3
        
        hirsch_path[old_cor[4]] = c;
#ifdef DEBUG
        fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
#endif
        //foward:
        hm->f[0].a = input_states[0];
        hm->f[0].ga = input_states[1];
        hm->f[0].gb = input_states[2];
        hm->b[0].a = 0.0;
        hm->b[0].ga = -FLOATINFTY;
        hm->b[0].gb = -FLOATINFTY;
        
        hm->starta = old_cor[0];
        hm->enda = old_cor[4]-1;
        
        hm->startb = old_cor[2];
        hm->endb = c-1;
#ifdef DEBUG
        fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        //backward:
        hm->starta = old_cor[4]+1;
        hm->enda = old_cor[1];
        hm->startb = c;
        hm->endb = old_cor[3];
        hm->f[0].a = -FLOATINFTY;
        hm->f[0].ga = -FLOATINFTY;
        hm->f[0].gb = 0.0;
        hm->b[0].a = input_states[3];
        hm->b[0].ga = input_states[4];
        hm->b[0].gb = input_states[5];

#ifdef DEBUG
        fprintf(stderr,"Following last: %d\n",c+1);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        break;
    case 5://ga -> a = 5
        hirsch_path[old_cor[4]+1] = c+1;
#ifdef DEBUG
        fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
#endif

        //foward:
        hm->f[0].a = input_states[0];
        hm->f[0].ga = input_states[1];
        hm->f[0].gb = input_states[2];
        hm->b[0].a = -FLOATINFTY;
        hm->b[0].ga = 0.0;
        hm->b[0].gb = -FLOATINFTY;
        
        hm->starta = old_cor[0];
        hm->enda = old_cor[4];
        
        hm->startb = old_cor[2];
        hm->endb = c-1;
#ifdef DEBUG
        fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);

        //backward:
        hm->starta = old_cor[4]+1;
        hm->enda = old_cor[1];
        hm->startb = c+1;
        hm->endb = old_cor[3];
        hm->f[0].a = 0.0;
        hm->f[0].ga = -FLOATINFTY;
        hm->f[0].gb = -FLOATINFTY;
        hm->b[0].a = input_states[3];
        hm->b[0].ga = input_states[4];
        hm->b[0].gb = input_states[5];

#ifdef DEBUG
        fprintf(stderr,"Following last: %d\n",c+1);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        break;
    case 6://gb->gb = 6;
        
        //foward:
        hm->f[0].a = input_states[0];
        hm->f[0].ga = input_states[1];
        hm->f[0].gb = input_states[2];
        hm->b[0].a = -FLOATINFTY;
        hm->b[0].ga = -FLOATINFTY;
        hm->b[0].gb = 0.0;
        
        hm->starta = old_cor[0];
        hm->enda = old_cor[4]-1;
        hm->startb = old_cor[2];
        hm->endb = c;
#ifdef DEBUG
        fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);

        //backward:
        hm->starta = old_cor[4]+1;
        hm->enda = old_cor[1];
        hm->startb = c;
        hm->endb = old_cor[3];
        hm->f[0].a = -FLOATINFTY;
        hm->f[0].ga = -FLOATINFTY;
        hm->f[0].gb = 0.0;
        hm->b[0].a = input_states[3];
        hm->b[0].ga = input_states[4];
        hm->b[0].gb = input_states[5];

#ifdef DEBUG
        fprintf(stderr,"Following last: %d\n",c+1);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);
        break;
    case 7://gb->a = 7;
        
        hirsch_path[old_cor[4]+1] = c+1;
#ifdef DEBUG
        fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
#endif
        //foward:
        hm->f[0].a = input_states[0];
        hm->f[0].ga = input_states[1];
        hm->f[0].gb = input_states[2];
        hm->b[0].a = -FLOATINFTY;
        hm->b[0].ga = -FLOATINFTY;
        hm->b[0].gb = 0.0;
        
        hm->starta = old_cor[0];
        hm->enda = old_cor[4]-1;
        hm->startb = old_cor[2];
        hm->endb = c;
#ifdef DEBUG
        fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path);

        //backward:
        hm->starta = old_cor[4]+1;
        hm->enda = old_cor[1];
        hm->startb = c+1;
        hm->endb = old_cor[3];
        hm->f[0].a = 0.0;
        hm->f[0].ga = -FLOATINFTY;
        hm->f[0].gb = -FLOATINFTY;
        hm->b[0].a = input_states[3];
        hm->b[0].ga = input_states[4];
        hm->b[0].gb = input_states[5];

#ifdef DEBUG
        fprintf(stderr,"Following last: %d\n",c+1);
#endif
        hirsch_path = hirsch_ss_dyn_new(subm,seq1,seq2,gpoArray1, gpeArray1, tgpeArray1, gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path);
        break;
}
    
return hirsch_path;
}/*}}}*/

struct states* foward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)/*{{{*/
{
    struct states* s = hm->f;
    float *subp = 0;
    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;
    register float xa = 0;
    register float xga = 0;
    register int i = 0;
    register int j = 0;


    s[startb].a = s[0].a;
    s[startb].ga = s[0].ga;
    s[startb].gb = s[0].gb;       /*initializing first row*/
    if(startb){
        for (j = startb+1; j < endb;j++){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);
            s[j].gb = -FLOATINFTY;
        }
    }else{
        for (j = startb+1; j < endb;j++){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;
            s[j].gb = -FLOATINFTY;
        }
    }
    s[endb].a = -FLOATINFTY;
    s[endb].ga = -FLOATINFTY;
    s[endb].gb = -FLOATINFTY;

    seq2--;
    for (i = starta;i < enda;i++){
        subp = subm[seq1[i]];

        pa = s[startb].a;
        pga = s[startb].ga;
        pgb = s[startb].gb;
        s[startb].a = -FLOATINFTY;
        s[startb].ga = -FLOATINFTY;

        xa = s[startb].a;
        xga = s[startb].ga;

        if(startb){
            s[startb].gb = MAX(pgb - gpe,pa - gpo);
        }else{
            s[startb].gb = MAX(pgb,pa) - tgpe;
        }
        for (j = startb+1; j < endb;j++){
            ca = s[j].a;
            pa = MAX3(pa,pga-gpo,pgb-gpo);
            pa += subp[seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
            s[j].ga = MAX(xga-gpe,xa-gpo);

            pgb = s[j].gb;
            s[j].gb = MAX(pgb-gpe ,ca-gpo);

            pa = ca;

            xa = s[j].a;
            xga = s[j].ga;

        }
        ca = s[j].a;
        pa = MAX3(pa,pga-gpo,pgb-gpo);
        pa += subp[seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
        if (endb != hm->len_b){
            s[j].gb = MAX(s[j].gb-gpe ,ca-gpo);
        }else{
            s[j].gb = MAX(s[j].gb,ca)-tgpe;
        }

    }
    return s;
}/*}}}*/
struct states* foward_hirsch_ss_dyn_new(float**subm,const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm)/*{{{*/
{
    struct states* s = hm->f;
    float *subp = 0;
    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;
    register float xa = 0;
    register float xga = 0;
    register int i = 0;
    register int j = 0;


    s[startb].a = s[0].a;
    s[startb].ga = s[0].ga;
    s[startb].gb = s[0].gb;       /*initializing first row*/
    if(startb){
        for (j = startb+1; j < endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);*/
            s[j].ga = MAX(s[j-1].ga - gpeArray1[starta],s[j-1].a - gpoArray1[starta]);
            s[j].gb = -FLOATINFTY;
        }
    }else{
        for (j = startb+1; j < endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpeArray1[starta];
            s[j].gb = -FLOATINFTY;
        }
    }
    s[endb].a = -FLOATINFTY;
    s[endb].ga = -FLOATINFTY;
    s[endb].gb = -FLOATINFTY;

    seq2--;
    for (i = starta;i < enda;i++){
        subp = subm[seq1[i]];

        pa = s[startb].a;
        pga = s[startb].ga;
        pgb = s[startb].gb;
        s[startb].a = -FLOATINFTY;
        s[startb].ga = -FLOATINFTY;

        xa = s[startb].a;
        xga = s[startb].ga;

        if(startb){
            /*s[startb].gb = MAX(pgb - gpe,pa - gpo);*/
            s[startb].gb = MAX(pgb - gpeArray2[startb],pa - gpoArray2[startb]);
        }else{
            /*s[startb].gb = MAX(pgb,pa) - tgpe;*/
            s[startb].gb = MAX(pgb,pa) - tgpeArray2[startb];
        }
        for (j = startb+1; j < endb;j++){
            ca = s[j].a; /*current a*/
            /*pa = MAX3(pa,pga-gpo,pgb-gpo);*/
            /*pa = MAX3(pa,pga-gpoArray1[i-1],pgb-gpoArray2[j-1]);*/
            pa = MAX3(pa,pga-gpoArray1[i-1],pgb-gpoArray2[j-1]);
#ifdef DEBUG_MEM
            if ((i-1)<-1 || (i-1) > hm->len_a ){
                fprintf(stderr,"invalid read, idx=%d, seqLength=%d. src at line %d file %s\n", i-1,hm->len_a,__LINE__, __FILE__ );
            }
            if ((j-1)<-1 || (j-1) > hm->len_b ){
                fprintf(stderr,"invalid read, idx=%d, seqLength=%d. src at line %d file %s\n", j-1,hm->len_b,__LINE__, __FILE__ );
            }
#endif
            pa += subp[seq2[j]];

            s[j].a = pa;

            pga = s[j].ga; /*xga and pga are different*/
            //s[j].ga = MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
            /*s[j].ga = MAX(xga-gpe,xa-gpo);*/
            s[j].ga = MAX(xga-gpeArray1[i],xa-gpoArray1[i]);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb-gpe ,ca-gpo);*/
            s[j].gb = MAX(pgb-gpeArray2[j] ,ca-gpoArray2[j]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
#ifdef DEBUG_SS
            fprintf(stderr,"%s: i=%d, j=%d, a=%f\tga=%f\tgb=%f\n",__FUNCTION__, i,j,s[j].a,s[j].ga,s[j].gb);
#endif 
        }
        ca = s[j].a;
        /*pa = MAX3(pa,pga-gpo,pgb-gpo);*/
        pa = MAX3(pa,pga-gpoArray1[i-1],pgb-gpoArray2[j-1]);

        pa += subp[seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
        if (endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb-gpe ,ca-gpo);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j] ,ca-gpoArray2[j]);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)-tgpe;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j];
        }

    }
    return s;
}/*}}}*/

struct states* backward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)/*{{{*/
{
    struct states* s = hm->b;
    float *subp = 0;
    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    s[endb].a = s[0].a ;
    s[endb].ga = s[0].ga;
    s[endb].gb = s[0].gb;


    //init of first row;

    //j = endb-startb;
    if(endb != hm->len_b){
        for(j = endb-1;j > startb;j--){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);	
            s[j].gb = -FLOATINFTY;
        }
    }else{
        for(j = endb-1;j > startb;j--){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;
            s[j].gb = -FLOATINFTY;
        }
    }


    s[startb].a = -FLOATINFTY;
    s[startb].ga = -FLOATINFTY;
    s[startb].gb = -FLOATINFTY;

    i = enda-starta;
    seq1+= starta;
    while(i--){
        subp = subm[seq1[i]];
        pa = s[endb].a;
        pga = s[endb].ga;
        pgb = s[endb].gb;
        s[endb].a = -FLOATINFTY;
        s[endb].ga = -FLOATINFTY;

        xa = s[endb].a;
        xga = s[endb].ga;

        if(endb != hm->len_b){
            s[endb].gb = MAX(pgb-gpe,pa-gpo);
        }else{
            s[endb].gb = MAX(pgb,pa)-tgpe;
        }

        for(j = endb-1;j > startb;j--){

            ca = s[j].a;

            pa = MAX3(pa,pga - gpo,pgb-gpo);

            pa += subp[seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

            s[j].ga = MAX(xga-gpe,xa-gpo);

            pgb = s[j].gb;
            s[j].gb = MAX(pgb-gpe,ca-gpo);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        ca = s[j].a;

        pa = MAX3(pa,pga - gpo,pgb-gpo);

        pa += subp[seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

        if(startb){
            s[j].gb = MAX(s[j].gb-gpe,ca-gpo);
        }else{
            s[j].gb = MAX(s[j].gb,ca)-tgpe;
        }


    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ss_dyn_new(float**subm,const int* seq1,const int* seq2,const float* gpoArray1, const float*gpeArray1, const float* tgpeArray1, const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm)/*{{{*/
{

    struct states* s = hm->b;
    float *subp = 0;
    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    s[endb].a = s[0].a ;
    s[endb].ga = s[0].ga;
    s[endb].gb = s[0].gb;


    //init of first row;

    //j = endb-startb;
    if(endb != hm->len_b){ /*initializing first row*/
        for(j = endb-1;j > startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);	*/
            s[j].ga = MAX(s[j+1].ga-gpeArray1[enda-1],s[j+1].a-gpoArray1[enda-1]);	
            s[j].gb = -FLOATINFTY;
        }
    }else{
        for(j = endb-1;j > startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;*/
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpeArray1[enda-1];
            s[j].gb = -FLOATINFTY;
        }
    }


    s[startb].a = -FLOATINFTY;
    s[startb].ga = -FLOATINFTY;
    s[startb].gb = -FLOATINFTY;

    i = enda-starta;
    seq1+= starta;
    while(i--){
        subp = subm[seq1[i]];
        pa = s[endb].a;
        pga = s[endb].ga;
        pgb = s[endb].gb;
        s[endb].a = -FLOATINFTY;
        s[endb].ga = -FLOATINFTY;

        xa = s[endb].a;
        xga = s[endb].ga;

        if(endb != hm->len_b){
            /*s[endb].gb = MAX(pgb-gpe,pa-gpo);*/
            s[endb].gb = MAX(pgb-gpeArray2[endb-1],pa-gpoArray2[endb-1]);
        }else{
            /*s[endb].gb = MAX(pgb,pa)-tgpe;*/
            s[endb].gb = MAX(pgb,pa)-tgpeArray2[endb-1];
        }

        for(j = endb-1;j > startb;j--){

            ca = s[j].a;

            /*pa = MAX3(pa,pga - gpo,pgb-gpo);*/
            pa = MAX3(pa,pga - gpoArray1[i+starta+1],pgb-gpoArray2[j+1]);

            pa += subp[seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

            /*s[j].ga = MAX(xga-gpe,xa-gpo);*/
            s[j].ga = MAX(xga-gpeArray1[i+starta],xa-gpoArray1[i+starta]);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb-gpe,ca-gpo);*/
            s[j].gb = MAX(pgb-gpeArray2[j],ca-gpoArray2[j]);

#ifdef DEBUG_SS
            fprintf(stderr,"%s: i=%d, j=%d, a=%f\tga=%f\tgb=%f\n",__FUNCTION__, i,j,s[j].a,s[j].ga,s[j].gb);
#endif 
            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga - gpo,pgb-gpo);*/
        pa = MAX3(pa,pga - gpoArray1[i+starta+1],pgb-gpoArray2[j+1]);
        /*pa = MAX3(pa,pga - gpoArray1[i+starta],pgb-gpoArray2[j]);*/
        /*pa = MAX3(pa,pga - gpeArray1[i+starta],pgb-gpeArray2[j]);*/

        pa += subp[seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);

        if(startb){
            /*s[j].gb = MAX(s[j].gb-gpe,ca-gpo);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j],ca-gpoArray2[j]);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)-tgpe;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j];
        }
    }		
    return s;
}/*}}}*/

/*==== consensus - sequence alignment ==== */
int* hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip)/*{{{*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


    if(hm->starta  >= hm->enda){
        return hirsch_path;
    }
    if(hm->startb  >= hm->endb){
        return hirsch_path;
    }

    hm->enda = mid;
    hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);

    /*int i;
      fprintf(stderr,"FOWARD\n");
      for (i = hm->startb; i <= hm->endb;i++){
      fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
      }*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);

    /*fprintf(stderr,"BaCKWARD\n");
      for (i = hm->startb; i <= hm->endb;i++){
      fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
      }*/

    hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
    return hirsch_path;
}/*}}}*/
int* hirsch_ps_dyn_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path,int sip)/*{{{*/
{ /*sip: nsip for prof1 */
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    hm->f = foward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);
    /*hm->f = foward_hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/

#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    hm->b = backward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);
    /*hm->b = backward_hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);*/
    hirsch_path = hirsch_align_two_ps_vector_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);
    /*hirsch_path = hirsch_align_two_ps_vector_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);*/
    return hirsch_path;
}/*}}}*/
int* hirsch_ps_dyn_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path,int sip)/*{{{*/
{ /*sip: nsip for prof1 */
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    hm->f = foward_hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);
    /*hm->f = foward_hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/

#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    hm->b = backward_hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);
    /*hm->b = backward_hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);*/
    hirsch_path = hirsch_align_two_ps_vector_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);
    /*hirsch_path = hirsch_align_two_ps_vector_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);*/
    return hirsch_path;
}/*}}}*/
int* hirsch_ps_dyn_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1,  int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->f = foward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->f = foward_hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1,  sip);

#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->b = backward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->b = backward_hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, sip);

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);*/
    /*hirsch_path = hirsch_align_two_ps_vector_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);*/
    hirsch_path = hirsch_align_two_ps_vector_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,nResArray1,  sip);
    return hirsch_path;
}/*}}}*/
int* hirsch_ps_dyn_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
    /*using wf1 and wf2 to make the output the same as the original kalign2
     * when constant gap penalty is used*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign
sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->f = foward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->f = foward_hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, weightArray1, sip);

#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->b = backward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->b = backward_hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, weightArray1, sip);

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);*/
    /*hirsch_path = hirsch_align_two_ps_vector_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);*/
    hirsch_path = hirsch_align_two_ps_vector_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,nResArray1, weightArray1, sip);
    return hirsch_path;
}/*}}}*/
int* hirsch_ps_dyn_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign
sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

wf1=1
wf2=nSeq

*/
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->f = foward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->f = foward_hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, weightArray1, sip);

#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->b = backward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->b = backward_hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, weightArray1, sip);

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);*/
    /*hirsch_path = hirsch_align_two_ps_vector_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);*/
    hirsch_path = hirsch_align_two_ps_vector_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,nResArray1, weightArray1, sip);
    return hirsch_path;
}/*}}}*/
int* hirsch_ps_dyn_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm, int* hirsch_path,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign
sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

wf1=1
wf2=1

*/
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->f = foward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->f = foward_hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, weightArray1, sip);

#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);*/
    /*hm->b = backward_hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,sip);*/
    hm->b = backward_hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,nResArray1, weightArray1, sip);

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);*/
    /*hirsch_path = hirsch_align_two_ps_vector_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,sip);*/
    hirsch_path = hirsch_align_two_ps_vector_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,input_states,old_cor,nResArray1, weightArray1, sip);
    return hirsch_path;
}/*}}}*/

int* hirsch_align_two_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)/*{{{*/
{
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;

    const float open = gpo * sip;


    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;


    prof1+= ((old_cor[4]+1)<<6);

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
        if(f[i].a+b[i].ga-open-sub > max){
            max = f[i].a+b[i].ga-open-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        if(f[i].a+b[i].gb+prof1[27]-sub > max){
            max = f[i].a+b[i].gb+prof1[27]-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        if(f[i].ga+b[i].a-open-sub > max){
            max = f[i].ga+b[i].a-open-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                max = f[i].gb+b[i].gb+prof1[29]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                max = f[i].gb+b[i].gb+prof1[28]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        if(f[i].gb+b[i].a+prof1[-37]-sub > max){
            max = f[i].gb+b[i].a+prof1[-37]-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }
    //i = hm->endb;
    i = old_cor[3];

    sub = abs(middle -i);
    sub /= 1000; 
    if(f[i].a+b[i].gb+prof1[27]-sub > max){
        max = f[i].a+b[i].gb+prof1[27]-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        if(f[i].gb+b[i].gb+prof1[29]-sub > max){
            max = f[i].gb+b[i].gb+prof1[29]-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        if(f[i].gb+b[i].gb+prof1[28]-sub > max){
            max = f[i].gb+b[i].gb+prof1[28]-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }



    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);			


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
            break;
    }

    return hirsch_path;
}/*}}}*/
int* hirsch_align_two_ps_vector_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)/*{{{*/
{/*return to the original ga gb calculation*/ 
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;

    const float open = gpo * sip;
    const float ext = gpe * sip;
    const float text = tgpe * sip;

    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;

    //prof1 += 64*(mid_a+1)
    prof1+= ((old_cor[4]+1)<<6);

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[i]*sip != open|| gpoArray2[i-1]*sip != (open)  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[i] =%f, gpoArray2[i-1]= %f, (open) =%f\n", old_cor[4]+1, i, gpoArray2[i]*sip, gpoArray2[i-1]*sip, open);
    }
    if (gpeArray2[i]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[i]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
        /*if(f[i].a+b[i].ga-open-sub > max)*/
        if(f[i].a+b[i].ga-gpoArray2[i]*sip-sub > max)
            /*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
        {
            max = f[i].a+b[i].ga - gpoArray2[i]*sip -sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        if(f[i].a+b[i].gb+prof1[27]-sub > max)
        /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
        {
            max = f[i].a+b[i].gb+prof1[27]-sub;
            //fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a-open-sub > max)*/
        if(f[i].ga+b[i].a-gpoArray2[i-1]*sip-sub > max)
            /*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
        {
            max = f[i].ga+b[i].a - gpoArray2[i-1]*sip-sub;
            //fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0)
        {
            if(f[i].gb+b[i].gb+prof1[29]-sub > max)
                /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
            {
                max = f[i].gb+b[i].gb+prof1[29]-sub;
                //fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            if(f[i].gb+b[i].gb+prof1[28]-sub > max)
                /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
            {
                max = f[i].gb+b[i].gb+prof1[28]-sub;
                // fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        if(f[i].gb+b[i].a+prof1[-37]-sub > max)
            /*if(f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub > max)*/
        {
            max = f[i].gb+b[i].a+prof1[-37]-sub;
            // fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }
    //i = hm->endb;
    i = old_cor[3];

    sub = abs(middle -i);
    sub /= 1000; 
    if(f[i].a+b[i].gb+prof1[27]-sub > max){
        max = f[i].a+b[i].gb+prof1[27]-sub;
        //fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        if(f[i].gb+b[i].gb+prof1[29]-sub > max){
            max = f[i].gb+b[i].gb+prof1[29]-sub;
            //fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        if(f[i].gb+b[i].gb+prof1[28]-sub > max){
            max = f[i].gb+b[i].gb+prof1[28]-sub;
            //fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }

    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);			


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new0(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
    }
    return hirsch_path;
}/*}}}*/
int* hirsch_align_two_ps_vector_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)/*{{{*/
{
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;

    const float open = gpo * sip;
    const float ext = gpe * sip;
    const float text = tgpe * sip;

    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;


    //prof1 += 64*(mid_a+1)
    prof1+= ((old_cor[4]+1)<<6);

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-open)) > FLOAT_ZERO|| fabs(prof1[-37] - (-open))>FLOAT_ZERO|| fabs(prof1[27] - (-gpoArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-ext))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-text))>FLOAT_ZERO||fabs( prof1[29] - (-tgpeArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        /*if(f[i].a+b[i].ga-open-sub > max)*/
        if(f[i].a+b[i].ga+prof1[27]-sub > max)
        {
            max = f[i].a+b[i].ga+prof1[27]-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
        if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max){
            max = f[i].a+b[i].gb-gpoArray2[i]*sip-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a-open-sub > max)*/
        if(f[i].ga+b[i].a+prof1[-37]-sub > max)
        {
            max = f[i].ga+b[i].a+prof1[-37]-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
            if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)
            {
                max = f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
            if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)
            {
                max = f[i].gb+b[i].gb-gpeArray2[i]*sip-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        /*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
        if(f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub > max)
        {
            max = f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }
    //i = hm->endb;
    i = old_cor[3];

    sub = abs(middle -i);
    sub /= 1000; 
    /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
    if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)
    {
        max = f[i].a+b[i].gb-gpoArray2[i]*sip-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
        if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)
        {
            max = f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
        if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)
        {
            max = f[i].gb+b[i].gb-gpeArray2[i]*sip-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[i]*sip))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-tgpeArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);			


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,sip);
            break;
    }

    return hirsch_path;
}/*}}}*/
int* hirsch_align_two_ps_vector_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, int sip)/*{{{*/
{/* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/ 
    struct states* f = hm->f;
    struct states* b = hm->b;
    register int i,j,c;
    int transition = -1;

    register int idx1 = 0; /*index for nResArray1*/
    
    const float open = gpo * sip;
    const float ext = gpe * sip;
    const float text = tgpe * sip;

    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;

    float wfi = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wfim1 = 0.0; /*weight for i-1*/
    float wfip1 = 0.0; /*weight for i+1*/

    //prof1 += 64*(mid_a+1)
    prof1+= ((old_cor[4]+1)<<6);
    idx1 += (old_cor[4]+1);


    wfi = ((double)sip/(double)nResArray1[idx1]);
    wfim1 = ((double)sip/(double)nResArray1[idx1-1]);
    wfip1 = ((double)sip/(double)nResArray1[idx1+1]);

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-open)) > FLOAT_ZERO|| fabs(prof1[-37] - (-open))>FLOAT_ZERO|| fabs(prof1[27] - (-gpoArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-ext))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-text))>FLOAT_ZERO||fabs( prof1[29] - (-tgpeArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        /*if(f[i].a+b[i].ga-open-sub > max)*/
        /*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
        if(f[i].a+b[i].ga+prof1[27]*wfi-sub > max)
        {
            max = f[i].a+b[i].ga+prof1[27]*wfi-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
        /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
        if(f[i].a+b[i].gb-gpoArray2[i]*nResArray1[idx1]-sub > max)
        {
            max = f[i].a+b[i].gb-gpoArray2[i]*nResArray1[idx1]-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a-open-sub > max)*/
        /*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
        if(f[i].ga+b[i].a+prof1[-37]*wfim1-sub > max)
        {
            max = f[i].ga+b[i].a+prof1[-37]*wfim1-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
            /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-tgpeArray2[i]*nResArray1[idx1]-sub > max)
            {
                max = f[i].gb+b[i].gb-tgpeArray2[i]*nResArray1[idx1]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
            /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-gpeArray2[i]*nResArray1[idx1]-sub > max)
            {
                max = f[i].gb+b[i].gb-gpeArray2[i]*nResArray1[idx1]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        /*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
        /*if(f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub > max)*/
        if(f[i].gb+b[i].a-gpoArray2[i-1]*nResArray1[idx1-1]-sub > max)
        {
            max = f[i].gb+b[i].a-gpoArray2[i-1]*nResArray1[idx1-1]-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }
    //i = hm->endb;
    i = old_cor[3];

    sub = abs(middle -i);
    sub /= 1000; 
    /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
    /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
    if(f[i].a+b[i].gb-gpoArray2[i]*nResArray1[idx1]-sub > max)
    {
        max = f[i].a+b[i].gb-gpoArray2[i]*nResArray1[idx1]-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
        /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-tgpeArray2[i]*nResArray1[idx1]-sub > max)
        {
            max = f[i].gb+b[i].gb-tgpeArray2[i]*nResArray1[idx1]-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
        /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-gpeArray2[i]*nResArray1[idx1]-sub > max)
        {
            max = f[i].gb+b[i].gb-gpeArray2[i]*nResArray1[idx1]-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[i]*sip))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-tgpeArray2[i]*sip))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,nResArray1,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,  sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,  sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2,gpeArray2,tgpeArray2,hm,hirsch_path,nResArray1, sip);


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,  sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,  sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,  sip);
            break;
    }
    return hirsch_path;
}/*}}}*/
int* hirsch_align_two_ps_vector_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{/* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/ 
    struct states* f = hm->f;
    struct states* b = hm->b;
    register int i,j,c;
    int transition = -1;

    register int idx1 = 0; /*index for nResArray1*/
    
    const float open = gpo * sip;
    const float ext = gpe * sip;
    const float text = tgpe * sip;

    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_im1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_im1 = 0.0;

    //prof1 += 64*(mid_a+1)
    prof1+= ((old_cor[4]+1)<<6);
    idx1 += (old_cor[4]+1);

    wf1_i = ((double)sip/(double)nResArray1[idx1])*(1.0/weightArray1[idx1]);
    wf1_im1 = ((double)sip/(double)nResArray1[idx1-1]*(1.0/weightArray1[idx1-1]));
    wf2_i = nResArray1[idx1]*weightArray1[idx1];
    wf2_im1 = nResArray1[idx1-1]*weightArray1[idx1-1];

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO|| fabs(prof1[-37]*wf1_im1 - (-open))>FLOAT_ZERO|| fabs(prof1[27] - (-gpoArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29]*wf1_i - (-text))>FLOAT_ZERO||fabs( prof1[29] - (-tgpeArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        /*if(f[i].a+b[i].ga-open-sub > max)*/
        /*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
        if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)
        {
            max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
        /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
        if(f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub > max)
        {
            max = f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a-open-sub > max)*/
        /*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
        if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)
        {
            max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
            /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub > max)
            {
                max = f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
            /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub > max)
            {
                max = f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        /*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
        /*if(f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub > max)*/
        if(f[i].gb+b[i].a-gpoArray2[i-1]*wf2_im1-sub > max)
        {
            max = f[i].gb+b[i].a-gpoArray2[i-1]*wf2_im1-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }

    wf2_i = nResArray1[idx1]*weightArray1[idx1];

    //i = hm->endb;
    i = old_cor[3];
    sub = abs(middle -i);
    sub /= 1000; 
    /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
    /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
    if(f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub > max)
    {
        max = f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
        /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub > max)
        {
            max = f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
        /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub > max)
        {
            max = f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[i]*wf2_i))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-tgpeArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,weightArray1, sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2,gpeArray2,tgpeArray2,hm,hirsch_path,nResArray1,weightArray1, sip);


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,weightArray1,  sip);
            break;
    }
    return hirsch_path;
}
/*}}}*/
int* hirsch_align_two_ps_vector_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{/* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

wf1 = 1
wf2 = sip
*/ 
    struct states* f = hm->f;
    struct states* b = hm->b;
    register int i,j,c;
    int transition = -1;

    register int idx1 = 0; /*index for nResArray1*/
    
    const float open = gpo * sip;
    const float ext = gpe * sip;
    const float text = tgpe * sip;

    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_im1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_im1 = 0.0;

    //prof1 += 64*(mid_a+1)
    prof1+= ((old_cor[4]+1)<<6);
    idx1 += (old_cor[4]+1);

    wf1_i = 1;
    wf1_im1 = 1;
    wf2_i = sip;
    wf2_im1 = sip;

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO|| fabs(prof1[-37]*wf1_im1 - (-open))>FLOAT_ZERO|| fabs(prof1[27] - (-gpoArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29]*wf1_i - (-text))>FLOAT_ZERO||fabs( prof1[29] - (-tgpeArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        /*if(f[i].a+b[i].ga-open-sub > max)*/
        /*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
        if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)
        {
            max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
        /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
        if(f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub > max)
        {
            max = f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a-open-sub > max)*/
        /*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
        if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)
        {
            max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
            /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub > max)
            {
                max = f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
            /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub > max)
            {
                max = f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        /*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
        /*if(f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub > max)*/
        if(f[i].gb+b[i].a-gpoArray2[i-1]*wf2_im1-sub > max)
        {
            max = f[i].gb+b[i].a-gpoArray2[i-1]*wf2_im1-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }

    wf2_i = sip;

    //i = hm->endb;
    i = old_cor[3];
    sub = abs(middle -i);
    sub /= 1000; 
    /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
    /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
    if(f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub > max)
    {
        max = f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
        /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub > max)
        {
            max = f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
        /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub > max)
        {
            max = f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[i]*wf2_i))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-tgpeArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,weightArray1, sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2,gpeArray2,tgpeArray2,hm,hirsch_path,nResArray1,weightArray1, sip);


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test1(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,weightArray1,  sip);
            break;
    }
    return hirsch_path;
}
/*}}}*/
int* hirsch_align_two_ps_vector_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{/* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

wf1 = 1
wf2 = 1
*/ 
    struct states* f = hm->f;
    struct states* b = hm->b;
    register int i,j,c;
    int transition = -1;

    register int idx1 = 0; /*index for nResArray1*/
    
    const float open = gpo * sip;
    const float ext = gpe * sip;
    const float text = tgpe * sip;

    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_im1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_im1 = 0.0;

    //prof1 += 64*(mid_a+1)
    prof1+= ((old_cor[4]+1)<<6);
    idx1 += (old_cor[4]+1);

    wf1_i = 1;
    wf1_im1 = 1;
    wf2_i = 1;
    wf2_im1 = 1;

    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++)
    for(i = old_cor[2]; i < old_cor[3];i++)
    {
        sub = abs(middle -i);
        sub /= 1000; 
        if(f[i].a+b[i].a-sub> max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO|| fabs(prof1[-37]*wf1_im1 - (-open))>FLOAT_ZERO|| fabs(prof1[27] - (-gpoArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29]*wf1_i - (-text))>FLOAT_ZERO||fabs( prof1[29] - (-tgpeArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        /*if(f[i].a+b[i].ga-open-sub > max)*/
        /*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
        if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)
        {
            max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
        /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
        if(f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub > max)
        {
            max = f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a-open-sub > max)*/
        /*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
        if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)
        {
            max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
            /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub > max)
            {
                max = f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
            /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
            if(f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub > max)
            {
                max = f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        /*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
        /*if(f[i].gb+b[i].a-gpoArray2[i-1]*sip-sub > max)*/
        if(f[i].gb+b[i].a-gpoArray2[i-1]*wf2_im1-sub > max)
        {
            max = f[i].gb+b[i].a-gpoArray2[i-1]*wf2_im1-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }

    wf2_i = sip;

    //i = hm->endb;
    i = old_cor[3];
    sub = abs(middle -i);
    sub /= 1000; 
    /*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
    /*if(f[i].a+b[i].gb-gpoArray2[i]*sip-sub > max)*/
    if(f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub > max)
    {
        max = f[i].a+b[i].gb-gpoArray2[i]*wf2_i-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        /*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
        /*if(f[i].gb+b[i].gb-tgpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub > max)
        {
            max = f[i].gb+b[i].gb-tgpeArray2[i]*wf2_i-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        /*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
        /*if(f[i].gb+b[i].gb-gpeArray2[i]*sip-sub > max)*/
        if(f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub > max)
        {
            max = f[i].gb+b[i].gb-gpeArray2[i]*wf2_i-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n",__LINE__, old_cor[4]+1, i , prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[i]*wf2_i))>FLOAT_ZERO ){
            fprintf(stderr, "line:%d, %s: prof1 and ext not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof1[29] - (-tgpeArray2[i]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr, "line:%d, %s: prof1 and text not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

    prof1-= ((old_cor[4]+1)<<6);

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2, hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,weightArray1, sip);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1,sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2,gpeArray2,tgpeArray2,hm,hirsch_path,nResArray1,weightArray1, sip);


            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1, weightArray1, sip);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_ps_dyn_new2_test2(prof1,seq2,gpoArray2, gpeArray2, tgpeArray2,hm,hirsch_path,nResArray1,weightArray1,  sip);
            break;
    }
    return hirsch_path;
}
/*}}}*/

struct states* foward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)/*{{{*/
{
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;



    prof1 += (hm->starta)<< 6;
    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	
            s[j].gb = -FLOATINFTY;
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
            s[j].gb = -FLOATINFTY;
        }

    }


    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;


        if(hm->startb){
            s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
        }else{
            s[hm->startb].gb = MAX(pgb,pa)+prof1[29];
        }
        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            pa = MAX3(pa,pga -open,pgb + prof1[-37]);

            pa += prof1[32 + seq2[j]];


            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            s[j].ga = MAX(xga-ext,xa-open);


            pgb = s[j].gb;

            s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

        }
        ca = s[j].a;

        pa = MAX3(pa,pga -open,pgb + prof1[-37]);

        pa += prof1[32 + seq2[j]];


        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
        }else{
            s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
        }

    }
    prof1 -= hm->enda << 6;
    return s;
}/*}}}*/
struct states* foward_hirsch_ps_dyn_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,int sip)/*{{{*/
{   /*return to the original ga gb calculation*/
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;
    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->starta)<< 6;
    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);	*/
            s[j].ga = MAX(s[j-1].ga- gpeArray2[j]*sip,s[j-1].a -gpoArray2[j]*sip);	
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f,  (open) =%f\n", hm->starta, j, gpoArray2[j]*sip, open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;*/
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) +prof1[29];*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) -tgpeArray2[j]*sip;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f,  (open) =%f\n", hm->starta, j, gpoArray2[j]*sip, open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        }
    }

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;


        if(hm->startb){
            s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
            /*s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*sip,pa-gpoArray2[hm->startb]*sip);*/
        }else{
            s[hm->startb].gb = MAX(pgb,pa)+prof1[29];
            /*s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*sip;*/
        }
        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb - gpoArray2[j-1]*sip);*/
            pa = MAX3(pa,pga -gpoArray2[j-1]*sip,pgb + prof1[-37]);

            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga-gpeArray2[j]*sip,xa-gpoArray2[j]*sip);

            pgb = s[j].gb;
            s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
            /*s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);*/

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open|| gpoArray2[j-1] *sip!= (open)  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f, gpoArray2[j-1]= %f, (open) =%f\n", i, j, gpoArray2[j]*sip, gpoArray2[j-1]*sip, open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga +prof1[-37] ,pgb -gpoArray2[j-1]*sip);*/
        pa = MAX3(pa,pga -gpoArray2[j-1]*sip,pgb + prof1[-37]);

        pa += prof1[32 + seq2[j]];


        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);*/
        }else{
            s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
        }
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open|| gpoArray2[j-1]*sip != (open)  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[i] =%f, gpoArray2[i-1]= %f, (open) =%f\n", i, j, gpoArray2[j]*sip, gpoArray2[j-1]*sip, open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
    }
    prof1 -= hm->enda << 6;
    return s;
}/*}}}*/
struct states* foward_hirsch_ps_dyn_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,int sip)/*{{{*/
{
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;
    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->starta)<< 6;
    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	*/
            s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);	
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27] - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28] - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29] - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) +prof1[29];
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27] - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28] - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29] - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;


        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*sip,pa-gpoArray2[hm->startb]*sip);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+prof1[29];*/
            s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*sip;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[hm->startb]*sip)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[hm->startb]*sip));
        }
        if (fabs(prof1[28] - (-gpeArray2[hm->startb]*sip)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[hm->startb]*sip)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
            pa = MAX3(pa,pga + prof1[-37],pgb - gpoArray2[j-1]*sip);

            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27] - (-open) )>FLOAT_ZERO|| fabs(prof1[-37] - (-gpoArray2[j-1]*sip))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-gpoArray2[j]*sip) =%f\n", __LINE__, i, j, prof1[27], prof1[-37], (-gpoArray2[j]*sip));
            }
            if (fabs(prof1[28] - (-gpeArray2[j]*sip))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-tgpeArray2[j]*sip))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
        pa = MAX3(pa,pga +prof1[-37] ,pgb -gpoArray2[j-1]*sip);

        pa += prof1[32 + seq2[j]];


        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-open)) > FLOAT_ZERO || fabs(prof1[-37] - (-gpoArray2[j-1]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[j]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }
    prof1 -= hm->enda << 6;
    return s;
}/*}}}*/
struct states* foward_hirsch_ps_dyn_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,unsigned int* nResArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/ 
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    register int idx1 = 0; /*index for nResArray1*/

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;
    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->starta)<< 6;
    idx1 += hm->starta;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;

    float wfi = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wfim1 = 0.0; /*weight for i-1*/
    float wfip1 = 0.0; /*weight for i+1*/

    wfi = ((double)sip/(double)nResArray1[idx1]);
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);	*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wfi,s[j-1].a+prof1[27]*wfi);	
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wfi - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wfi - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29]*wfi - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) +prof1[29]*wfi;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wfi - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wfi - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29]*wfi - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
        idx1 += 1;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;


        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*sip,pa-gpoArray2[hm->startb]*sip);*/
            s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*nResArray1[idx1],pa-gpoArray2[hm->startb]*nResArray1[idx1]);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*sip;*/
            s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*nResArray1[idx1];
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[hm->startb]*nResArray1[idx1])) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[hm->startb]*sip));
        }
        if (fabs(prof1[28] - (-gpeArray2[hm->startb]*nResArray1[idx1])) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[hm->startb]*nResArray1[idx1])) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
        wfi = ((double)sip/(double)nResArray1[idx1]);
        wfim1 = ((double)sip/(double)nResArray1[idx1-1]);
        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb - gpoArray2[j-1]*sip);*/
            pa = MAX3(pa,pga + prof1[-37]*wfim1,pgb - gpoArray2[j-1]*nResArray1[idx1-1]);

            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wfi,xa+prof1[27]*wfi);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*nResArray1[idx1],ca-gpoArray2[j]*nResArray1[idx1]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wfi - (-open) )>FLOAT_ZERO|| fabs(prof1[-37] - (-gpoArray2[j-1]*nResArray1[idx1-1]))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-gpoArray2[j]*sip) =%f\n", __LINE__, i, j, prof1[27], prof1[-37], (-gpoArray2[j]*sip));
            }
            if (fabs(prof1[28] - (-gpeArray2[j]*nResArray1[idx1]))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-tgpeArray2[j]*nResArray1[idx1]))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga +prof1[-37] ,pgb -gpoArray2[j-1]*sip);*/
        pa = MAX3(pa,pga +prof1[-37]*wfim1 ,pgb -gpoArray2[j-1]*nResArray1[idx1-1]);

        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*nResArray1[idx1] ,ca-gpoArray2[j]*nResArray1[idx1]);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*nResArray1[idx1];
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wfi - (-open)) > FLOAT_ZERO || fabs(prof1[-37] - (-gpoArray2[j-1]*nResArray1[idx1-1]))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[j]*nResArray1[idx1]))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*nResArray1[idx1]))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }
    prof1 -= hm->enda << 6;
	idx1 -= hm->enda;
    return s;
}/*}}}*/
struct states* foward_hirsch_ps_dyn_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/ 
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    register int idx1 = 0; /*index for nResArray1*/

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;
    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->starta)<< 6;
    idx1 += hm->starta;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_im1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_im1 = 0.0;

    wf1_i = ((double)sip/(double)nResArray1[idx1])*(1.0/weightArray1[idx1]);
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);	*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);	
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wf1_i - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) +prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[29]*wf1_i - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
        idx1 += 1;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

        wf2_i = nResArray1[idx1]*weightArray1[idx1];
        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*sip,pa-gpoArray2[hm->startb]*sip);*/
            s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*wf2_i,pa-gpoArray2[hm->startb]*wf2_i);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*sip;*/
            s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[hm->startb]*sip));
        }
        if (fabs(prof1[28] - (-gpeArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/

        wf1_i = ((double)sip/(double)nResArray1[idx1])*(1.0/weightArray1[idx1]);
        wf1_im1 = ((double)sip/(double)nResArray1[idx1-1]*(1.0/weightArray1[idx1-1]));
        wf2_i = nResArray1[idx1]*weightArray1[idx1];
        wf2_im1 = nResArray1[idx1-1]*weightArray1[idx1-1];

        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb - gpoArray2[j-1]*sip);*/
            pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb - gpoArray2[j-1]*wf2_im1);

            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,   xa+prof1[27]*wf1_i);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*wf2_i,ca-gpoArray2[j]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef CHECK_PS
            fprintf(stderr,"prof[27](%g)/nRes(%d)=%g, weightArray=%g\n",prof1[27],nResArray1[idx1], prof1[27]/nResArray1[idx1], weightArray1[idx1]);
#endif
#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wf1_i - (-open) )>FLOAT_ZERO|| fabs(prof1[-37] - (-gpoArray2[j-1]*wf2_im1))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-gpoArray2[j]*sip) =%f\n", __LINE__, i, j, prof1[27], prof1[-37], (-gpoArray2[j]*sip));
            }
            if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        wf1_im1 = ((double)sip/(double)nResArray1[idx1-1]*(1.0/weightArray1[idx1-1]));
        wf2_i = nResArray1[idx1]*weightArray1[idx1];
        wf2_im1 = nResArray1[idx1-1]*weightArray1[idx1-1];

        /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga +prof1[-37] ,pgb -gpoArray2[j-1]*sip);*/
        pa = MAX3(pa,pga +prof1[-37]*wf1_im1 ,pgb -gpoArray2[j-1]*wf2_im1);

        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*wf2_i ,ca-gpoArray2[j]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO || fabs(prof1[-37] - (-gpoArray2[j-1]*wf2_im1))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }
    prof1 -= hm->enda << 6;
	idx1 -= hm->enda;
    return s;
}/*}}}*/
struct states* foward_hirsch_ps_dyn_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

wf1 = 1
wf2 = sip
*/ 
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    register int idx1 = 0; /*index for nResArray1*/

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;
    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->starta)<< 6;
    idx1 += hm->starta;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_im1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_im1 = 0.0;

    wf1_i = 1;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);	*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);	
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wf1_i - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) +prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[29]*wf1_i - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
        idx1 += 1;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

        wf2_i = sip;
        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*sip,pa-gpoArray2[hm->startb]*sip);*/
            s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*wf2_i,pa-gpoArray2[hm->startb]*wf2_i);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*sip;*/
            s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[hm->startb]*sip));
        }
        if (fabs(prof1[28] - (-gpeArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/

        wf1_i = 1;
        wf1_im1 = 1;
        wf2_i = sip;
        wf2_im1 = sip;

        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb - gpoArray2[j-1]*sip);*/
            pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb - gpoArray2[j-1]*wf2_im1);

            pa += prof1[32 + seq2[j]];
#ifdef CHECK_PA
            fprintf(stderr,"forward_ps: j=%d, seq2[j]=%d, pa += prof1[32 + seq2[j]] =%g, pa/sip(%d)=%g\n", j, seq2[j], prof1[32+seq2[j]], sip,prof1[32+seq2[j]]/sip );
#endif 

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,   xa+prof1[27]*wf1_i);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*wf2_i,ca-gpoArray2[j]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef CHECK_PS
            fprintf(stderr,"prof[27](%g)/nRes(%d)=%g, weightArray=%g\n",prof1[27],nResArray1[idx1], prof1[27]/nResArray1[idx1], weightArray1[idx1]);
#endif
#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wf1_i - (-open) )>FLOAT_ZERO|| fabs(prof1[-37] - (-gpoArray2[j-1]*wf2_im1))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-gpoArray2[j]*sip) =%f\n", __LINE__, i, j, prof1[27], prof1[-37], (-gpoArray2[j]*sip));
            }
            if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        wf1_im1 = 1;
        wf2_i = sip;
        wf2_im1 = sip;

        /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga +prof1[-37] ,pgb -gpoArray2[j-1]*sip);*/
        pa = MAX3(pa,pga +prof1[-37]*wf1_im1 ,pgb -gpoArray2[j-1]*wf2_im1);

        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*wf2_i ,ca-gpoArray2[j]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO || fabs(prof1[-37] - (-gpoArray2[j-1]*wf2_im1))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }
    prof1 -= hm->enda << 6;
	idx1 -= hm->enda;
    return s;
}/*}}}*/
struct states* foward_hirsch_ps_dyn_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2, struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

wf1 = 1
wf2 = sip
*/ 
    struct states* s = hm->f;

    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    register int idx1 = 0; /*index for nResArray1*/

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;
    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->starta)<< 6;
    idx1 += hm->starta;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_im1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_im1 = 0.0;

    wf1_i = 1;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);	*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);	
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wf1_i - (-ext)) > FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a) +prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[29]*wf1_i - (-text))> FLOAT_ZERO){
        fprintf(stderr,"line:%d, i=%d, j=%d, ps, tgpe not equal\n",__LINE__,  i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;
    seq2--;

    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
        idx1 += 1;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

        wf2_i = 1;
        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*sip,pa-gpoArray2[hm->startb]*sip);*/
            s[hm->startb].gb = MAX(pgb-gpeArray2[hm->startb]*wf2_i,pa-gpoArray2[hm->startb]*wf2_i);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*sip;*/
            s[hm->startb].gb = MAX(pgb,pa)-tgpeArray2[hm->startb]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[hm->startb]*sip));
        }
        if (fabs(prof1[28] - (-gpeArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[hm->startb]*wf2_i)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/

        wf1_i = 1;
        wf1_im1 = 1;
        wf2_i = 1;
        wf2_im1 = 1;

        for (j = hm->startb+1; j < hm->endb;j++){
            ca = s[j].a;

            /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb - gpoArray2[j-1]*sip);*/
            pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb - gpoArray2[j-1]*wf2_im1);

            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,   xa+prof1[27]*wf1_i);

            pgb = s[j].gb;
            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*wf2_i,ca-gpoArray2[j]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef CHECK_PS
            fprintf(stderr,"prof[27](%g)/nRes(%d)=%g, weightArray=%g\n",prof1[27],nResArray1[idx1], prof1[27]/nResArray1[idx1], weightArray1[idx1]);
#endif
#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wf1_i - (-open) )>FLOAT_ZERO|| fabs(prof1[-37] - (-gpoArray2[j-1]*wf2_im1))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-gpoArray2[j]*sip) =%f\n", __LINE__, i, j, prof1[27], prof1[-37], (-gpoArray2[j]*sip));
            }
            if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        wf1_im1 = 1;
        wf2_i = 1;
        wf2_im1 = 1;

        /*pa = MAX3(pa,pga -open,pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga +prof1[-37] ,pgb -gpoArray2[j-1]*sip);*/
        pa = MAX3(pa,pga +prof1[-37]*wf1_im1 ,pgb -gpoArray2[j-1]*wf2_im1);

        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip ,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*wf2_i ,ca-gpoArray2[j]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]*wf1_i - (-open)) > FLOAT_ZERO || fabs(prof1[-37] - (-gpoArray2[j-1]*wf2_im1))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[-37]=%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[-37], (-open));
        }
        if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }
    prof1 -= hm->enda << 6;
	idx1 -= hm->enda;
    return s;
}/*}}}*/

struct states* backward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)/*{{{*/
{
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    /*const float text = tgpe * sip;*/

    /*const int starta = hm->starta;*/
    /*const int enda = hm->enda;*/
    /*const int startb =hm->startb;*/
    /*const int endb = hm->endb;*/

    prof1 += (hm->enda+1) << 6;

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            s[j].gb = -FLOATINFTY;
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29];
            s[j].gb = -FLOATINFTY;
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;


        if(hm->endb != hm->len_b){
            s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
        }else{
            s[hm->endb].gb = MAX(pgb,pa) +prof1[29];
        }

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            pa = MAX3(pa,pga - open,pgb +prof1[91]);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;
            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            s[j].ga = MAX(xga-ext,xa-open);

            pgb = s[j].gb;
            s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;


        }
        ca = s[j].a;

        pa = MAX3(pa,pga - open,pgb +prof1[91]);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;


        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){
            s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
        }else{
            s[j].gb = MAX(s[j].gb,ca)+prof1[29];
        }

    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ps_dyn_new0(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int sip)/*{{{*/
{ /*return to the original ga gb calculation*/ 
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;

    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    prof1 += (hm->enda+1) << 6;

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga-gpeArray2[j]*sip,s[j+1].a-gpoArray2[j]*sip);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f, gpoArray2[j-1]= %f, (open) =%f\n", starta, j , gpoArray2[j]*sip, gpoArray2[j-1*sip], open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;*/
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]; */
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpeArray2[j]*sip;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f, gpoArray2[j-1]= %f, (open) =%f\n", starta, j , gpoArray2[j]*sip, gpoArray2[j-1]*sip, open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;


        if(hm->endb != hm->len_b){
            s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
            /*s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*sip,pa-gpoArray2[endb-1]*sip);*/
        }else{
            s[hm->endb].gb = MAX(pgb,pa) +prof1[29];
            /*s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*sip;*/
        }

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb -gpoArray2[j+1]*sip);*/
            pa = MAX3(pa,pga - gpoArray2[j+1]*sip,pgb +prof1[91]);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga-gpeArray2[j]*sip,xa-gpoArray2[j]*sip);

            pgb = s[j].gb;

            s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
            /*s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);*/

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open|| gpoArray2[j-1]*sip != (open)  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f, gpoArray2[j-1]= %f, (open) =%f\n", i, j, gpoArray2[j]*sip, gpoArray2[j-1]*sip, open);
    }
    if (gpeArray2[j]*sip != (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j] *sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
        /*pa = MAX3(pa,pga +prof1[91],pgb -gpoArray2[j+1]*sip);*/
        pa = MAX3(pa,pga - gpoArray2[j+1]*sip,pgb +prof1[91]);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){
            s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip, ca-gpeArray2[j]*sip);*/
        }else{
            s[j].gb = MAX(s[j].gb,ca)+prof1[29];
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
        }
#ifdef DEBUG_PS/*{{{*/
    if (gpoArray2[j]*sip != open|| gpoArray2[j+1]*sip != (open)  ){
        fprintf(stderr,"i=%d, j=%d, ps, gpo not equal, gpoArray2[j] =%f, gpoArray2[j+1]= %f, (open) =%f\n", i, j, gpoArray2[j]*sip, gpoArray2[j+1]*sip, open);
    }
    if (gpeArray2[j] *sip!= (ext) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
    if (tgpeArray2[j]*sip != (text) ){
        fprintf(stderr, "%s: gpe not equal\n", __FUNCTION__);
    }
#endif/*}}}*/

    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ps_dyn_new(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,int sip)/*{{{*/
{
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;

    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    prof1 += (hm->enda+1) << 6;

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27] - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28] - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;*/
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]; 
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27] - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28] - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;


        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*sip,pa-gpoArray2[endb-1]*sip);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa) +prof1[29];*/
            s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*sip;
        }

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
            pa = MAX3(pa,pga + prof1[91],pgb -gpoArray2[j+1]*sip);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]- (-open) )>FLOAT_ZERO|| fabs(prof1[91] - (-gpoArray2[j+1]*sip))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open));
            }
            if (fabs(prof1[28] - (-ext))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
        pa = MAX3(pa,pga +prof1[91],pgb -gpoArray2[j+1]*sip);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip, ca-gpeArray2[j]*sip);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;
        }

#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-open)) >FLOAT_ZERO || fabs(prof1[91] - (-gpoArray2[j+1]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[91] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[91], (-open));
        }
        if (fabs(prof1[28] - (-ext))> FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ps_dyn_new1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/ 
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;

    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    register int idx1 = 0; /*index for nResArray1*/

    prof1 += (hm->enda+1) << 6;
    idx1 += (hm->enda+1);

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    float wfi = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wfim1 = 0.0; /*weight for i-1*/
    float wfip1 = 0.0; /*weight for i+1*/

    wfi = ((double)sip/(double)nResArray1[idx1]);
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wfi,s[j+1].a+prof1[27]*wfi);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27] - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28] - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;*/
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]; */
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wfi; 
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27] - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28] - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
    if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){

        prof1 -= 64;
        idx1 -= 1;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;


        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*sip,pa-gpoArray2[endb-1]*sip);*/
            s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*nResArray1[idx1],pa-gpoArray2[endb-1]*nResArray1[idx1]);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa) +prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*sip;*/
            s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*nResArray1[idx1];
        }

        wfi = ((double)sip/(double)nResArray1[idx1]);
        wfip1 = ((double)sip/(double)nResArray1[idx1+1]);

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb -gpoArray2[j+1]*sip);*/
            pa = MAX3(pa,pga + prof1[91]*wfip1,pgb -gpoArray2[j+1]*nResArray1[idx1+1]);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wfi,xa+prof1[27]*wfi);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*nResArray1[idx1],ca-gpoArray2[j]*nResArray1[idx1]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]- (-open) )>FLOAT_ZERO|| fabs(prof1[91] - (-gpoArray2[j+1]*sip))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open));
            }
            if (fabs(prof1[28] - (-ext))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, j);
            }
            if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, j);
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
        /*pa = MAX3(pa,pga +prof1[91],pgb -gpoArray2[j+1]*sip);*/
        pa = MAX3(pa,pga +prof1[91]*wfip1,pgb -gpoArray2[j+1]*nResArray1[idx1+1]);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip, ca-gpeArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*nResArray1[idx1], ca-gpeArray2[j]*nResArray1[idx1]);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*nResArray1[idx1];
        }

#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-open)) >FLOAT_ZERO || fabs(prof1[91] - (-gpoArray2[j+1]*sip))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof[91] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], prof1[91], (-open));
        }
        if (fabs(prof1[28] - (-ext))> FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-text))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ps_dyn_new2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence
*/ 
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;

    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    register int idx1 = 0; /*index for nResArray1*/

    prof1 += (hm->enda+1) << 6;
    idx1 += (hm->enda+1);

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_ip1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_ip1 = 0.0;

    wf1_i = ((double)sip/(double)nResArray1[idx1])*(1.0/weightArray1[idx1]);

    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wf1_i - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;*/
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]; */
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[29] *wf1_i- (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){

        prof1 -= 64;
        idx1 -= 1;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

        wf2_i = nResArray1[idx1]*weightArray1[idx1];
        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*sip,pa-gpoArray2[endb-1]*sip);*/
            s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*wf2_i,pa-gpoArray2[endb-1]*wf2_i);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa) +prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*sip;*/
            s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-gpoArray2[endb-a]) =%f, wf2_i=%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[endb-1]), wf2_i);
        }
        if (fabs(prof1[28] - (-gpeArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/

        wf1_i = ((double)sip/(double)nResArray1[idx1])*(1.0/weightArray1[idx1]);
        wf1_ip1 = ((double)sip/(double)nResArray1[idx1+1]*(1.0/weightArray1[idx1+1]));
        wf2_i = nResArray1[idx1]*weightArray1[idx1];
        wf2_ip1 = nResArray1[idx1+1]*weightArray1[idx1+1];

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb -gpoArray2[j+1]*sip);*/
            pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb -gpoArray2[j+1]*wf2_ip1);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*wf2_i,ca-gpoArray2[j]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wf1_i- (-open) )>FLOAT_ZERO|| fabs(prof1[27]-(-gpoArray2[j]*wf2_i)) >FLOAT_ZERO || fabs(prof1[91]*wf1_ip1 - (-open))>FLOAT_ZERO || fabs(prof1[91]-(-gpoArray2[j+1]*wf2_ip1)) > FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open), wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO || fabs(prof1[28]- (-gpeArray2[j]*wf2_i))> FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal, prof1[28]=%f, gpeArray2[j]=%f, wf1_i=%f,wf2_i=%f\n", __LINE__, i, j, prof1[28], (-gpeArray2[j], wf1_i, wf2_i));
            }
            if (fabs(prof1[29]*wf1_i - (-text))>FLOAT_ZERO || fabs(prof1[29]- (-tgpeArray2[j]*wf2_i))> FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal, prof1[29]=%f, tgpeArray2[j]=%f, wf1_i=%f,wf2_i=%f\n", __LINE__, i, j, prof1[29], (-tgpeArray2[j], wf1_i, wf2_i));
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        wf1_ip1 = ((double)sip/(double)nResArray1[idx1+1]*(1.0/weightArray1[idx1+1]));
        wf2_i = nResArray1[idx1]*weightArray1[idx1];
        wf2_ip1 = nResArray1[idx1+1]*weightArray1[idx1+1];

        /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
        /*pa = MAX3(pa,pga +prof1[91],pgb -gpoArray2[j+1]*sip);*/
        pa = MAX3(pa,pga +prof1[91]*wf1_ip1,pgb -gpoArray2[j+1]*wf2_ip1);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){ /*the difference between the kalingP and kaign is caused by this code as of 2010-11-03 14:45:34 Wednesday Week 44, fixed now*/
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip, ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*wf2_i, ca-gpoArray2[j]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]-(-gpoArray2[j]*wf2_i)) >FLOAT_ZERO || fabs(prof1[91]*wf1_ip1 - (-open))>FLOAT_ZERO || fabs(prof1[91]-(-gpoArray2[j+1]*wf2_ip1)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open), wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
        if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))> FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ps_dyn_new2_test1(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

2011-02-24 
set 
wf1 = 1;
wf2 = sip;
*/ 
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;

    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    register int idx1 = 0; /*index for nResArray1*/

    prof1 += (hm->enda+1) << 6;
    idx1 += (hm->enda+1);

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_ip1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_ip1 = 0.0;

    wf1_i = 1;

    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wf1_i - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;*/
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]; */
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[29] *wf1_i- (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){

        prof1 -= 64;
        idx1 -= 1;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

        wf2_i = sip;
        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*sip,pa-gpoArray2[endb-1]*sip);*/
            s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*wf2_i,pa-gpoArray2[endb-1]*wf2_i);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa) +prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*sip;*/
            s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-gpoArray2[endb-a]) =%f, wf2_i=%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[endb-1]), wf2_i);
        }
        if (fabs(prof1[28] - (-gpeArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/

        wf1_i = 1;
        wf1_ip1 = 1;
        wf2_i = sip;
        wf2_ip1 = sip;

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb -gpoArray2[j+1]*sip);*/
            pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb -gpoArray2[j+1]*wf2_ip1);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*wf2_i,ca-gpoArray2[j]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wf1_i- (-open) )>FLOAT_ZERO|| fabs(prof1[27]-(-gpoArray2[j]*wf2_i)) >FLOAT_ZERO || fabs(prof1[91]*wf1_ip1 - (-open))>FLOAT_ZERO || fabs(prof1[91]-(-gpoArray2[j+1]*wf2_ip1)) > FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open), wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO || fabs(prof1[28]- (-gpeArray2[j]*wf2_i))> FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal, prof1[28]=%f, gpeArray2[j]=%f, wf1_i=%f,wf2_i=%f\n", __LINE__, i, j, prof1[28], (-gpeArray2[j], wf1_i, wf2_i));
            }
            if (fabs(prof1[29]*wf1_i - (-text))>FLOAT_ZERO || fabs(prof1[29]- (-tgpeArray2[j]*wf2_i))> FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal, prof1[29]=%f, tgpeArray2[j]=%f, wf1_i=%f,wf2_i=%f\n", __LINE__, i, j, prof1[29], (-tgpeArray2[j], wf1_i, wf2_i));
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        wf1_ip1 = 1;
        wf2_i = sip;
        wf2_ip1 = sip;

        /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
        /*pa = MAX3(pa,pga +prof1[91],pgb -gpoArray2[j+1]*sip);*/
        pa = MAX3(pa,pga +prof1[91]*wf1_ip1,pgb -gpoArray2[j+1]*wf2_ip1);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){ /*the difference between the kalingP and kaign is caused by this code as of 2010-11-03 14:45:34 Wednesday Week 44, fixed now*/
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip, ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*wf2_i, ca-gpoArray2[j]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]-(-gpoArray2[j]*wf2_i)) >FLOAT_ZERO || fabs(prof1[91]*wf1_ip1 - (-open))>FLOAT_ZERO || fabs(prof1[91]-(-gpoArray2[j+1]*wf2_ip1)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open), wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
        if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))> FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_ps_dyn_new2_test2(const float* prof1,const int* seq2,const float* gpoArray2, const float* gpeArray2, const float* tgpeArray2,struct hirsch_mem* hm,unsigned int* nResArray1, float* weightArray1, int sip)/*{{{*/
{ /* 
gap penalties is still effective, but when default (constant) gap penalties are used, the alignment is the same as the original kalign

sip:        nsip for prof1 
nResArray1: number of residues at each position of the consensus sequence

2011-02-24 
set 
wf1 = 1;
wf2 = 1;
*/ 
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;

    const float open = gpo * sip;
    const float ext = gpe *sip;
    const float text = tgpe * sip;

    const int starta = hm->starta;
    const int enda = hm->enda;
    const int startb =hm->startb;
    const int endb = hm->endb;

    register int idx1 = 0; /*index for nResArray1*/

    prof1 += (hm->enda+1) << 6;
    idx1 += (hm->enda+1);

    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;

    float wf1_i = 0.0; /*weight factor to make the gap penalties are the same as the original kalign when constant gap penalties are used, 2010-10-20*/
    float wf1_ip1 = 0.0; /*weight for i-1*/
    float wf2_i = 0.0;
    float wf2_ip1 = 0.0;

    wf1_i = 1;

    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[27]*wf1_i - (-open))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-open) =%f\n", __LINE__, i, hm->startb, prof1[27], (-open));
    }
    if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;*/
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]; */
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PS/*{{{*/
    if (fabs(prof1[29] *wf1_i- (-text))>FLOAT_ZERO){
        fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
    }
#endif/*}}}*/
        }
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){

        prof1 -= 64;
        idx1 -= 1;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

        wf2_i = 1;
        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*sip,pa-gpoArray2[endb-1]*sip);*/
            s[hm->endb].gb = MAX(pgb-gpeArray2[endb-1]*wf2_i,pa-gpoArray2[endb-1]*wf2_i);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa) +prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*sip;*/
            s[hm->endb].gb = MAX(pgb,pa) -tgpeArray2[endb-1]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27] - (-gpoArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, (-gpoArray2[endb-a]) =%f, wf2_i=%f\n", __LINE__, i, hm->startb, prof1[27], (-gpoArray2[endb-1]), wf2_i);
        }
        if (fabs(prof1[28] - (-gpeArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[endb-1]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/

        wf1_i = 1;
        wf1_ip1 = 1;
        wf2_i = 1;
        wf2_ip1 = 1;

        for(j = hm->endb-1;j > hm->startb;j--){
            ca = s[j].a;

            /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb -gpoArray2[j+1]*sip);*/
            pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb -gpoArray2[j+1]*wf2_ip1);
            pa += prof1[32 + seq2[j]];

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
            /*s[j].ga = MAX(xga-ext,xa-open);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb-gpeArray2[j]*sip,ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(pgb-gpeArray2[j]*wf2_i,ca-gpoArray2[j]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;

#ifdef DEBUG_PS/*{{{*/
            if (fabs(prof1[27]*wf1_i- (-open) )>FLOAT_ZERO|| fabs(prof1[27]-(-gpoArray2[j]*wf2_i)) >FLOAT_ZERO || fabs(prof1[91]*wf1_ip1 - (-open))>FLOAT_ZERO || fabs(prof1[91]-(-gpoArray2[j+1]*wf2_ip1)) > FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open), wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof1[28]*wf1_i - (-ext))>FLOAT_ZERO || fabs(prof1[28]- (-gpeArray2[j]*wf2_i))> FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal, prof1[28]=%f, gpeArray2[j]=%f, wf1_i=%f,wf2_i=%f\n", __LINE__, i, j, prof1[28], (-gpeArray2[j], wf1_i, wf2_i));
            }
            if (fabs(prof1[29]*wf1_i - (-text))>FLOAT_ZERO || fabs(prof1[29]- (-tgpeArray2[j]*wf2_i))> FLOAT_ZERO){
                fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal, prof1[29]=%f, tgpeArray2[j]=%f, wf1_i=%f,wf2_i=%f\n", __LINE__, i, j, prof1[29], (-tgpeArray2[j], wf1_i, wf2_i));
            }
#endif/*}}}*/
        }
        ca = s[j].a;

        wf1_ip1 = 1;
        wf2_i = 1;
        wf2_ip1 = 1;

        /*pa = MAX3(pa,pga - open,pgb +prof1[91]);*/
        /*pa = MAX3(pa,pga +prof1[91],pgb -gpoArray2[j+1]*sip);*/
        pa = MAX3(pa,pga +prof1[91]*wf1_ip1,pgb -gpoArray2[j+1]*wf2_ip1);
        pa += prof1[32 + seq2[j]];

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
        if(hm->startb){ /*the difference between the kalingP and kaign is caused by this code as of 2010-11-03 14:45:34 Wednesday Week 44, fixed now*/
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb-gpeArray2[j]*sip, ca-gpoArray2[j]*sip);*/
            s[j].gb = MAX(s[j].gb-gpeArray2[j]*wf2_i, ca-gpoArray2[j]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*sip;*/
            s[j].gb = MAX(s[j].gb,ca)-tgpeArray2[j]*wf2_i;
        }
#ifdef DEBUG_PS/*{{{*/
        if (fabs(prof1[27]-(-gpoArray2[j]*wf2_i)) >FLOAT_ZERO || fabs(prof1[91]*wf1_ip1 - (-open))>FLOAT_ZERO || fabs(prof1[91]-(-gpoArray2[j+1]*wf2_ip1)) > FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpo not equal, prof1[27] =%f, prof1[91] =%f, (-open) =%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n", __LINE__, i, j, prof1[27], prof1[91], (-open), wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
        if (fabs(prof1[28] - (-gpeArray2[j]*wf2_i))> FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, gpe not equal\n", __LINE__, i, hm->startb);
        }
        if (fabs(prof1[29] - (-tgpeArray2[j]*wf2_i))>FLOAT_ZERO){
            fprintf(stderr,"line=%d, i=%d, j=%d, ps, tgpe not equal\n", __LINE__, i, hm->startb);
        }
#endif/*}}}*/
    }		
    return s;
}/*}}}*/

/*==== consensus - consensus alignment ==== */
int* hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


    //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);


    if(hm->starta  >= hm->enda){
        return hirsch_path;
    }
    if(hm->startb  >= hm->endb){
        return hirsch_path;
    }

    hm->enda = mid;
    hm->f = foward_hirsch_pp_dyn(prof1,prof2,hm);
    /*int i;
      fprintf(stderr,"FOWARD\n");
      for (i = hm->startb; i <= hm->endb;i++){
      fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
      }*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    hm->b = backward_hirsch_pp_dyn(prof1,prof2,hm);
    /*fprintf(stderr,"BaCKWARD\n");

      for (i = hm->startb; i <= hm->endb;i++){
      fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
      }*/

    hirsch_path = hirsch_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
    return hirsch_path;
}/*}}}*/
int* hirsch_pp_dyn_new(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};
    //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);

    if(hm->starta  >= hm->enda){
        return hirsch_path;
    }
    if(hm->startb  >= hm->endb){
        return hirsch_path;
    }

    hm->enda = mid;
    hm->f = foward_hirsch_pp_dyn_new(prof1,prof2,hm);
#ifdef DEBUG
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif

    hm->starta = mid;
    hm->enda = old_cor[1];
    hm->b = backward_hirsch_pp_dyn_new(prof1,prof2,hm);

#ifdef DEBUG
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif

    hirsch_path = hirsch_align_two_pp_vector_new(prof1,prof2,hm,hirsch_path,input_states,old_cor);
    return hirsch_path;
}/*}}}*/
int* hirsch_pp_dyn_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};
    //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_pp_dyn_new(prof1,prof2,hm);*/
    hm->f = foward_hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2, hm);
    /*hm->f = foward_hirsch_pp_dyn_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2, hm);*/
#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_pp_dyn_new(prof1,prof2,hm);*/
    hm->b = backward_hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm);
    /*hm->b = backward_hirsch_pp_dyn_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm);*/

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_pp_vector_new(prof1,prof2,hm,hirsch_path,input_states,old_cor);*/
    hirsch_path = hirsch_align_two_pp_vector_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);
    /*hirsch_path = hirsch_align_two_pp_vector_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);*/
    return hirsch_path;
}/*}}}*/
int* hirsch_pp_dyn_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
    /*using wf1 and wf2 to make the output the same as the original kalign2
     * when constant gap penalty is used*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};
    //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_pp_dyn_new(prof1,prof2,hm);*/
    hm->f = foward_hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2, hm);
    /*hm->f = foward_hirsch_pp_dyn_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2, hm);*/
#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_pp_dyn_new(prof1,prof2,hm);*/
    hm->b = backward_hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm);
    /*hm->b = backward_hirsch_pp_dyn_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm);*/

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_pp_vector_new(prof1,prof2,hm,hirsch_path,input_states,old_cor);*/
    hirsch_path = hirsch_align_two_pp_vector_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);
    /*hirsch_path = hirsch_align_two_pp_vector_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);*/
    return hirsch_path;
}/*}}}*/
int* hirsch_pp_dyn_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
    /*set wf1=1 and  wf2=1*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};
    //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    /*hm->f = foward_hirsch_pp_dyn_new(prof1,prof2,hm);*/
    hm->f = foward_hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2, hm);
    /*hm->f = foward_hirsch_pp_dyn_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2, hm);*/
#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/

    hm->starta = mid;
    hm->enda = old_cor[1];
    /*hm->b = backward_hirsch_pp_dyn_new(prof1,prof2,hm);*/
    hm->b = backward_hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm);
    /*hm->b = backward_hirsch_pp_dyn_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm);*/

#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/

    /*hirsch_path = hirsch_align_two_pp_vector_new(prof1,prof2,hm,hirsch_path,input_states,old_cor);*/
    hirsch_path = hirsch_align_two_pp_vector_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);
    /*hirsch_path = hirsch_align_two_pp_vector_new2_0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);*/
    return hirsch_path;
}/*}}}*/
int* hirsch_pp_dyn_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm, int* hirsch_path)/*{{{*/
    /*set wf1=1/nSeq2 and  wf2=1/nSeq1*/
{
    int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
    float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
    int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};
    //fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);

    if(hm->starta  >= hm->enda){ return hirsch_path; }
    if(hm->startb  >= hm->endb){ return hirsch_path; }

    hm->enda = mid;
    hm->f = foward_hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2, nsip1, nsip2, weightArray1, weightArray2, hm);
#ifdef DEBUG/*{{{*/
    int i;
    fprintf(stderr,"%s: FOWARD\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->f[i].a,i, hm->f[i].ga,i, hm->f[i].gb);
    }
#endif/*}}}*/
    hm->starta = mid;
    hm->enda = old_cor[1];
    hm->b = backward_hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2, nsip1, nsip2, weightArray1, weightArray2,hm);
#ifdef DEBUG/*{{{*/
    fprintf(stderr,"%s: Backward\n", __FUNCTION__);
    for (i = hm->startb; i <= hm->endb;i++){
        fprintf(stderr,"%s: a[%d]=%f	ga[%d]=%f	gb[%d]=%f\n",__FUNCTION__, i, hm->b[i].a,i, hm->b[i].ga,i, hm->b[i].gb);
    }
#endif/*}}}*/
    hirsch_path = hirsch_align_two_pp_vector_new2_test2(prof1,prof2,nResArray1, nResArray2, nsip1, nsip2, weightArray1, weightArray2,hm,hirsch_path,input_states,old_cor);
    return hirsch_path;
}/*}}}*/

int* hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;


    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;	
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;


    prof1+= ((old_cor[4]+1) << 6);
    //prof2 += 64 * (hm->startb);
    //i = hm->startb;
    prof2 += old_cor[2] << 6;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++){
    for(i = old_cor[2]; i < old_cor[3];i++){
        sub = abs(middle -i);
        sub /= 1000; 
        prof2 += 64;
        //fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
        if(f[i].a+b[i].a-sub > max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
        if(f[i].a+b[i].ga+prof2[27]-sub > max){
            max = f[i].a+b[i].ga+prof2[27]-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        if(f[i].a+b[i].gb+prof1[27] -sub> max){
            max = f[i].a+b[i].gb+prof1[27]-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        if(f[i].ga+b[i].a+prof2[-37]-sub > max){
            max = f[i].ga+b[i].a+prof2[-37]-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            if(f[i].gb+b[i].gb+prof1[29]-sub > max){
                max = f[i].gb+b[i].gb+prof1[29]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            if(f[i].gb+b[i].gb+prof1[28]-sub > max){
                max = f[i].gb+b[i].gb+prof1[28]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        if(f[i].gb+b[i].a+prof1[-37]-sub > max){
            max = f[i].gb+b[i].a+prof1[-37]-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
    }
    //i = hm->endb;
    i = old_cor[3];
    sub = abs(middle -i);
    sub /= 1000; 
    if(f[i].a+b[i].gb+prof1[27]-sub > max){
        max = f[i].a+b[i].gb+prof1[27]-sub;
        //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
        transition = 3;
        c = i;
    }
    if(hm->endb == hm->len_b){
        if(f[i].gb+b[i].gb+prof1[29]-sub > max){
            max = f[i].gb+b[i].gb+prof1[29]-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }	
    }else{
        if(f[i].gb+b[i].gb+prof1[28]-sub > max){
            max = f[i].gb+b[i].gb+prof1[28]-sub;
            //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
            transition = 6;
            c = i;
        }
    }



    prof1-= (old_cor[4]+1)<<6;
    //prof2 -= hm->endb << 6;
    prof2 -= old_cor[3] << 6;

    //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
    //if(transition == -1){
    //	exit(0);
    //}

    j = hirsch_path[0];
    switch(transition){
        case 1: //a -> a = 1

            hirsch_path[old_cor[4]] = c;
            hirsch_path[old_cor[4]+1] = c+1;

            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;
            //fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
            break;
        case 2:// a -> ga = 2

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;


            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4];
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = 0.0;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
            break;
        case 3:// a -> gb = 3

            hirsch_path[old_cor[4]] = c;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = 0.0;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
            break;
        case 5://ga -> a = 5
            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = 0.0;
            hm->b[0].gb = -FLOATINFTY;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4];

            hm->startb = old_cor[2];
            hm->endb = c-1;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
            break;
        case 6://gb->gb = 6;

            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c;
            hm->endb = old_cor[3];
            hm->f[0].a = -FLOATINFTY;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = 0.0;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
            break;
        case 7://gb->a = 7;

            hirsch_path[old_cor[4]+1] = c+1;
            //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
            //foward:
            hm->f[0].a = input_states[0];
            hm->f[0].ga = input_states[1];
            hm->f[0].gb = input_states[2];
            hm->b[0].a = -FLOATINFTY;
            hm->b[0].ga = -FLOATINFTY;
            hm->b[0].gb = 0.0;

            hm->starta = old_cor[0];
            hm->enda = old_cor[4]-1;
            hm->startb = old_cor[2];
            hm->endb = c;
            //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

            //backward:
            hm->starta = old_cor[4]+1;
            hm->enda = old_cor[1];
            hm->startb = c+1;
            hm->endb = old_cor[3];
            hm->f[0].a = 0.0;
            hm->f[0].ga = -FLOATINFTY;
            hm->f[0].gb = -FLOATINFTY;
            hm->b[0].a = input_states[3];
            hm->b[0].ga = input_states[4];
            hm->b[0].gb = input_states[5];

            //fprintf(stderr,"Following last: %d\n",c+1);
            hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
            break;
    }

    return hirsch_path;
}/*}}}*/
int* hirsch_align_two_pp_vector_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{ /*for test purpose, the usage of prof1 prof2 return to original*/
	struct states* f = hm->f;
	struct states* b = hm->b;
	int i,j,c;
	int transition = -1;


	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;

	//int max = -INFTY;
	float max = -INFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


	//prof1 += 64 *(mid+1)
	prof1+= ((old_cor[4]+1) << 6);
	idx1 += (old_cor[4]+1);
	//prof2 += 64 * (hm->startb);
	prof2 += old_cor[2] << 6;
	idx2 += (old_cor[2]);
	//i = hm->startb;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++)
	for(i = old_cor[2]; i < old_cor[3];i++)
	{
		sub = abs(middle -i);
		sub /= 1000; 

		prof2 += 64;
		idx2 += 1;

        wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
        wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
        wf1_im1 = (double)(nResArray2[idx2-1])/(double)(nResArray1[idx1-1])*(weightArray2[idx2-1]/weightArray1[idx1-1]);
        wf2_im1 = (double)(nResArray1[idx1-1])/(double)(nResArray2[idx2-1])*(weightArray1[idx1-1]/weightArray2[idx2-1]);

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1 - prof2[-37])> FLOAT_ZERO * fz_weight || fabs(prof2[-37]*wf2_im1-prof1[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

		//fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
			//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		if(f[i].a+b[i].ga+prof2[27]-sub > max)
		/*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
		/*if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)*/
		{
			max = f[i].a+b[i].ga+prof2[27]-sub;
			/*max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;*/
			//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb+prof1[27] -sub> max)
		/*if(f[i].a+b[i].gb+prof2[27] -sub> max)*/
		/*if(f[i].a+b[i].gb+prof2[27]*wf2_i -sub> max)*/
		{
			max = f[i].a+b[i].gb+prof1[27]-sub;
			/*max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;*/
			//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a+prof2[-37]-sub > max)
		/*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
		/*if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)*/
		{
			max = f[i].ga+b[i].a+prof2[-37]-sub;
			/*max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;*/
			//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb+prof1[29]-sub > max)
			/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)*/
			{
				max = f[i].gb+b[i].gb+prof1[29]-sub;
				/*max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;*/
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb+prof1[28]-sub > max)
			/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)*/
			{
				max = f[i].gb+b[i].gb+prof1[28]-sub;
				/*max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;*/
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[-37]-sub > max)
		/*if(f[i].gb+b[i].a+prof2[-37]-sub > max)*/
		/*if(f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub > max)*/
		{
			max = f[i].gb+b[i].a+prof1[-37]-sub;
			/*max = f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub;*/
			//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}

    wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);

#ifdef DEBUG_PP/*{{{*/
    fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
    if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
        fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
    }
    if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
    }
    if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
    }
#endif/*}}}*/
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[27]-sub > max)
	/*if(f[i].a+b[i].gb+prof2[27]-sub > max)*/
	/*if(f[i].a+b[i].gb+prof2[27]*wf2_i-sub > max)*/
	{
		max = f[i].a+b[i].gb+prof1[27]-sub;
		/*max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;*/
		//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb+prof1[29]-sub > max)
		/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)*/
		{
			max = f[i].gb+b[i].gb+prof1[29]-sub;
			/*max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;*/
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb+prof1[28]-sub > max)
		/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)*/
		{
			max = f[i].gb+b[i].gb+prof1[28]-sub;
			/*max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;*/
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}

	prof1-= (old_cor[4]+1)<<6;
	idx1 -= (old_cor[4]+1);
	//prof2 -= hm->endb << 6;
	prof2 -= old_cor[3] << 6;
	idx2 -= old_cor[3];

	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	//if(transition == -1){
	//	exit(0);
	//}

	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1

			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;

			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			//fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,  nResArray1, nResArray2, weightArray1, weightArray2, hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2, nResArray1, nResArray2, weightArray1, weightArray2, hm,hirsch_path);
			break;
		case 2:// a -> ga = 2

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;


			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4];

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;

			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new0(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
	}

	return hirsch_path;
}/*}}}*/
int* hirsch_align_two_pp_vector_new(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{
    struct states* f = hm->f;
    struct states* b = hm->b;
    int i,j,c;
    int transition = -1;


    //code:
    // a -> a = 1
    // a -> ga = 2
    // a -> gb = 3
    // ga ->ga = 4
    // ga -> a = 5
    //gb->gb = 6;
    //gb->a = 7;

    //int max = -INFTY;
    float max = -INFTY;	
    //float middle =  (hm->endb - hm->startb)/2 + hm->startb;
    float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
    float sub = 0.0;


    //prof1 += 64 *(mid+1)
    prof1+= ((old_cor[4]+1) << 6);
    //prof2 += 64 * (hm->startb);
    prof2 += old_cor[2] << 6;
    //i = hm->startb;
    i = old_cor[2];
    c = -1;
    //for(i = hm->startb; i < hm->endb;i++){
    for(i = old_cor[2]; i < old_cor[3];i++){
        sub = abs(middle -i);
        sub /= 1000; 
        prof2 += 64;
        //fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
        if(f[i].a+b[i].a-sub > max){
            max = f[i].a+b[i].a-sub;
            //		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
            transition = 1;
            c = i;
        }
        /*if(f[i].a+b[i].ga+prof2[27]-sub > max){*/
        if(f[i].a+b[i].ga+prof1[27]-sub > max){
            max = f[i].a+b[i].ga+prof1[27]-sub;
            //		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
            transition = 2;
            c = i;
        }
        /*if(f[i].a+b[i].gb+prof1[27] -sub> max){*/
        if(f[i].a+b[i].gb+prof2[27] -sub> max){
            max = f[i].a+b[i].gb+prof2[27]-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        /*if(f[i].ga+b[i].a+prof2[-37]-sub > max){*/
        if(f[i].ga+b[i].a+prof1[-37]-sub > max){
            max = f[i].ga+b[i].a+prof1[-37]-sub;
            //		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
            transition = 5;
            c = i;
        }


        if(hm->startb == 0){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max){*/
            if(f[i].gb+b[i].gb+prof2[29]-sub > max){
                max = f[i].gb+b[i].gb+prof2[29]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max){*/
            if(f[i].gb+b[i].gb+prof2[28]-sub > max){
                max = f[i].gb+b[i].gb+prof2[28]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }
        /*if(f[i].gb+b[i].a+prof1[-37]-sub > max){*/
        if(f[i].gb+b[i].a+prof2[-37]-sub > max){
            max = f[i].gb+b[i].a+prof2[-37]-sub;
            //		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
            transition = 7;
            c = i;
        }
        }
        //i = hm->endb;
        i = old_cor[3];
        sub = abs(middle -i);
        sub /= 1000; 
        /*if(f[i].a+b[i].gb+prof1[27]-sub > max){*/
        if(f[i].a+b[i].gb+prof2[27]-sub > max){
            max = f[i].a+b[i].gb+prof2[27]-sub;
            //		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
            transition = 3;
            c = i;
        }
        if(hm->endb == hm->len_b){
            /*if(f[i].gb+b[i].gb+prof1[29]-sub > max){*/
            if(f[i].gb+b[i].gb+prof2[29]-sub > max){
                max = f[i].gb+b[i].gb+prof2[29]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }	
        }else{
            /*if(f[i].gb+b[i].gb+prof1[28]-sub > max){*/
            if(f[i].gb+b[i].gb+prof2[28]-sub > max){
                max = f[i].gb+b[i].gb+prof2[28]-sub;
                //			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
                transition = 6;
                c = i;
            }
        }



        prof1-= (old_cor[4]+1)<<6;
        //prof2 -= hm->endb << 6;
        prof2 -= old_cor[3] << 6;

        //fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
        //if(transition == -1){
        //	exit(0);
        //}

        j = hirsch_path[0];
        switch(transition){
            case 1: //a -> a = 1

                hirsch_path[old_cor[4]] = c;
                hirsch_path[old_cor[4]+1] = c+1;

                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLOATINFTY;
                hm->b[0].gb = -FLOATINFTY;
                //fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLOATINFTY;
                hm->f[0].gb = -FLOATINFTY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);
                break;
            case 2:// a -> ga = 2

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLOATINFTY;
                hm->b[0].gb = -FLOATINFTY;


                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4];
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLOATINFTY;
                hm->f[0].ga = 0.0;
                hm->f[0].gb = -FLOATINFTY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);
                break;
            case 3:// a -> gb = 3

                hirsch_path[old_cor[4]] = c;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = 0.0;
                hm->b[0].ga = -FLOATINFTY;
                hm->b[0].gb = -FLOATINFTY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLOATINFTY;
                hm->f[0].ga = -FLOATINFTY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);
                break;
            case 5://ga -> a = 5
                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLOATINFTY;
                hm->b[0].ga = 0.0;
                hm->b[0].gb = -FLOATINFTY;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4];

                hm->startb = old_cor[2];
                hm->endb = c-1;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLOATINFTY;
                hm->f[0].gb = -FLOATINFTY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);
                break;
            case 6://gb->gb = 6;

                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLOATINFTY;
                hm->b[0].ga = -FLOATINFTY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c;
                hm->endb = old_cor[3];
                hm->f[0].a = -FLOATINFTY;
                hm->f[0].ga = -FLOATINFTY;
                hm->f[0].gb = 0.0;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);
                break;
            case 7://gb->a = 7;

                hirsch_path[old_cor[4]+1] = c+1;
                //		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
                //foward:
                hm->f[0].a = input_states[0];
                hm->f[0].ga = input_states[1];
                hm->f[0].gb = input_states[2];
                hm->b[0].a = -FLOATINFTY;
                hm->b[0].ga = -FLOATINFTY;
                hm->b[0].gb = 0.0;

                hm->starta = old_cor[0];
                hm->enda = old_cor[4]-1;
                hm->startb = old_cor[2];
                hm->endb = c;
                //fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);

                //backward:
                hm->starta = old_cor[4]+1;
                hm->enda = old_cor[1];
                hm->startb = c+1;
                hm->endb = old_cor[3];
                hm->f[0].a = 0.0;
                hm->f[0].ga = -FLOATINFTY;
                hm->f[0].gb = -FLOATINFTY;
                hm->b[0].a = input_states[3];
                hm->b[0].ga = input_states[4];
                hm->b[0].gb = input_states[5];

                //fprintf(stderr,"Following last: %d\n",c+1);
                hirsch_path = hirsch_pp_dyn_new(prof1,prof2,hm,hirsch_path);
                break;
        }

        return hirsch_path;
        }/*}}}*/
int* hirsch_align_two_pp_vector_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{
	struct states* f = hm->f;
	struct states* b = hm->b;
	int i,j,c;
	int transition = -1;


	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;

	//int max = -INFTY;
	float max = -INFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


	//prof1 += 64 *(mid+1)
	prof1+= ((old_cor[4]+1) << 6);
	idx1 += (old_cor[4]+1);
	//prof2 += 64 * (hm->startb);
	prof2 += old_cor[2] << 6;
	idx2 += (old_cor[2]);
	//i = hm->startb;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++)
	for(i = old_cor[2]; i < old_cor[3];i++)
	{
		sub = abs(middle -i);
		sub /= 1000; 

		prof2 += 64;
		idx2 += 1;

        wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
        wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
        wf1_im1 = (double)(nResArray2[idx2-1])/(double)(nResArray1[idx1-1])*(weightArray2[idx2-1]/weightArray1[idx1-1]);
        wf2_im1 = (double)(nResArray1[idx1-1])/(double)(nResArray2[idx2-1])*(weightArray1[idx1-1]/weightArray2[idx2-1]);

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1 - prof2[-37])> FLOAT_ZERO * fz_weight || fabs(prof2[-37]*wf2_im1-prof1[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

		//fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
			//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		/*if(f[i].a+b[i].ga+prof2[27]-sub > max)*/
		/*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
		if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)
		{
			max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;
			//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		/*if(f[i].a+b[i].gb+prof1[27] -sub> max)*/
		/*if(f[i].a+b[i].gb+prof2[27] -sub> max)*/
		if(f[i].a+b[i].gb+prof2[27]*wf2_i -sub> max)
		{
			max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;
			//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		/*if(f[i].ga+b[i].a+prof2[-37]-sub > max)*/
		/*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
		if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)
		{
			max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;
			//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			/*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
			if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)
			{
				max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			/*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
			if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)
			{
				max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		/*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
		/*if(f[i].gb+b[i].a+prof2[-37]-sub > max)*/
		if(f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub > max)
		{
			max = f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub;
			//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}

    wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);

#ifdef DEBUG_PP/*{{{*/
    fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
    if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
        fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
    }
    if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
    }
    if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
    }
#endif/*}}}*/
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	/*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
	/*if(f[i].a+b[i].gb+prof2[27]-sub > max)*/
	if(f[i].a+b[i].gb+prof2[27]*wf2_i-sub > max)
	{
		max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;
		//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		/*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
		if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)
		{
			max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		/*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
		if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)
		{
			max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}

	prof1-= (old_cor[4]+1)<<6;
	idx1 -= (old_cor[4]+1);
	//prof2 -= hm->endb << 6;
	prof2 -= old_cor[3] << 6;
	idx2 -= old_cor[3];

	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	//if(transition == -1){
	//	exit(0);
	//}

	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1

			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;

			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			//fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,  nResArray1, nResArray2, weightArray1, weightArray2, hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2, nResArray1, nResArray2, weightArray1, weightArray2, hm,hirsch_path);
			break;
		case 2:// a -> ga = 2

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;


			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4];

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;

			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
	}

	return hirsch_path;
}/*}}}*/
int* hirsch_align_two_pp_vector_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{
	struct states* f = hm->f;
	struct states* b = hm->b;
	int i,j,c;
	int transition = -1;


	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;

	//int max = -INFTY;
	float max = -INFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


	//prof1 += 64 *(mid+1)
	prof1+= ((old_cor[4]+1) << 6);
	idx1 += (old_cor[4]+1);
	//prof2 += 64 * (hm->startb);
	prof2 += old_cor[2] << 6;
	idx2 += (old_cor[2]);
	//i = hm->startb;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++)
	for(i = old_cor[2]; i < old_cor[3];i++)
	{
		sub = abs(middle -i);
		sub /= 1000; 

		prof2 += 64;
		idx2 += 1;

        wf1_i = 1;
        wf2_i = 1;
        wf1_im1 = 1;
        wf2_im1 =1;

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1 - prof2[-37])> FLOAT_ZERO * fz_weight || fabs(prof2[-37]*wf2_im1-prof1[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

		//fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
			//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		/*if(f[i].a+b[i].ga+prof2[27]-sub > max)*/
		/*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
		if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)
		{
			max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;
			//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		/*if(f[i].a+b[i].gb+prof1[27] -sub> max)*/
		/*if(f[i].a+b[i].gb+prof2[27] -sub> max)*/
		if(f[i].a+b[i].gb+prof2[27]*wf2_i -sub> max)
		{
			max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;
			//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		/*if(f[i].ga+b[i].a+prof2[-37]-sub > max)*/
		/*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
		if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)
		{
			max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;
			//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			/*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
			if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)
			{
				max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			/*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
			if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)
			{
				max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		/*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
		/*if(f[i].gb+b[i].a+prof2[-37]-sub > max)*/
		if(f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub > max)
		{
			max = f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub;
			//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}

    wf2_i = 1;

#ifdef DEBUG_PP/*{{{*/
    fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
    if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
        fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
    }
    if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
    }
    if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
    }
#endif/*}}}*/
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	/*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
	/*if(f[i].a+b[i].gb+prof2[27]-sub > max)*/
	if(f[i].a+b[i].gb+prof2[27]*wf2_i-sub > max)
	{
		max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;
		//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		/*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
		if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)
		{
			max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		/*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
		if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)
		{
			max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}

	prof1-= (old_cor[4]+1)<<6;
	idx1 -= (old_cor[4]+1);
	//prof2 -= hm->endb << 6;
	prof2 -= old_cor[3] << 6;
	idx2 -= old_cor[3];

	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	//if(transition == -1){
	//	exit(0);
	//}

	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1

			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;

			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			//fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,  nResArray1, nResArray2, weightArray1, weightArray2, hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2, nResArray1, nResArray2, weightArray1, weightArray2, hm,hirsch_path);
			break;
		case 2:// a -> ga = 2

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;


			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4];

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;

			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test1(prof1,prof2,nResArray1, nResArray2, weightArray1, weightArray2,hm,hirsch_path);
			break;
	}

	return hirsch_path;
}/*}}}*/
int* hirsch_align_two_pp_vector_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])/*{{{*/
{
	struct states* f = hm->f;
	struct states* b = hm->b;
	int i,j,c;
	int transition = -1;


	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;

	//int max = -INFTY;
	float max = -INFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


	//prof1 += 64 *(mid+1)
	prof1+= ((old_cor[4]+1) << 6);
	idx1 += (old_cor[4]+1);
	//prof2 += 64 * (hm->startb);
	prof2 += old_cor[2] << 6;
	idx2 += (old_cor[2]);
	//i = hm->startb;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++)
	for(i = old_cor[2]; i < old_cor[3];i++)
	{
		sub = abs(middle -i);
		sub /= 1000; 

		prof2 += 64;
		idx2 += 1;

        wf1_i = 1.0/nsip2;
        wf2_i = 1.0/nsip1;
        wf1_im1 = 1.0/nsip2;
        wf2_im1 =1.0/nsip1;

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1 - prof2[-37])> FLOAT_ZERO * fz_weight || fabs(prof2[-37]*wf2_im1-prof1[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

		//fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
			//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		/*if(f[i].a+b[i].ga+prof2[27]-sub > max)*/
		/*if(f[i].a+b[i].ga+prof1[27]-sub > max)*/
		if(f[i].a+b[i].ga+prof1[27]*wf1_i-sub > max)
		{
			max = f[i].a+b[i].ga+prof1[27]*wf1_i-sub;
			//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		/*if(f[i].a+b[i].gb+prof1[27] -sub> max)*/
		/*if(f[i].a+b[i].gb+prof2[27] -sub> max)*/
		if(f[i].a+b[i].gb+prof2[27]*wf2_i -sub> max)
		{
			max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;
			//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		/*if(f[i].ga+b[i].a+prof2[-37]-sub > max)*/
		/*if(f[i].ga+b[i].a+prof1[-37]-sub > max)*/
		if(f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub > max)
		{
			max = f[i].ga+b[i].a+prof1[-37]*wf1_im1-sub;
			//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			/*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
			if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)
			{
				max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			/*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
			/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
			if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)
			{
				max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;
				//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		/*if(f[i].gb+b[i].a+prof1[-37]-sub > max)*/
		/*if(f[i].gb+b[i].a+prof2[-37]-sub > max)*/
		if(f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub > max)
		{
			max = f[i].gb+b[i].a+prof2[-37]*wf2_im1-sub;
			//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}

    wf2_i = 1.0/nsip1;

#ifdef DEBUG_PP/*{{{*/
    fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
    if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
        fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
    }
    if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
    }
    if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
        fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
    }
#endif/*}}}*/
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	/*if(f[i].a+b[i].gb+prof1[27]-sub > max)*/
	/*if(f[i].a+b[i].gb+prof2[27]-sub > max)*/
	if(f[i].a+b[i].gb+prof2[27]*wf2_i-sub > max)
	{
		max = f[i].a+b[i].gb+prof2[27]*wf2_i-sub;
		//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		/*if(f[i].gb+b[i].gb+prof1[29]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[29]-sub > max)*/
		if(f[i].gb+b[i].gb+prof2[29]*wf2_i-sub > max)
		{
			max = f[i].gb+b[i].gb+prof2[29]*wf2_i-sub;
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		/*if(f[i].gb+b[i].gb+prof1[28]-sub > max)*/
		/*if(f[i].gb+b[i].gb+prof2[28]-sub > max)*/
		if(f[i].gb+b[i].gb+prof2[28]*wf2_i-sub > max)
		{
			max = f[i].gb+b[i].gb+prof2[28]*wf2_i-sub;
			//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}

	prof1-= (old_cor[4]+1)<<6;
	idx1 -= (old_cor[4]+1);
	//prof2 -= hm->endb << 6;
	prof2 -= old_cor[3] << 6;
	idx2 -= old_cor[3];

	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	//if(transition == -1){
	//	exit(0);
	//}

	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1

			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;

			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			//fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,  nResArray1, nResArray2, nsip1, nsip2, weightArray1, weightArray2, hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2, nResArray1, nResArray2,nsip1, nsip2, weightArray1, weightArray2, hm,hirsch_path);
			break;
		case 2:// a -> ga = 2

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;


			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2, weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2, nsip1, nsip2, weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3

			hirsch_path[old_cor[4]] = c;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4];

			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;

			hirsch_path[old_cor[4]+1] = c+1;
			//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;

			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];

			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn_new2_test2(prof1,prof2,nResArray1, nResArray2,nsip1, nsip2,  weightArray1, weightArray2,hm,hirsch_path);
			break;
	}

	return hirsch_path;
}/*}}}*/

struct states* foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)/*{{{*/
{
    unsigned int freq[26];

    struct states* s = hm->f;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

    prof1 += (hm->starta) << 6;
    prof2 +=  (hm->startb) << 6;
    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            s[j].gb = -FLOATINFTY;
        }
        prof2+=64;
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];
            s[j].gb = -FLOATINFTY;
        }	
        prof2+=64;	
    }

    prof2 -= (hm->endb-hm->startb) << 6;

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;


    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;


        if(hm->startb){
            s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
        }else{
            s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];
        }
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2 += 64;
            ca = s[j].a;

            pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);

            pgb = s[j].gb;

            s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);

            pa = ca;


            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 += 64;
        ca = s[j].a;

        pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);

        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;

        if (hm->endb != hm->len_b){
            s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
        }else{
            s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
        }
        prof2 -= (hm->endb-hm->startb) << 6;

    }
    prof1 -=  (hm->enda) << 6;
    return s;
}/*}}}*/
struct states* foward_hirsch_pp_dyn_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{ /*for test purpose, the usage of prof1 prof2 return to original*/ 
    unsigned int freq[26];

    struct states* s = hm->f;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;


	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


    prof1 += (hm->starta) << 6;
    prof2 +=  (hm->startb) << 6;
	idx1 += hm->starta;
	idx2 += hm->startb;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;


    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);

            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);*/
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight= (nResArray1[idx1]+nResArray2[idx2])/2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }
        prof2+=64;
		idx2 += 1;
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29];*/
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29]*wf1_i;*/
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]=%f, prof1[27]*wf1_i =%f, nRes2[%d]=%d, nRes1[%d]=%d, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], prof1[27]*wf1_i, idx2, nResArray2[idx2],idx1, nResArray1[idx1],wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }	
        prof2+=64;	
		idx2 += 1;
    }

    prof2 -= (hm->endb-hm->startb) << 6;
	idx2 -= (hm->endb-hm->startb);

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;


    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
		idx1 += 1;
        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

		wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);

        if(hm->startb){
            s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
            /*s[hm->startb].gb = MAX(pgb+prof2[28],pa+prof2[27]);*/
            /*s[hm->startb].gb = MAX(pgb+prof2[28]*wf2_i,pa+prof2[27]*wf2_i);*/
        }else{
            s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof2[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof2[29]*wf2_i;*/
        }

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO*fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] *wf2_i=%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO *fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO*fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2 += 64;
			idx2 += 1;
            ca = s[j].a;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
			wf1_im1 = (double)(nResArray2[idx2-1])/(double)(nResArray1[idx1-1])*(weightArray2[idx2-1]/weightArray1[idx1-1]);
            wf2_im1 = (double)(nResArray1[idx1-1])/(double)(nResArray2[idx2-1])*(weightArray1[idx1-1]/weightArray2[idx2-1]);

            pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);
            /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);*/

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            /*s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);*/

            pgb = s[j].gb;

            s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);
            /*s[j].gb = MAX(pgb+prof2[28] ,ca+prof2[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);*/

            pa = ca;


            xa = s[j].a;
            xga = s[j].ga;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 += 64;
		idx2 += 1;
        ca = s[j].a;

        wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
        wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
        wf1_im1 = (double)(nResArray2[idx2-1])/(double)(nResArray1[idx1-1])*(weightArray2[idx2-1]/weightArray1[idx1-1]);
        wf2_im1 = (double)(nResArray1[idx1-1])/(double)(nResArray2[idx2-1])*(weightArray1[idx1-1]/weightArray2[idx2-1]);

        pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);
        /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
        /*pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);*/

        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;

        if (hm->endb != hm->len_b){
            s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
            /*s[j].gb = MAX(s[j].gb+prof2[28] ,ca+prof2[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);*/
        }else{
            s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
            /*s[j].gb = MAX(s[j].gb,ca)+ prof2[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+ prof2[29]*wf2_i;*/
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        prof2 -= (hm->endb-hm->startb) << 6;
		idx2 -= (hm->endb-hm->startb);

    }
    prof1 -=  (hm->enda) << 6;
	idx1 -= (hm->enda);
    return s;
}/*}}}*/
struct states* foward_hirsch_pp_dyn_new(const float* prof1,const float* prof2,struct hirsch_mem* hm)/*{{{*/
{
    unsigned int freq[26];

    struct states* s = hm->f;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

    prof1 += (hm->starta) << 6;
    prof2 +=  (hm->startb) << 6;
    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;
    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);*/
            s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);
            s[j].gb = -FLOATINFTY;
        }
        prof2+=64;
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29];
            s[j].gb = -FLOATINFTY;
        }	
        prof2+=64;	
    }

    prof2 -= (hm->endb-hm->startb) << 6;

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;


    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;


        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            s[hm->startb].gb = MAX(pgb+prof2[28],pa+prof2[27]);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];*/
            s[hm->startb].gb = MAX(pgb,pa)+ prof2[29];
        }
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2 += 64;
            ca = s[j].a;

            /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
            pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);*/
            s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);*/
            s[j].gb = MAX(pgb+prof2[28] ,ca+prof2[27]);

            pa = ca;


            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 += 64;
        ca = s[j].a;

        /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
        pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);

        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28] ,ca+prof2[27]);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            s[j].gb = MAX(s[j].gb,ca)+ prof2[29];
        }
        prof2 -= (hm->endb-hm->startb) << 6;

    }
    prof1 -=  (hm->enda) << 6;
    return s;
}/*}}}*/
struct states* foward_hirsch_pp_dyn_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{
    unsigned int freq[26];

    struct states* s = hm->f;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;


	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


    prof1 += (hm->starta) << 6;
    prof2 +=  (hm->startb) << 6;
	idx1 += hm->starta;
	idx2 += hm->startb;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;


    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);

            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight= (nResArray1[idx1]+nResArray2[idx2])/2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }
        prof2+=64;
		idx2 += 1;
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];*/
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29];*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]=%f, prof1[27]*wf1_i =%f, nRes2[%d]=%d, nRes1[%d]=%d, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], prof1[27]*wf1_i, idx2, nResArray2[idx2],idx1, nResArray1[idx1],wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }	
        prof2+=64;	
		idx2 += 1;
    }

    prof2 -= (hm->endb-hm->startb) << 6;
	idx2 -= (hm->endb-hm->startb);

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;


    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
		idx1 += 1;
        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

		wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);

        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb+prof2[28],pa+prof2[27]);*/
            s[hm->startb].gb = MAX(pgb+prof2[28]*wf2_i,pa+prof2[27]*wf2_i);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof2[29];*/
            s[hm->startb].gb = MAX(pgb,pa)+ prof2[29]*wf2_i;
        }

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO*fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] *wf2_i=%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO *fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO*fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2 += 64;
			idx2 += 1;
            ca = s[j].a;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
			wf1_im1 = (double)(nResArray2[idx2-1])/(double)(nResArray1[idx1-1])*(weightArray2[idx2-1]/weightArray1[idx1-1]);
            wf2_im1 = (double)(nResArray1[idx1-1])/(double)(nResArray2[idx2-1])*(weightArray1[idx1-1]/weightArray2[idx2-1]);
#ifdef CHECK_PP/*{{{*/
            fprintf(stderr,"\nprof1[27](%g)/nRes1(%d)/weightArray1(%g)=%g \t prof2[27](%g)/nRes2(%d)/weightArray2(%g)=%g\n",prof1[27],nResArray1[idx1], weightArray1[idx1], prof1[27]/nResArray1[idx1]/weightArray1[idx1], prof2[27], nResArray2[idx2],weightArray2[idx2],  prof2[27]/nResArray2[idx2]/weightArray2[idx2]);
#endif/*}}}*/

            /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
            pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28] ,ca+prof2[27]);*/
            s[j].gb = MAX(pgb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);

            pa = ca;


            xa = s[j].a;
            xga = s[j].ga;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 += 64;
		idx2 += 1;
        ca = s[j].a;

        wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
        wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
        wf1_im1 = (double)(nResArray2[idx2-1])/(double)(nResArray1[idx1-1])*(weightArray2[idx2-1]/weightArray1[idx1-1]);
        wf2_im1 = (double)(nResArray1[idx1-1])/(double)(nResArray2[idx2-1])*(weightArray1[idx1-1]/weightArray2[idx2-1]);

        /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
        pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);

        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28] ,ca+prof2[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+ prof2[29];*/
            s[j].gb = MAX(s[j].gb,ca)+ prof2[29]*wf2_i;
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        prof2 -= (hm->endb-hm->startb) << 6;
		idx2 -= (hm->endb-hm->startb);

    }
    prof1 -=  (hm->enda) << 6;
	idx1 -= (hm->enda);
    return s;
}/*}}}*/
struct states* foward_hirsch_pp_dyn_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{
    /*set wf1=wf2=1*/
    unsigned int freq[26];

    struct states* s = hm->f;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;


	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


    prof1 += (hm->starta) << 6;
    prof2 +=  (hm->startb) << 6;
	idx1 += hm->starta;
	idx2 += hm->startb;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;


    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = 1;

            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight= (nResArray1[idx1]+nResArray2[idx2])/2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }
        prof2+=64;
		idx2 += 1;
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = 1;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];*/
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29];*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]=%f, prof1[27]*wf1_i =%f, nRes2[%d]=%d, nRes1[%d]=%d, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], prof1[27]*wf1_i, idx2, nResArray2[idx2],idx1, nResArray1[idx1],wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }	
        prof2+=64;	
		idx2 += 1;
    }

    prof2 -= (hm->endb-hm->startb) << 6;
	idx2 -= (hm->endb-hm->startb);

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;


    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
		idx1 += 1;
        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

		wf2_i = 1;

        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb+prof2[28],pa+prof2[27]);*/
            s[hm->startb].gb = MAX(pgb+prof2[28]*wf2_i,pa+prof2[27]*wf2_i);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof2[29];*/
            s[hm->startb].gb = MAX(pgb,pa)+ prof2[29]*wf2_i;
        }

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO*fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] *wf2_i=%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO *fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO*fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2 += 64;
			idx2 += 1;
            ca = s[j].a;

			wf1_i = 1;
            wf2_i = 1;
			wf1_im1 = 1;
            wf2_im1 = 1;
#ifdef CHECK_PP/*{{{*/
            fprintf(stderr,"\nprof1[27](%g)/nRes1(%d)/weightArray1(%g)=%g \t prof2[27](%g)/nRes2(%d)/weightArray2(%g)=%g\n",prof1[27],nResArray1[idx1], weightArray1[idx1], prof1[27]/nResArray1[idx1]/weightArray1[idx1], prof2[27], nResArray2[idx2],weightArray2[idx2],  prof2[27]/nResArray2[idx2]/weightArray2[idx2]);
#endif/*}}}*/

            /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
            pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
#ifdef CHECK_PA
            fprintf(stderr,"forward_pp: j=%d, freq[0]=%d, freq[%d]=%d, prof1[freq[c]]=%g, prof2[freq[c]]=%g, nsip1=%g, nsip2=%g\n", j, freq[0], c, freq[c], prof1[freq[c]], prof2[freq[c]], weightArray1[idx1], weightArray2[idx2] );
#endif 
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28] ,ca+prof2[27]);*/
            s[j].gb = MAX(pgb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);

            pa = ca;


            xa = s[j].a;
            xga = s[j].ga;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 += 64;
		idx2 += 1;
        ca = s[j].a;

        wf1_i = 1;
        wf2_i = 1;
        wf1_im1 = 1;
        wf2_im1 = 1;

        /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
        pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);

        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28] ,ca+prof2[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+ prof2[29];*/
            s[j].gb = MAX(s[j].gb,ca)+ prof2[29]*wf2_i;
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        prof2 -= (hm->endb-hm->startb) << 6;
		idx2 -= (hm->endb-hm->startb);

    }
    prof1 -=  (hm->enda) << 6;
	idx1 -= (hm->enda);
    return s;
}/*}}}*/
struct states* foward_hirsch_pp_dyn_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2,unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{
    /*set wf1=wf2=1*/
    unsigned int freq[26];

    struct states* s = hm->f;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;


	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_im1 = 0.0;
	float wf2_im1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif


    prof1 += (hm->starta) << 6;
    prof2 +=  (hm->startb) << 6;
	idx1 += hm->starta;
	idx2 += hm->startb;

    s[hm->startb].a = s[0].a;
    s[hm->startb].ga = s[0].ga;
    s[hm->startb].gb = s[0].gb;


    if(hm->startb){
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = 1.0/nsip2;

            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);*/
            /*s[j].ga = MAX(s[j-1].ga+prof1[28],s[j-1].a+prof1[27]);*/
            s[j].ga = MAX(s[j-1].ga+prof1[28]*wf1_i,s[j-1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight= (nResArray1[idx1]+nResArray2[idx2])/2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }
        prof2+=64;
		idx2 += 1;
    }else{
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2+=64;
			idx2 += 1;

			wf1_i = 1.0/nsip2;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];*/
            /*s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29];*/
            s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO*fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]=%f, prof1[27]*wf1_i =%f, nRes2[%d]=%d, nRes1[%d]=%d, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], prof1[27]*wf1_i, idx2, nResArray2[idx2],idx1, nResArray1[idx1],wf1_i);
            }
            if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO *fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO*fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }	
        prof2+=64;	
		idx2 += 1;
    }

    prof2 -= (hm->endb-hm->startb) << 6;
	idx2 -= (hm->endb-hm->startb);

    s[hm->endb].a = -FLOATINFTY;
    s[hm->endb].ga = -FLOATINFTY;
    s[hm->endb].gb = -FLOATINFTY;


    for (i = hm->starta;i < hm->enda;i++){
        prof1 += 64;
		idx1 += 1;
        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->startb].a;
        pga = s[hm->startb].ga;
        pgb = s[hm->startb].gb;
        s[hm->startb].a = -FLOATINFTY;
        s[hm->startb].ga = -FLOATINFTY;

        xa = s[hm->startb].a;
        xga = s[hm->startb].ga;

		wf2_i = 1.0/nsip1;

        if(hm->startb){
            /*s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);*/
            /*s[hm->startb].gb = MAX(pgb+prof2[28],pa+prof2[27]);*/
            s[hm->startb].gb = MAX(pgb+prof2[28]*wf2_i,pa+prof2[27]*wf2_i);
        }else{
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];*/
            /*s[hm->startb].gb = MAX(pgb,pa)+ prof2[29];*/
            s[hm->startb].gb = MAX(pgb,pa)+ prof2[29]*wf2_i;
        }

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO*fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] *wf2_i=%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO *fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO*fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        for (j = hm->startb+1; j < hm->endb;j++){
            prof2 += 64;
			idx2 += 1;
            ca = s[j].a;

			wf1_i = 1.0/nsip2;
            wf2_i = 1.0/nsip1;
			wf1_im1 = 1.0/nsip2;
            wf2_im1 = 1.0/nsip1;
#ifdef CHECK_PP/*{{{*/
            fprintf(stderr,"\nprof1[27](%g)/nRes1(%d)/weightArray1(%g)=%g \t prof2[27](%g)/nRes2(%d)/weightArray2(%g)=%g\n",prof1[27],nResArray1[idx1], weightArray1[idx1], prof1[27]/nResArray1[idx1]/weightArray1[idx1], prof2[27], nResArray2[idx2],weightArray2[idx2],  prof2[27]/nResArray2[idx2]/weightArray2[idx2]);
#endif/*}}}*/

            /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
            /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
            pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
#ifdef CHECK_PA
            fprintf(stderr,"forward_pp: j=%d, freq[0]=%d, freq[%d]=%d, prof1[freq[c]]=%g, prof2[freq[c]]=%g, nsip1=%g, nsip2=%g\n", j, freq[0], c, freq[c], prof1[freq[c]], prof2[freq[c]], weightArray1[idx1], weightArray2[idx2] );
#endif 
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);*/
            /*s[j].ga = MAX(xga+prof1[28],xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i,xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28] ,ca+prof2[27]);*/
            s[j].gb = MAX(pgb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);

            pa = ca;


            xa = s[j].a;
            xga = s[j].ga;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 += 64;
		idx2 += 1;
        ca = s[j].a;

        wf1_i = 1.0/nsip2;
        wf2_i = 1.0/nsip1;
        wf1_im1 = 1.0/nsip2;
        wf2_im1 = 1.0/nsip1;

        /*pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);*/
        /*pa = MAX3(pa,pga + prof1[-37],pgb + prof2[-37]);*/
        pa = MAX3(pa,pga + prof1[-37]*wf1_im1,pgb + prof2[-37]*wf2_im1);

        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;

        s[j].a = pa;

        s[j].ga = -FLOATINFTY;

        if (hm->endb != hm->len_b){
            /*s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28] ,ca+prof2[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i ,ca+prof2[27]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+ prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+ prof2[29];*/
            s[j].gb = MAX(s[j].gb,ca)+ prof2[29]*wf2_i;
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight || fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight|| fabs(prof2[-37]*wf2_im1 -prof1[-37] ) > FLOAT_ZERO * fz_weight || fabs(prof1[-37]*wf1_im1-prof2[-37])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[-37]=%f, prof2[-37]=%f, wf1_i=%f, wf2_i=%f, wf1_im1=%f, wf2_im1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[-37],prof2[-37], wf1_i, wf2_i, wf1_im1, wf2_im1);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        prof2 -= (hm->endb-hm->startb) << 6;
		idx2 -= (hm->endb-hm->startb);

    }
    prof1 -=  (hm->enda) << 6;
	idx1 -= (hm->enda);
    return s;
}/*}}}*/

struct states* backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)/*{{{*/
{
    unsigned int freq[26];
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

    prof1 += (hm->enda+1) << 6;
    prof2 += (hm->endb+1) << 6;
    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);
            s[j].gb = -FLOATINFTY;
        }
        prof2 -= 64;
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];
            s[j].gb = -FLOATINFTY;
        }
        prof2 -= 64;
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;

        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

        if(hm->endb != hm->len_b){
            s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);
        }else{
            s[hm->endb].gb = MAX(pgb,pa)+prof1[29];
        }

        prof2 += (hm->endb-hm->startb) << 6;
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
            ca = s[j].a;

            pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
            s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);

            pgb = s[j].gb;

            s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 -= 64;
        ca = s[j].a;

        pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;
        s[j].a = pa;

        //pga = s[j].ga;
        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

        //pgb = s[j].gb;
        if(hm->startb){
            s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
        }else{
            s[j].gb = MAX(s[j].gb,ca)+prof1[29];
        }

        //pa = ca;
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_pp_dyn_new0(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{ /*for test purpose, the usage of prof1 prof2 return to original*/ 
    unsigned int freq[26];
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_ip1 = 0.0;
	float wf2_ip1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif

    prof1 += (hm->enda+1) << 6;
    prof2 += (hm->endb+1) << 6;
	idx1 += (hm->enda+1);
	idx2 += (hm->endb+1);
    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);*/
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            s[j].a = -FLOATINFTY;
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;*/
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27] =%f, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
		idx1 -= 1;

        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

		wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        if(hm->endb != hm->len_b){
            s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);
            /*s[hm->endb].gb = MAX(pgb+prof2[28] ,pa+prof2[27]);*/
            /*s[hm->endb].gb = MAX(pgb+prof2[28]*wf2_i ,pa+prof2[27]*wf2_i);*/
        }else{
            s[hm->endb].gb = MAX(pgb,pa)+prof1[29];
            /*s[hm->endb].gb = MAX(pgb,pa)+prof2[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa)+prof2[29]*wf2_i;*/
        }

        prof2 += (hm->endb-hm->startb) << 6;
		idx2  += (hm->endb-hm->startb);
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
            ca = s[j].a;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
			wf1_ip1 = (double)(nResArray2[idx2+1])/(double)(nResArray1[idx1+1])*(weightArray2[idx2+1]/weightArray1[idx1+1]);
            wf2_ip1 = (double)(nResArray1[idx1+1])/(double)(nResArray2[idx2+1])*(weightArray1[idx1+1]/weightArray2[idx2+1]);
            pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
            /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
            /*pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91])*wf2_ip1;*/

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
            s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);
            /*s[j].ga = MAX(xga+prof1[28], xa+prof1[27]);*/
            /*s[j].ga = MAX(xga+prof1[28]*wf1_i, xa+prof1[27]*wf1_i);*/

            pgb = s[j].gb;

            s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);
            /*s[j].gb = MAX(pgb+prof2[28], ca+prof2[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);*/

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[91]=%f, prof2[91]=%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[91],prof2[91], wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
        ca = s[j].a;

        wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
        wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
        wf1_ip1 = (double)(nResArray2[idx2+1])/(double)(nResArray1[idx1+1])*(weightArray2[idx2+1]/weightArray1[idx1+1]);
        wf2_ip1 = (double)(nResArray1[idx1+1])/(double)(nResArray2[idx2+1])*(weightArray1[idx1+1]/weightArray2[idx2+1]);
        pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
        /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
        /*pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);*/
        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;
        s[j].a = pa;

        //pga = s[j].ga;
        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

        //pgb = s[j].gb;
        if(hm->startb){
            s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
            /*s[j].gb = MAX(s[j].gb+prof2[28], ca+prof2[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);*/
        }else{
            s[j].gb = MAX(s[j].gb,ca)+prof1[29];
            /*s[j].gb = MAX(s[j].gb,ca)+prof2[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+prof2[29]*wf2_i;*/
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

        //pa = ca;
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_pp_dyn_new(const float* prof1,const float* prof2,struct hirsch_mem* hm)/*{{{*/
{
    unsigned int freq[26];
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

    prof1 += (hm->enda+1) << 6;
    prof2 += (hm->endb+1) << 6;
    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);
            s[j].gb = -FLOATINFTY;
        }
        prof2 -= 64;
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];*/
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29];
            s[j].gb = -FLOATINFTY;
        }
        prof2 -= 64;
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;

        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);*/
            s[hm->endb].gb = MAX(pgb+prof2[28] ,pa+prof2[27]);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa)+prof1[29];*/
            s[hm->endb].gb = MAX(pgb,pa)+prof2[29];
        }

        prof2 += (hm->endb-hm->startb) << 6;
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
            ca = s[j].a;

            /*pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);*/
            pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);*/
            s[j].ga = MAX(xga+prof1[28], xa+prof1[27]);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);*/
            s[j].gb = MAX(pgb+prof2[28], ca+prof2[27]);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 -= 64;
        ca = s[j].a;

        /*pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);*/
        pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);
        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;
        s[j].a = pa;

        //pga = s[j].ga;
        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

        //pgb = s[j].gb;
        if(hm->startb){
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28], ca+prof2[27]);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            s[j].gb = MAX(s[j].gb,ca)+prof2[29];
        }

        //pa = ca;
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_pp_dyn_new2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{
    unsigned int freq[26];
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_ip1 = 0.0;
	float wf2_ip1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif

    prof1 += (hm->enda+1) << 6;
    prof2 += (hm->endb+1) << 6;
	idx1 += (hm->enda+1);
	idx2 += (hm->endb+1);
    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];*/
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27] =%f, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
		idx1 -= 1;

        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

		wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb+prof2[28] ,pa+prof2[27]);*/
            s[hm->endb].gb = MAX(pgb+prof2[28]*wf2_i ,pa+prof2[27]*wf2_i);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa)+prof2[29];*/
            s[hm->endb].gb = MAX(pgb,pa)+prof2[29]*wf2_i;
        }

        prof2 += (hm->endb-hm->startb) << 6;
		idx2  += (hm->endb-hm->startb);
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
            ca = s[j].a;

			wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
            wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
			wf1_ip1 = (double)(nResArray2[idx2+1])/(double)(nResArray1[idx1+1])*(weightArray2[idx2+1]/weightArray1[idx1+1]);
            wf2_ip1 = (double)(nResArray1[idx1+1])/(double)(nResArray2[idx2+1])*(weightArray1[idx1+1]/weightArray2[idx2+1]);
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[91]=%f, prof2[91]=%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[91],prof2[91], wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/

            /*pa = MAX3(pa,pga + prof2[91],        pgb + prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
            pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);*/
            /*s[j].ga = MAX(xga+prof1[28], xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i, xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28], ca+prof2[27]);*/
            s[j].gb = MAX(pgb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 -= 64;
		idx2 -= 1;
        ca = s[j].a;

        wf1_i = (double)(nResArray2[idx2])/(double)(nResArray1[idx1])*(weightArray2[idx2]/weightArray1[idx1]);
        wf2_i = (double)(nResArray1[idx1])/(double)(nResArray2[idx2])*(weightArray1[idx1]/weightArray2[idx2]);
        wf1_ip1 = (double)(nResArray2[idx2+1])/(double)(nResArray1[idx1+1])*(weightArray2[idx2+1]/weightArray1[idx1+1]);
        wf2_ip1 = (double)(nResArray1[idx1+1])/(double)(nResArray2[idx2+1])*(weightArray1[idx1+1]/weightArray2[idx2+1]);

        /*pa = MAX3(pa,pga + prof2[91],        pgb + prof1[91]);*/
        /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
        pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);
        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;
        s[j].a = pa;

        //pga = s[j].ga;
        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

        //pgb = s[j].gb;
        if(hm->startb){
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28], ca+prof2[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+prof2[29];*/
            s[j].gb = MAX(s[j].gb,ca)+prof2[29]*wf2_i;
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

        //pa = ca;
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_pp_dyn_new2_test1(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{
    /*set wf1=wf2=1*/
    unsigned int freq[26];
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_ip1 = 0.0;
	float wf2_ip1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif

    prof1 += (hm->enda+1) << 6;
    prof2 += (hm->endb+1) << 6;
	idx1 += (hm->enda+1);
	idx2 += (hm->endb+1);
    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = 1;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = 1;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];*/
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27] =%f, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
		idx1 -= 1;

        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

		wf2_i = 1;

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb+prof2[28] ,pa+prof2[27]);*/
            s[hm->endb].gb = MAX(pgb+prof2[28]*wf2_i ,pa+prof2[27]*wf2_i);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa)+prof2[29];*/
            s[hm->endb].gb = MAX(pgb,pa)+prof2[29]*wf2_i;
        }

        prof2 += (hm->endb-hm->startb) << 6;
		idx2  += (hm->endb-hm->startb);
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
            ca = s[j].a;

			wf1_i = 1;
            wf2_i = 1;
			wf1_ip1 = 1;
            wf2_ip1 = 1;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[91]=%f, prof2[91]=%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[91],prof2[91], wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/

            /*pa = MAX3(pa,pga + prof2[91],        pgb + prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
            pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);*/
            /*s[j].ga = MAX(xga+prof1[28], xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i, xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28], ca+prof2[27]);*/
            s[j].gb = MAX(pgb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 -= 64;
		idx2 -= 1;
        ca = s[j].a;

        wf1_i = 1;
        wf2_i = 1;
        wf1_ip1 = 1;
        wf2_ip1 = 1;

        /*pa = MAX3(pa,pga + prof2[91],        pgb + prof1[91]);*/
        /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
        pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);
        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;
        s[j].a = pa;

        //pga = s[j].ga;
        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

        //pgb = s[j].gb;
        if(hm->startb){
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28], ca+prof2[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+prof2[29];*/
            s[j].gb = MAX(s[j].gb,ca)+prof2[29]*wf2_i;
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

        //pa = ca;
    }		
    return s;
}/*}}}*/
struct states* backward_hirsch_pp_dyn_new2_test2(const float* prof1,const float* prof2,unsigned int* nResArray1, unsigned int* nResArray2, unsigned int nsip1, unsigned int nsip2, float* weightArray1, float* weightArray2, struct hirsch_mem* hm)/*{{{*/
{
    /*set wf1=wf2=1*/
    unsigned int freq[26];
    struct states* s = hm->b;
    register float pa = 0;
    register float pga = 0;
    register float pgb = 0;
    register float ca = 0;

    register float xa = 0;
    register float xga = 0;

    register int i = 0;
    register int j = 0;
    register int c = 0;

	register int idx1 = 0;
	register int idx2 = 0;
	float wf1_i = 0.0;
	float wf2_i = 0.0;
	float wf1_ip1 = 0.0;
	float wf2_ip1 = 0.0;
#ifdef DEBUG_PP
    float fz_weight = 0.0;
#endif

    prof1 += (hm->enda+1) << 6;
    prof2 += (hm->endb+1) << 6;
	idx1 += (hm->enda+1);
	idx2 += (hm->endb+1);
    s[hm->endb].a = s[0].a;
    s[hm->endb].ga = s[0].ga;
    s[hm->endb].gb = s[0].gb;
    if(hm->endb != hm->len_b){
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = 1.0/nsip2;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);*/
            /*s[j].ga = MAX(s[j+1].ga+prof1[28],s[j+1].a+prof1[27]);*/
            s[j].ga = MAX(s[j+1].ga+prof1[28]*wf1_i,s[j+1].a+prof1[27]*wf1_i);
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27]*wf1_i =%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27]*wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }else{
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
			wf1_i = 1.0/nsip2;
            s[j].a = -FLOATINFTY;
            /*s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];*/
            s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof1[29]*wf1_i;
            s[j].gb = -FLOATINFTY;
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27] -prof1[27]*wf1_i ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27] =%f, prof1[27] =%f, wf1_i=%f\n",__LINE__, idx1, idx2 , prof2[27], prof1[27], wf1_i);
        }
        if (fabs(prof2[28] - prof1[28]*wf1_i)>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29] - prof1[29]*wf1_i)>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        }
        prof2 -= 64;
		idx2 -= 1;
    }

    s[hm->startb].a = -FLOATINFTY;
    s[hm->startb].ga = -FLOATINFTY;
    s[hm->startb].gb = -FLOATINFTY;

    i = hm->enda-hm->starta;
    while(i--){
        prof1 -= 64;
		idx1 -= 1;

        c = 1;
        for (j = 0;j < 26; j++){
            if(prof1[j]){
                freq[c] = j;
                c++;
            }
        }
        freq[0] = c;

        pa = s[hm->endb].a;
        pga = s[hm->endb].ga;
        pgb = s[hm->endb].gb;
        s[hm->endb].a = -FLOATINFTY;
        s[hm->endb].ga = -FLOATINFTY;

        xa = s[hm->endb].a;
        xga = s[hm->endb].ga;

		wf2_i = 1.0/nsip1;

#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/
        if(hm->endb != hm->len_b){
            /*s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);*/
            /*s[hm->endb].gb = MAX(pgb+prof2[28] ,pa+prof2[27]);*/
            s[hm->endb].gb = MAX(pgb+prof2[28]*wf2_i ,pa+prof2[27]*wf2_i);
        }else{
            /*s[hm->endb].gb = MAX(pgb,pa)+prof1[29];*/
            /*s[hm->endb].gb = MAX(pgb,pa)+prof2[29];*/
            s[hm->endb].gb = MAX(pgb,pa)+prof2[29]*wf2_i;
        }

        prof2 += (hm->endb-hm->startb) << 6;
		idx2  += (hm->endb-hm->startb);
        for(j = hm->endb-1;j > hm->startb;j--){
            prof2 -= 64;
			idx2 -= 1;
            ca = s[j].a;

			wf1_i = 1.0/nsip2;
            wf2_i = 1.0/nsip1;
			wf1_ip1 = 1.0/nsip2;
            wf2_ip1 = 1.0/nsip1;
#ifdef DEBUG_PP/*{{{*/
            fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
            if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
                fprintf(stderr,"\nline=%d, i=%d, j=%d, pp, gpo not equal, prof1[27]=%f, prof2[27]=%f, prof1[91]=%f, prof2[91]=%f, wf1_i=%f, wf2_i=%f, wf1_ip1=%f, wf2_ip1=%f\n",__LINE__, idx1, idx2 , prof1[27], prof2[27], prof1[91],prof2[91], wf1_i, wf2_i, wf1_ip1, wf2_ip1);
            }
            if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
            }
            if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
                fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
            }
#endif/*}}}*/

            /*pa = MAX3(pa,pga + prof2[91],        pgb + prof1[91]);*/
            /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
            pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);

            prof2 += 32;
            for (c = 1;c < freq[0];c++){
                pa += prof1[freq[c]]*prof2[freq[c]];
            }
            prof2 -= 32;

            s[j].a = pa;

            pga = s[j].ga;

            //s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
            /*s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);*/
            /*s[j].ga = MAX(xga+prof1[28], xa+prof1[27]);*/
            s[j].ga = MAX(xga+prof1[28]*wf1_i, xa+prof1[27]*wf1_i);

            pgb = s[j].gb;

            /*s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(pgb+prof2[28], ca+prof2[27]);*/
            s[j].gb = MAX(pgb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);

            pa = ca;
            xa = s[j].a;
            xga = s[j].ga;
        }
        prof2 -= 64;
		idx2 -= 1;
        ca = s[j].a;

        wf1_i = 1.0/nsip2;
        wf2_i = 1.0/nsip1;
        wf1_ip1 = 1.0/nsip2;
        wf2_ip1 = 1.0/nsip1;

        /*pa = MAX3(pa,pga + prof2[91],        pgb + prof1[91]);*/
        /*pa = MAX3(pa,pga + prof1[91],pgb + prof2[91]);*/
        pa = MAX3(pa,pga + prof1[91]*wf1_ip1,pgb + prof2[91]*wf2_ip1);
        prof2 += 32;
        for (c = 1;c < freq[0];c++){
            pa += prof1[freq[c]]*prof2[freq[c]];
        }
        prof2 -= 32;
        s[j].a = pa;

        //pga = s[j].ga;
        s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

        //pgb = s[j].gb;
        if(hm->startb){
            /*s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);*/
            /*s[j].gb = MAX(s[j].gb+prof2[28], ca+prof2[27]);*/
            s[j].gb = MAX(s[j].gb+prof2[28]*wf2_i, ca+prof2[27]*wf2_i);
        }else{
            /*s[j].gb = MAX(s[j].gb,ca)+prof1[29];*/
            /*s[j].gb = MAX(s[j].gb,ca)+prof2[29];*/
            s[j].gb = MAX(s[j].gb,ca)+prof2[29]*wf2_i;
        }
#ifdef DEBUG_PP/*{{{*/
        fz_weight = (nResArray1[idx1]+nResArray2[idx2]) / 2.0;
        if (fabs(prof2[27]*wf2_i -prof1[27] ) > FLOAT_ZERO * fz_weight || fabs(prof1[27]*wf1_i-prof2[27]) > FLOAT_ZERO * fz_weight || fabs(prof1[91]*wf1_ip1 - prof2[91])> FLOAT_ZERO * fz_weight || fabs(prof2[91]*wf2_ip1-prof1[91])> FLOAT_ZERO * fz_weight){
            fprintf(stderr,"line=%d, i=%d, j=%d, pp, gpo not equal, prof2[27]*wf2_i =%f, prof1[27] =%f\n",__LINE__, idx1, idx2 , prof2[27]*wf2_i, prof1[27]);
        }
        if (fabs(prof2[28]*wf2_i - prof1[28])>FLOAT_ZERO * fz_weight ){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n", __LINE__, __FUNCTION__);
        }
        if (fabs(prof2[29]*wf2_i - prof1[29])>FLOAT_ZERO * fz_weight){
            fprintf(stderr, "line:%d, %s: prof2 and prof1 not equal\n",__LINE__,  __FUNCTION__);
        }
#endif/*}}}*/

        //pa = ca;
    }		
    return s;
}/*}}}*/

/*other functions*/
int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b)/*{{{*/
{
    int* np = 0;

    int i;
    np =malloc(sizeof(int)*(len_a+2));
    for(i =0; i < len_a+2;i++){
        np[i] = -1;
    }

    for(i = 1; i <= len_b;i++){
        if(hirsch_path[i] != -1){
            np[hirsch_path[i]] = i;
        }
    }

    free(hirsch_path);
    return np;
}/*}}}*/

int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b)/*{{{*/
{
    int i,j;
    int a = 0;
    int b = 0;

    int* np = 0;
    np =malloc(sizeof(int)*(len_a+len_b+2));
    for(i =0; i < len_a+len_b+2;i++){
        np[i] = 0;
    }

    j = 1;
    b = -1;
    if(hirsch_path[1] == -1){
        np[j] = 2;
        j++;
    }else{
        if(hirsch_path[1] != 1){
            for ( a = 0;a < hirsch_path[1] -1;a++){
                np[j] = 1;
                j++;
            }
            np[j] = 0;
            j++;
        }else{
            np[j] = 0;
            j++;
        }
    }
    b = hirsch_path[1];

    /*for ( i= 0;i <= len_a;i++){
        fprintf(stderr,"%d,",hirsch_path[i]);
        } 
        fprintf(stderr,"\n");*/

    for(i = 2; i <= len_a;i++){

        if(hirsch_path[i] == -1){
            np[j] = 2;
            j++;
        }else{
            if(hirsch_path[i]-1 != b && b != -1){
                for ( a = 0;a < hirsch_path[i] - b-1;a++){
                    np[j] = 1;
                    j++;
                }
                np[j] = 0;
                j++;
            }else{
                np[j] = 0;
                j++;
            }
        }
        b = hirsch_path[i];
    }





    if(hirsch_path[len_a] < len_b && hirsch_path[len_a] != -1){
        //	fprintf(stderr,"WARNING:%d	%d\n",hirsch_path[len_a],len_b);
        for ( a = 0;a < len_b - hirsch_path[len_a];a++){
            np[j] = 1;
            j++;
        }
    } 
    np[0] = j-1;
    np[j] = 3;
    np = realloc(np,sizeof(int)* (np[0]+2));
    //for ( i= 0;i <= np[0];i++){
    //	fprintf(stderr,"%d,",np[i]);
    //} 
    //fprintf(stderr,"\n");

    free(hirsch_path);

    //add gap info..
    i = 2;
    while(np[i] != 3){
        if ((np[i-1] &3) && !(np[i] & 3)){
            if(np[i-1] & 8){
                np[i-1] += 8;
            }else{
                np[i-1] |= 16;
            }
        }else if (!(np[i-1] & 3) &&(np[i] &3)){
            np[i] |= 4;
        }else if ((np[i-1] & 1) && (np[i] & 1)){
            np[i] |= 8;
        }else if ((np[i-1] & 2) && (np[i] & 2)){
            np[i] |= 8;
        }
        i++;
    }
    //add terminal gap...
    i = 1;
    while(np[i] != 0){
        np[i] |= 32;
        i++;
    }
    j = i;
    i = np[0];
    while(np[i] != 0){
        np[i] |= 32;
        i--;
    }
    //for ( i= 0;i <= np[0];i++){
    //	fprintf(stderr,"%d,",np[i]);
    //} 
    //fprintf(stderr,"\n");
    return np;
}/*}}}*/

/*{{{*//*deleted code*/
/*
int* foward_pp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b)
{
unsigned int freq[26];

struct states* s = 0;
char** trace = 0;
char* tracep = 0;
register int pa = 0;
register int pga = 0;
register int pgb = 0;
register int ca = 0;
register int i = 0;
register int j = 0;
register int c = 0;

s = dp->s;

trace = dp->tb;

trace[0][0] = 32;


s[0].a = 0;
s[0].ga = -INFTY;
s[0].gb = -INFTY;
//init of first row;
tracep = trace[0];

for (j = 1; j < len_b;j++){
    s[j].a = -INFTY;
    
    s[j].ga = s[j-1].a+prof2[29];
    if (s[j-1].ga+prof2[29] > s[j].ga){
        s[j].ga = s[j-1].ga+prof2[29];
    }
    s[j].gb = -INFTY;
    tracep[j] = 8;
}

s[len_b].a = -INFTY;
s[len_b].ga = -INFTY;
s[len_b].gb = -INFTY;

for ( i = 1;i <len_a;i++){
    prof1 += 64;

    c = 1;
    for (j = 26; j--;){
        if(prof1[j]){
            freq[c] = j;
            c++;	
        }
    }
    freq[0] = c;
    
    tracep = trace[i];
    pa = s[0].a;
    pga = s[0].ga;
    pgb = s[0].gb;
    s[0].a = -INFTY;
    s[0].ga = -INFTY;
    
    s[0].gb = pa+prof1[29];
    if(pgb+prof1[29] > s[0].gb){
        s[0].gb = pgb+prof1[29];
    }

    tracep[0] = 16;

    for (j = 1; j < len_b;j++){
        prof2 += 64;
        ca = s[j].a;

        c = 1;
        if((pga += prof2[-37]) > pa){
            pa = pga;
            c = 2;
        }
        if((pgb += prof1[-37]) > pa){
            pa = pgb;
            c = 4;
        }
        
        prof2 += 32;
        for (pga = freq[0];--pga;){
            pgb = freq[pga];
            pa += prof1[pgb]*prof2[pgb];
        }
        prof2 -= 32;

        s[j].a = pa;
        
        pga = s[j].ga;
        
        s[j].ga = s[j-1].a+prof2[27];
        if (s[j-1].ga+prof2[28] > s[j].ga){
            s[j].ga = s[j-1].ga+prof2[28];
            c |= 8;
        }
        
        pgb = s[j].gb;
        
        s[j].gb = ca+prof1[27];
        if(pgb+prof1[28] > s[j].gb){
            s[j].gb = pgb+prof1[28];
            c |= 16;
        }
        tracep[j] = c;
        pa = ca;

    }
    

    prof2 += 64;
    //LAST CELL (0)
    ca = s[len_b].a;

    c = 1;
    if((pga+=prof2[-37]) > pa){
        pa = pga;
        c = 2;
    }
    if((pgb+=prof1[-37]) > pa){
        pa = pgb;
        c = 4;
    }
    
    prof2 += 32;
    for (pga = freq[0];--pga;){
        pgb = freq[pga];
        pa += prof1[pgb]*prof2[pgb];
    }
    prof2 -= 32;
    
    s[len_b].a = pa;
    
    s[len_b].ga = -INFTY;
    
    pgb = s[len_b].gb;
    s[len_b].gb = ca+prof1[27]+prof1[29];
    if(pgb+prof1[29] > s[len_b].gb){
        s[len_b].gb = pgb+prof1[29];
        c |= 16;
    }
    tracep[len_b] = c;	
    prof2 -= len_b << 6;
    
}
prof1 += 64;

c = 1;
for (j = 26; j--;){
    if(prof1[j]){
        freq[c] = j;
        c++;	
    }
}
freq[0] = c;

tracep = trace[len_a];

pa = s[0].a;
pga = s[0].ga;
pgb = s[0].gb;
s[0].a = -INFTY;
s[0].ga = -INFTY;

s[0].gb = pa+prof1[29];
if(pgb+prof1[29] > s[0].gb){
    s[0].gb = pgb+prof1[29];
}
tracep[0] = 16;

for (j = 1;j< len_b;j++){	

    prof2 += 64;
    ca = s[j].a;

    c = 1;

    if((pga+=prof2[-37]) > pa){
        pa = pga;
        c = 2;
    }

    if((pgb+=prof1[-37]) > pa){
        pa = pgb;
        c = 4;
    }
    
    prof2+=32;
    
    for (pga = freq[0];--pga;){
        pgb = freq[pga];
        pa += prof1[pgb]*prof2[pgb];
    }
    prof2-=32;
    
    s[j].a = pa;
    pga = s[j].ga;
    s[j].ga = s[j-1].a+prof2[27]+prof2[29];
    if (s[j-1].ga+prof2[29] > s[j].ga){
        s[j].ga = s[j-1].ga+prof2[29];
        c |= 8;
    }	
    pgb = s[j].gb;
    s[j].gb = -INFTY;	
    
    tracep[j] = c;
    pa = ca;
}
prof2 += 64;

ca = s[len_b].a;

c = 1;

if((pga+=prof2[-37]) > pa){
    pa = pga;
    c = 2;
}
if((pgb+=prof1[-37]) > pa){
    pa = pgb;
    c = 4;
}
prof2+=32;
for (pga = freq[0];--pga;){	
    pgb = freq[pga];
    pa += prof1[pgb]*prof2[pgb];
}
prof2-=32;

s[len_b].a = pa;

s[len_b].ga = s[len_b-1].a+prof2[27]+prof2[29];
if (s[len_b-1].ga+prof2[29] > s[len_b].ga){
    s[len_b].ga = s[len_b-1].ga+prof2[29];
    c |= 8;
}

pgb = s[len_b].gb;
s[len_b].gb = ca+prof1[27]+prof1[29];
if(pgb +prof1[29]> s[len_b].gb){
    s[len_b].gb = pgb+prof1[29];
    c |= 16;
}	
tracep[len_b] = c;

pgb = s[len_b].gb;
c = 2;
if(s[len_b].ga > pgb){
    pgb = s[len_b].ga;
    c = 1;
}
if(s[len_b].a >= pgb){
    pgb = s[len_b].a;
    c = 0;
}

ca = c;

i = len_a;
j = len_b;
c = 1;
while(trace[i][j] < 32){
//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
    switch(ca){
        case 0:
            if (trace[i][j] & 2){
                ca = 1;
                if(i-1!= 0){
                    path[c+1] |= 16;
//					fprintf(stderr,"GAP_CLOSE\n");
                }else{
                    path[c+1] |= 32+16;
                }
            }else if (trace[i][j] & 4){
                ca = 2;
                if(j-1!= 0){
                    path[c+1] |= 16;
//					fprintf(stderr,"GAP_CLOSE\n");
                }else{
                    path[c+1] |= 32+16;
                }
            }

            //path[c] = 0;
            i--;
            j--;
        break;
        case 1:
            if(trace[i][j] & 8){
                ca = 1;
                if(i!=0 && i!= len_a){
//				/	fprintf(stderr,"GAP_EXT\n");
                    if(!(path[c]&16)){
                        path[c] |= 8;
                    }
                }else{
                    if(!(path[c]&16)){
                        path[c] |= 32+8;
                    }
                }
            }else{
                ca = 0;
                if(i!=0 && i!= len_a){
//					fprintf(stderr,"GAP_OPEN\n");
                    path[c] |= 4;
                }else{
                    path[c] |= 32+4;
                }
            }
            path[c] |= 1;
            j--;
        break;
        case  2:
            if(trace[i][j] & 16){
                ca = 2;
                if(j !=0 && j != len_b){
//					fprintf(stderr,"GAP_EXT\n");
                    if(!(path[c]&16)){
                        path[c] |= 8;
                    }
                }else{
                    if(!(path[c]&16)){
                        path[c] |= 32+8;
                    }
                }
            }else{
                ca = 0;
                if(j !=0 && j != len_b){
//					fprintf(stderr,"GAP_OPEN\n");
                    path[c] |= 4;
                }else{
                    path[c] |= 32+4;
                }
                
            }
            path[c] |= 2;
            i--;
        break;
    }
    c++;
}



path[0] = c-1;
path[c] = 3;
path[c+1] = pgb;

j = path[0];
for(i =0 ;i < path[0]/2;i++){
    c = path[i+1];
    path[i+1] = path[j-i];
    path[j -i] = c;
}
return path;
}

int* backward_pp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b)
{
unsigned int freq[26];

struct states* s = 0;
char** trace = 0;
char* tracep = 0;
register int pa = 0;
register int pga = 0;
register int pgb = 0;
register int ca = 0;
register int i = 0;
register int j = 0;
register int c = 0;

prof1+= 64;
prof2 += 64;


s = dp->s;

trace = dp->tb;

trace[len_a][len_b] = 32;

prof1 +=  len_a << 6;

s[len_b].a = 0;
s[len_b].ga = -INFTY;
s[len_b].gb = -INFTY;
//init of first row;
tracep = trace[len_a];

j = len_b;
while(--j){
    s[j].a = -INFTY;
    
    s[j].ga = s[j+1].a+prof2[29];
    if (s[j+1].ga+prof2[29] > s[j].ga){
        s[j].ga = s[j+1].ga+prof2[29];
    }
    s[j].gb = -INFTY;
    tracep[j] = 8;
}

s[0].a = -INFTY;
s[0].ga = -INFTY;
s[0].gb = -INFTY;
i = len_a;
while(--i){
    prof1 -= 64;

    c = 1;
    for (j = 26; j--;){
        if(prof1[j]){
            freq[c] = j;
            c++;	
        }
    }
    freq[0] = c;
    
    tracep = trace[i];
    pa = s[len_b].a;
    pga = s[len_b].ga;
    pgb = s[len_b].gb;
    s[len_b].a = -INFTY;
    s[len_b].ga = -INFTY;
    
    s[len_b].gb = pa+prof1[29];
    if(pgb+prof1[29] > s[len_b].gb){
        s[len_b].gb = pgb+prof1[29];
    }

    tracep[len_b] = 16;
    
    j = len_b;
    prof2 += len_b << 6;
    while(--j){
        prof2 -= 64;
        ca = s[j].a;

        c = 1;
        if((pga += prof2[91]) > pa){
            pa = pga;
            c = 2;
        }
        if((pgb += prof1[91]) > pa){
            pa = pgb;
            c = 4;
        }
        
        prof2 += 32;
        for (pga = freq[0];--pga;){
            pgb = freq[pga];
            pa += prof1[pgb]*prof2[pgb];
        }
        prof2 -= 32;

        s[j].a = pa;
        
        pga = s[j].ga;
        
        s[j].ga = s[j+1].a+prof2[27];
        if (s[j+1].ga+prof2[28] > s[j].ga){
            s[j].ga = s[j+1].ga+prof2[28];
            c |= 8;
        }
        
        pgb = s[j].gb;
        
        s[j].gb = ca+prof1[27];
        if(pgb+prof1[28] > s[j].gb){
            s[j].gb = pgb+prof1[28];
            c |= 16;
        }
        tracep[j] = c;
        pa = ca;

    }

    prof2 -= 64;
    //LAST CELL (0)
    ca = s[0].a;

    c = 1;
    if((pga+=prof2[91]) > pa){
        pa = pga;
        c = 2;
    }
    if((pgb+=prof1[91]) > pa){
        pa = pgb;
        c = 4;
    }
    
    prof2 += 32;
    for (pga = freq[0];--pga;){
        pgb = freq[pga];
        pa += prof1[pgb]*prof2[pgb];
    }
    prof2 -= 32;
    
    s[0].a = pa;
    
    s[0].ga = -INFTY;
    
    pgb = s[0].gb;
    s[0].gb = ca+prof1[27]+prof1[29];
    if(pgb+prof1[29] > s[0].gb){
        s[0].gb = pgb+prof1[29];
        c |= 16;
    }
    tracep[0] = c;	
    
}
prof1 -= 64;

c = 1;
for (j = 26; j--;){
    if(prof1[j]){
        freq[c] = j;
        c++;	
    }
}
freq[0] = c;

tracep = trace[0];
j = len_b;
prof2 += len_b << 6;
pa = s[j].a;
pga = s[j].ga;
pgb = s[j].gb;
s[j].a = -INFTY;
s[j].ga = -INFTY;

s[len_b].gb = pa+prof1[29];
if(pgb+prof1[29] > s[len_b].gb){
    s[len_b].gb = pgb+prof1[29];
}



while(--j){
    prof2 -= 64;
    ca = s[j].a;

    c = 1;

    if((pga+=prof2[91]) > pa){
        pa = pga;
        c = 2;
    }

    if((pgb+=prof1[91]) > pa){
        pa = pgb;
        c = 4;
    }
    
    prof2+=32;
    
    for (pga = freq[0];--pga;){
        pgb = freq[pga];
        pa += prof1[pgb]*prof2[pgb];
    }
    prof2-=32;
    
    s[j].a = pa;
    pga = s[j].ga;
    s[j].ga = s[j+1].a+prof2[27]+prof2[29];
    if (s[j+1].ga+prof2[29] > s[j].ga){
        s[j].ga = s[j+1].ga+prof2[29];
        c |= 8;
    }	
    pgb = s[j].gb;
    s[j].gb = -INFTY;	
    
    tracep[j] = c;
    pa = ca;
}
prof2 -= 64;

ca = s[0].a;

c = 1;

if((pga+=prof2[91]) > pa){
    pa = pga;
    c = 2;
}
if((pgb+=prof1[91]) > pa){
    pa = pgb;
    c = 4;
}
prof2+=32;
for (pga = freq[0];--pga;){	
    pgb = freq[pga];
    pa += prof1[pgb]*prof2[pgb];
}
prof2-=32;

s[0].a = pa;

s[0].ga = s[1].a+prof2[27]+prof2[29];
if (s[1].ga+prof2[29] > s[0].ga){
    s[0].ga = s[1].ga+prof2[29];
    c |= 8;
}

pgb = s[0].gb;
s[0].gb = ca+prof1[27]+prof1[29];
if(pgb +prof1[29]> s[0].gb){
    s[0].gb = pgb+prof1[29];
    c |= 16;
}	
tracep[0] = c;

pgb = s[0].gb;
c = 2;
if(s[0].ga > pgb){
    pgb = s[0].ga;
    c = 1;
}
if(s[0].a >= pgb){
    pgb = s[0].a;
    c = 0;
}

//fprintf(stderr,"SCORE:%d\n",ca);
ca = c;

i = 0;
j = 0;
c = 1;
while(trace[i][j] < 32){
//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
    switch(ca){
        case 0:
            if (trace[i][j] & 2){
                ca = 1;
                if(i+1!= len_a){
                    path[c+1] |= 16;
//					fprintf(stderr,"GAP_CLOSE\n");
                }else{
                    path[c+1] |= 32+16;
                }
            }else if (trace[i][j] & 4){
                ca = 2;
                if(j+1!= len_b){
                    path[c+1] |= 16;
//					fprintf(stderr,"GAP_CLOSE\n");
                }else{
                    path[c+1] |= 32+16;
                }
            }

            //path[c] = 0;
            i++;
            j++;
        break;
        case 1:
            if(trace[i][j] & 8){
                ca = 1;
                if(i!=0 && i!= len_a){
//				/	fprintf(stderr,"GAP_EXT\n");
                    if(!(path[c]&16)){
                        path[c] |= 8;
                    }
                }else{
                    if(!(path[c]&16)){
                        path[c] |= 32+8;
                    }
                }
            }else{
                ca = 0;
                if(i!=0 && i!= len_a){
//					fprintf(stderr,"GAP_OPEN\n");
                    path[c] |= 4;
                }else{
                    path[c] |= 32+4;
                }
            }
            path[c] |= 1;
            j++;
        break;
        case  2:
            if(trace[i][j] & 16){
                ca = 2;
                if(j !=0 && j != len_b){
//					fprintf(stderr,"GAP_EXT\n");
                    if(!(path[c]&16)){
                        path[c] |= 8;
                    }
                }else{
                    if(!(path[c]&16)){
                        path[c] |= 32+8;
                    }
                }
            }else{
                ca = 0;
                if(j!=0 && j != len_b){
//					fprintf(stderr,"GAP_OPEN\n");
                    path[c] |= 4;
                }else{
                    path[c] |= 32+4;
                }
                
            }
            path[c] |= 2;
            i++;
        break;
    }
    c++;
}
path[0] = c-1;
path[c] = 3;
path[c+1] = pgb;


return path;
}

*//*}}}*/
