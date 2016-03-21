/*
	kalignP_main.c 
	
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



unsigned int numseq = 0;
unsigned int numprofiles = 0;  /*numprofiles = numseq << 1 - 1*/
/*for amino acid sequence, by default*/
/*gpo = (float)54.94941;*/
/*gpe = (float)8.52492;*/
/*tgpe = (float)4.42410;*/
float gpo = (float)0.0;
float gpe = (float)0.0;
float tgpe = (float)0.0;
float FLOAT_ZERO = 1e-2;

int main(int argc,char **argv)
{
    int i;
    int j;
    int* tree = 0;
    int a, b, c;

    struct alignment* aln = 0;
    struct parameters* param = 0;
    struct aln_tree_node* tree2 = 0;

    param = malloc(sizeof(struct parameters));

    param =  interface(param,argc,argv);
//     fprintf(stdout,"seq_type=%d\n", param->seq_type); /*debug*/
    fprintf(stdout,"custom_sub_matrix_file=%s\n", param->custom_sub_matrix_file); /*debug*/
	
	/*aln = detect_and_read_sequences(aln,param);*/
	aln = detect_and_read_sequences_new(aln,param);



#ifdef DEBUG/*{{{*/
    /*FILE *fpLog = fopen("logseq.tmp", "w");*/
    FILE *fpLog = stderr;
    fprintf(fpLog,"\nGap penalties read in\n\n");
    for (i = 0 ; i < numseq;i ++){
        fprintf(fpLog,">%s\n", aln->sn[i]);
        fprintf(fpLog,"%s\n", aln->seq[i]);
        if(aln->gpo[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{gpo: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->gpo[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"gpo[0]=%g, gpo[%d]=%g\n", aln->gpo[i][0],aln->sl[i]+1, aln->gpo[i][aln->sl[i]+1]);
        }
        if(aln->gpe[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{gpe: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->gpe[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"gpe[0]=%g, gpe[%d]=%g\n", aln->gpe[i][0],aln->sl[i]+1, aln->gpe[i][aln->sl[i]+1]);
        }
        if(aln->tgpe[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{tgpe: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->tgpe[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"tgpe[0]=%g, tgpe[%d]=%g\n", aln->tgpe[i][0],aln->sl[i]+1, aln->tgpe[i][aln->sl[i]+1]);
        }
    }
    if (fpLog != stdout && fpLog != stderr){
        fclose(fpLog);
    }
    fprintf(stderr,"readin sequence output to logseq.tmp\n");
#endif/*}}}*/

	if(param->ntree > numseq){
		param->ntree = numseq;
	}

	//DETECT DNA
    if(param->dna == -1){
        for (i = 0; i < numseq;i++){
            param->dna = byg_detect(aln->s[i],aln->sl[i]);
            if(param->dna){
                break;
            }
        }
    }
	//param->dna = 0;
	//fprintf(stderr,"DNA:%d\n",param->dna);
	//exit(0);
	
	if(param->dna == 1){
		//brief sanity check...
		for (i = 0; i < numseq;i++){
			if(aln->sl[i] < 6){
				fprintf(stderr,"Dna/Rna alignments are only supported for sequences longer than 6.");
				free(param);
				free_aln(aln);
				exit(0);
			}
		}
		aln =  make_dna(aln);
	}

	
	if(param->reformat){
		for (i = 0 ;i < numseq;i++){
			aln->nsip[i] = i;
			for (j = 0; j < aln->sl[i];j++){
				aln->s[i][j] = 0;
			}
		}
		param->format = "fasta";//param->reformat;
		output(aln,param);
		exit(1);
	}
	
	
	
	//fast distance calculation;
	float** submatrix = 0;
// 	submatrix = read_matrix(submatrix,param); // sets gap penalties as well.....
	submatrix = read_matrix_new(submatrix,param); /* modified 2014-11-27 by Nanjiang, accept also customer substitution matrix*/

/*debug print the matrix*/
    fprintf(stdout,"\n\n Substitution matrix:\n");
    for (i=0;i<23;i++){
        for (j=0;j<23;j++){
            fprintf(stdout,"%3.0f ", submatrix[i][j]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n\n\n");
	
	if(!param->quiet){
		parameter_message(param);
	}
	
	if(byg_start(param->alignment_type,"profPROFprofilePROFILE") != -1){
		profile_alignment_main(aln,param,submatrix);
	}
	
	float** dm = 0;
    /*debug*/
    /*print para*/
    /*debug*/
	if(param->ntree > 1)
    {
		if(byg_start(param->distance,"pairclustalPAIRCLUSTAL") != -1)
        {
			if(byg_start(param->tree,"njNJ") != -1)
            {
				dm = protein_pairwise_alignment_distance(aln,dm,param,submatrix,1);
			}else
            {
				dm = protein_pairwise_alignment_distance(aln,dm,param,submatrix,0);
			}
		}
        else if(byg_start("wu",param->alignment_type) != -1)
        {
			dm =  protein_wu_distance2(aln,dm,param);
		//	param->feature_type = "wumanber";
		}
        else if(param->dna == 1)
        {
			if(byg_start(param->tree,"njNJ") != -1)
            {
				dm =  dna_distance(aln,dm,param,1);
			}else{
				dm =  dna_distance(aln,dm,param,0);
			}
		}else{
			if(byg_start(param->tree,"njNJ") != -1)
            {
				dm =  protein_wu_distance(aln,dm,param,1);
			}
            else
            {
				dm =  protein_wu_distance(aln,dm,param,0); /*this function will be run by default*/
			}
		}
		/*
		for (i = 0; i< numseq;i++){
			for (j = 0; j< numseq;j++){
				fprintf(stderr,"%f	",dm[i][j]);
			}
			fprintf(stderr,"\n");
		}*/
#ifdef DEBUG
        /*print the distance matrix*/
        fprintf(stderr,"\nThe dm matrix of the size %d x %d:\n",numseq, numseq);
        int ii1, jj1;
        for (ii1=0;ii1 < numseq; ii1++ ) {
            for (jj1=0;jj1 < numseq; jj1++ ) {
                fprintf(stderr,"%6g ", dm[ii1][jj1]);
            }
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"\n");
#endif

		if(byg_start(param->tree,"njNJ") != -1){
			tree2 = real_nj(dm,param->ntree);
		}else{
			tree2 = real_upgma(dm,param->ntree);
		}
		if(param->print_tree){
			print_tree(tree2,aln,param->print_tree);
		}
	}

	tree = malloc(sizeof(int)*(numseq*3+1));
	for ( i = 1; i < (numseq*3)+1;i++){
		tree[i] = 0;
	}
	tree[0] = 1; 
	
	if(param->ntree < 2){
		tree[0] = 0;
		tree[1] = 1;
		
		c = numseq;
		tree[2] = c;
		a = 2;
		for ( i = 3; i < (numseq-1)*3;i+=3){
			tree[i] = c;
			tree[i+1] = a;
			c++;
			tree[i+2] = c;
			a++;
		}
	}else if(param->ntree > 2){
		ntreeify(tree2,param->ntree);
	}else{
		tree = readtree(tree2,tree);
		for (i = 0; i < (numseq*3);i++){
			tree[i] = tree[i+1];
		}
		free(tree2->links);
		free(tree2->internal_lables);
		free(tree2);
	}

	

    /*add to gpo gpe tgpe array*//*{{{*/
    for (i = 0 ; i < numseq;i ++){
        j=0;
        if(!aln->gpo[i]){
            aln->gpo[i] = malloc(sizeof(float) * (aln->sl[i]+2));
            for (j = 0; j < aln->sl[i]+2; j ++) {
                aln->gpo[i][j] = 1.0;
            }
        } else{
            aln->gpo[i][0] = aln->gpo[i][aln->sl[i]+1] = 1.0; /*reset the value of gpo [0] and gpo[length+1] to default, note that when the function detect_and_read_sequences_new is run, gpo, gpe and tgpe haven't been set, fixed 2010-10-29*/
        }

        if(!aln->gpe[i]){
            aln->gpe[i] = malloc(sizeof(float) *(aln->sl[i]+2));
            for (j = 0; j < aln->sl[i]+2; j ++) {
                aln->gpe[i][j] = 1.0;
            }
        }  else{
            aln->gpe[i][0] = aln->gpe[i][aln->sl[i]+1] = 1.0; /*reset the value of gpo [0] and gpo[length+1] to default, note that when the function detect_and_read_sequences_new is run, gpo, gpe and tgpe haven't been set*/
        }  
        if(!aln->tgpe[i]){
            aln->tgpe[i] = malloc(sizeof(float) *(aln->sl[i]+2));
            for (j = 0; j < aln->sl[i]+2; j ++) {
                aln->tgpe[i][j] = 1.0;
            }
        }  else{
            aln->tgpe[i][0] = aln->tgpe[i][aln->sl[i]+1] = 1.0; /*reset the value of gpo [0] and gpo[length+1] to default, note that when the function detect_and_read_sequences_new is run, gpo, gpe and tgpe haven't been set*/
        } 
    }/*}}}*/
#ifdef DEBUG/*{{{*/
    fpLog = stderr;
    fprintf(fpLog,"\nGap penalties after filled in with defaults\n\n");
    /*fprintf(fpLog,"\nAfter added gap arrays\n");*/
    for (i = 0 ; i < numseq;i ++){
        fprintf(fpLog,">%s\n", aln->sn[i]);
        fprintf(fpLog,"%s\n", aln->seq[i]);
        if(aln->gpo[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{gpo: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->gpo[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"gpo[0]=%g, gpo[%d]=%g\n", aln->gpo[i][0],aln->sl[i]+1, aln->gpo[i][aln->sl[i]+1]);
        }
        if(aln->gpe[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{gpe: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->gpe[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"gpe[0]=%g, gpe[%d]=%g\n", aln->gpe[i][0],aln->sl[i]+1, aln->gpe[i][aln->sl[i]+1]);
        }
        if(aln->tgpe[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{tgpe: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->tgpe[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"tgpe[0]=%g, tgpe[%d]=%g\n", aln->tgpe[i][0],aln->sl[i]+1, aln->tgpe[i][aln->sl[i]+1]);
        }
    }
    if (fpLog != stdout && fpLog != stderr){
        fclose(fpLog);
    }
    fprintf(stderr,"readin sequence output to logseq.tmp2\n");
#endif/*}}}*/
    for (i = 0 ; i < numseq;i ++){
        j=0;
        for (j = 0; j < aln->sl[i]+2; j ++) {
            aln->gpo[i][j] *= gpo;
            aln->gpe[i][j] *= gpe;
            aln->tgpe[i][j] *= tgpe;
        }
    }

#ifdef DEBUG/*{{{*/
    fpLog = stderr;
    fprintf(fpLog,"\n\n\nGap penalties multipled with gpo (%f), gpe (%f) and tgpe (%f)\n\n", gpo, gpe, tgpe);
    /*fprintf(fpLog,"\nAfter added gap arrays\n");*/
    for (i = 0 ; i < numseq;i ++){
        fprintf(fpLog,">%s\n", aln->sn[i]);
        fprintf(fpLog,"%s\n", aln->seq[i]);
        if(aln->gpo[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{gpo: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->gpo[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"gpo[0]=%g, gpo[%d]=%g\n", aln->gpo[i][0],aln->sl[i]+1, aln->gpo[i][aln->sl[i]+1]);
        }
        if(aln->gpe[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{gpe: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->gpe[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"gpe[0]=%g, gpe[%d]=%g\n", aln->gpe[i][0],aln->sl[i]+1, aln->gpe[i][aln->sl[i]+1]);
        }
        if(aln->tgpe[i]){
            j = 0; 
            fprintf(fpLog,"%s", "{tgpe: ");
            for (j = 1;j <= aln->sl[i]; j ++ ) {
                fprintf(fpLog,"%g ", aln->tgpe[i][j]);
            }
            fprintf(fpLog,"%s\n", "}");
            fprintf(fpLog,"tgpe[0]=%g, tgpe[%d]=%g\n", aln->tgpe[i][0],aln->sl[i]+1, aln->tgpe[i][aln->sl[i]+1]);
        }
    }
    if (fpLog != stdout && fpLog != stderr){
        fclose(fpLog);
    }
    fprintf(stderr,"readin sequence output to logseq.tmp2\n");
#endif/*}}}*/

	//get matrices... 
	struct feature_matrix* fm = 0;
	
	struct ntree_data* ntree_data = 0;
	
	int** map = 0;
	if(param->ntree > 2){
		ntree_data = malloc(sizeof(struct ntree_data));
		ntree_data->realtree = tree2;
		ntree_data->aln = aln;
		ntree_data->profile = 0;
		ntree_data->map = 0;
		ntree_data->ntree = param->ntree;
		ntree_data->submatrix = submatrix;
		ntree_data->tree = tree; 
		
		ntree_data = ntree_alignment(ntree_data);
		map = ntree_data->map;
		tree = ntree_data->tree;
		for (i = 0; i < (numseq*3);i++){
			tree[i] = tree[i+1];
		}
		free(ntree_data);
	}else if (param->feature_type){
		fm = get_feature_matrix(fm,aln,param);
		if(!fm){
			for (i = 32;i--;){
				free(submatrix[i]);
			}
			free(submatrix);
			free_param(param);
			free(map);
			free(tree);
			exit(0);
		}
		
		map = feature_hirschberg_alignment(aln,tree,submatrix,map,fm);
		//exit(0);
		//map =  feature_alignment(aln,tree,submatrix, map,fm);

	}else if (byg_start("pairwise",param->alignment_type) != -1){
		if(param->dna == 1){
			map = dna_alignment_against_a(aln,tree,submatrix, map,param->gap_inc);
		}else{
			map = hirschberg_alignment_against_a(aln,tree,submatrix, map,param->smooth_window,param->gap_inc);
		}
		//map =  default_alignment(aln,tree,submatrix, map);
	}else if (byg_start("fast",param->alignment_type) != -1){
		map =  default_alignment(aln,tree,submatrix, map);
	}else if(param->dna == 1){
		map =  dna_alignment(aln,tree,submatrix, map,param->gap_inc);
	/*}else if (byg_start("test",param->alignment_type) != -1){
		map =  test_alignment(aln,tree,submatrix, map,param->internal_gap_weight,param->smooth_window,param->gap_inc);
	}else if (param->aa){
		map =  aa_alignment(aln,tree,submatrix, map,param->aa);
	}else if (param->alter_gaps){
		map = alter_gaps_alignment(aln,tree,submatrix,map,param->alter_gaps,param->alter_range,param->alter_weight);
	}else if (byg_start("altergaps",param->alignment_type) != -1){
		map = alter_gaps_alignment(aln,tree,submatrix,map,param->alter_gaps,param->alter_range,param->alter_weight);
	}else if(byg_start("simple",param->alignment_type) != -1){
		map =  simple_hirschberg_alignment(aln,tree,submatrix, map);*/
	}else if(byg_start("advanced",param->alignment_type) != -1){
		map =  advanced_hirschberg_alignment(aln,tree,submatrix, map,param->smooth_window,param->gap_inc,param->internal_gap_weight);
	}else{ /*by default, this will be run*/
		map =  hirschberg_alignment(aln,tree,submatrix, map,param->smooth_window,param->gap_inc);
	}

	//clear up sequence array to be reused as gap array....
    /*After that, aln->s[i][j] will have another meaning
     * aln->s[i][j] == 0 means at position j of sequence i, there is no gap, 
     * aln->s[i][j] == n (n>0) means there are n gaps before position j in sequence i*/
	int *p = 0;
	for (i = 0; i < numseq;i++){
		p = aln->s[i];
		for (a = 0; a < aln->sl[i];a++){
			p[a] = 0;
		}
	}
	//clear up
	for (i = 0; i < (numseq-1)*3;i +=3){
		a = tree[i];
		b = tree[i+1];
		aln = make_seq(aln,a,b,map[tree[i+2]]);
	}
#ifdef DEBUG
    fprintf(stderr,"\nAfter make_seq\n");
	for (i = 0; i < numseq;i++){
		fprintf(stderr,"aln->sn[%d]=%s	aln->nsip[%d]=%d\n",i, aln->sn[i],i, aln->nsip[i]);
	}
#endif
	for (i = 0; i < numseq;i++){
		aln->nsip[i] = 0;
	}
	aln =  sort_sequences(aln,tree,param->sort);
#ifdef DEBUG
    fprintf(stderr,"\nAfter sort_sequences\n");
	for (i = 0; i < numseq;i++){
		fprintf(stderr,"i=%d	aln->nsip[%d]=%d	aln->sip[%d][0]=%d\n",i,i, aln->nsip[i],i, aln->sip[i][0]);
	}
#endif

#ifdef DEBUG
    fprintf(stderr,"\nAligned sequence in digit format\n");
	for (i = 0; i < numseq;i++){
		fprintf(stderr,">%s\n",aln->sn[aln->nsip[i]]);
        j = 0;
        for (j = 0; j <= aln->sl[aln->nsip[i]]; j ++){
            fprintf(stderr,"%2d ",aln->s[aln->nsip[i]][j]);
        }
        fprintf(stderr,"\n");
	}
#endif
	output(aln,param);
/*	if(!param->format){
		fasta_output(aln,param->outfile);
	}else{
		if (byg_start("msf",param->format) != -1){
			msf_output(aln,param->outfile);
		}else if (byg_start("clustal",param->format) != -1){
			clustal_output(aln,param->outfile);
		}else if (byg_start("macsim",param->format) != -1){
			macsim_output(aln,param->outfile,param->infile[0]);
		}
	}
	free_param(param);*/
	
	free(map);
	free(tree);
	return 0;
}



