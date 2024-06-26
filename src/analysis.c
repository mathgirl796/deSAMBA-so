/*
 * PBsimTest.c
 *
 *  Created on: 2018-4-11
 *      Author: fenghe
 *      FUNCTION: analysis the result of deSPI3rd and minimap2
 *      USED for pbSim data
 *      input:
 *      	(1) pbSim MAF file, it records the original taxonomy of each reads
 *      	(2) deSPI result file, it records the result of deSPI3rd
 *      	(3) taxonomy file: that is "node.dmp"
 *      output:
 *      	each record separated by '\n', each item separated by '\t'
 *      	(1) read name
 *      	(2) read true taxonomy
 *      	(3) read true taxonomy rank
 *      	(4) read classified taxonomy
 *      	(5) read classified taxonomy rank
 *      	(6) distance between "true taxonomy" and "classified taxonomy"
 *      	(7) LCA tid
 *      	(8) LCA tid rank
 *
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>

#include "lib/desc.h"
#include "lib/utils.h"
#include "lib/kvec.h"

#define READ_NAME_LEN 100

typedef struct {
	uint32_t 	tid;
	char 		read_name[READ_NAME_LEN];
} MAF;//8byte

typedef struct {
	char 		read_name[READ_NAME_LEN];
	char 		isClassify;
	uint32_t 	tid;
	uint32_t 	read_length;
	uint8_t 	MAPQ;
	uint32_t 	score;
}RST;//24byte

typedef struct {
	char 		read_name[READ_NAME_LEN];
	uint32_t 	true_tid;
	char 		true_rank[20];
	uint32_t 	cly_tid;
	char 		cly_rank[20];
	int 		distance;
	uint32_t 	LCA_TID;
	char 		LCA_rank[20];
}MAF_CMP;

typedef struct {
	uint32_t p_tid;
	char rank[20];
}TAXONOMY_rank;

int max_tid_global = 0;
//------------------------------------------TAX TREE-----------------------------------------//
#define MAX_BUFF_LEN 10000000//1000K
//store taxonomy tree in list, with rank info; return max TID
static int taxonTree_rank(const char *taxonomyNodesPath,TAXONOMY_rank ** taxonomyTree_)
{
	TAXONOMY_rank * taxonomyTree;
	FILE *fp = xopen(taxonomyNodesPath,"r");
	char *line = (char*)xmalloc(1024);
	size_t max_l = 1024;
	//step1: get line number
	char *token;
	uint32_t max_tid = 0;
	while(1){
		if(getline(&line,&max_l,fp) <= 0)
			break;
		//get tid
		token = strtok(line,"\t|");
		max_tid = strtoul(token,NULL,10);
	}
	//step2: reset file
	fclose(fp);
	fp = xopen(taxonomyNodesPath,"r");
	//step 3: store
	max_tid += 1000000;
	max_tid_global = max_tid;
	taxonomyTree = (TAXONOMY_rank*)malloc(sizeof(TAXONOMY_rank)*(max_tid + 1));
	for(uint32_t i = 0; i <= max_tid; i++)
		taxonomyTree[i].p_tid = MAX_uint32_t;
	while(1){
		if(getline(&line,&max_l,fp) <= 0)
			break;
		//get tid
		token = strtok(line,"\t|");
		uint32_t tid = strtoul(token,NULL,10);
		//get p_tid
		token = strtok(NULL,"\t|");
		taxonomyTree[tid].p_tid = strtoul(token,NULL,10);
		//get rank
		token = strtok(NULL,"\t|");
		strcpy(taxonomyTree[tid].rank,token);
	}
	//end mark
	taxonomyTree[1].p_tid = 0;//set the root node to be 0(instead of 1)
	strcpy(taxonomyTree[1].rank,"root");
	strcpy(taxonomyTree[0].rank,"CLY_FAIL");

	//step4 : close file
	fclose(fp);
	*taxonomyTree_ = taxonomyTree;
	return max_tid;
}
//
////return distance, and store the LCA in LCA
//static int TID_distance(uint32_t tid1, uint32_t tid2, TAXONOMY_rank * taxonomy, uint32_t *LCA)
//{
//	//get distance
//	int distance1 = 0;//distance between tid1 and LCA
//	int distance2 = 0;//distance between tid1 and LCA
//	while(1)
//	{
//		if(tid1 == 0)//stop when search to the end, it`s a fatal wrong
//			break;
//		//reset distance2
//		distance2 = 0;
//		uint32_t current_tid2 = tid2;
//		while(1)
//		{
//			if(current_tid2 == 0)//stop when search to the end
//				break;
//			if(tid1 == current_tid2){//return when found LCA
//				*LCA = tid1;
//				return distance1 + distance2;
//			}
//			//get parent node
//			current_tid2 = taxonomy[current_tid2].p_tid;
//			//add distance
//			distance2++;
//		}
//		//get parent node
//		tid1 = taxonomy[tid1].p_tid;
//		//add distance
//		distance1++;
//	}
//	if(tid2 == 0) {*LCA = 0; return distance1;}
//	xassert(0,"[TID_distance]:FATAL WRONG, no LCA are found");
//	return 0;
//}

//------------------------------------------RST-----------------------------------------//

//get one line of RST file, store result in RST
//return -1 when reach the end of the file
//buff pool are needed, "char *buff = (char*)malloc(50000);"
static int getOneRST(FILE * RST_file, RST * rst)
{
	size_t max_l = 1024;
	char static_BUFF[1024];
	char *buff = static_BUFF;
	char *tokens;
	if(getline(&buff,&max_l,RST_file) <= 0)
		return -1;
	tokens = strtok(buff,"\t");
	strcpy(rst->read_name,tokens);
	tokens = strtok(NULL,"\t");
	rst->isClassify = tokens[0];
	tokens = strtok(NULL,"\t");
	rst->tid = strtoul(tokens,NULL, 10);
	tokens = strtok(NULL,"\t");
	rst->read_length = strtoul(tokens,NULL, 10);
	tokens = strtok(NULL,"\t");
	if(tokens == 0)
		rst->MAPQ = 0;
	else
		rst->MAPQ = strtoul(tokens,NULL, 10);

	tokens = strtok(NULL,"\t");
	if(tokens == 0)
		rst->score = 0;
	else
		rst->score = strtoul(tokens,NULL, 10);
	return 0;
}
//------------------------------------------SAM-----------------------------------------//
static int getOneSAM(FILE * SAM_file, char *buff, RST * rst)
{
	size_t max_l = MAX_BUFF_LEN;
	char *tokens;
	int read_L = 0;
	read_L = getline(&buff,&max_l,SAM_file);
	if(read_L <= 0)
		return -1;
	//if(buff[0] != 'S' && buff[0] != 'D' && buff[0] != 'E')
	//{
	//	fprintf(stderr, "Not Normal SAM.\n");
	//}
	//get read name
	tokens = strtok(buff,"\t");
	strcpy(rst->read_name,tokens);
	//ignore flag
	tokens = strtok(NULL,"\t");
	//get refNAME
	rst->read_length = 0;
	rst->score = 0;
	tokens = strtok(NULL,"\t");
	if(tokens[0] == '*') {
		rst->isClassify = 'U';
		rst->tid = 0;
		rst->MAPQ = 0;
	}
	else{
		rst->isClassify = 'C';
		char * ref_tokens = tokens;
		//ignore the POS part
		tokens = strtok(NULL,"\t");
		//get MAQ part
		tokens = strtok(NULL,"\t");
		rst->MAPQ = strtoul(tokens,NULL,10);
		//get CIGAR
		tokens = strtok(NULL,"\t");
		char * CIGAR = tokens;
		//get *
		tokens = strtok(NULL,"\t");
		//get 0
		tokens = strtok(NULL,"\t");
		//get 0
		tokens = strtok(NULL,"\t");
		//get SEQ
		tokens = strtok(NULL,"\t");
		//get QUAL
		tokens = strtok(NULL,"\t");
		//with other label
		//get AS
		tokens = strtok(NULL,":");
		//if(tokens != NULL && ((tokens[0] == 'N' && tokens[1] == 'M')) )//code for minimap2
		if(tokens != NULL && ((tokens[0] == 'A' && tokens[1] == 'S') || (tokens[0] == 'N' && tokens[1] == 'M')) )//AS:i code for deSAMBA; NM:i code for minimap2
		{
			//get i
			tokens = strtok(NULL,":");
			//get score
			tokens = strtok(NULL,"\t");
			rst->score = strtoul(tokens,NULL,10);
			tokens = strtok(NULL,":");
			//minimap2: ms
			if(tokens != NULL && ((tokens[0] == 'm' && tokens[1] == 's') ))//ms:i code for deSAMBA; NM:i code for minimap2
			{
				tokens = strtok(NULL,":");//ignore 'i'
				//get score
				tokens = strtok(NULL,"\t");
				rst->score = strtoul(tokens,NULL,10);
			}
			tokens = strtok(NULL,":");
			//get score
			tokens = strtok(NULL,"\t");

		}
		//for the ref name part
		{
			//ignore 'tid|'
			ref_tokens = strtok(ref_tokens,"|");
			//get tid
			ref_tokens = strtok(NULL,"|");
			rst->tid = strtoul(ref_tokens,NULL,10);
			//get read length
		}
		//for the read length part
		{
			int read_len = 0;
			int type_len = 0;
			while(1)
			{
				char c_char = *CIGAR++;
				if(c_char == 0)
					break;
				if(c_char <= '9' && c_char >= '0')
					type_len = (type_len * 10) + (c_char - '0');
				else
				{
					if(c_char == 'M' || c_char == 'I' || c_char == 'S' || c_char == 'X')
						read_len += type_len;
					type_len = 0;
				}
			}
			rst->read_length = read_len;
		}
	}
	return 0;
}

//compare two string, but when one char are number[0~9], it bigger than any other char
//used to sort the string of read name
//static int strcmp_name(const void *a_, const void *b_)
//{
//	char * s1 = (char*)a_; char * s2 = (char*)b_;
//	for(;*s1 == *s2 && *s1 != 0;s1++,s2++);
//	int s1_is_num = 0, s2_is_num = 0;
//	if(*s1 <= '9' && *s1 >= '0')
//		s1_is_num = 1;
//	if(*s2 <= '9' && *s2 >= '0')
//		s2_is_num = 1;
//	if(s1_is_num == 1 && s2_is_num == 0)//a is number and b not
//		return 1;
//	if(s1_is_num == 0 && s2_is_num == 1)//a is number and b not
//		return - 1;
//	if(s1_is_num == 0 && s2_is_num == 0)//both not
//		return *s1 - *s2;
//	//both is number
//	//copy to new memory
//	char s1_[1024];
//	char s2_[1024];
//	strcpy(s1_ + 1,s1);
//	strcpy(s2_ + 1,s2);
//	s1_[0] = '1';
//	s2_[0] = '1';
//	int int_a = atoi(s1_);
//	int int_b = atoi(s2_);
//	return int_a - int_b;
//}

//all char that are more than '9' are use less than '9'
//static int cmp_rst(const void *a_, const void *b_){
//	RST *a = (RST *) a_, *b= (RST *) b_;
//	return strcmp_name(a->read_name, b->read_name);
//}

//static int cmp_maf(const void *a_, const void *b_){
//	MAF *a = (MAF *) a_, *b= (MAF *) b_;
//	return strcmp_name(a->read_name, b->read_name);
//}

// ignore the @SQ and @PG line
static void skip_sam_head(FILE * SAM_file, char *buff)
{
	size_t max_l = MAX_BUFF_LEN;
	//step2: ignore the @SQ and @PG line
	while(1){
		int read_L = getline(&buff,&max_l,SAM_file);
		xassert(read_L >= 0,"Read SAM file FAILED\n");
		if(buff[0] != '@'){
			//reset this line
			err_fseek(SAM_file, - read_L, SEEK_CUR);
			break;
		}
	}
}
//
////store minimap SAM in RST format, only the result with the best MAPQ are stored
////USAGE: deSPI3rd cmp dump_minimap [SAM file name] [RST file name]
//static void dump_sam(char * SAM_file_name)
//{
//	fprintf(stderr,"Dump_sam begin\n");
//	//SAM to RST
//	//step1: OPEN files
//	char   RST_file_name[1024];
//	strcpy(RST_file_name, SAM_file_name); strcat(RST_file_name, ".dump.tmp");
//	FILE * SAM_file = xopen(SAM_file_name,"r");
//	FILE * RST_file = xopen(RST_file_name,"wb");
//	char *buff = (char*)malloc(MAX_BUFF_LEN);
//	//step2: ignore the @SQ and @PG line
//	skip_sam_head(SAM_file, buff);
//
//	//step3: store SAM into RST
//	uint32_t record_num = 0;
//	RST temp_rst;
//	while(1){
//		//get SAM
//		if(getOneSAM(SAM_file, buff, &temp_rst) < 0)
//			break;
//		//write MAF
//		fwrite(&temp_rst,sizeof(RST),1,RST_file);
//		record_num++;
//	}
//	//step4: close files
//	fclose(SAM_file);
//	fclose(RST_file);
//	//step5: reload file and sort
//	FILE * rst_file = xopen(RST_file_name, "rb");
//	RST * rst = (RST *)xmalloc(sizeof(RST)*record_num);
//	xread(rst,sizeof(RST),record_num,rst_file);
//	qsort(rst,record_num,sizeof(RST),cmp_rst);
//	fclose(rst_file);
//	//step6: delete results that is not the best MAPQ
//	RST * uni_rst = (RST *)xmalloc(sizeof(RST)*record_num);
//	int uni_record_num = 0;
//	uint32_t rst_index = 0;
//	while(rst_index < record_num){
//		uint8_t best_MAPQ = rst[rst_index].MAPQ;
//		char * read_name = rst[rst_index].read_name;
//		int begin_index = rst_index;
//		//get best MAPQ
//		while(1){
//			rst_index ++;
//			if(strcmp(read_name,rst[rst_index].read_name) == 0)
//				best_MAPQ = MAX(best_MAPQ,rst[rst_index].MAPQ);
//			else
//				break;
//		}
//		//store the rst that has best MAPQ
//		for(uint32_t i = begin_index; i < rst_index; i++){
//			//if(rst[i].MAPQ == best_MAPQ){
//				memcpy(uni_rst + uni_record_num, rst + i, sizeof(RST));
//				uni_record_num++;
//			//}
//		}
//	}
//	//dump
//	xrm(RST_file_name);
//	RST * c_uni_rst = uni_rst;
//	for(int i = 0; i < uni_record_num; i++, c_uni_rst++)
//		printf("%s\t%c\t%d\t%d\t%d\n",
//				c_uni_rst->read_name,
//				c_uni_rst->isClassify,
//				c_uni_rst->tid,
//				c_uni_rst->read_length,
//				c_uni_rst->MAPQ);
//	free(uni_rst);
//	free(buff);
//	free(rst);
//}
#define MIN_read_N 0
#define MAX_read_N MAX_uint32_t
//store minimap SAM in RST format, only the result with the best MAPQ are stored
//USAGE: deSPI3rd cmp dump_minimap [SAM file name] [RST file name]
static void dump_des_sam_file(char * SAM_file_name, char*dump_file_name)
{
	//fprintf(stderr,"Dump_sam begin\n");
	//SAM to RST
	//step1: OPEN files
	FILE * SAM_file = xopen(SAM_file_name,"r");
	FILE * dump_file = xopen(dump_file_name,"w");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step2: ignore the @SQ and @PG line
	skip_sam_head(SAM_file, buff);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	RST temp_rst;
	while(1){
		//get SAM
		if(getOneSAM(SAM_file, buff, &temp_rst) < 0)
			break;
		//write MAF
		record_num++;
		if(record_num > MAX_read_N)
			break;
		if(record_num < MIN_read_N)
			continue;
		fprintf(dump_file, "%s\t%c\t%d\t%d\t%d\t%d\n",
				temp_rst.read_name,
				temp_rst.isClassify,
				temp_rst.tid,
				temp_rst.read_length,
				temp_rst.MAPQ,
				temp_rst.score);
	}
	free(buff);
	fclose(SAM_file);
	fclose(dump_file);
}


//------------------------------------------SAM-----------------------------------------//
int getOneMATEMAP(FILE * SAM_file, char *buff, RST * rst, int *exchange_list)
{
	size_t max_l = MAX_BUFF_LEN;
	char *tokens;
	int read_L = 0;
	read_L = getline(&buff,&max_l,SAM_file);
	if(read_L <= 0)
		return -1;
	//change \0 to blank
	for(int i = 0; i < read_L; i++)
	{
		if(buff[i] < 10)
		{
			buff[i] = ' ';
		}
	}

	//if(buff[0] != 'S' && buff[0] != 'D' && buff[0] != 'E')
	//{
	//	fprintf(stderr, "Not Normal SAM.\n");
	//}
	//get read name
	tokens = strtok(buff," ");
	strcpy(rst->read_name,tokens);
	//get read length
	tokens = strtok(NULL," ");
	rst->read_length = strtoul(tokens,NULL,10);
	//ignore 0
	tokens = strtok(NULL," ");
	//ignore read len
	tokens = strtok(NULL," ");
	//ignore +/-
	tokens = strtok(NULL," ");
	//ignore cXXXXX
	tokens = strtok(NULL,"d");
	//ignore cXXXXX
	tokens = strtok(NULL,"|");
	//get tid
	if(tokens[0] == 'x')
		rst->tid = exchange_list[strtoul(tokens + 1,NULL,10)];
	else
		rst->tid = strtoul(tokens + 0,NULL,10);
	//ignore other |NZ_XXXX
	tokens = strtok(NULL," ");
	//ignore ref part
	tokens = strtok(NULL," ");
	tokens = strtok(NULL," ");
	tokens = strtok(NULL," ");
	tokens = strtok(NULL," ");
	tokens = strtok(NULL," ");
	tokens = strtok(NULL," ");
	tokens = strtok(NULL," ");
	//get score
	float score = strtof(tokens,NULL);
	rst->score = score*10000;
	rst->isClassify = 'C';
	return 0;
}

//store minimap SAM in RST format, only the result with the best MAPQ are stored
//USAGE: deSPI3rd cmp dump_minimap [SAM file name] [RST file name]
void dump_matemaps_file(char * SAM_file_name, char*dump_file_name, char * exchange_file)
{
	//load exchange file
	int *exchange_list = xcalloc(10000, sizeof(int));
	FILE * EXC_file = xopen(exchange_file,"r");
	int subspecies_ID = 0;
	int species_ID = 0;
	while(fscanf(EXC_file, "x%d\t%d\n", &subspecies_ID, &species_ID) > 0)
		exchange_list[subspecies_ID] = species_ID;
	fclose(EXC_file);
	//SAM to RST
	//step1: OPEN files
	FILE * SAM_file = xopen(SAM_file_name,"r");
	FILE * dump_file = xopen(dump_file_name,"w");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step2: ignore the @SQ and @PG line
	skip_sam_head(SAM_file, buff);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	RST temp_rst;
	RST *max_record_L = xmalloc(sizeof(RST)*30000);
	int number_max_record = 0;
	char old_name[128] = "";
	int max_score = 0;
	while(1){
		//get SAM
		if(getOneMATEMAP(SAM_file, buff, &temp_rst, exchange_list) < 0)
			break;
		if(strcmp(old_name, temp_rst.read_name) == 0)
		{
			if(max_score < temp_rst.score)
			{
				memcpy(max_record_L, &temp_rst,sizeof(RST));
				max_score = temp_rst.score;
				number_max_record = 1;
			}
			else if(max_score == temp_rst.score)
			{
				memcpy(max_record_L + number_max_record, &temp_rst,sizeof(RST));
				number_max_record++;
			}
		}
		else if(record_num != 0)
		{
			if(number_max_record > 10000)
				fprintf(stderr, "read: %s, number %d\n", max_record_L[0].read_name, number_max_record);
			for(int i = 0; i < number_max_record; i++)
			{
				fprintf(dump_file, "%s\t%c\t%d\t%d\t%d\t%d\n",
						max_record_L[i].read_name,
						max_record_L[i].isClassify,
						max_record_L[i].tid,
						max_record_L[i].read_length,
						max_record_L[i].MAPQ,
						max_record_L[i].score);
			}
			max_score = 0;
			strcpy(old_name, temp_rst.read_name);
			memcpy(max_record_L, &temp_rst,sizeof(RST));
			number_max_record = 1;
		}
		else
			strcpy(old_name, temp_rst.read_name);
		//write MAF
		record_num++;
	}
	free(buff);
	fclose(SAM_file);
	fclose(dump_file);
}

//------------------------------------------SAM-----------------------------------------//
static int getOnePAF(FILE * PAF_file, char *buff, RST * rst)
{
	size_t max_l = MAX_BUFF_LEN;
	char *tokens;
	int read_L = 0;
	read_L = getline(&buff,&max_l,PAF_file);
	if(read_L <= 0)
		return -1;
	//if(buff[0] != 'S' && buff[0] != 'D' && buff[0] != 'E')
	//{
	//	fprintf(stderr, "Not Normal SAM.\n");
	//}
	//get read name
	tokens = strtok(buff,"\t");
	strcpy(rst->read_name,tokens);
	//ignore four part
	tokens = strtok(NULL,"\t");
	tokens = strtok(NULL,"\t");
	tokens = strtok(NULL,"\t");
	tokens = strtok(NULL,"\t");
	tokens = strtok(NULL,"\t");
	//get refNAME
	char * ref_tokens = tokens;
	ref_tokens = strtok(ref_tokens,"|");
	//get tid
	ref_tokens = strtok(NULL,"|");
	rst->tid = strtoul(ref_tokens,NULL,10);
	rst->MAPQ = 0;
	rst->read_length = 0;
	return 0;
}

//store minimap SAM in RST format, only the result with the best MAPQ are stored
//USAGE: deSPI3rd cmp dump_minimap [SAM file name] [RST file name]
static void dump_des_PAF_file(char * PAF_file_name, char*dump_file_name)
{
	//fprintf(stderr,"Dump_sam begin\n");
	//SAM to RST
	//step1: OPEN files
	FILE * PAF_file = xopen(PAF_file_name,"r");
	FILE * dump_file = xopen(dump_file_name,"w");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	RST temp_rst;
	while(1){
		//get SAM
		if(getOnePAF(PAF_file, buff, &temp_rst) < 0)
			break;
		//write MAF
		record_num++;
		temp_rst.isClassify = 'C';
		fprintf(dump_file, "%s\t%c\t%d\t%d\t%d\n",
				temp_rst.read_name,
				temp_rst.isClassify,
				temp_rst.tid,
				temp_rst.read_length,
				temp_rst.MAPQ);
	}
	free(buff);
	fclose(PAF_file);
	fclose(dump_file);
}


//------------------------------------------MAF-----------------------------------------//
//get one line of MAF file, store result in maf
//return -1 when reach the end of the file
//buff pool are needed, "char *buff = (char*)malloc(50000);"
//static int getOneMAF(FILE * MAF_file, char *buff, MAF * maf)
//{
//	size_t max_l = 50000;
//	char *tokens;
//	for(int flag = 0;flag < 4;flag++){
//		if(getline(&buff,&max_l,MAF_file) <= 0)
//			return -1;
//		switch(flag) {
//		case 0: break;//ignore that line
//		case 1:	tokens = strtok(buff + 2,"|"); tokens = strtok(NULL,"|");
//			maf->tid = strtoul(tokens,NULL, 10); break;//get tid
//		case 2:
//			tokens = strtok(buff + 2," ");
//			strcpy(maf->read_name,tokens); break;//get name
//		case 3: break;//ignore that line
//		}
//	}
//	return 0;
//}

//store pbsim MAF file in maf format
//USAGE: deSPI3rd cmp dump_pb_maf [PB MAF file name] [MAF file name]
//static void dump_maf(char * pacbio_maf_file_name)
//{
//	//MAF to maf.srt
//	//open maf file
//	char   maf_file_name[1024];
//	strcpy(maf_file_name, pacbio_maf_file_name); strcat(maf_file_name, ".dump");
//	FILE * pacbio_maf_file = xopen(pacbio_maf_file_name, "r");
//	FILE * maf_file_srt = xopen(maf_file_name, "wb");
//	uint32_t maf_number = 0;
//	char *buff = (char*)xmalloc(50000);
//	MAF maf;
//	while(1){
//		if(getOneMAF(pacbio_maf_file, buff, &maf) < 0) break;//get MAF
//		fwrite(&maf,sizeof(MAF),1,maf_file_srt);			//write MAF
//		maf_number++;
//	}
//	//close file
//	fclose(pacbio_maf_file);
//	fclose(maf_file_srt);
//	//reload and sort
//	maf_file_srt = xopen(maf_file_name, "rb");
//	MAF *maf_v = (MAF*)malloc(sizeof(MAF)*maf_number);
//	xread(maf_v,sizeof(MAF),maf_number,maf_file_srt);
//	qsort(maf_v,maf_number,sizeof(MAF),cmp_maf);
//	fclose(maf_file_srt);
//	//dump
//	xrm(maf_file_name);
//	MAF *c_maf = maf_v;
//	for(uint32_t i = 0; i < maf_number; i++, c_maf++)
//		printf("%s\t%d\n", c_maf->read_name, c_maf->tid);
//	free(maf_v);
//}
//
//static void inline cmp_maf_rst(MAF * maf, RST * rst, MAF_CMP * cmp, TAXONOMY_rank * taxonomy)
//{
//	//step1: assert the same read name and store read name
//	xassert(strcmp_name(maf->read_name,rst->read_name) == 0,": Wrong read name\n");
//	//store read name
//	strcpy(cmp->read_name,maf->read_name);
//	//step2: store true tid and rank
//	cmp->true_tid = maf->tid;
//	strcpy(cmp->true_rank,taxonomy[maf->tid].rank);
//	//step3: store cly tid and rank
//	cmp->cly_tid = rst->tid;
//	strcpy(cmp->cly_rank,taxonomy[rst->tid].rank);
//	//step4: get tid distance and LCA
//	uint32_t LCA;
//	cmp->distance = TID_distance(cmp->true_tid,cmp->cly_tid, taxonomy, &LCA);
//	cmp->LCA_TID = LCA;
//	strcpy(cmp->LCA_rank,taxonomy[LCA].rank);
//	//step5: print result
//	printf("%10s\t%d\t%12s\t%d\t%12s\t%d\t%u\t%12s\n",
//		cmp->read_name,
//		cmp->true_tid,
//		cmp->true_rank,
//		cmp->cly_tid,
//		cmp->cly_rank,
//		cmp->distance,
//		cmp->LCA_TID,
//		cmp->LCA_rank);
//}

//static int getOneMAF_DUMP(FILE * MAF_file,MAF * maf)
//{
//	size_t max_l = 1024;
//	char static_BUFF[1024];
//	char *buff = static_BUFF;
//	char *tokens;
//	if(getline(&buff,&max_l,MAF_file) <= 0)
//		return -1;
//	tokens = strtok(buff,"\t");
//	strcpy(maf->read_name,tokens);
//	tokens = strtok(NULL,"\t");
//	maf->tid = strtoul(tokens,NULL, 10);
//	return 0;
//}

//compare minimap rst file with maf file
//USAGE: deSPI3rd cmp cmp_minimap [rst file name] [maf file name] [tax_file_name]
//void rst_ana(char * rst_file_name, char * maf_file_name, char * tax_file_name)
//{
//	//step1: load taxonomyTree
//	TAXONOMY_rank * taxonomyTree = NULL;
//	taxonTree_rank(tax_file_name, &taxonomyTree);
//
//	//step2: compare
//	FILE * rst_file_srt = xopen(rst_file_name, "rb");
//	FILE * maf_file_srt = xopen(maf_file_name, "rb");
//	MAF_CMP cmp;
//	int total_distance = 0;
//	int total_unclassified_reads = 0;
//	uint32_t maf_number = 0, rst_number = 0;
//	RST rst; int eof_rst = getOneRST(rst_file_srt, &rst);
//	MAF maf; int eof_maf = getOneMAF_DUMP(maf_file_srt, &maf);
//	for(; eof_rst >= 0 && eof_maf >= 0;){
//		int cmp_rst = strcmp_name(maf.read_name, rst.read_name);
//		if(cmp_rst == 0){
//			cmp_maf_rst(&maf, &rst, &cmp, taxonomyTree);
//			if(rst.isClassify == 'U') 	total_unclassified_reads++;
//			else						total_distance += cmp.distance;
//			eof_rst = getOneRST(rst_file_srt, &rst); rst_number++;
//		}
//		else if(cmp_rst < 0){
//			eof_maf = getOneMAF_DUMP(maf_file_srt, &maf); maf_number++;
//		}
//		else {
//			fprintf(stderr,"UNKNOW read\n");
//			eof_rst = getOneRST(rst_file_srt, &rst); rst_number++;
//		}
//	}
//	//[1] total simulate reads number
//	//[2] total results number
//	//[3] total distance
//	//[4] total unclassified reads
//	char print_string[10000];
//	sprintf(print_string,
//			"[1] total simulate reads number: [%d]\n"
//			"[2] total results number: [%d]\n"
//			"[3] total distance: [%d]\n"
//			"[4] total unclassified reads: [%d]\n"
//			,maf_number,rst_number,total_distance,total_unclassified_reads);
//	printf("%s",print_string);
//	fprintf(stderr,"%s",print_string);
//	//step3: end PRO
//	free(taxonomyTree);
//	fclose(rst_file_srt);
//	fclose(maf_file_srt);
//}

//------------------------------------------CEN-----------------------------------------//

static int getOnecenSAM(FILE * SAM_file, char *buff, RST * rst)
{
	size_t max_l = MAX_BUFF_LEN;
	char *tokens;
	int read_L = 0;
	read_L = getline(&buff,&max_l,SAM_file);
	if(read_L <= 0)
		return -1;
	//get read name
	tokens = strtok(buff,"\t");
	strcpy(rst->read_name,tokens);
	//ignore flag
	tokens = strtok(NULL,"\t");
	//get tax_ID
	tokens = strtok(NULL,"\t");
	rst->tid = strtoul(tokens,NULL,10);
	rst->MAPQ = 0;
	rst->read_length = 0;
	rst->score = 1;
	if(rst->tid == 0) {
		rst->isClassify = 'U';
	}
	else{
		rst->isClassify = 'C';
		//ignore 0
		tokens = strtok(NULL,"\t");
		//ignore 0
		tokens = strtok(NULL,"\t");
		//ignore *0
		tokens = strtok(NULL,"\t");
		//ignore cid|***
		tokens = strtok(NULL,"\t");
		//ignore 0
		tokens = strtok(NULL,"\t");
		//get read length
		tokens = strtok(NULL,"\t");
		rst->read_length = strtoul(tokens,NULL,10);
	}
	return 0;
}

//store minimap SAM in RST format, only the result with the best MAPQ are stored
//USAGE: deSPI3rd cmp dump_minimap [SAM file name] [RST file name]
static void dump_CEN_file(char * SAM_file_name, char*dump_file_name)
{
	//SAM to RST
	//step1: OPEN files
	FILE * SAM_file = xopen(SAM_file_name,"r");
	FILE * dump_file = xopen(dump_file_name,"w");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	RST temp_rst;
	while(1){
		//get SAM
		if(getOnecenSAM(SAM_file, buff, &temp_rst) < 0)
			break;
		//write MAF
		record_num++;
		fprintf(dump_file, "%s\t%c\t%d\t%d\t%d\t%d\n",
				temp_rst.read_name,
				temp_rst.isClassify,
				temp_rst.tid,
				temp_rst.read_length,
				temp_rst.MAPQ,
				temp_rst.score);
	}
	free(buff);
	fclose(SAM_file);
	fclose(dump_file);
}
//------------------------------------------KAI-----------------------------------------//
static int getOnekaiSAM(FILE * SAM_file, char *buff, RST * rst)
{
	size_t max_l = MAX_BUFF_LEN;
	char *tokens;
	int read_L = 0;
	read_L = getline(&buff,&max_l,SAM_file);
	if(read_L <= 0)
		return -1;
	//get U/C
	rst->isClassify = buff[0];
	//get read name
	tokens = strtok(buff + 2,"\t");
	strcpy(rst->read_name,tokens);
	if(rst->isClassify == 'C')
	{
		tokens = strtok(NULL,"\t");
		tokens = strtok(NULL,"\t");
		tokens = strtok(NULL,",");
		rst->tid = strtoul(tokens,NULL,10);
	}
	else
	{
		rst->tid = 0;
	}
	//get tax_ID

	rst->MAPQ = 0;
	rst->read_length = 0;
	return 0;
}

static void dump_KAI_file(char * SAM_file_name, char*dump_file_name)
{
	//step1: OPEN files
	FILE * SAM_file = xopen(SAM_file_name,"r");
	FILE * dump_file = xopen(dump_file_name,"w");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	RST temp_rst;
	while(1){
		//get SAM
		if(getOnekaiSAM(SAM_file, buff, &temp_rst) < 0)
			break;
		//write MAF
		record_num++;
		fprintf(dump_file, "%s\t%c\t%d\t%d\t%d\n",
				temp_rst.read_name,
				temp_rst.isClassify,
				temp_rst.tid,
				temp_rst.read_length,
				temp_rst.MAPQ);
	}
	free(buff);
	fclose(SAM_file);
	fclose(dump_file);
}

//------------------------------------------DES-----------------------------------------//
//static void print_rst(RST * rst){printf("%s\t%c\t%d\t%d\t%d\n", rst->read_name, rst->isClassify, rst->tid, rst->read_length, rst->MAPQ);}
////
//static void dump_des(char * deSPI_rst_file_name)
//{
//	//DES to RST
//	//step1: OPEN files
//	//TO
//	FILE * des_file = xopen(deSPI_rst_file_name,"r");
//	char *buff = (char*)malloc(MAX_BUFF_LEN);
//	size_t max_l = MAX_BUFF_LEN; char *tokens;  RST rst;
//
//	while(1){
//		int read_L = getline(&buff,&max_l,des_file);
//		if(read_L <= 0) break;
//		///head line:
//		//get read name
//		tokens = strtok(buff,"\t");
//		strcpy(rst.read_name,tokens);
//		//is classify ?
//		tokens = strtok(NULL,"\t");
//		rst.isClassify = tokens[0];
//		//ignore speed
//		tokens = strtok(NULL,"\t");
//		//get read length
//		tokens = strtok(NULL,"\t");
//		rst.read_length = strtol(tokens,NULL,10);
//		rst.MAPQ = 0;
//		tokens = strtok(NULL,"[");
//		tokens = strtok(NULL,"]");
//		int rst_number = strtol(tokens,NULL,10);
//		tokens = strtok(NULL,"[");
//		tokens = strtok(NULL,"]");
//		int anc_number = strtol(tokens,NULL,10);
//
//		if(rst.isClassify == 'U'){
//			rst.tid = 0;
//			print_rst(&rst);
//		}
//
//		for(int i = 0; i < rst_number; i++)	{
//			read_L = getline(&buff,&max_l,des_file);
//			xassert(read_L > 0, "");
//			if(buff[0] == '\n')
//				break;
//			//cnt and PRI
//			tokens = strtok(buff + 4," ");
//			//F/R
//			tokens = strtok(NULL," ");
//			//tid|
//			tokens = strtok(NULL,"|");
//			tokens = strtok(NULL,"|");
//			rst.tid = strtoul(tokens,NULL,10);
//			print_rst(&rst);
//		}
//		if(rst.isClassify == 'U')
//			for(int i = 0; i < anc_number; i++)
//			{
//				read_L = getline(&buff,&max_l,des_file);
//				xassert(read_L > 0, "");
//				if(buff[0] == '\n')
//					break;
//			}
//		if(buff[0] == '\n')
//			continue;
//		//read "\n" line
//		read_L = getline(&buff,&max_l,des_file);
//		if(read_L > 0)
//			xassert(buff[0] == '\n', "");
//	}
//	free(buff);
//	fclose(des_file);
//}

//-----------------------------------------ANA TAX-----------------------------------------//
//rank equal one of "species"/"genus"/"superkingdom"/"family"/"order"/"class"/"phylum"
static uint32_t get_tax_by_rank(TAXONOMY_rank * taxonomyTree, uint32_t tax, char *rank)
{
	uint32_t c_tax = tax;
	//get species name
	uint32_t rst_tax_id = 0;
	while(1){
		if(strcmp(taxonomyTree[c_tax].rank, rank) == 0)
		{
			rst_tax_id = c_tax;
			break;
		}
		c_tax = taxonomyTree[c_tax].p_tid;
		if(c_tax <= 1 || c_tax == 0xffffffff)
			break;
	}
	//if(c_tax == 0xffffffff)
	//	fprintf(stderr, "");
	return rst_tax_id;
}

//-----------------------------------------ANA TAX-----------------------------------------//
//compare whether tax A is the ancestors node of tax B
static bool compare_tax(TAXONOMY_rank * taxonomyTree, uint32_t tax_A, uint32_t tax_B)
{
	uint32_t c_tax = tax_B;
	//get species name
	while(1){
		if(c_tax == tax_A)
			return true;
		c_tax = taxonomyTree[c_tax].p_tid;
		if(c_tax <= 1 || c_tax == 0xffffffff)
			break;
	}
	//if(c_tax == 0xffffffff)
	//	fprintf(stderr, "");
	return false;
}

//compare
//USAGE: deSPI3rd cmp cmp_minimap [rst file name] [maf file name] [tax_file_name]
//rank equal one of "species"/"genus"/"superkingdom"/"family"/"order"/"class"/"phylum"
//when rank is null, the program will determine the rank automatically
#define ANA_SHOW_TITLE false
#define SHOW_DETAIL 1
static void ana_tax(char * rst_file_name, uint32_t right_tax, char * tax_file_name, char *rank)
{
	//PART 1: load result file
	//char * rst_file_name = argv[1];
	fprintf(stderr, "%s\t", rst_file_name);
	FILE * rst_file_srt = xopen(rst_file_name, "rb");

	//PART 2:get right tax
	//uint32_t right_tax = strtoul(argv[2],0,10);
	bool no_rank = false;
	if(strcmp(rank, "null") == 0)
		no_rank = true;

	//PART 3:load taxonomyTree
	//char * tax_file_name = argv[3];
	TAXONOMY_rank * taxonomyTree = NULL;
	taxonTree_rank(tax_file_name, &taxonomyTree);

	//step4: analysis
	int wrong_alignment = 0;
	int total_read_number = 0;
	int umapped_read_number = 0;
	char old_read_name[READ_NAME_LEN] = {0};
	bool right_alignment = 0;
	int right_aligment_first = 0;
	int right_aligment_second =0;
	RST rst; //int eof_rst = getOneRST(rst_file_srt, &rst);
	uint32_t c_tax;


	if(getOneRST(rst_file_srt, &rst) < 0)
		return;
	while(1)
	{
		//get the first result
		total_read_number++;
		if(SHOW_DETAIL)
			printf("\n%s ",rst.read_name);
		if(rst.isClassify == 'U')
		{
			umapped_read_number++;
			if(SHOW_DETAIL)
				printf("UM");
			int eof_ = getOneRST(rst_file_srt, &rst);
			if(eof_ < 0)
				break;
			else
				continue;
		}

		bool right_classify = false;
		if(no_rank)
			right_classify = compare_tax(taxonomyTree, right_tax, rst.tid);
		else
		{
			c_tax = get_tax_by_rank(taxonomyTree, rst.tid, rank);
			if(right_tax == c_tax)
				right_classify = true;
		}
		if(right_classify == true){
			right_alignment = true;
			right_aligment_first++;
			if(SHOW_DETAIL)
				printf("PRI");
		}
		else
			right_alignment = false;
		strcpy(old_read_name, rst.read_name);
		int eof_ = 0;
		//get other results
		while(1)
		{
			eof_ = getOneRST(rst_file_srt, &rst);
			if(eof_ < 0)
				break;
			if(strcmp(old_read_name, rst.read_name) == 0)
			{
				if(right_alignment == true)
					continue;

				bool right_classify = false;
				if(no_rank)
					right_classify = compare_tax(taxonomyTree, right_tax, rst.tid);
				else
				{
					c_tax = get_tax_by_rank(taxonomyTree, rst.tid, rank);
					if(right_tax == c_tax)
						right_classify = true;
				}

				if(right_classify)
				{
					right_alignment = true;
					right_aligment_second++;
					if(SHOW_DETAIL)
						printf("SEC");
				}
			}
			else
				break;
		}
		if(eof_ < 0)
			break;

		if(right_alignment == false)
			wrong_alignment++;
	}
	//title:
	if(ANA_SHOW_TITLE)
	{
		fprintf(stderr,"[1] total reads number\n");
		fprintf(stderr,"[2] total unmapped read number\n");
		fprintf(stderr,"[3] total right PRIMARY\n");
		fprintf(stderr,"[4] total right ALL[PRIMARY + SEN + SUP]\n");

		fprintf(stderr,"[5] unmapped rate\n");
		fprintf(stderr,"[6] SEN PRI(RIGHT/ALL)[PRI]\n");
		fprintf(stderr,"[7] ACC PRI(RIGHT/(RIGHT+WRONG))[PRI]\n");

		fprintf(stderr,"[8] SEN ALL(RIGHT/ALL)[PRI + SEN + SUP]\n");
		fprintf(stderr,"[9] ACC ALL(RIGHT/(RIGHT+WRONG))[PRI + SEN + SUP]\n");
	}

	//number
	fprintf(stderr,"%d\t",total_read_number);
	fprintf(stderr,"%d\t",umapped_read_number);
	fprintf(stderr,"%d\t",right_aligment_first);
	fprintf(stderr,"%d\t",right_aligment_second + right_aligment_first);

	fprintf(stderr,"%f%%\t", (float)umapped_read_number/total_read_number *100);
	fprintf(stderr,"%f%%\t", (float)right_aligment_first/total_read_number *100);
	fprintf(stderr,"%f%%\t", (float)(right_aligment_first)/(total_read_number - umapped_read_number)*100);

	fprintf(stderr,"%f%%\t",(float)(right_aligment_second + right_aligment_first)/total_read_number * 100);
	fprintf(stderr,"%f%%\n",(float)(right_aligment_second + right_aligment_first)/(total_read_number - umapped_read_number)*100);

	//step3: end PRO
	free(taxonomyTree);
	fclose(rst_file_srt);
	//fclose(maf_file_srt);
}

//-----------------------------------------ANA META-----------------------------------------//
//compare
//USAGE: deSPI3rd cmp cmp_minimap [rst file name] [maf file name] [tax_file_name]
//rank equal one of "species"/"genus"/"superkingdom"/"family"/"order"/"class"/"phylum"
typedef struct{
	uint32_t	tid;
	uint32_t	next;
}CN_CHILD;

typedef struct{
	char 		tax_name[32];
	uint32_t	tid;
	char 		rank[20];
	uint64_t 		weight;
	uint32_t	child_list_begin;
	uint64_t 		total_mapQ;
}CLY_NODE;


// 2024.3.29需求：按照网信格式打印desamba的结果
int ANA_PRINT_USE_LIST = 0;
static void ana_meta_loop_fprint(FILE * file, TAXONOMY_rank * taxonomyTree, CLY_NODE *list, uint32_t node_ID, CN_CHILD *child_list, int level, uint64_t total_weight, bool is_base)
{
	CLY_NODE *node = list + node_ID;
	if (node->weight == 0) return; // 剪枝
	float rate =  (float)node->weight/total_weight;
	// if (DEBUG) {
	// 	for (int i = 0; i < level; i ++) fprintf(stderr, "  ");
	// }
	// if (DEBUG) fprintf(stderr, "hello ana_meta_loop_fprint: %d\t%s\t%ld/%ld\n", node_ID, taxonomyTree[node_ID].name, node->weight, total_weight);
	if(node->child_list_begin != 0)
	{
		uint32_t child = node->child_list_begin;
		while(1)
		{
			ana_meta_loop_fprint(file, taxonomyTree, list, child_list[child].tid, child_list, level + 1, total_weight, is_base);
			if(child_list[child].next == 0)
				break;
			else
				child = child_list[child].next;
		}
	}
	else // 叶节点
	{
		char species_type[20] = {0};
		uint32_t leaf_ID = node_ID;
		if (leaf_ID == 0 || leaf_ID == 1) { // root or CLY_FAIL
			strcpy(species_type, "no_match");
		}
		else {
			while (node_ID != MAX_uint32_t) {
				if (node_ID == 9606) {
					strcpy(species_type, "human");
					break;
				}
				else if (node_ID == 33208 || node_ID == 33090) {
					strcpy(species_type, "animal_and_plant");
					break;
				}
				else {
					node_ID = taxonomyTree[node_ID].p_tid;
				}
			}
			if (strlen(species_type) == 0)
			{
				strcpy(species_type, "microbe");
			}
		}
		fprintf(file, "%s\t%d|%s\tnull\t%f\n", species_type, leaf_ID, taxonomyTree[leaf_ID].rank, rate);
		if (DEBUG) {
			for (int i = 0; i < level; i ++) fprintf(stderr, "  ");
			fprintf(stderr, "DEBUG: %s\t%d|%s\tnull\t%f\n", species_type, leaf_ID, taxonomyTree[leaf_ID].rank, rate);
		} 
	}
}

static void ana_meta_loop_print(TAXONOMY_rank * taxonomyTree, CLY_NODE *list, uint32_t node_ID, CN_CHILD *child_list, int level, uint64_t total_read_number, bool is_base)
{
	CLY_NODE *node = list + node_ID;
	float rate =  (float)node->weight/total_read_number*100;
	float map_q = (float)node->total_mapQ/node->weight * rate;
	if(rate < 0.01)
		return;
	for(int i = 0; i < level; i++)
		printf("|");
	if(is_base)
		printf("%s TID:%d %s %f%%, mapQ:%f\n",taxonomyTree[node_ID].rank, node_ID, node->tax_name, rate, map_q);
	else
		printf("%s TID:%d %s %f%%\n",taxonomyTree[node_ID].rank, node_ID, node->tax_name, rate);
	if(node->child_list_begin != 0)
	{
		uint32_t child = node->child_list_begin;
		while(1)
		{
			ana_meta_loop_print(taxonomyTree, list, child_list[child].tid, child_list, level + 1, total_read_number, is_base);
			if(child_list[child].next == 0)
				break;
			else
				child = child_list[child].next;
		}
	}
}

typedef struct{
	uint32_t tid;
	int count;
}COUNT_SORT;
//all char that are more than '9' are use less than '9'
static int cmp_count_sort(const void *a_, const void *b_){
	COUNT_SORT *a = (COUNT_SORT*) a_, *b= (COUNT_SORT *) b_;
	return a->count < b->count;
}

//analysis
static uint32_t ana_get_tid(RST *rst, int max_tid, FILE * rst_file_srt, int *eof_, TAXONOMY_rank * taxonomyTree, int *read_len, float* coverage)
{
	char old_read_name[READ_NAME_LEN];
	uint32_t tid = 0;
	uint32_t score = 0;
	*eof_ = 0;
	*read_len = rst->read_length;
	//when UNMAPPED
	if(rst->isClassify != 'C')//UNMAPPED
	{
		if(getOneRST(rst_file_srt, rst) < 0)
			*eof_ = -1;
		return 0;
	}
	//get PRIMARY result
	//copy name
	strcpy(old_read_name, rst->read_name);
	//get tid and score
	if(rst->tid <=  max_tid)//store tid
	{
		tid = rst->tid;
		score = rst->score;
		if(rst->read_length > 0)
			*coverage = (float)score/rst->read_length;
		else
			*coverage = 0;
	}
	//get other results
	while(1)//get and ignore other results
	{
		*eof_ = getOneRST(rst_file_srt, rst);
		//return when reach end of file
		if(*eof_ < 0)
			return 0;
		//stop when search to other read
		if(strcmp(old_read_name, rst->read_name) != 0)
			break;
		if(score == 0)
			break;//with out score information
		//search other results:
		if(rst->score != score)
			continue;
		if(rst->tid >  max_tid)//UNKNOW TID
			continue;
		//when score are the same. detect whether rst->tid is child node of tid
		uint32_t p_tid = rst->tid;
		while(1)
		{
			if(p_tid == tid)//same tid or is parent
			{
				tid = rst->tid;
				break;
			}
			if(p_tid < 1 || p_tid == 4294967295)//over root node, without result
				break;
			p_tid = taxonomyTree[p_tid].p_tid;//update p_tid
		}
	}
	return tid;
}

static void ana_meta(char * rst_file_name, char * tax_file_name)
{
	printf("Current read %s\t", rst_file_name);
	//PART 1: load result file
	//char * rst_file_name = argv[1];
	printf("%s\t", rst_file_name);
	FILE * rst_file_srt = xopen(rst_file_name, "rb");

	//PART 2:load taxonomyTree
	//char * tax_file_name = argv[2];
	TAXONOMY_rank * taxonomyTree = NULL;
	int max_tid = taxonTree_rank(tax_file_name, &taxonomyTree);

	//step3: analysis
	uint32_t *node_count = xcalloc(max_tid, sizeof(uint32_t));
	//CLY_NODE_V cly_v = {0,0,0};
	//CLY_NODE root={"root", 1, "root", 0, NULL, NULL};
	// = {"root", 0. "root", 0, NULL, NULL}
	int total_read_number = 0;
	RST rst;
	int eof_ = 0;
	float coverage = 0;

	//count all PRIMARY tax
	if(getOneRST(rst_file_srt, &rst) < 0)
		return;
	while(1)
	{
		total_read_number++;
		int read_len = 0;
		//get all results for one read
		uint32_t final_tid = ana_get_tid(&rst, max_tid, rst_file_srt, &eof_, taxonomyTree, &read_len, &coverage);
		if(final_tid > 0)
		{
			node_count[final_tid]++;
		}
		if(eof_ < 0)//end if file
			break;
//
//
//
//
//		//get the first result
//		total_read_number++;
//		if(rst.isClassify == 'C')//PRIMARY classified
//		{
//			if(rst.tid <=  max_tid)//store tid
//				node_count[rst.tid]++;
//		}
//		else
//		{
//			if(getOneRST(rst_file_srt, &rst) < 0)
//				break;
//			else
//				continue;
//		}
//		strcpy(old_read_name, rst.read_name);
//		int eof_ = 0;
//		while(1)//get and ignore other results
//		{
//			eof_ = getOneRST(rst_file_srt, &rst);
//			if(eof_ < 0 || (strcmp(old_read_name, rst.read_name) != 0))
//				break;
//		}
//		if(eof_ < 0)
//			break;
	}

	CLY_NODE *node_table = xcalloc(max_tid, sizeof(CLY_NODE));
	CN_CHILD *child_list = xcalloc(max_tid*2, sizeof(CN_CHILD));
	uint32_t child_count = 1;

	COUNT_SORT *sort = xmalloc(max_tid* sizeof(COUNT_SORT));
	int rst_num = 0;
	for(uint32_t i = 0; i <= max_tid; i++)
	{
		if(node_count[i] == 0)
			continue;
		sort[rst_num].tid = i;
		sort[rst_num++].count = node_count[i];
	}
	qsort(sort,rst_num,sizeof(COUNT_SORT),cmp_count_sort);

	//calculate middle nodes
	for(int i = 0; i < rst_num; i++)
	{
		uint32_t c_tid = sort[i].tid;
		node_table[sort[i].tid].weight += node_count[sort[i].tid];
		while(1)
		{
			uint32_t p_tid = taxonomyTree[c_tid].p_tid;
			if(p_tid < 1 || p_tid == 4294967295)//over root node
				break;
			node_table[p_tid].weight += node_count[sort[i].tid];//add weight for p_tid
			//store p_tid->c_tid
			if(node_table[p_tid].child_list_begin == 0)
			{
				node_table[p_tid].child_list_begin = child_count++;
				child_list[child_count - 1].tid = c_tid;
			}
			else
			{
				int list_begin = node_table[p_tid].child_list_begin;
				while(child_list[list_begin].tid != c_tid && child_list[list_begin].next != 0)
					list_begin = child_list[list_begin].next;
				if(child_list[list_begin].tid != c_tid && child_list[list_begin].next == 0)//store new node
				{
					child_list[list_begin].next = child_count++;
					child_list[child_count - 1].tid = c_tid;
				}
			}
			//end, reset c_tid
			c_tid = p_tid;
		}
	}
	//out put
	printf("Data:\n");
	if (ANA_PRINT_USE_LIST == 1)
		ana_meta_loop_fprint(stdout, taxonomyTree, node_table, 1, child_list, 0, total_read_number, false);
	else
		ana_meta_loop_print(taxonomyTree, node_table, 1, child_list, 0, total_read_number, false);
	//printf("</Data>\n");
	//number
	printf("total_read_number :%d\t",total_read_number);

	//step3: end PRO
	free(taxonomyTree);
	fclose(rst_file_srt);
	//fclose(maf_file_srt);
}

typedef struct{
	uint32_t tid;
	uint64_t base;
}BASE_SORT;
//all char that are more than '9' are use less than '9'
static int cmp_base_sort(const void *a_, const void *b_){
	BASE_SORT *a = (BASE_SORT*) a_, *b= (BASE_SORT *) b_;
	return a->base < b->base;
}

#define MIN_SCORE 10
static void ana_meta_base(char * rst_file_name, char * tax_file_name)
{
	printf("Current read %s\t", rst_file_name);
	//PART 1: load result file
	//char * rst_file_name = argv[1];
	printf("%s\t", rst_file_name);
	FILE * rst_file_srt = xopen(rst_file_name, "rb");

	//PART 2:load taxonomyTree
	//char * tax_file_name = argv[2];
	TAXONOMY_rank * taxonomyTree = NULL;
	int max_tid = taxonTree_rank(tax_file_name, &taxonomyTree);

	//step3: analysis
	uint64_t *node_base = xcalloc(max_tid, sizeof(uint64_t));
	//CLY_NODE_V cly_v = {0,0,0};
	//CLY_NODE root={"root", 1, "root", 0, NULL, NULL};

	// = {"root", 0. "root", 0, NULL, NULL}
	int total_read_number = 0;
	uint64_t total_base_num = 0;
	uint64_t low_identity_read_num = 0;
	uint64_t low_identity_read_base = 0;
	RST rst;
	int eof_ = 0;//end of file
	float coverage = 0;
	//count all PRIMARY tax
	if(getOneRST(rst_file_srt, &rst) < 0)
		return;
	while(1)
	{
		total_read_number++;
		int read_len = 0;
		//get all results for one read
		uint32_t final_tid = ana_get_tid(&rst, max_tid, rst_file_srt, &eof_, taxonomyTree, &read_len, &coverage);
		if(final_tid > 0)
		{
			if(coverage * read_len > MIN_SCORE)
			{
				total_base_num += read_len;
				node_base[final_tid] += read_len;
				if(coverage < 0.08)
				{
					low_identity_read_base += read_len;
					low_identity_read_num++;
				}
			}
		}
		if(eof_ < 0)//end if file
			break;
//
//		if(rst.isClassify == 'C')//PRIMARY classified
//		{
//			if(rst.tid <=  max_tid)//store tid
//			{
//				total_base_num += rst.read_length;
//				node_base[rst.tid] += rst.read_length;
//			}
//		}
//		else
//		{
//			if(getOneRST(rst_file_srt, &rst) < 0)
//				break;
//			else
//				continue;
//		}
//		strcpy(old_read_name, rst.read_name);
//		int eof_ = 0;
//		while(1)//get and ignore other results
//		{
//			eof_ = getOneRST(rst_file_srt, &rst);
//			if(eof_ < 0 || (strcmp(old_read_name, rst.read_name) != 0))
//				break;
//		}
//		if(eof_ < 0)
//			break;
	}

	CLY_NODE *node_table = xcalloc(max_tid, sizeof(CLY_NODE));
	CN_CHILD *child_list = xcalloc(max_tid*2, sizeof(CN_CHILD));
	uint32_t child_count = 1;

	BASE_SORT *sort = xmalloc(max_tid* sizeof(BASE_SORT));
	int rst_num = 0;
	for(uint32_t i = 0; i <= max_tid; i++)
	{
		if(node_base[i] == 0)
			continue;
		sort[rst_num].tid = i;
		sort[rst_num++].base = node_base[i];
	}
	qsort(sort,rst_num,sizeof(BASE_SORT),cmp_base_sort);

	//calculate middle nodes
	for(int i = 0; i < rst_num; i++)
	{
		uint32_t c_tid = sort[i].tid;
		node_table[sort[i].tid].weight += node_base[sort[i].tid];
		while(1)
		{
			uint32_t p_tid = taxonomyTree[c_tid].p_tid;
			if(p_tid < 1 || p_tid == 4294967295)//over root node
				break;
			node_table[p_tid].weight += node_base[sort[i].tid];//add weight for p_tid
			//store p_tid->c_tid
			if(node_table[p_tid].child_list_begin == 0)
			{
				node_table[p_tid].child_list_begin = child_count++;
				child_list[child_count - 1].tid = c_tid;
			}
			else
			{
				int list_begin = node_table[p_tid].child_list_begin;
				while(child_list[list_begin].tid != c_tid && child_list[list_begin].next != 0)
					list_begin = child_list[list_begin].next;
				if(child_list[list_begin].tid != c_tid && child_list[list_begin].next == 0)//store new node
				{
					child_list[list_begin].next = child_count++;
					child_list[child_count - 1].tid = c_tid;
				}
			}
			//end, reset c_tid
			c_tid = p_tid;
		}
	}
	//out put
	printf("Analysis based on base number:\n");
	if (ANA_PRINT_USE_LIST == 1)
		ana_meta_loop_fprint(stdout, taxonomyTree, node_table, 1, child_list, 0, total_base_num, false);
	else
		ana_meta_loop_print(taxonomyTree, node_table, 1, child_list, 0, total_base_num, false);
	//printf("</Data>\n");
	//number
	printf("total_mapped_base_number :%ld\n",total_base_num);
	printf("low identity read (identity <= 75%%) number :%ld\t",low_identity_read_num);
	printf("total base %ld\t", low_identity_read_base);
	//step3: end PRO
	free(taxonomyTree);
	fclose(rst_file_srt);
	//fclose(maf_file_srt);
}

typedef struct{
	uint32_t tid;
	uint64_t base;
	uint64_t map_q;
}NODE_BASE_Q;
//all char that are more than '9' are use less than '9'
static int cmp_base_Q(const void *a_, const void *b_){
	NODE_BASE_Q *a = (NODE_BASE_Q*) a_, *b= (NODE_BASE_Q *) b_;
	return a->base < b->base;
}
static void ana_meta_base_M2(char * rst_file_name, char * tax_file_name)
{
	printf("Current read %s\t", rst_file_name);
	//PART 1: load result file
	//char * rst_file_name = argv[1];
	printf("%s\t", rst_file_name);
	FILE * rst_file_srt = xopen(rst_file_name, "rb");

	//PART 2:load taxonomyTree
	//char * tax_file_name = argv[2];
	TAXONOMY_rank * taxonomyTree = NULL;
	int max_tid = taxonTree_rank(tax_file_name, &taxonomyTree);

	//step3: analysis
	NODE_BASE_Q *node_base = xcalloc(max_tid, sizeof(NODE_BASE_Q));
	//CLY_NODE_V cly_v = {0,0,0};
	//CLY_NODE root={"root", 1, "root", 0, NULL, NULL};

	// = {"root", 0. "root", 0, NULL, NULL}
	int total_read_number = 0;
	uint64_t total_base_num = 0;
	uint64_t low_identity_read_num = 0;
	uint64_t low_identity_read_base = 0;
	RST rst;
	int eof_ = 0;//end of file
	float coverage = 0;
	//count all PRIMARY tax
	if(getOneRST(rst_file_srt, &rst) < 0)
		return;
	while(1)
	{
		total_read_number++;
		int read_len = 0;
		//get all results for one read
		int map_Q = rst.MAPQ;
		uint32_t final_tid = ana_get_tid(&rst, max_tid, rst_file_srt, &eof_, taxonomyTree, &read_len, &coverage);
		if(final_tid > 0)
		{
			if(coverage * read_len > MIN_SCORE)
			{
				total_base_num += read_len;
				node_base[final_tid].base += read_len;
				node_base[final_tid].map_q += read_len*map_Q;
				if(coverage < 0.08)
				{
					low_identity_read_base += read_len;
					low_identity_read_num++;
				}
			}
		}
		if(eof_ < 0)//end if file
			break;
	}

	CLY_NODE *node_table = xcalloc(max_tid, sizeof(CLY_NODE));
	CN_CHILD *child_list = xcalloc(max_tid*2, sizeof(CN_CHILD));
	uint32_t child_count = 1;

	NODE_BASE_Q *sort = xmalloc(max_tid* sizeof(NODE_BASE_Q));
	int rst_num = 0;
	for(uint32_t i = 0; i <= max_tid; i++)
	{
		if(node_base[i].base == 0)
			continue;
		sort[rst_num].tid = i;
		sort[rst_num].base = node_base[i].base;
		sort[rst_num++].map_q = node_base[i].map_q;
	}
	qsort(sort,rst_num,sizeof(NODE_BASE_Q), cmp_base_Q);

	//calculate middle nodes
	for(int i = 0; i < rst_num; i++)
	{
		uint32_t c_tid = sort[i].tid;
		node_table[sort[i].tid].weight += node_base[sort[i].tid].base;
		node_table[sort[i].tid].total_mapQ += node_base[sort[i].tid].map_q;
		while(1)
		{
			uint32_t p_tid = taxonomyTree[c_tid].p_tid;
			if(p_tid < 1 || p_tid == 4294967295)//over root node
				break;
			node_table[p_tid].weight += node_base[sort[i].tid].base;//add weight for p_tid
			node_table[p_tid].total_mapQ += node_base[sort[i].tid].map_q;//add weight for p_tid
			//store p_tid->c_tid
			if(node_table[p_tid].child_list_begin == 0)
			{
				node_table[p_tid].child_list_begin = child_count++;
				child_list[child_count - 1].tid = c_tid;
			}
			else
			{
				int list_begin = node_table[p_tid].child_list_begin;
				while(child_list[list_begin].tid != c_tid && child_list[list_begin].next != 0)
					list_begin = child_list[list_begin].next;
				if(child_list[list_begin].tid != c_tid && child_list[list_begin].next == 0)//store new node
				{
					child_list[list_begin].next = child_count++;
					child_list[child_count - 1].tid = c_tid;
				}
			}
			//end, reset c_tid
			c_tid = p_tid;
		}
	}
	//out put
	printf("Analysis based on base number:\n");
	if (ANA_PRINT_USE_LIST == 1)
		ana_meta_loop_fprint(stdout, taxonomyTree, node_table, 1, child_list, 0, total_base_num, true);
	else
		ana_meta_loop_print(taxonomyTree, node_table, 1, child_list, 0, total_base_num, true);
	//printf("</Data>\n");
	//number
	printf("total_mapped_base_number :%ld\n",total_base_num);
	printf("low identity read (identity <= 75%%) number :%ld\t",low_identity_read_num);
	printf("total base %ld\t", low_identity_read_base);
	//step3: end PRO
	free(taxonomyTree);
	fclose(rst_file_srt);
	//fclose(maf_file_srt);
}


static void ana_meta_metamaps_base(char * sam_file_name, char * tax_file_name, char * exchange_file)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_matemaps_file(sam_file_name, temp_file_name, exchange_file);
	//analysis
	ana_meta_base(temp_file_name, tax_file_name);
	//xrm(temp_file_name);
}

//-----------------------------------------uni_v_analysis-----------------------------------------//
//compare minimap rst file with maf file
//USAGE: deSPI3rd cmp cmp_minimap [rst file name] [maf file name] [tax_file_name]
#define uni_v_analysis_NUM 1000
static void uni_v_analysis(char * univ_file_name)
{
	typedef struct {
		uint32_t ref_list;
		uint32_t length;
		//uint64_t offset;
	}UNITIG;
	FILE * uni_v_file = xopen(univ_file_name, "rb");
	uint64_t uni_number = 0;

	xassert(fread(&uni_number, 8, 1, uni_v_file) > 0,"");
	UNITIG *data = xmalloc(uni_number*8);
	xassert(fread(data, 8, uni_number, uni_v_file) == uni_number, "");
	uint32_t count[uni_v_analysis_NUM] = {0};
	uint64_t over_100 = 0;
	for(int i = 0; i < uni_number; i++)
	{
		if(data[i].length < uni_v_analysis_NUM)
			count[data[i].length]++;
		else
			over_100 += data[i].length;
	}

	for(int i = 0; i < uni_v_analysis_NUM; i++)
	{
		printf("%d %d\n", i, count[i]);
	}
	printf("over_%d %ld\n",uni_v_analysis_NUM, over_100);

	fclose(uni_v_file);
}

//-----------------------------------------U/C count-----------------------------------------//
//compare minimap rst file with maf file
//USAGE: deSPI3rd cmp cmp_minimap [rst file name] [maf file name] [tax_file_name]
static void rst_stat(char * rst_file_name)
{
	FILE * rst_file_srt = xopen(rst_file_name, "rb");
	uint32_t rst_number = 0, classified = 0, un_classified = 0;
	RST rst; int eof_rst = getOneRST(rst_file_srt, &rst);
	for(; eof_rst >= 0;){
		eof_rst = getOneRST(rst_file_srt, &rst); rst_number++;
		if(rst.isClassify == 'U')
			un_classified++;
		else
			classified++;
	}
	fprintf(stderr,"total:%d, U:%d,C:%d",rst_number,un_classified,classified);
	fclose(rst_file_srt);
}

//----------------------------------------File name list-----------------------------------------//
//compare minimap rst file with maf file
//USAGE: deSPI3rd cmp cmp_minimap [rst file name] [maf file name] [tax_file_name]
static void file_name(char * ref_file_name)
{
	gzFile fp = xzopen(ref_file_name, "r");
	kseq_t temp;
	temp.f = ks_init(fp);
	while(kseq_read(&temp) > 0)
	{
		//char *tokens;
		strtok(temp.name.s + 4,"|");
		printf("%s\t", temp.name.s);
		printf("%s\n", temp.name.s+4);
		//tokens = strtok(NULL,"|");
		//tokens = strtok(NULL,"|");
		//printf("%s\n", tokens);
		//tokens = strtok(temp.name.s,"|");
		//tokens = strtok(NULL,0);
		//printf("%s\n", tokens);
	}
	gzclose(fp);
}
//-----------------------------------------main function DES-----------------------------------------//
static void ana_meta_des(char * sam_file_name, char * tax_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_des_sam_file(sam_file_name, temp_file_name);
	//analysis
	ana_meta(temp_file_name, tax_file_name);
	xrm(temp_file_name);
}

static void ana_meta_des_base(char * sam_file_name, char * tax_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_des_sam_file(sam_file_name, temp_file_name);
	//analysis
	ana_meta_base_M2(temp_file_name, tax_file_name);
	xrm(temp_file_name);
}

static void ana_meta_cen_base(char * sam_file_name, char * tax_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_CEN_file(sam_file_name, temp_file_name);
	//analysis
	ana_meta_base(temp_file_name, tax_file_name);
	xrm(temp_file_name);
}

static void ana_meta_cen(char * sam_file_name, char * tax_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_CEN_file(sam_file_name, temp_file_name);
	//analysis
	ana_meta(temp_file_name, tax_file_name);
	xrm(temp_file_name);
}

static void ana_meta_kai(char * sam_file_name, char * tax_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_KAI_file(sam_file_name, temp_file_name);
	//analysis
	ana_meta(temp_file_name, tax_file_name);
	xrm(temp_file_name);
}

static void dump_KAI_file_with_length(char * SAM_file_name, char*dump_file_name, int *read_length_list)
{
	//step1: OPEN files
	FILE * SAM_file = xopen(SAM_file_name,"r");
	FILE * dump_file = xopen(dump_file_name,"w");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	RST temp_rst;
	while(1){
		//get SAM
		if(getOnekaiSAM(SAM_file, buff, &temp_rst) < 0)
			break;
		//write MAF
		record_num++;
		long int read_ID = strtol(temp_rst.read_name + 11, NULL, 10);
		fprintf(dump_file, "%s\t%c\t%d\t%d\t%d\n",
				temp_rst.read_name,
				temp_rst.isClassify,
				temp_rst.tid,
				read_length_list[read_ID],
				temp_rst.MAPQ);
	}
	free(buff);
	fclose(SAM_file);
	fclose(dump_file);
}

//read_length_fn: read ID with length, format: XXX XXX
static void ana_meta_kai_base(char * sam_file_name, char * tax_file_name, char *read_length_fn)
{
	//load length file
	FILE * read_length_file = xopen(read_length_fn,"r");
	int *read_length_list = xcalloc(5000000, sizeof(int));//at most 5M reads support
	int read_ID = 0;
	int read_len = 0;
	while(fscanf(read_length_file, "%d %d\n", &read_ID, &read_len) > 0)
		read_length_list[read_ID] = read_len;
	fclose(read_length_file);
	//dump kai file with length
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//temp file
	//dump des sam
	dump_KAI_file_with_length(sam_file_name, temp_file_name, read_length_list);
	//analysis
	ana_meta_base(temp_file_name, tax_file_name);
	xrm(temp_file_name);
}

static void ana_tax_des(char * sam_file_name, uint32_t right_tax, char * tax_file_name ,char *rank)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//dump des sam
	dump_des_sam_file(sam_file_name, temp_file_name);
	//analysis
	ana_tax(temp_file_name, right_tax, tax_file_name, rank);
	xrm(temp_file_name);
}


static void ana_tax_PAF(char * paf_file_name, uint32_t right_tax, char * tax_file_name ,char *rank)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, paf_file_name);
	strcat(temp_file_name, ".temp");
	//dump des sam
	dump_des_PAF_file(paf_file_name, temp_file_name);
	//analysis
	ana_tax(temp_file_name, right_tax, tax_file_name, rank);
	xrm(temp_file_name);
}

static void ana_tax_CEN(char * cen_file_name, uint32_t right_tax, char * tax_file_name ,char *rank)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, cen_file_name);
	strcat(temp_file_name, ".temp");
	//dump des sam
	dump_CEN_file(cen_file_name, temp_file_name);
	//analysis
	ana_tax(temp_file_name, right_tax, tax_file_name, rank);
	xrm(temp_file_name);
}

static void ana_tax_KAI(char * kai_file_name, uint32_t right_tax, char * tax_file_name ,char *rank)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, kai_file_name);
	strcat(temp_file_name, ".temp");
	//dump des sam
	dump_KAI_file(kai_file_name, temp_file_name);
	//analysis
	ana_tax(temp_file_name, right_tax, tax_file_name, rank);
	xrm(temp_file_name);
}

/*******************************************
BLASTn tabular output format 6
1. 	 qseqid 	 query (e.g., gene) sequence id
2. 	 sseqid 	 subject (e.g., reference geome) sequence id
3. 	 pident 	 percentage of identical matches
4. 	 length 	 alignment length
5. 	 mismatch 	 number of mismatches
6. 	 gapopen 	 number of gap openings
7. 	 qstart 	 start of alignment in query
8. 	 qend 	 end of alignment in query
9. 	 sstart 	 start of alignment in subject
10. 	 send 	 end of alignment in subject
11. 	 evalue 	 expect value
12. 	 bitscore 	 bit score
**************************************************************************/

typedef struct blast_RST{
	char 		read_name[50];
	uint32_t 	mapping_length;
	float 		indentify;
	int 		read_st;
	int 		read_ed;
}blast_RST;//24byte

//------------------------------------------KAI-----------------------------------------//
static int getOneBlast(FILE * blast_file, char *buff, blast_RST * rst)
{
	size_t max_l = MAX_BUFF_LEN;
	int read_L = 0;
	read_L = getline(&buff,&max_l,blast_file);
	if(read_L <= 0)
		return -1;
	sscanf(buff,
			"%s\t" 	//1
			"%*s\t"	//2
			"%f"	//3
			"%d"	//4
			"%*d"	//5
			"%*d"	//6
			"%d"	//7
			"%d"	//8
			"%*d"	//9
			"%*d"	//10
			"%*f"	//11
			"%*d",	//12
			rst->read_name,
			&(rst->indentify),
			&(rst->mapping_length),
			&(rst->read_st),
			&(rst->read_ed)
			);
	return 0;
}

static void ana_BLASTN(char * blastn_file_name)
{
	//step1: OPEN files
	FILE * BLASTN_file = xopen(blastn_file_name,"r");
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	//step3: store SAM into RST
	uint32_t record_num = 0;
	blast_RST rst = {0};
	blast_RST old_rst = {0};
	uint64_t total_length = 0;
	int st_l[1000];
	int ed_l[1000];
	int region_n = 0;
	while(1){

		//get SAM
		if(getOneBlast(BLASTN_file, buff, &rst) < 0)
			break;
		if(strcmp(rst.read_name, old_rst.read_name) != 0)//new reads
		{
			st_l[0] = rst.read_st;
			ed_l[0] = rst.read_ed;
			region_n = 1;
			//add length
			total_length += rst.mapping_length;
		}
		else//old reads
		{
			int i = 0;
			for(; i < region_n; i++)//search all regions
				if(rst.read_st <= ed_l[i] && rst.read_ed >= st_l[i])
					break;
			if(i == region_n)
			{
				//store new block
				{
					st_l[region_n] = rst.read_st;
					ed_l[region_n] = rst.read_ed;
					region_n++;
				}
				//add length
				total_length += rst.mapping_length;
			}
			continue;
		}
		record_num++;
		old_rst = rst;
	}
	free(buff);
	fclose(BLASTN_file);
	fprintf(stderr, "%s\t %d\t %ld\n", blastn_file_name, record_num, total_length);
}


//---------------------------------ANALYSIS with filter----------------------------------------------//
//return 'P' when pass, return 'F' when fail
static char get_filter_result(FILE * filter_file, char * read_name)
{
	size_t max_l = 1024;
	char static_BUFF[1024];
	char *buff = static_BUFF;
	char *tokens;
	int reset = false;
	while(1)
	{
		if(getline(&buff,&max_l,filter_file) <= 0)
		{
			fprintf(stderr, "With out filter info! file reset");
			fseek(filter_file, 0, SEEK_SET);
			if(reset == true)
				xassert(0, "With out filter info!");
			reset = true;
			if(getline(&buff,&max_l,filter_file) <= 0)
				xassert(0, "Filter info file no data!");
		}
		tokens = strtok(buff," ");
		if(strcmp(read_name, tokens) == 0)
		{
			tokens = strtok(NULL,"\n");
			return tokens[0];
		}
	}
	return 'F';
}

static void ana_tax_with_filter(char * rst_file_name, uint32_t right_tax, char * tax_file_name, char *rank, char * filter_file_name)
{
	//PART 1: load result file
	//char * rst_file_name = argv[1];
	fprintf(stderr, "%s\t", rst_file_name);
	FILE * rst_file_srt = xopen(rst_file_name, "rb");

	FILE * filter_file = xopen(filter_file_name, "rb");
	//PART 2:get right tax
	//uint32_t right_tax = strtoul(argv[2],0,10);

	//PART 3:load taxonomyTree
	//char * tax_file_name = argv[3];
	TAXONOMY_rank * taxonomyTree = NULL;
	taxonTree_rank(tax_file_name, &taxonomyTree);

	//step4: analysis
	int wrong_alignment = 0;
	int total_read_number = 0;
	int umapped_read_number = 0;
	char old_read_name[READ_NAME_LEN] = {0};
	bool right_alignment = 0;
	int right_aligment_first = 0;
	int right_aligment_second =0;
	RST rst; //int eof_rst = getOneRST(rst_file_srt, &rst);

	uint32_t c_tax;
	if(getOneRST(rst_file_srt, &rst) < 0)
		return;
	while(1)
	{
		char is_filtered = get_filter_result(filter_file, rst.read_name);//'P';

		//get the first result
		if(is_filtered == 'P')
			total_read_number++;
		if(SHOW_DETAIL)
			printf("\n%s ",rst.read_name);
		if(rst.isClassify == 'U')
		{
			if(is_filtered == 'P')
				umapped_read_number++;
			if(SHOW_DETAIL)
				printf("UM");
			int eof_ = getOneRST(rst_file_srt, &rst);
			if(eof_ < 0)
				break;
			else
				continue;
		}

		c_tax = get_tax_by_rank(taxonomyTree, rst.tid, rank);
		if(right_tax == c_tax){
			right_alignment = true;
			if(is_filtered == 'P')
				right_aligment_first++;
			if(SHOW_DETAIL)
				printf("PRI");
		}
		else
			right_alignment = false;
		strcpy(old_read_name, rst.read_name);
		int eof_ = 0;
		//get other results
		while(1)
		{
			eof_ = getOneRST(rst_file_srt, &rst);
			if(eof_ < 0)
				break;
			if(strcmp(old_read_name, rst.read_name) == 0)
			{
				if(right_alignment == true)
					continue;
				c_tax = get_tax_by_rank(taxonomyTree, rst.tid, rank);
				if(right_tax == c_tax){
					right_alignment = true;
					if(is_filtered == 'P')
						right_aligment_second++;
					if(SHOW_DETAIL)
						printf("SEC");
				}
			}
			else
				break;
		}
		if(eof_ < 0)
			break;

		if(right_alignment == false)
			if(is_filtered == 'P')
				wrong_alignment++;
	}
	//title:
	if(ANA_SHOW_TITLE)
	{
		fprintf(stderr,"[1] total reads number\n");
		fprintf(stderr,"[2] total unmapped read number\n");
		fprintf(stderr,"[3] total right PRIMARY\n");
		fprintf(stderr,"[4] total right ALL[PRIMARY + SEN + SUP]\n");

		fprintf(stderr,"[5] unmapped rate\n");
		fprintf(stderr,"[6] SEN PRI(RIGHT/ALL)[PRI]\n");
		fprintf(stderr,"[7] ACC PRI(RIGHT/(RIGHT+WRONG))[PRI]\n");

		fprintf(stderr,"[8] SEN ALL(RIGHT/ALL)[PRI + SEN + SUP]\n");
		fprintf(stderr,"[9] ACC ALL(RIGHT/(RIGHT+WRONG))[PRI + SEN + SUP]\n");
	}

	//number
	fprintf(stderr,"%d\t",total_read_number);
	fprintf(stderr,"%d\t",umapped_read_number);
	fprintf(stderr,"%d\t",right_aligment_first);
	fprintf(stderr,"%d\t",right_aligment_second + right_aligment_first);

	fprintf(stderr,"%f%%\t", (float)umapped_read_number/total_read_number *100);
	fprintf(stderr,"%f%%\t", (float)right_aligment_first/total_read_number *100);
	fprintf(stderr,"%f%%\t", (float)(right_aligment_first)/(total_read_number - umapped_read_number)*100);

	fprintf(stderr,"%f%%\t",(float)(right_aligment_second + right_aligment_first)/total_read_number * 100);
	fprintf(stderr,"%f%%\n",(float)(right_aligment_second + right_aligment_first)/(total_read_number - umapped_read_number)*100);

	//step3: end PRO
	free(taxonomyTree);
	fclose(rst_file_srt);
	fclose(filter_file);
}

//for sam files

static void ana_tax_SAM_filter(char * sam_file_name, uint32_t right_tax, char * tax_file_name ,char *rank, char * filter_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, sam_file_name);
	strcat(temp_file_name, ".temp");
	//dump des sam
	dump_des_sam_file(sam_file_name, temp_file_name);
	//analysis
	ana_tax_with_filter(temp_file_name, right_tax, tax_file_name, rank, filter_file_name);
	xrm(temp_file_name);
}

//for paf files
static void ana_tax_PAF_filter(char * paf_file_name, uint32_t right_tax, char * tax_file_name ,char *rank, char * filter_file_name)
{
	char temp_file_name[1024];
	strcpy(temp_file_name, paf_file_name);
	strcat(temp_file_name, ".temp");
	//dump des sam
	dump_des_PAF_file(paf_file_name, temp_file_name);
	//analysis
	ana_tax_with_filter(temp_file_name, right_tax, tax_file_name, rank, filter_file_name);
	xrm(temp_file_name);
}

//for dump files (cen and kai)
static void ana_tax_DUMP_filter(char * dump_file_name, uint32_t right_tax, char * tax_file_name ,char *rank, char * filter_file_name)
{
	ana_tax_with_filter(dump_file_name, right_tax, tax_file_name, rank, filter_file_name);
}

void file_cmp_bin(char * file_name1, char * file_name2)
{
	FILE * f1 = xopen(file_name1, "rb");
	FILE * f2 = xopen(file_name2, "rb");

	while(1)
	{
		uint8_t d1, d2, r1, r2;
		r1 = fread(&d1, 1, 1, f1);
		r2 = fread(&d2, 1, 1, f2);
		xassert(r1 == r2, "");
		xassert(d1 == d2, "");
		if(r1 == 0)
			return;
	}
	fclose(f1);
	fclose(f2);
}

//-----------------------------------------main-----------------------------------------//

//-----------------------------------------mark SAM file-----------------------------------------//
static int mark_SAM(char * sam_file_name, char * tax_file_name ,char *rank)
{
	//PART 1: load result file
	fprintf(stderr, "%s\t", sam_file_name);
	FILE * rst_file_srt = xopen(sam_file_name, "rb");

	//PART 2:load taxonomyTree
	TAXONOMY_rank * taxonomyTree = NULL;
	taxonTree_rank(tax_file_name, &taxonomyTree);

	//PART3: line buff
	size_t max_l = MAX_BUFF_LEN;
	char *ori_str = (char*)malloc(MAX_BUFF_LEN);
	char *buff = (char*)malloc(MAX_BUFF_LEN);
	char *tokens;
	uint32_t tid = 0;
	while(getline(&buff,&max_l,rst_file_srt) > 0)
	{
		if(buff[0] == '@')
			continue;
		strcpy(ori_str, buff);
		//get SAM tid
		tokens = strtok(buff,"\t");//read name
		tokens = strtok(NULL,"\t");//flag
		tokens = strtok(NULL,"\t");//refNAME
		if(tokens[0] == '*')
			tid = 0;
		else{
			char * ref_tokens = tokens;
			//ignore the POS part
			tokens = strtok(NULL,"\t");
			//get MAQ part
			tokens = strtok(NULL,"\t");
			//for the ref name part
			{
				//ignore 'tid|'
				ref_tokens = strtok(ref_tokens,"|");
				//get tid
				ref_tokens = strtok(NULL,"|");
				tid = strtoul(ref_tokens,NULL,10);
			}
		}
		ori_str[100] = 0;
		if(tid == 0)
			printf("0\t");
		else
			printf("%d\t", get_tax_by_rank(taxonomyTree, tid, rank));
		printf("%s\n", ori_str);
	}
	//step3: end PRO
	free(taxonomyTree);
	fclose(rst_file_srt);
	return 0;
}

static void count_base(char * fastq_file_name)
{
	gzFile fp = xzopen(fastq_file_name, "r");
	kstream_t *_fp = ks_init(fp);
	kseq_t seq = {0};
	seq.f = _fp;
	uint64_t total_length = 0;
	uint64_t read_number = 0;
	while( kseq_read(&seq) >= 0)
	{
		read_number ++;
		total_length += seq.seq.l;
	}
	gzclose(fp);
	fprintf(stderr, "%s read number: %ld base number %ld ( %f Mbp)\n", fastq_file_name, read_number, total_length, (float)total_length/1000000);
}

static void get_read_by_NAME(char * fastq_file_name , char * read_name)
{
	gzFile fp = xzopen(fastq_file_name, "r");
	kstream_t *_fp = ks_init(fp);
	kseq_t seq = {0};
	seq.f = _fp;
	while( kseq_read(&seq) >= 0)
	{
		if(strcmp(seq.name.s, read_name) == 0)
			printf(	"@%s %s\n"
					"%s\n"
					"+\n"
					"%s\n",
					seq.name.s,
					seq.comment.s,
					seq.seq.s,
					seq.qual.s);
	}
	gzclose(fp);
}

static void reverse_read(char * fastq_file_name)
{
	gzFile fp = xzopen(fastq_file_name, "r");
	kstream_t *_fp = ks_init(fp);
	kseq_t seq = {0};
	seq.f = _fp;
	uint64_t total_length = 0;
	uint64_t read_number = 0;
	while( kseq_read(&seq) >= 0)
	{
		read_number ++;
		total_length += seq.seq.l;
		for(int j = 0; j < seq.seq.l; j++)
		{
			char new_char = 'X';
			switch(seq.seq.s[seq.seq.l - j - 1])
			{
			case 'A': new_char = 'T'; break;
			case 'C': new_char = 'G'; break;
			case 'G': new_char = 'C'; break;
			case 'T': new_char = 'A'; break;
			}
			fprintf(stderr, "%c", new_char);
		}
		fprintf(stderr, "\n\n\n");
	}
	gzclose(fp);
	fprintf(stderr, "%s read number: %ld base number %ld ( %f Mbp)\n", fastq_file_name, read_number, total_length, (float)total_length/1000000);
}

static void split_fastq(char * fastq_file_name, int begin, int step)
{
	gzFile fp = xzopen(fastq_file_name, "r");
	kstream_t *_fp = ks_init(fp);
	kseq_t seq = {0};
	seq.f = _fp;
	uint64_t total_length = 0;
	uint64_t read_number = 0;
	while( kseq_read(&seq) >= 0)
	{
		if((read_number >= begin) && ((read_number - begin) % step == 0))
		{
			printf(	"@%s %s\n"
					"%s\n"
					"+\n"
					"%s\n",
					seq.name.s,
					seq.comment.s,
					seq.seq.s,
					seq.qual.s);
			total_length += seq.seq.l;
		}
		read_number ++;
	}
	gzclose(fp);
	fprintf(stderr, "%s read number: %ld base number %ld ( %f Mbp)\n", fastq_file_name, read_number, total_length, (float)total_length/1000000);
}

static void get_centrifuge_map_file(char * fasta_file_name)
{
	gzFile fp = xzopen(fasta_file_name, "r");
	kstream_t *_fp = ks_init(fp);
	kseq_t seq = {0};
	seq.f = _fp;
	int tid = 0;
	while( kseq_read(&seq) >= 0)
	{
		for(int i = 4; ;i++)
			if(seq.name.s[i] == '|')
			{
				seq.name.s[i] = 0;
				break;
			}

		printf("%s    ", seq.name.s);
		char * ref_tokens = seq.name.s;
		//for the ref name part
		//ignore 'tid|'
		ref_tokens = strtok(ref_tokens,"|");
		//get tid
		ref_tokens = strtok(NULL,"|");
		tid = strtoul(ref_tokens,NULL,10);
		printf("%d\n", tid);
	}
	gzclose(fp);
}

static int is_low_complex(char * st, int len)
{
	int number_A = 0;
	int number_C = 0;
	int number_G = 0;
	int number_T = 0;
	int MAX = len*0.7;

	for(int i = 0; i < len; i ++)
	{
		switch(st[i])
		{
		case 'A':
		case 'a':
			number_A++; break;
		case 'C':
		case 'c':
			number_C++; break;
		case 'G':
		case 'g':
			number_G++; break;
		case 'T':
		case 't':
			number_T++; break;
		}
	}

	if(number_A >= MAX ||number_C >= MAX ||number_G >= MAX ||number_T >= MAX)
	{
		return true;
	}
	return false;

}

#define ONLY_READ_FILTER_INFO 1
#define READ_FILTER_MIN_LEN 1000

static void pacbio_filter(char * fastq_file_name)
{
	gzFile fp = xzopen(fastq_file_name, "r");
	kstream_t *_fp = ks_init(fp);
	kseq_t seq = {0};
	seq.f = _fp;
	int read_number = 0;
	int filterd_read = 0;
	while( kseq_read(&seq) >= 0)
	{
		read_number ++;
		int pass_filter = true;
		if(seq.seq.l < READ_FILTER_MIN_LEN)
			pass_filter = false;
		else
		{
			char* str = seq.seq.s;
			int max_normal_length = 0;
			int abnormal_length = 0;
			for(int i = 0; i < seq.seq.l - 28; i++)
				if(is_low_complex(str + i, 27))
					abnormal_length++;
			max_normal_length = seq.seq.l - abnormal_length;
			if(max_normal_length < READ_FILTER_MIN_LEN)
				pass_filter = false;
		}

		if(ONLY_READ_FILTER_INFO)//output reads info
		{
			if(pass_filter == false)
			{
				filterd_read++;
				printf("%s F\n", seq.name.s);
			}
			else
				printf("%s P\n", seq.name.s);
		}
		else//output passed reads
		{
			if(pass_filter == false)
				filterd_read++;
			else
				printf("@%s\n%s\n+%s\n%s\n", seq.name.s, seq.seq.s, seq.name.s, seq.qual.s);
		}
	}
	gzclose(fp);
	fprintf(stderr, "file name: %s total number: %d filtered number: %d\n", fastq_file_name, read_number, filterd_read);
}

static void fastq_to_fasta(char *fastq_file)
{
	gzFile fp = xzopen(fastq_file, "r");
	kstream_t *_fp = ks_init(fp);

	kseq_t temp = {0};
	temp.f = _fp;
	while( kseq_read(&temp) >= 0)
	{
		printf(">%s %s\n", temp.name.s, temp.comment.s);
		printf("%s\n", temp.seq.s);
	}
}

static void fastq_to_name(char *fastq_file)
{
	gzFile fp = xzopen(fastq_file, "r");
	kstream_t *_fp = ks_init(fp);

	kseq_t temp = {0};
	temp.f = _fp;
	while( kseq_read(&temp) >= 0)
		printf("%s %s\n", temp.name.s, temp.comment.s);
}

#define ANALYSIS_MAIN_COMMAND "analysis"
static int cmp_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME);
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact:   %s\n\n", CONTACT);
	fprintf(stderr, "  Usage:     %s %s <command> [file]\n\n", PACKAGE_NAME, ANALYSIS_MAIN_COMMAND);
	fprintf(stderr, "  Command list: \n");
	fprintf(stderr, "    %s ana_meta    	 [SAM_file.sam] [node.dmp]\n", ANALYSIS_MAIN_COMMAND);
	fprintf(stderr, "    %s ana_meta_base    [SAM_file.sam] [node.dmp]\n", ANALYSIS_MAIN_COMMAND);
	fprintf(stderr, "    %s count_base    	 [FASTQ_file.fq]\n", ANALYSIS_MAIN_COMMAND);
	fprintf(stderr, "    %s split_fastq    	 [FASTQ_file.fq] [start_number] [step_length]\n", ANALYSIS_MAIN_COMMAND);
	fprintf(stderr, "    %s fastq_to_fasta   [FASTQ_file.fq] \n", ANALYSIS_MAIN_COMMAND);
	//fprintf(stderr, "    %s ana_species [SAM_file.sam] [species taxonomy] [node.dmp]\n", ANALYSIS_MAIN_COMMAND);
	//fprintf(stderr, "    %s ana_genus   [SAM_file.sam] [genus taxonomy]   [node.dmp]\n", ANALYSIS_MAIN_COMMAND);
	//fprintf(stderr, "    %s ana_sam     [SAM_file.sam] [XXX taxonomy] [node.dmp] [rank]\n", ANALYSIS_MAIN_COMMAND);
	fprintf(stderr, "  Basic:\n");
	fprintf(stderr, "    [SAM_file.sam]  FILE  Classify file generated from \"classify\" command\n");
	fprintf(stderr, "    [XXX taxonomy]  INT   taxonomy you want to test in the result file\n");
	fprintf(stderr, "    [rank]          STR   rank you want to test in the result file\n");
	fprintf(stderr, "    [node.dmp]      FILE  node.dmp file download from: \n");
	fprintf(stderr, "                          ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n");
    fprintf(stderr, "\n");
	return 0;
}

#define ASSERT_PARA(parameter_number) if(argc != parameter_number) { \
										fprintf(stderr,"\nERROR: File name parameters is not enough, %d files needed, read usage for more information\n\n", parameter_number - 2);\
										return cmp_usage();}
//you need to first dump files first, than compare
int simDataTest(int argc, char *argv[])
{
	if (0 == strcmp(argv[argc - 1], "print_list")) {ANA_PRINT_USE_LIST = 1; fprintf(stderr, "ANA_PRINT_USE_LIST = 1\n");} 
	if		(argc <= 1)						  			{cmp_usage();}
	else if	(0 == strcmp(argv[1], "ana_meta"))			{ana_meta_des(	argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "ana_meta_base"))		{ana_meta_des_base(	argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "ana_meta_cen_base"))	{ana_meta_cen_base(	argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "ana_meta_kai"))		{ana_meta_kai(	argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "ana_meta_kai_base"))	{ana_meta_kai_base(	argv[2], argv[3],  argv[4]);}
	else if	(0 == strcmp(argv[1], "ana_meta_cen"))		{ana_meta_cen(	argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "ana_species"))		{ana_tax_des(	argv[2], strtoul(argv[3],0,10), argv[4], "species");}
	else if	(0 == strcmp(argv[1], "ana_genus"))			{ana_tax_des(	argv[2], strtoul(argv[3],0,10), argv[4], "genus");}
	else if	(0 == strcmp(argv[1], "mark_genus"))		{mark_SAM(		argv[2], argv[3], "genus");}
	else if	(0 == strcmp(argv[1], "ana_meta_rst"))		{ana_meta(		argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "ana_matemaps_base"))	{ana_meta_metamaps_base(	argv[2], argv[3], argv[4]);}
	//else if	(0 == strcmp(argv[1], "ana_species_rst"))	{ana_tax(		argv[2], strtoul(argv[3],0,10), argv[4], "species");}
	//else if	(0 == strcmp(argv[1], "ana_genus_rst"))		{ana_tax(		argv[2], strtoul(argv[3],0,10), argv[4], "genus");}
	//else if	(0 == strcmp(argv[1], "ana_rank_rst"))		{ana_tax(		argv[2], strtoul(argv[3],0,10), argv[4], argv[5]);}
	else if	(0 == strcmp(argv[1], "ana_sam"))			{ana_tax_des(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5]);}
	else if	(0 == strcmp(argv[1], "ana_paf"))			{ana_tax_PAF(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5]);}
	else if	(0 == strcmp(argv[1], "ana_cen"))			{ana_tax_CEN(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5]);}//centrifuge
	else if	(0 == strcmp(argv[1], "ana_kai"))			{ana_tax_KAI(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5]);}//kaiju
	else if	(0 == strcmp(argv[1], "ana_BLASTN"))		{ana_BLASTN(	argv[2]);}//kaiju
	else if	(0 == strcmp(argv[1], "ana_dump_filter"))	{ana_tax_DUMP_filter(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5], argv[6]);}
	else if	(0 == strcmp(argv[1], "ana_sam_filter"))	{ana_tax_SAM_filter(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5], argv[6]);}
	else if	(0 == strcmp(argv[1], "ana_paf_filter"))	{ana_tax_PAF_filter(	argv[2], strtoul(argv[3],0,10), argv[4], argv[5], argv[6]);}
	else if	(0 == strcmp(argv[1], "count_base"))		{count_base(	argv[2]);}
	else if	(0 == strcmp(argv[1], "get_read_by_NAME"))	{get_read_by_NAME(	argv[2], argv[3]);}
	else if	(0 == strcmp(argv[1], "reverse_read"))		{reverse_read(	argv[2]);}
	else if	(0 == strcmp(argv[1], "cen_map"))			{get_centrifuge_map_file(	argv[2]);}
	else if	(0 == strcmp(argv[1], "split_fastq"))		{split_fastq(	argv[2], strtoul(argv[3],0,10), strtoul(argv[4],0,10));}
	else if	(0 == strcmp(argv[1], "pacbio_filter"))		{pacbio_filter(	argv[2]);}
	else if	(0 == strcmp(argv[1], "fastq_to_fasta"))	{fastq_to_fasta(	argv[2]);}
	else if	(0 == strcmp(argv[1], "fastq_to_name"))		{fastq_to_name(	argv[2]);}
	//else if	(0 == strcmp(argv[1], "dump_sam"))			{ASSERT_PARA(3);dump_sam(argv[2]);}
	//else if	(0 == strcmp(argv[1], "dump_maf"))			{ASSERT_PARA(3);dump_maf(argv[2]);}
	//else if	(0 == strcmp(argv[1], "dump_des"))			{ASSERT_PARA(3);dump_des(argv[2]);}
	else if	(0 == strcmp(argv[1], "ana_univ"))			{ASSERT_PARA(3);uni_v_analysis(argv[2]);}
	//else if	(0 == strcmp(argv[1], "cmp_rst")) 			{ASSERT_PARA(5);rst_ana(argv[2],argv[3],argv[4]);}
	else if	(0 == strcmp(argv[1], "rst_stat"))			{ASSERT_PARA(3);rst_stat(argv[2]);}
	else if	(0 == strcmp(argv[1], "file_name"))			{ASSERT_PARA(3);file_name(argv[2]);}
	else if	(0 == strcmp(argv[1], "file_cmp"))			{file_cmp_bin(argv[2], argv[3]);}
	else					   							{fprintf(stderr, "command [%s] unsupported!\n\n", argv[1]); cmp_usage();}
	return 0;
}
