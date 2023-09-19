#ifndef DUANRAN_H
#define DUANRAN_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "lib/utils.h"
#include "lib/desc.h"
#include "cly.h"

#define MAX_BUFF_LEN 10000000//1000K
#define READ_NAME_LEN 100
#define MIN_read_N 0
#define MAX_read_N MAX_uint64_t

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

typedef struct
{
	kseq_t *seqs;
	cly_r *results;
	long int readNum;
	void *share_data_pointer;
} CLASSIFY_THREAD_DATA;

typedef struct
{
	// shared
	DA_IDX *idx;
	kstream_t *_fp;
	MAP_opt *o;
	// for each pipeline thread
	Classify_buff_pool *buff;
	CLASSIFY_THREAD_DATA *data;
} CLASSIFY_SHARE_DATA;

typedef struct{
	uint32_t tid;
	long int count;
}COUNT_SORT;

typedef struct {
	char 		read_name[READ_NAME_LEN];
	char 		isClassify;
	uint32_t 	tid;
	uint32_t 	read_length;
	uint8_t 	MAPQ;
	uint32_t 	score;
	char*		seq;
}RST;//24byte

typedef struct RM_buffer
{
	FILE *read_classify_input_tmpfile;
	FILE *read_classify_output_tmpfile;
	MAP_opt map_opt;
	CLASSIFY_THREAD_DATA *classify_thread_data;
	Classify_buff_pool *classify_buff_pool;
	CLASSIFY_SHARE_DATA classify_share_data;

	FILE *meta_analysis_input_tmpfile;
	FILE *meta_analysis_dump_tmpfile;
	FILE *meta_analysis_output_tmpfile;
	char *getOneSAM_buff;
	uint64_t *node_count;
	CLY_NODE *node_table;
	CN_CHILD *child_list;
	COUNT_SORT *sort;

} RM_buffer;

#endif