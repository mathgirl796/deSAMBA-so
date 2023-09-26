/*
 * classify_main.c
 *
 *  Created on: 2018-5-14
 *      Author: fenghe
 */
#include <string.h>
#include "lib/desc.h"
#include <getopt.h>
#include "lib/utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "lib/kthread.h"
#include "lib/kvec.h"
#include "zlib.h"

#include "bwt.h"
#include "cly.h"
#include "cly_mt.h"
#define N_NEEDED 5000
#define MAX_read_size 10000000 // 10M
#define MAX_HUMAN_SNAPSHOT_LEN (64 * 1024)
uint64_t total_sequences = 0;



long int read_reads(kstream_t *_fp, kseq_t *_seqs, long int n_needed)
{
	kseq_t *temp = _seqs;
	long int i, rst = 0, total_length = 0;
	for (i = 0; i < n_needed && total_length < MAX_read_size; ++i)
	{
		temp[i].f = _fp;
		rst = kseq_read(temp + i);
		if (rst < 0)
			break;
		total_length += temp[i].seq.l;
	}
	total_sequences += i;
	return i;
}

static char primary_string[3][4] = {"PRI", "SEC", "SUP"};

static inline void print_hit(chain_item *c, REF_INFO *r_i, int rst_cnt, FILE *outfile)
{
	fprintf(outfile,
			"%3d "
			"%s "
			"%s "
			"%20s "

			"ts:%-10d "
			"te:%-10d "
			"qs:%-10d "
			"qe:%-10d "

			// debug value:
			//"%6d "
			//"%5d "
			"%-5d\t"
			"%d\t"
			//"%d\t"
			//"%d\t"
			//"%s\t"
			//"ln:%-10d\t"
			//"%d\t"
			"\n",

			rst_cnt,
			primary_string[c->primary - 1],
			(c->direction) ? "F" : "R",
			r_i[c->ref_ID].ref_name,

			c->t_st,
			c->t_ed,
			c->q_st,
			c->q_ed,
			// debug value:
			// c->ref_ID,
			// c->chain_id,
			c->sum_score,
			c->indel
			//(c->t_ed - c->t_st) - (c->q_ed - c->q_st)
			// c->anchor_number,
			//(c->with_top_anchor)?"WT":"NT",
			// c->q_ed - c->q_st,
			// c->pri_index
	);
}

void print_anchor(Anchor *anchor_b, REF_INFO *r_i, FILE *outfile)
{
	// basic part
	fprintf(outfile,
			"%c "	  // dir
			"%20s\t"  // ref ID
			"%5d\t"	  // index in read
			"%10u\t", // ref-offset
			(anchor_b->direction) ? 'F' : 'R',
			r_i[anchor_b->ref_ID].ref_name,
			anchor_b->index_in_read,
			anchor_b->ref_offset);

	// debug_part
	if (DEBUG)
		fprintf(outfile,
				"%4d\t" // seed ID
				"%4d\t" // chain ID
				"%c "	// useless
				"%d\t"
				"%d\t"
				"%d\t"
				"%d\t"
				"%c ",
				anchor_b->seed_ID,
				anchor_b->chain_id,
				(anchor_b->anchor_useless) ? '*' : ' ',
				anchor_b->a_m.left_len,
				anchor_b->a_m.left_ED,
				anchor_b->a_m.rigt_len,
				anchor_b->a_m.rigt_ED,
				(anchor_b->duplicate) ? 'D' : 'N');
	// am part
	fprintf(outfile,
			"%d\t"
			"%d\t"
			"\n",
			anchor_b->a_m.score,
			anchor_b->a_m.mtch_len); //,
}

static int inline cmp_anchor(const void *a_, const void *b_)
{
	Anchor *a = (Anchor *)a_, *b = (Anchor *)b_;
	if (a->chain_id != b->chain_id)
		return a->chain_id - b->chain_id;
	return a->index_in_read - b->index_in_read;
}

#define SHOW_ANCHOR_WHEN_NO_RST 0
void output_one_result_des(DA_IDX *idx, cly_r *p_rst, MAP_opt *o)
{
	REF_INFO *r_i = idx->r_i_v.a;
	// NAME:
	fprintf(o->outfile,
			"%s\t"
			"%s\t"
			"%s\t"
			"%ld\t"
			"n_rst:[%ld]\t"
			"n_anc:[%ld]\t"
			"\n",
			p_rst->read->name.s,
			(p_rst->hit.n) ? "CLASSIFY" : "UNCLASSIFY",
			(p_rst->fast_classify) ? "FAST" : "SLOW",
			p_rst->read->seq.l,
			p_rst->hit.n,
			p_rst->anchor_v.n);
	// RESULT
	chain_item *c_s = p_rst->hit.a, *c_e = c_s + p_rst->hit.n;
	int rst_cnt = 0;
	for (chain_item *c = c_s; c < c_e; c++) // primary + sup
		if (c->pri_index == 0)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	for (chain_item *c = c_s; c < c_e; c++) // secondary
		if (c->pri_index > 0 && c->pri_index <= o->max_sec_N)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	// ANCHORs
	// sort anchor
	if (o->show_anchor == true || (rst_cnt == 0 && SHOW_ANCHOR_WHEN_NO_RST)) // when result == 0, show anchors
	{
		qsort(p_rst->anchor_v.a, p_rst->anchor_v.n, sizeof(Anchor), cmp_anchor);
		Anchor *anchor_b = p_rst->anchor_v.a, *anchor_e = anchor_b + p_rst->anchor_v.n;
		for (; anchor_b < anchor_e; anchor_b++)
		{
			if (anchor_b->a_m.score == 0)
				continue;
			print_anchor(anchor_b, r_i, o->outfile);
		}
	}
	fprintf(o->outfile, "\n");
}

void output_one_result_full(DA_IDX *idx, cly_r *p_rst, MAP_opt *o)
{
	REF_INFO *r_i = idx->r_i_v.a;
	// NAME:
	fprintf(o->outfile,
			"%s\t"
			"%s\t"
			"%s\t"
			"%ld\t"
			"n_rst:[%ld]\t"
			"n_anc:[%ld]\t"
			"\n",
			p_rst->read->name.s,
			(p_rst->hit.n) ? "CLASSIFY" : "UNCLASSIFY",
			(p_rst->fast_classify) ? "FAST" : "SLOW",
			p_rst->read->seq.l,
			p_rst->hit.n,
			p_rst->anchor_v.n);
	// RESULT
	chain_item *c_s = p_rst->hit.a, *c_e = c_s + p_rst->hit.n;
	int rst_cnt = 0;
	for (chain_item *c = c_s; c < c_e; c++) // primary + sup
		if (c->pri_index == 0)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	for (chain_item *c = c_s; c < c_e; c++) // secondary
		if (c->pri_index > 0)
			print_hit(c, r_i, rst_cnt++, o->outfile);
	// ANCHORs
	if ((o->show_anchor == true && rst_cnt == 0)) // when rst_cnt == 0; always show anchors
	{
		qsort(p_rst->anchor_v.a, p_rst->anchor_v.n, sizeof(Anchor), cmp_anchor); // sort anchor
		Anchor *anchor_b = p_rst->anchor_v.a, *anchor_e = anchor_b + p_rst->anchor_v.n;
		for (; anchor_b < anchor_e; anchor_b++)
		{
			if (anchor_b->a_m.score == 0)
				continue;
			print_anchor(anchor_b, r_i, o->outfile);
		}
	}
	fprintf(o->outfile, "\n");
}

void output_one_result_sam(DA_IDX *idx, cly_r *p_rst, int output_seq, MAP_opt *o)
{
	REF_INFO *r_i = idx->r_i_v.a;
	char star[2] = "*", *seq_s = (output_seq) ? p_rst->read->seq.s : star, *qual_s = (output_seq) ? p_rst->read->qual.s : star;
	// unmapped
	if (p_rst->hit.n == 0)
	{ // NAME:
		fprintf(o->outfile, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t"
							"%s\t"
							"%s\t"
							"\n",
				p_rst->read->name.s,
				seq_s,
				qual_s);
		return;
	}
	// mapped
	// primary
	uint32_t read_l = p_rst->read->seq.l;
	chain_item *c_s = p_rst->hit.a, *c_e = c_s + p_rst->hit.n;
	int flag = c_s->direction ? 0 : 0x10; // direction
	int mapQ_PRI = 0;
	if (p_rst->hit.n == 1 || (c_s->sum_score - c_s[1].sum_score > 5))
		mapQ_PRI = 30;
	else
		mapQ_PRI = (c_s->sum_score - c_s[1].sum_score) << 2;
	fprintf(o->outfile,
			"%s\t"
			"%d\t"
			"%s\t"
			"%d\t"
			"%d\t"
			"%dS%dM%dS\t" // CIGAR
			"*\t0\t0\t"
			"%s\t"
			"%s\t"
			"AS:i:%d\t" // number of mapping 9-mers, used as score
			//"di:i:%d\t"//sum number of big deletion and insertions(only an es30timated value)
			"\n",
			p_rst->read->name.s,
			flag,
			r_i[c_s->ref_ID].ref_name,
			c_s->t_st,
			mapQ_PRI,
			c_s->q_st, c_s->q_ed - c_s->q_st, read_l - c_s->q_ed,
			seq_s,
			qual_s,
			c_s->sum_score
			// c_s->indel
	);
	// supplementary and secondary
	int rst_cnt = 0;
	for (int loop = 0; loop <= 1; loop++)
	{
		for (chain_item *c = c_s + 1; c < c_e; c++) // sup + sec
		{
			int show_rst = false;
			int flag = c->direction ? 0 : 0x10; // direction
			int mapQ = 0;
			if (loop == 0 && c->pri_index == 0) // supplementary
			{
				show_rst = true;
				flag += 0x800; // supplementary alignment
				mapQ = MIN(30, mapQ_PRI);
			}
			else if (loop == 1 && c->pri_index > 0 && c->pri_index <= o->max_sec_N) // secondary
			{
				show_rst = true;
				flag += 0x100; // secondary alignment
			}
			if (show_rst == true)
			{
				rst_cnt++;
				fprintf(o->outfile,
						"%s\t"
						"%d\t"
						"%s\t"
						"%d\t"
						"%d\t"
						"%d%c%dM%d%c\t" // CIGAR
						"*\t0\t0\t"
						"*\t"
						"*\t"
						"AS:i:%d\t" // number of mapping 9-mers, used as score
						//"di:i:%d\t"//sum number of big deletion and insertions(only an estimated value)
						"\n",
						p_rst->read->name.s,
						flag,
						r_i[c->ref_ID].ref_name,
						c->t_st,
						mapQ,
						c->q_st, (loop == 0) ? 'H' : 'S', c->q_ed - c->q_st, read_l - c->q_ed, (loop == 0) ? 'H' : 'S',
						c->sum_score
						// c->indel
				);
			}
		}
	}
}

#define OUTPUT_MODE_SAM 1
#define OUTPUT_MODE_SAM_FULL 2
#define OUTPUT_MODE_DES 3
#define OUTPUT_MODE_DES_FULL 4
void output_results(DA_IDX *idx, cly_r *results, long int n_results, MAP_opt *o)
{
	cly_r *p_rst = results, *p_e_rst = p_rst + n_results;
	if (o->out_format == OUTPUT_MODE_SAM)
		for (; p_rst < p_e_rst; ++p_rst)				 // out put SAM results ,but with out read sequence
			output_one_result_sam(idx, p_rst, false, o); // the last 0 occupy the MAPQ part
	else if (o->out_format == OUTPUT_MODE_SAM_FULL)
		for (; p_rst < p_e_rst; ++p_rst) // out put SAM results
			output_one_result_sam(idx, p_rst, true, o);
	else if (o->out_format == OUTPUT_MODE_DES)
		for (; p_rst < p_e_rst; ++p_rst)		  // out put DES results
			output_one_result_des(idx, p_rst, o); // the last 0 occupy the MAPQ part
	else if (o->out_format == OUTPUT_MODE_DES_FULL)
		for (; p_rst < p_e_rst; ++p_rst)		   // out put FULL results
			output_one_result_full(idx, p_rst, o); // the last 0 occupy the MAPQ part
}

extern void classify_seq(kseq_t *read, DA_IDX *idx, cly_r *results, Classify_buff_pool *buff);

static void inline worker_for(void *_data, long data_index, int thread_index)
{ // kt_for() callback
	CLASSIFY_THREAD_DATA *d = (CLASSIFY_THREAD_DATA *)_data;
	CLASSIFY_SHARE_DATA *s = (CLASSIFY_SHARE_DATA *)(d->share_data_pointer);
	classify_seq(d->seqs + data_index,
				 s->idx, d->results + data_index,
				 s->buff + thread_index);
}

void *classify_pipeline(void *shared, int step, int tid, void *_data)
{
	CLASSIFY_SHARE_DATA *s = (CLASSIFY_SHARE_DATA *)shared;
	// step0: read read data from files; step1: process; step2: output result
	if (step == 0)
	{
		if ((s->data[tid].readNum = read_reads(s->_fp, s->data[tid].seqs, N_NEEDED)))
			return (void *)1;
	}
	else if (step == 1)
	{
		kt_for(s->o->thread_num, worker_for, s->data + tid, s->data[tid].readNum);
		return (void *)1;
	}
	else if (step == 2)
	{
		output_results(s->idx, s->data[tid].results, s->data[tid].readNum, s->o);
		return (void *)1;
	}
	return 0;
}

// function usage
/**
 *  int Q_MEM[Q_MEM_MAX];
 *  int Q_LV[MAX_LV_WRONG][MAX_LV_R_LEN];
 *  double P_E = 0.15;
 *  uint64_t L_REF = (0x1<<30);
 *  L_REF *= 10;
 *
 *  calculate_MAPQ_TABLE(Q_MEM, Q_LV, P_E, L_REF);
 *  int (*Q_CONFLICT)[Q_CONF_MAX][Q_CONF_MAX][3] = (int(*)[Q_CONF_MAX][Q_CONF_MAX][3])malloc(sizeof(int)*Q_CONF_MAX*Q_CONF_MAX*Q_CONF_MAX*3);
 *  calculate_ANCHOR_CONFLICT_MAPQ_TABLE(Q_CONFLICT);
 *
 */
void calculate_MAPQ_TABLE(
	int *Q_MEM,
	int (*Q_LV)[MAX_LV_R_LEN],
	double P_E, uint64_t L_REF)
{
	// step0: scores
	double REF_SIZE_PUNALTY = -10 * log(L_REF) / log(10);
	double MATCH_SCORE = -10 * log(0.25 / (1 - P_E)) / log(10);
	double MISMATCH_PUNALTY = -10 * log(0.75 / (P_E)) / log(10);
	// step1: get Q_MEM
	for (int i = 0; i < Q_MEM_MAX; i++)
		Q_MEM[i] = REF_SIZE_PUNALTY + i * MATCH_SCORE + 0.5;

	// step2: get Q_LV
	for (int j = 0; j < MAX_LV_R_LEN; j++)
	{
		for (int i = 0; i < MAX_LV_WRONG; i++)
		{
			Q_LV[i][j] = (j - i) * MATCH_SCORE + i * MISMATCH_PUNALTY + 0.5;
			if (j < 5)
				Q_LV[i][j] += 15; // gain score when ref_l is too small
			Q_LV[i][j] = MAX(Q_LV[i][j], -8);
		}
	}
}

void report_stats(struct timeval start)
{
	double seconds = realduration(start);
	fprintf(stderr, "%ld sequences processed in %.3fs (%.1f Kseq/m).\n",
			total_sequences,
			seconds,
			total_sequences / 1.0e3 / (seconds / 60));
}

static void classify_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:   %s\n", PACKAGE_NAME);
	fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact:   %s\n\n", CONTACT);
	fprintf(stderr, "  Usage:     %s  classify  [Options] <IndexDir> [ReadFiles.fa][...]>\n", PACKAGE_NAME);
	fprintf(stderr, "  Basic:   \n");
	fprintf(stderr, "    <IndexDir>      FOLDER   the directory contains %s index\n", PACKAGE_NAME);
	fprintf(stderr, "    [ReadFiles.fa]  FILES    reads files, FASTQ(A) format, separated by space\n");
	fprintf(stderr, "  Options:\n");
	fprintf(stderr, "    -h,             help\n");
	// fprintf(stderr, "    -u,             run in \"Strain mode\": this mode will increase strain\n");
	// fprintf(stderr, "                    level classification sensitivity and accuracy, but the speed\n");
	// fprintf(stderr, "                    may be slower.\n");
	fprintf(stderr, "    -t, INT         number of threads[4]\n");
	fprintf(stderr, "    -l, INT         minimum matching length, ignored for NGS reads [170]\n");
	fprintf(stderr, "    -r, INT         max Output number of secondary alignments[5]\n");
	fprintf(stderr, "    -o, FILE        output results into file [stdout]\n");
	fprintf(stderr, "    -s, INT         MIN score[64]\n");
	fprintf(stderr, "    -f, STR         output format, one of:\n");
	fprintf(stderr, "                    - SAM: SAM-like results without SEQ and QUAL and header, default\n");
	fprintf(stderr, "                    - SAM_FULL: SAM-like results with SEQ and QUAL\n");
	fprintf(stderr, "                    - DES: smallest format\n");
	fprintf(stderr, "                    - DES_FULL: all results are showed, ignore '-r' opinion\n");
	fprintf(stderr, "\n");
	// fprintf(stderr, "Options:   -e,           Error rate[0.85]\n");
	// fprintf(stderr, "Options:   -a,           Output all results\n");
}

#define PIPELINE_T_NUM 3 // one for reading; one for classifying; one for writing
#define STEP_NUM PIPELINE_T_NUM
kvec_T(kstring_t, kstring_V)

	int classify_main(int argc, char *argv[])
{					   // only use one thread, random means random seed
	double P_E = 0.15; // ERROR rate
	int c = 0;
	MAP_opt o = {170, 4, 5, OUTPUT_MODE_SAM, false, stdout, 64}; //, false};
	while ((c = getopt(argc, argv, "ht:l:r:f:o:s:")) >= 0)
	{
		if (c == 'h')
		{
			classify_usage();
			return 0;
		}
		else if (c == 't')
			o.thread_num = (int)atoi(optarg);
		// else if(c == 'u') o.strain_mode = true;
		else if (c == 'l')
			o.L_min_matching = (int)atoi(optarg);
		else if (c == 'r')
			o.max_sec_N = (int)atoi(optarg);
		else if (c == 'o')
			o.outfile = xopen(optarg, "w");
		else if (c == 's')
			o.min_score = (int)atoi(optarg);
		else if (c == 'f')
		{
			if (strcmp(optarg, "SAM") == 0)
				o.out_format = OUTPUT_MODE_SAM;
			else if (strcmp(optarg, "SAM_FULL") == 0)
				o.out_format = OUTPUT_MODE_SAM_FULL;
			else if (strcmp(optarg, "DES") == 0)
				o.out_format = OUTPUT_MODE_DES;
			else if (strcmp(optarg, "DES_FULL") == 0)
				o.out_format = OUTPUT_MODE_DES_FULL;
		}
	}
	if (optind + 2 > argc)
	{
		classify_usage();
		return 0;
	}
	char *index_dir = argv[optind++];
	kstring_V reads_files = {0};
	while (optind < argc)
	{
		kstring_t *reads_file;
		kv_pushp(kstring_t, reads_files, &reads_file);
		kstring_initp(reads_file);
		kputs(argv[optind++], reads_file);
	}
	fprintf(stderr, "loading index\t"); // load index
	DA_IDX idx = {0};
	// 第一个函数；索引载入
	load_idx(&idx, index_dir);
	idx.filter_min_length = o.L_min_matching;
	idx.filter_min_score = o.min_score;
	idx.filter_min_score_LV3 = o.min_score + 10;
	// idx.strain_mode = o.strain_mode;
	idx.mapQ.Q_MEM = xmalloc_t(int, Q_MEM_MAX);
	idx.mapQ.Q_LV = (int(*)[MAX_LV_R_LEN])xmalloc(MAX_LV_WRONG * MAX_LV_R_LEN * sizeof(int));
	calculate_MAPQ_TABLE(idx.mapQ.Q_MEM, idx.mapQ.Q_LV, P_E, idx.ref_bin.n * 4);

	struct timeval start;
	gettimeofday(&start, NULL);
	fprintf(stderr, "Start classify\n");
	double cpu_time = cputime();

	CLASSIFY_THREAD_DATA data[PIPELINE_T_NUM];
	o.thread_num = MAX(o.thread_num, 1);
	Classify_buff_pool *buff = xcalloc_t(Classify_buff_pool, o.thread_num);
	// init: buff
	for (int i = 0; i < o.thread_num; i++)
	{
		buff[i].sa_hash[0] = xmalloc(sizeof(sparse_align_HASH) * 0x100000); // 1M for each
		buff[i].sa_hash[1] = xmalloc(sizeof(sparse_align_HASH) * 0x100000); // 1M for each
	}
	CLASSIFY_SHARE_DATA share = {&idx, NULL, &o, buff, data};

	for (int i = 0; i < PIPELINE_T_NUM; i++)
	{
		data[i].seqs = xcalloc_t(kseq_t, N_NEEDED);
		data[i].results = xcalloc_t(cly_r, N_NEEDED);
		data[i].share_data_pointer = &share;
	}
	for (uint32_t i = 0; i < reads_files.n; ++i)
	{
		gzFile fp = xzopen(reads_files.a[i].s, "r");
		share._fp = ks_init(fp);
		fprintf(stderr, "Processing file: [%s].\n", reads_files.a[i].s);
		// 第二个函数： read载入以及序列比对
		kt_pipeline(PIPELINE_T_NUM, classify_pipeline, &share, STEP_NUM);
		gzclose(fp);
	}
	report_stats(start);
	fprintf(stderr, "Classify CPU: %.3f sec\n", cputime() - cpu_time);

	// 第三个函数： 数据分析：
	return 0;
}

void print_sensei()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "                                 D G       . D         ; :                       \n");
	fprintf(stderr, "                             E       D K E E E   E E .       E                   \n");
	fprintf(stderr, "                               G E E E E . ; E E   E           D                 \n");
	fprintf(stderr, "                       E E E E E E E E f     E . E E E E E E E                   \n");
	fprintf(stderr, "                       E E E E K .     E K     K       K       K                 \n");
	fprintf(stderr, "                       E E                         ; E E E E D                   \n");
	fprintf(stderr, "                   E E E                                 E E E                   \n");
	fprintf(stderr, "                     E                                         E                 \n");
	fprintf(stderr, "                   E                                             E               \n");
	fprintf(stderr, "               E D                                               K E             \n");
	fprintf(stderr, "             E   E                                                 E             \n");
	fprintf(stderr, "             E E                                                   E             \n");
	fprintf(stderr, "           E                                                       E             \n");
	fprintf(stderr, "           E E                                                     E E           \n");
	fprintf(stderr, "         D   K                                                   K E E           \n");
	fprintf(stderr, "         E D t       E D E L                                       E .           \n");
	fprintf(stderr, "         E E       . E       E                                     E E E     E   \n");
	fprintf(stderr, "     E   E E       E         E                                     E E D       D \n");
	fprintf(stderr, "   E     E E       E                   E                           i E E E     E \n");
	fprintf(stderr, "   j     E E       E                   E E             E   K .       E E E     E \n");
	fprintf(stderr, " E       E K       D                   E E           E       E       E E E     D \n");
	fprintf(stderr, " E     E E E                                                 E       K E E   E   \n");
	fprintf(stderr, "   E   E . E                                                         f E E E     \n");
	fprintf(stderr, "     E E E E                                                         ; K E       \n");
	fprintf(stderr, "           E                               E                         E E E       \n");
	fprintf(stderr, "           j                   .           E                         E E E       \n");
	fprintf(stderr, "             K                 E           E                         .   E       \n");
	fprintf(stderr, "             E                 i         E                         E     E       \n");
	fprintf(stderr, "                                 E     . D                       E               \n");
	fprintf(stderr, "               E                   E E                         . D               \n");
	fprintf(stderr, "                 E                                           . D                 \n");
	fprintf(stderr, "                 E                                           E                   \n");
	fprintf(stderr, "                   E                                       E                     \n");
	fprintf(stderr, "                     E                                   E                       \n");
	fprintf(stderr, "                       E                               E                         \n");
	fprintf(stderr, "                           K                       E                             \n");
	fprintf(stderr, "                               E E             t E                               \n");
	fprintf(stderr, "                                     . K E E E                                   \n");
	fprintf(stderr, "\n");
}


static int cmp_count_sort(const void *a_, const void *b_){
	COUNT_SORT *a = (COUNT_SORT*) a_, *b= (COUNT_SORT *) b_;
	return a->count < b->count;
}

static int max_tid_global = 0;
static int taxonTree_rank(const char *taxonomyPath,TAXONOMY_rank ** taxonomyTree_)
{
	TAXONOMY_rank * taxonomyTree;
	char buf[2048];
	sprintf(buf, "%s/nodes.dmp", taxonomyPath);
	FILE *fp = xopen(buf,"r");
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
	sprintf(buf, "%s/nodes.dmp", taxonomyPath);
	fp = xopen(buf,"r");
	//step 3: store
	max_tid += 1000000;
	max_tid_global = max_tid;
	taxonomyTree = (TAXONOMY_rank*)malloc(sizeof(TAXONOMY_rank)*(max_tid + 1));
	for(uint32_t i = 0; i <= max_tid; i++) {
		taxonomyTree[i].p_tid = MAX_uint32_t;
		taxonomyTree[i].name[0] = '\0';
	}
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
	// set vitural root [MAX_uint32_t] for normal tax tree
	taxonomyTree[1].p_tid = MAX_uint32_t;//this tree's root is [1]
	// set vitural root [MAX_uint32_t] for another tax tree, which only has a root [0] for unmapped reads
	taxonomyTree[0].p_tid = MAX_uint32_t;
	strcpy(taxonomyTree[0].rank,"no rank");
	strcpy(taxonomyTree[0].name,"CLY_FAIL");

	//step4 : close file
	fclose(fp);

	//step5 : read names
	sprintf(buf, "%s/names.dmp", taxonomyPath);
	fp = xopen(buf,"r");
	while(1){
		if(getline(&line,&max_l,fp) <= 0)
			break;
		//get tid
		token = strtok(line,"|\t");
		uint32_t tid = strtoul(token,NULL,10);
		//get name
		token = strtok(NULL,"\t|");
		char * name = token;
		//skip column
		token = strtok(NULL,"|");
		//get name type
		token = strtok(NULL,"|");
		if (strncmp("\tscien", token, 6) == 0) {
			strncpy(taxonomyTree[tid].name, name, MAX_taxon_name_N);
			// fprintf(stderr, "[%d][%s][%s]\n", tid, token, name);
		}
	}
	fclose(fp);

	// return
	free(line);
	*taxonomyTree_ = taxonomyTree;
	return max_tid;
}

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
	char* seq_tokens;
	if(tokens[0] == '*') {
		rst->isClassify = 'U';
		rst->tid = 0;
		rst->MAPQ = 0;
		//ignore the POS part
		tokens = strtok(NULL,"\t");
		//ignore MAQ part
		tokens = strtok(NULL,"\t");
		//ignore CIGAR
		tokens = strtok(NULL,"\t");
		//ignore *
		tokens = strtok(NULL,"\t");
		//ignore 0
		tokens = strtok(NULL,"\t");
		//ignore 0
		tokens = strtok(NULL,"\t");
		//get SEQ
		tokens = strtok(NULL,"\t");
		seq_tokens = tokens;
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
		// char * CIGAR = tokens;
		//get *
		tokens = strtok(NULL,"\t");
		//get 0
		tokens = strtok(NULL,"\t");
		//get 0
		tokens = strtok(NULL,"\t");
		//get SEQ
		tokens = strtok(NULL,"\t");
		seq_tokens = tokens;
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
		// //for the read length part
		// {
		// 	int read_len = 0;
		// 	int type_len = 0;
		// 	while(1)
		// 	{
		// 		char c_char = *CIGAR++;
		// 		if(c_char == 0)
		// 			break;
		// 		if(c_char <= '9' && c_char >= '0')
		// 			type_len = (type_len * 10) + (c_char - '0');
		// 		else
		// 		{
		// 			if(c_char == 'M' || c_char == 'I' || c_char == 'S' || c_char == 'X')
		// 				read_len += type_len;
		// 			type_len = 0;
		// 		}
		// 	}
		// 	rst->read_length = read_len;
		// }
	}
	rst->read_length = strlen(seq_tokens);
	rst->seq = xmalloc(rst->read_length + 1);
	strcpy(rst->seq, seq_tokens);
	rst->seq[rst->read_length] = 0;
	return 0;
}

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
		fprintf(file, "%s\t%s|%s\tnull\t%f\n", species_type, taxonomyTree[leaf_ID].name, taxonomyTree[leaf_ID].rank, rate);
		if (DEBUG) {
			for (int i = 0; i < level; i ++) fprintf(stderr, "  ");
			fprintf(stderr, "DEBUG: %s\t%s|%s\tnull\t%f\n", species_type, taxonomyTree[leaf_ID].name, taxonomyTree[leaf_ID].rank, rate);
		} 
	}
}


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
			return tid;
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

void init_buff_core(void *buff_, void *idx, int thread_num)
{

	RM_buffer *buff = (RM_buffer *)buff_;
	{ // read_classify
		buff->read_classify_input_tmpfile = tmpfile();
		buff->read_classify_output_tmpfile = tmpfile();
		buff->map_opt.L_min_matching = 170;
		buff->map_opt.thread_num = thread_num;
		buff->map_opt.max_sec_N = 5;
		buff->map_opt.out_format = OUTPUT_MODE_SAM_FULL;
		buff->map_opt.show_anchor = false;
		buff->map_opt.outfile = buff->read_classify_output_tmpfile;
		buff->map_opt.min_score = 64;
		buff->classify_thread_data = xcalloc_t(CLASSIFY_THREAD_DATA, PIPELINE_T_NUM);
		for (int i = 0; i < PIPELINE_T_NUM; i++)
		{
			buff->classify_thread_data[i].seqs = xcalloc_t(kseq_t, N_NEEDED);
			buff->classify_thread_data[i].results = xcalloc_t(cly_r, N_NEEDED);
			buff->classify_thread_data[i].share_data_pointer = &(buff->classify_share_data);
		}
		buff->classify_buff_pool = xcalloc_t(Classify_buff_pool, buff->map_opt.thread_num);
		for (int i = 0; i < buff->map_opt.thread_num; i++)
		{
			buff->classify_buff_pool[i].sa_hash[0] = xmalloc(sizeof(sparse_align_HASH) * 0x100000); // 1M for each
			buff->classify_buff_pool[i].sa_hash[1] = xmalloc(sizeof(sparse_align_HASH) * 0x100000); // 1M for each
		}
		buff->classify_share_data.idx = idx;
		buff->classify_share_data._fp = NULL;
		buff->classify_share_data.o = &(buff->map_opt);
		buff->classify_share_data.buff = buff->classify_buff_pool;
		buff->classify_share_data.data = buff->classify_thread_data;
	}
	{ // meta_analysis
		buff->meta_analysis_input_tmpfile = tmpfile();
		buff->meta_analysis_dump_tmpfile = tmpfile();
		buff->meta_analysis_output_tmpfile = tmpfile();
		buff->getOneSAM_buff = (char *)malloc(MAX_BUFF_LEN);
		buff->node_count = xcalloc(((DA_IDX *)idx)->max_tid, sizeof(uint64_t));
		buff->node_table = xcalloc(((DA_IDX *)idx)->max_tid, sizeof(CLY_NODE));
		buff->child_list = xcalloc(((DA_IDX *)idx)->max_tid * 2, sizeof(CN_CHILD));
		buff->sort = xmalloc(((DA_IDX *)idx)->max_tid * sizeof(COUNT_SORT));
	}
}

// 释放init_buff函数中分配的缓存；
void free_buff_core(void *buff_)
{
	RM_buffer *buff = (RM_buffer *)buff_;
	{ // read_classify
		fclose(buff->read_classify_input_tmpfile);
		fclose(buff->read_classify_output_tmpfile);
		for (int i = 0; i < PIPELINE_T_NUM; i++)
		{
			free(buff->classify_thread_data[i].seqs);
			free(buff->classify_thread_data[i].results);
		}
		free(buff->classify_thread_data);
		for (int i = 0; i < buff->map_opt.thread_num; i++)
		{
			free(buff->classify_buff_pool[i].sa_hash[0]);
			free(buff->classify_buff_pool[i].sa_hash[1]);
		}
		free(buff->classify_buff_pool);
	}
	{ // meta_analysis
		fclose(buff->meta_analysis_input_tmpfile);
		fclose(buff->meta_analysis_dump_tmpfile);
		fclose(buff->meta_analysis_output_tmpfile);
		free(buff->getOneSAM_buff);
		free(buff->node_count);
		free(buff->node_table);
		free(buff->child_list);
		free(buff->sort);
	}
	free(buff);
}

void read_classify_core(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, void *buff)
{
	RM_buffer *buff_ = buff;

	// 把输入的fastq字符串或是fastq文件地址统统转换成 FILE* 变量，并且套成gzFile
	FILE *fp = NULL;
	if (input_n == -1)
	{
		fp = xopen(input, "r");
	}
	else
	{
		fp = buff_->read_classify_input_tmpfile;
		rewind(fp);
		fprintf(fp, "%s\n", input);
		rewind(fp);
	}
	gzFile gzfp = gzdopen(fileno(fp), "r");

	// 调用deSAMBA执行分类运算
	buff_->classify_share_data._fp = ks_init(gzfp);
	kt_pipeline(PIPELINE_T_NUM, classify_pipeline, &(buff_->classify_share_data), STEP_NUM);

	// 把结果输出到文件
	*output_n = ftell(buff_->read_classify_output_tmpfile);
	*output = xmalloc(*output_n + 1);
	memset(*output, 0, *output_n + 1);
	rewind(buff_->read_classify_output_tmpfile);
	xread(*output, 1, *output_n, buff_->read_classify_output_tmpfile);

	// 释放必要的内存和文件
	ks_destroy(buff_->classify_share_data._fp);
	gzclose(gzfp);
	if (input_n == -1)
	{
		fclose(fp);
	}
}

void meta_analysis_core(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, void *buff, int flag)
{
	RM_buffer *buff_ = buff;

	// 重置缓冲区中的一些数据
	memset(buff_->node_count, 0, ((DA_IDX *)idx)->max_tid * sizeof(uint64_t));
	memset(buff_->node_table, 0, ((DA_IDX *)idx)->max_tid * sizeof(CLY_NODE));
	memset(buff_->child_list, 0, ((DA_IDX *)idx)->max_tid * sizeof(CN_CHILD) * 2);
	memset(buff_->sort, 0, ((DA_IDX *)idx)->max_tid * sizeof(COUNT_SORT));

	// 把输入的格式转换为deSAMBA需要的格式
	rewind(buff_->meta_analysis_input_tmpfile);
	fwrite(input, 1, input_n, buff_->meta_analysis_input_tmpfile);
	buff_->meta_analysis_input_tmpfile = freopen(NULL, "r", buff_->meta_analysis_input_tmpfile); // for getline
	skip_sam_head(buff_->meta_analysis_input_tmpfile, buff_->getOneSAM_buff);
	uint64_t record_num = 0;
	RST temp_rst;

	if (DEBUG)
		fprintf(stderr, "begin read sam\n");
	rewind(buff_->meta_analysis_dump_tmpfile);

	rewind(buff_->meta_analysis_output_tmpfile); // 收集人类序列

	while (1)
	{
		int getOneSAM_rst = getOneSAM(buff_->meta_analysis_input_tmpfile, buff_->getOneSAM_buff, &temp_rst);
		if (getOneSAM_rst < 0)
			break;
		record_num++;
		if (record_num >= MAX_read_N)
			break;
		if (record_num < MIN_read_N)
			continue;
		if (temp_rst.tid == 9606 || temp_rst.tid == 63221 || temp_rst.tid == 741158) {
			fprintf(buff_->meta_analysis_output_tmpfile, "%s", temp_rst.seq);
		}
		free(temp_rst.seq);
		fprintf(buff_->meta_analysis_dump_tmpfile, "%s\t%c\t%d\t%d\t%d\t%d\n",
				temp_rst.read_name,
				temp_rst.isClassify,
				temp_rst.tid,
				temp_rst.read_length,
				temp_rst.MAPQ,
				temp_rst.score);
		// if (DEBUG) fprintf(stderr, "getOneSAM:\t%s\t%c\t%d\t%d\t%d\t%d\n", temp_rst.read_name, temp_rst.isClassify, temp_rst.tid, temp_rst.read_length, temp_rst.MAPQ, temp_rst.score);
	}
	if (DEBUG) fprintf(stderr, "end read [%ld] sam, find [%ld] human bases\n", record_num, ftell(buff_->meta_analysis_output_tmpfile));
	fprintf(buff_->meta_analysis_output_tmpfile, "\n");
	rewind(buff_->meta_analysis_dump_tmpfile);

	// 调用deSAMBA进行数据分析
	long int total_read_number = 0;
	long int total_weight = 0;
	RST rst;
	int eof_ = 0;
	float coverage = 0;

	if (DEBUG)
		fprintf(stderr, "begin read RST and map unique tid to total count in node_count\n");
	int getOneRST_result = getOneRST(buff_->meta_analysis_dump_tmpfile, &rst);
	if (DEBUG)
		fprintf(stderr, "*** first getOneRST_result: %d (0 means continue, -1 means EOF)\n", getOneRST_result);
	if (getOneRST_result < 0)
	{
		return;
	}
	while (1)
	{
		// if (DEBUG) fprintf(stderr, "RST:\t%s\t%c\t%d\t%d\t%d\t%d\t", rst.read_name, rst.isClassify, rst.tid, rst.read_length, rst.MAPQ, rst.score);
		uint32_t current_read_length = rst.read_length;
		uint32_t current_weight = ((flag & 0x1) == 0) ? 1 : current_read_length; // 选择按read条数统计碱基数量统计
		total_read_number++;
		total_weight += current_weight;
		int read_len = 0;
		uint32_t final_tid = ana_get_tid(&rst, ((DA_IDX *)idx)->max_tid, buff_->meta_analysis_dump_tmpfile, &eof_, ((DA_IDX *)idx)->taxonomyTree, &read_len, &coverage);
		buff_->node_count[final_tid] += current_weight; // 因为此时rst已经是下一条read的了
		// if (DEBUG) fprintf(stderr, "node_count[%d]=%ld, total_weight=%ld\n", final_tid, buff_->node_count[final_tid], total_weight);
		if (eof_ < 0)
			break;
	}
	if (DEBUG)
		fprintf(stderr, "end read [%ld] RST\n", total_read_number);

	uint32_t child_count = 1;

	// sort node_count
	if (DEBUG)
		fprintf(stderr, "begin sort node_count by count\n");
	int rst_num = 0;
	for (uint32_t i = 0; i <= ((DA_IDX *)idx)->max_tid; i++)
	{
		if (buff_->node_count[i] == 0)
			continue;
		buff_->sort[rst_num].tid = i;
		buff_->sort[rst_num++].count = buff_->node_count[i];
		if (DEBUG) fprintf(stderr, "taxid: %d\tcount: %ld\n", i, buff_->node_count[i]);
	}
	qsort(buff_->sort, rst_num, sizeof(COUNT_SORT), cmp_count_sort);
	if (DEBUG)
		fprintf(stderr, "end sort node_count by count, total [%d] kinds of taxid\n", rst_num);

	// calculate middle nodes
	if (DEBUG)
		fprintf(stderr, "begin calculate middle nodes\n");
	for (int i = 0; i < rst_num; i++)
	{
		uint32_t c_tid = buff_->sort[i].tid;
		while (1)
		{
			uint32_t p_tid = ((DA_IDX *)idx)->taxonomyTree[c_tid].p_tid;
			// if (DEBUG) fprintf(stderr, "c_tid:%u, p_tid:%u\n", c_tid, p_tid);
			buff_->node_table[c_tid].weight += buff_->node_count[buff_->sort[i].tid]; // add weight for p_tid
			// if (DEBUG) fprintf(stderr, "update node_table[%d].weight to %lu\n", c_tid, buff_->node_table[c_tid].weight);
			if (p_tid == MAX_uint32_t) {
				break;
			}
			if (buff_->node_table[p_tid].child_list_begin == 0)
			{
				buff_->node_table[p_tid].child_list_begin = child_count++;
				buff_->child_list[child_count - 1].tid = c_tid;
			}
			else
			{
				int list_begin = buff_->node_table[p_tid].child_list_begin;
				while (buff_->child_list[list_begin].tid != c_tid && buff_->child_list[list_begin].next != 0)
					list_begin = buff_->child_list[list_begin].next;
				if (buff_->child_list[list_begin].tid != c_tid && buff_->child_list[list_begin].next == 0) // store new node
				{
					buff_->child_list[list_begin].next = child_count++;
					buff_->child_list[child_count - 1].tid = c_tid;
				}
			}
			c_tid = ((DA_IDX *)idx)->taxonomyTree[c_tid].p_tid;
		}
	}
	if (DEBUG)
		fprintf(stderr, "end calculate middle nodes\n");


	// 格式化运算结果
	ana_meta_loop_fprint(buff_->meta_analysis_output_tmpfile, ((DA_IDX *)idx)->taxonomyTree, buff_->node_table, 0, buff_->child_list, 0, total_weight, false);
	ana_meta_loop_fprint(buff_->meta_analysis_output_tmpfile, ((DA_IDX *)idx)->taxonomyTree, buff_->node_table, 1, buff_->child_list, 0, total_weight, false);
	*output_n = ftell(buff_->meta_analysis_output_tmpfile);
	*output = xmalloc(*output_n + 1);
	memset(*output, 0, *output_n + 1);
	rewind(buff_->meta_analysis_output_tmpfile);
	xread(*output, 1, *output_n, buff_->meta_analysis_output_tmpfile);
}

// 索引载入函数。其中，void **idx 是一个 void * 类型的指针的指针。 dirPath是字符串，存储索引的路径；
void load_index(void **idx, const char *dirPath)
{
	if (DEBUG)
		print_sensei();
	// decode DA_IDX type
	if (DEBUG)
		fprintf(stderr, "BEGIN loading index INIT\n");
	DA_IDX *index = xcalloc(1, sizeof(DA_IDX));
	memset(index, 0, sizeof(DA_IDX));
	if (DEBUG)
		fprintf(stderr, "END loading index INIT\n");
	// default parameter
	double P_E = 0.15;
	// load index
	if (DEBUG)
		fprintf(stderr, "BEGIN loading index core\n");
	load_idx(index, dirPath);
	if (DEBUG)
		fprintf(stderr, "END loading index core\n");
	(*index).filter_min_length = 170;
	(*index).filter_min_score = 64;
	(*index).filter_min_score_LV3 = 64 + 10;
	(*index).mapQ.Q_MEM = xmalloc_t(int, Q_MEM_MAX);
	(*index).mapQ.Q_LV = (int(*)[MAX_LV_R_LEN])xmalloc(MAX_LV_WRONG * MAX_LV_R_LEN * sizeof(int));
	calculate_MAPQ_TABLE((*index).mapQ.Q_MEM, (*index).mapQ.Q_LV, P_E, (*index).ref_bin.n * 4);
	// load taxonomy
	if (DEBUG)
		fprintf(stderr, "BEGIN loading taxonTree_rank\n");
	(*index).max_tid = taxonTree_rank(dirPath, &((*index).taxonomyTree));
	if (DEBUG)
		fprintf(stderr, "END loading taxonTree_rank\n");
	// init thread2buff
	(*index).thread2buff = kh_init(i2p);
	pthread_mutex_init(&((*index).thread2buff_mutex), NULL);
	// return
	*idx = index;
}

/**
 * @param thread_num how many thread to use in a single read_classify. 
*/
void* find_and_init_buff_for_thread_mutex(int thread_id, void *idx, int thread_num) {
	pthread_mutex_lock(&(((DA_IDX *)idx)->thread2buff_mutex)); ////////////////// 上锁
	// 查找thread_id是否在idx的thread2buff中
	int ret, is_missing;
	khiter_t k;
	k = kh_get(i2p, ((DA_IDX *)idx)->thread2buff, thread_id);
	is_missing = (k == kh_end(((DA_IDX *)idx)->thread2buff));
	// 如果thread_id存在，但是线程数与需要的数量不相等，则将其删除，后面重新分配缓冲区
	if (!is_missing && thread_num != -1) {
		int current_thread_num = ((RM_buffer*)kh_value(((DA_IDX *)idx)->thread2buff, k))->map_opt.thread_num;
		if (current_thread_num != thread_num) {
			if (DEBUG) fprintf(stderr, "delete buff of thread %d with %d sub-threads because of sub-thread-num change\n", thread_id, current_thread_num);
			free_buff_core(kh_value(((DA_IDX *)idx)->thread2buff, k));
			kh_del(i2p, ((DA_IDX *)idx)->thread2buff, k);
			is_missing = 1;
		}
	}
	// 如果thread_id不存在，则将其插入字典中
	if (is_missing) {
		if (DEBUG) fprintf(stderr, "create buff for thread %d with %d sub-threads, creating buff for it\n", thread_id, thread_num);
		k = kh_put(i2p, ((DA_IDX *)idx)->thread2buff, thread_id, &ret);
		// 如果缓冲区不存在，则创建缓冲区
		kh_value(((DA_IDX *)idx)->thread2buff, k) = xcalloc(1, sizeof(RM_buffer));
		init_buff_core(kh_value(((DA_IDX *)idx)->thread2buff, k), idx, thread_num);
		if (DEBUG) fprintf(stderr, "done create buff for thread %d, buff is at address %p\n", thread_id, kh_value(((DA_IDX *)idx)->thread2buff, k));
	}
	pthread_mutex_unlock(&(((DA_IDX *)idx)->thread2buff_mutex)); ////////////////// 解锁
	return kh_value(((DA_IDX *)idx)->thread2buff, k);
}

void read_classify(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int thread_num){
	if (input_n == 0) {
		*output_n = 0;
		return;
	}
	void *buff = find_and_init_buff_for_thread_mutex(thread_id, idx, thread_num); // TODO: add to desamba.h interface
	read_classify_core(idx, input, input_n, output, output_n, buff);
}

typedef struct {
	char type[256];
	char species[256];
	char tech[256];
	float rate;
} MetaRST;
int cmp_MetaRST(const void *a,const void *b) {
	if (((MetaRST*)a)->rate > ((MetaRST*)b)->rate) return -1;
	else if (((MetaRST*)a)->rate < ((MetaRST*)b)->rate) return 1;
	else return 0;
}
void meta_analysis(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int flag, uint64_t max_snap_shot_len, char **human_snapshot, uint64_t *human_snapshot_n)
{
	*human_snapshot = xmalloc(max_snap_shot_len + 1);
	memset(*human_snapshot, 0, max_snap_shot_len + 1);

	if (input_n == 0) {
		*output_n = 0;
		*human_snapshot_n = 0;
		return;
	}

	void *buff = find_and_init_buff_for_thread_mutex(thread_id, idx, -1);
	meta_analysis_core(idx, input, input_n, output, output_n, buff, flag);

	char *cursor = *output;
	if (cursor[0] == '\n') { // 第一行：快照.
		strncpy(*human_snapshot, "", max_snap_shot_len);
		cursor += 1;
	}
	else {
		char *rst_line = strtok(cursor, "\n");
		// if (DEBUG) fprintf(stderr, "parse human_base: %s\n", rst_line);
		strncpy(*human_snapshot, rst_line, max_snap_shot_len);
		// if (DEBUG) fprintf(stderr, "      human_base: %s\n", *human_snapshot);
		*human_snapshot_n = strlen(*human_snapshot);
		cursor += strlen(rst_line) + 1;
		// cursor += strlen(*human_snapshot) + 1;
	}

	if (DEBUG) fprintf(stderr, "total human base num: %ld\n", cursor - *output - 1);
	if (DEBUG) fprintf(stderr, "original meta_analysis file (item part):\n%s\n", cursor);

	kvec_t(MetaRST) results;
	kv_init(results);
	float no_match_rate = 0;
	while(1) { // 其他行：转换格式
		char *rst_line = strtok(cursor, "\n");
		if (rst_line == NULL) break;
		if (DEBUG) fprintf(stderr, "parse line: %s\n", rst_line);
		cursor += strlen(rst_line) + 1;

		MetaRST rst;
		sscanf(rst_line, "%[^\t]\t%[^\t]\t%[^\t]\t%f", rst.type, rst.species, rst.tech, &(rst.rate));
		// if (DEBUG) fprintf(stderr, "      line: %s\t%s\t%s\t%f\n", rst.type, rst.species, rst.tech, rst.rate);
		if (strcmp("no_match", rst.type) == 0) {
			no_match_rate += rst.rate;
		}
		else {
			kv_push(MetaRST, results, rst);
		}
	}
	if (DEBUG) fprintf(stderr, "After sort:\n");
	if (no_match_rate > 0.95) {
		sprintf(*output, "no_match\tnull|null\tnull\t0\n");
		if (DEBUG) fprintf(stderr, "%s", "no_match\tnull|null\tnull\t0\n");
		*output_n = strlen(*output);
		return;
	}
	else {
		// 归一化
		for (long i = 0; i < kv_size(results); ++i) {
			kv_A(results, i).rate = kv_A(results, i).rate / (1 - no_match_rate);
		}
		// 排序
		qsort(results.a, results.n, sizeof(MetaRST), cmp_MetaRST);
		*output_n = 0;
		for (long i = 0; i < kv_size(results); ++i) {
			MetaRST rst = kv_A(results, i);
			if ((i < 3) || (strcmp("human", rst.type) == 0 && rst.rate > 0.05)) {
				sprintf(*output + *output_n, "%s\t%s\t%s\t%f\n", rst.type, rst.species, rst.tech, rst.rate);
				if (DEBUG) fprintf(stderr, "%s\t%s\t%s\t%f\n", rst.type, rst.species, rst.tech, rst.rate);
				*output_n = strlen(*output);
			}
		}
	}

	kv_destroy(results);
}