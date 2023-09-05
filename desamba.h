#ifndef DESAMBA_H
#define DESAMBA_H
#include <stdint.h>

/**
 * @brief 把dirPath中的索引加载到*idx中。所需内存会在函数内进行分配
 * @param idx Pointer to a wild pointer, memory will be allocated inside function
 * @param dirPath Path to directory containing index files, including deSAMBA.*, nodes.dmp and names.dmp
*/
void load_index(void **idx, const char *dirPath);

/**
 * @brief 将来自文件或来自内存的格式正确的fastq数据进行分类，对每一条fastq数据，输出一条或多条sam格式字符串到内存中，输出字符串内存会在函数内部申请，需要调用者手动释放
 * @param idx Pointer to a index loaded by load_index
 * @param input Input fastq format string | Input file path
 * @param input_n length of input string | -1 (indicate file-path mode)
 * @param output Pointer to a wild pointer, memory will be allocated inside function
 * @param output_n length of output string
 * @param thread_id a unique thread_id indicating which thread calls this function
 * @param thread_num 本函数内部使用的线程个数，可以为每个thread_id设置一个thread_num。对于同一个thread_id，在两次调用之间修改thread_num会导致缓冲区重新分配的开销
 * 举个例子：你在网络框架中设置了10个成员的线程池，但你一共有100个cpu core，那么你可以为每个线程分配10个cpu core来调用read_classify，实现方法为设置thread_num为10
*/
void read_classify(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int thread_num);

/**
 * @brief 对read_classify的结果以分类树上的节点（通常是叶节点）为单位进行统计
 * 输出格式：[1]\t[2]\t[3]\t[4]
 * 可能为空串，也可能有多个结果，以\n分隔
 * 输出字符串内存会在函数内部申请，需要调用者手动释放
 * [1]表示该物种大致类别，取值范围为[no_match, human, animal_and_plant, microbe]。no_match表示在read_classify中比对到分类树根节点或者比对失败
 * [2]以a|b的形式表示，其中a为该节点的名字，例如"Homo sapiens","Zaire ebolavirus"；b为该节点的分类级别，例如"species","no rank"。
 * 特别的，a为"root"或"CLY_FAIL"时，b将是"no_rank"，而且这两种情况完全对应[1]为no_match的情况
 * [3]为保留字段
 * [4]为该节点上映射到的fastq数据的比例
 * 该函数保证返回结果若不为空，则返回的各行中[4]的总和一定为1(浮点数精度导致的误差不计)
 * @param idx Pointer to a index loaded by load_index
 * @param input Read_classify output
 * @param input_n Read_classify output length
 * @param output Pointer to a wild pointer, memory will be allocated inside function
 * @param output_n length of output string
 * @param thread_id a unique thread_id indicating which thread calls this function
 * @param flag 0表示使用各节点数据条数进行统计，1表示使用各节点数据总碱基数进行统计
*/
void meta_analysis(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int flag);

#endif