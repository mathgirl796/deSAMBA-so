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
 * @param output Pointer to a wild pointer, memory will be allocated inside function,需要调用者释放内存
 * @param output_n length of output string
 * @param thread_id a unique thread_id indicating which thread calls this function
 * @param thread_num 本函数内部使用的线程个数，可以为每个thread_id设置一个thread_num。对于同一个thread_id，在两次调用之间修改thread_num会导致缓冲区重新分配的时间开销
 * 举个例子：你在网络框架中设置了10个成员的线程池，但你一共有100个cpu core，那么你可以为每个线程分配10个cpu core来调用read_classify，实现方法为设置thread_num为10
*/
void read_classify(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int thread_num);

/**
 * @brief 对read_classify的结果以分类树上的节点（通常是叶节点）为单位进行统计
 * 输出格式：[1]\t[2]\t[3]\t[4]\n
 * 当成功分类数据少于5%时，输出"no_match   null|null       null    0"
 * 否则输出去除分类失败数据正则化后比例排名前三的类别，输出的[4]为正则化后该类别所占比例
 * 若分类得到的人类数据大于5%，且不排在前三的话，额外输出一条人类条目，即一共输出四条
 * 
 * @param idx Pointer to a index loaded by load_index
 * @param input Read_classify output
 * @param input_n Read_classify output length
 * @param output Pointer to a wild pointer, memory will be allocated inside function，需要调用者释放内存
 * @param output_n length of output string
 * @param thread_id a unique thread_id indicating which thread calls this function
 * @param flag META_USE_READ_NUM表示使用各节点数据条数进行统计，META_USE_BASE_NUM表示使用各节点数据总碱基数进行统计
 * @param human_snapshot Pointer to a wild pointer, memory will be allocated inside function，人类基因快照，需要调用者释放内存
 * @param max_snap_shot_len 指定human_snapshot最大碱基数，程序会尽可能填充指定的碱基数量
*/
#define META_USE_READ_NUM 0
#define META_USE_BASE_NUM 1
void meta_analysis(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int flag, uint64_t max_snap_shot_len, char **human_snapshot, uint64_t *human_snapshot_n);

#endif