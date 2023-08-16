#ifndef DESAMBA_H
#define DESAMBA_H
#include <stdint.h>

/**
 * @brief Load index from directory
 * @param idx Pointer to a wild pointer, memory will be allocated inside function
 * @param dirPath Path to directory containing index files, including deSAMBA.*, nodes.dmp and names.dmp
*/
void load_index(void **idx, const char *dirPath);

/**
 * @brief classify reads
 * @param idx Pointer to a index loaded by load_index
 * @param input Input fastq format string | Input file path
 * @param input_n length of input string | -1 (indicate file-path mode)
 * @param output Pointer to a wild pointer, memory will be allocated inside function
 * @param output_n length of output string
 * @param thread_id a unique thread_id indicating which thread calls this function
*/
void read_classify(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id);

/**
 * @brief convert read_classify output to a more readable format
 * 输出格式：种\t亚种\t测序方式（暂时全部为null）\t准确率
 * 可能为空串，也可能有多个结果，以\n分隔
 * @param idx Pointer to a index loaded by load_index
 * @param input Read_classify output
 * @param input_n Read_classify output length
 * @param output Pointer to a wild pointer, memory will be allocated inside function
 * @param output_n length of output string
 * @param thread_id a unique thread_id indicating which thread calls this function
*/
void meta_analysis(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id);

#endif