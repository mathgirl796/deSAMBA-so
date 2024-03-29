#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "desamba.h"

int main(int argc, char *argv[])
{

    void *idx = NULL;
    char *soPath = argv[1];
    char *dirPath = argv[2];
    char fqPath[1000];

    char *sam = "no sam";
    uint64_t sam_n;
    char *ana = "no ana";
    uint64_t ana_n;
    char *human_snapshot;
    uint64_t human_snapshot_n;

    // 加载库函数
    void *handle = dlopen(soPath, RTLD_LAZY);
    void (*load_index_func)(void **idx, const char *dirPath) = dlsym(handle, "load_index");
    void (*read_classify_func)(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int thread_num) = dlsym(handle, "read_classify");
    void (*meta_analysis_func)(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int flag, uint64_t max_s, char **human_snapshot, uint64_t *hss_n) = dlsym(handle, "meta_analysis");


    if (0) // 测试文件路径输入
    {
        load_index_func(&idx, dirPath);
        while (1)
        {
            puts("input fqPath:");
            gets(fqPath, 1024);
            if (access(fqPath, 0) != 0) {
                fprintf(stderr, "error: [%s] not exist!\n", fqPath);
                continue;
            }
            read_classify_func(idx, fqPath, -1, &sam, &sam_n, 0, 4); // 模拟1个线程，以四个子线程运行
            meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, 0, META_USE_BASE_NUM, 12, &human_snapshot, &human_snapshot_n);
            puts(ana);
            free(sam);
            free(ana);
        }
    }

    if (1) // 测试fastq字符串输入
    {
        load_index_func(&idx, dirPath);
        FILE* fin;
        int fsize;
        char* fq;
        int fread_result;
        while (1)
        {
            puts("input fqPath:");
            gets(fqPath, 1024);
            if (access(fqPath, 0) != 0) {
                fprintf(stderr, "error: [%s] not exist!\n", fqPath);
                continue;
            }
            fin = fopen(fqPath, "r");
            fseek(fin, 0, SEEK_END);
            fsize = ftell(fin);
            rewind(fin);
	        fprintf(stderr, "ftell: fsize = %d\n", fsize);
            fq = (char*)malloc(fsize+1);
	        fprintf(stderr, "malloc: fq at 0x%p\n", fq);
            fread_result = fread(fq, 1, fsize, fin);
	        fprintf(stderr, "fread: fread_result = %d\n", fread_result);
            read_classify_func(idx, fq, fsize+1, &sam, &sam_n, 0, 4); // 模拟1个线程，以四个子线程运行
            if (sam_n == 0) {
	            fprintf(stderr, "main: sam_n = 0, classify failed\n");
                free(fq); fclose(fin);
                continue;
            }
            meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, 0, META_USE_BASE_NUM, 12, &human_snapshot, &human_snapshot_n);
            if (ana_n == 0) {
	            fprintf(stderr, "main: ana_n = 0, analysis failed\n");
                free(fq); free(sam); fclose(fin);
                continue;
            }
            free(sam);
            free(ana);
            free(fq);
            fclose(fin);
        }
    }
}