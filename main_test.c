#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <dlfcn.h>

int main(int argc, char *argv[]) {

    if (argc < 4) {
        fprintf(stderr, "usage: %s soPath idxPath fqPath", argv[0]);
    }

    void *idx = NULL;
    char *soPath = argv[1];
    char *dirPath = argv[2];
    char *fqPath = argv[3];
    
    char *sam = "no sam";
    uint64_t sam_n;
    char *ana = "no ana";
    uint64_t ana_n;
    void *buff = NULL;

    // 加载库函数
    void *handle = dlopen(soPath, RTLD_LAZY);
    void (*load_index_func)(void **idx, const char *dirPath) = dlsym(handle, "load_index");
    void (*read_classify_func)(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id) = dlsym(handle, "read_classify");
    void (*meta_analysis_func)(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id) = dlsym(handle, "meta_analysis");

    // 测试字符串输入api
    // load_index_func(&idx, dirPath);
    // FILE* fq_file = fopen(fqPath, "r");
    // char *input = malloc(1000 * 1000 * 10); // 10MB
    // memset(input, 0, 1000 * 1000 * 10);
    // uint64_t input_n;
    // input_n = fread(input, 1, 1000 * 1000 * 10, fq_file);
    // read_classify_func(idx, input, input_n, &sam, &sam_n, &buff); 
    // meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, &buff);
    // read_classify_func(idx, input, input_n, &sam, &sam_n, &buff); 
    // meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, &buff);
    // read_classify_func(idx, input, input_n, &sam, &sam_n, &buff); 
    // meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, &buff);
    // printf("%s", ana);

    // 测试仅对于sam文件的meta_analysis
    // load_index_func(&idx, dirPath);
    // FILE* sam_file = fopen("demo/beijing.sam", "r");
    // char *input = malloc(1000 * 1000 * 10); // 10MB
    // memset(input, 0, 1000 * 1000 * 10);
    // uint64_t input_n = fread(input, 1, 1000 * 1000 * 10, sam_file);
    // printf("%s", input);
    // meta_analysis_func(idx, input, input_n, &ana, &ana_n);
    // printf("%s", ana);

    // 测试文件输入api
    // load_index_func(&idx, dirPath);
    // read_classify_func(idx, fqPath, -1, &sam, &sam_n, 0); 
    // meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, 0);
    // read_classify_func(idx, fqPath, -1, &sam, &sam_n, 2); 
    // meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, 2);
    // read_classify_func(idx, fqPath, -1, &sam, &sam_n, 1); 
    // meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, 1);
    // printf("%s", ana);

    // 测试内存泄漏
    load_index_func(&idx, dirPath);
    FILE* fq_file = fopen(fqPath, "r");
    char *input = malloc(1000 * 1000 * 10); // 10MB
    memset(input, 0, 1000 * 1000 * 10);
    uint64_t input_n;
    input_n = fread(input, 1, 1000 * 1000 * 10, fq_file);
    for (int i = 0; i < 60000000; ++i) {
        read_classify_func(idx, input, input_n, &sam, &sam_n, i % 3); // 模拟3个线程池
        meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, i % 3);
        printf("%d\t%s", i, ana);
        free(sam);
        free(ana);
    }
}