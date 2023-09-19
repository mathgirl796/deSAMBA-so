#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "desamba.h"

int main(int argc, char *argv[])
{

    if (argc < 4)
    {
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
    char *human_snapshot;
    uint64_t human_snapshot_n;

    // 加载库函数
    void *handle = dlopen(soPath, RTLD_LAZY);
    void (*load_index_func)(void **idx, const char *dirPath) = dlsym(handle, "load_index");
    void (*read_classify_func)(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int thread_num) = dlsym(handle, "read_classify");
    void (*meta_analysis_func)(void *idx, char *input, uint64_t input_n, char **output, uint64_t *output_n, int thread_id, int flag, uint64_t max_s, char **human_snapshot, uint64_t *hss_n) = dlsym(handle, "meta_analysis");

    // 测试文件输入api
    if (1) {
        load_index_func(&idx, dirPath);
        read_classify_func(idx, fqPath, -1, &sam, &sam_n, 0, 4);
        meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, 0, META_USE_BASE_NUM, 64, &human_snapshot, &human_snapshot_n);
        fprintf(stderr, "<SOS>%s<EOS>\n", ana);
        fprintf(stderr, "<SOS>%lu-%s<EOS>\n", strlen(human_snapshot), human_snapshot);
    }

    // 测试内存泄漏
    if (0)
    {
        for (int i = 0; i < 60000000; ++i)
        {
            read_classify_func(idx, fqPath, -1, &sam, &sam_n, i % 3, 4); // 模拟3个线程池
            meta_analysis_func(idx, sam, sam_n, &ana, &ana_n, i % 3, META_USE_BASE_NUM, 12, &human_snapshot, &human_snapshot_n);
            printf("%d\t%s", i, ana);
            free(sam);
            free(ana);
            free(human_snapshot);
        }
    }
}