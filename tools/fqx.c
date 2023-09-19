#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>

void processFile(const char *input_file_path, const char *output_file_path, int threshold)
{
    FILE *input_file = fopen(input_file_path, "rb");
    gzFile output_file = gzopen(output_file_path, "wb");
    gzbuffer(output_file, 1000000);

    if (input_file == NULL || output_file == NULL)
    {
        printf("无法打开输入或输出文件。\n");
        return;
    }

    char buf[threshold];
    int length = 0;
    int id = 0;
    int ch;
    while (!feof(input_file))
    {
        ch = fgetc(input_file);
        int is_acgt = ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't' ;
            
        if (is_acgt)
        { 
            if (length < threshold) // 新acgt串，长度尚未达标
            { 
                buf[length] = ch;
                length += 1;
                if (length == threshold) // 达标后立刻写入头部
                {
                    gzprintf(output_file, "@%d\n", id);
                    for (int i = 0; i < length; ++i) {
                        gzprintf(output_file, "%c", buf[i]);
                    }
                }
            }
            else // acgt串长度达标，继续写入
            { 
                gzprintf(output_file, "%c", ch);
                length += 1;
            }
        }
        if (!is_acgt || feof(input_file))
        { 
            if (length < threshold) // 长度不达标，丢弃
            {
                length = 0;
            }
            else // 长度达标，一定已有头部写入，此处补充质量分数
            {
                gzprintf(output_file, "\n+\n");
                for (int i = 0; i < length; ++i) {
                    gzprintf(output_file, "%c", 'I');
                }
                gzprintf(output_file, "\n");
                length = 0;
                id += 1;
            }
        }
    }
    fclose(input_file);
    gzclose(output_file);
}

// 测试
int main(int argc, char ** argv) {
    processFile(argv[1], argv[2], 10);
}