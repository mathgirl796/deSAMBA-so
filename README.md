# 功能
Add desamba.h and makefile support for .so file, to allow keeping index in memory and classify and analysis countless times.

If you really want to use it, please read the source code. If you can't read well, please contact me.

# 编译so库，结果在bin目录下
```bash
make rebuild -C src
```

# 构建索引
```bash
!!!BUILD jellyfish-1.1.12 by YOURSELF!!!
mkdir bin
mv path/to/jellyfish bin/jellyfish
bash build-index demo/human.fa demo_index/human
cp nodes.dmp names.dmp demo_index/human
```

# 测试程序
```bash
gcc main_test.c -o main_test && ./main_test bin/libdesamba.so demo_index/human demo/human.fq && rm ./main_test
```

# 测试deSAMBA
```bash
bin/deSAMBA classify demo_index/human demo/human.fq > /tmp/deSAMBA.sam
bin/deSAMBA analysis ana_meta /tmp/deSAMBA.sam demo_index/nodes.dmp
```
