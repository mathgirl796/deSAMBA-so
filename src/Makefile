SOURCE  := $(wildcard lib/*.c) $(wildcard *.c) $(wildcard *.cpp)
OBJS    := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
  
TARGET  := deSAMBA
TARGET_SO := libdesamba.so
  
#compile and lib parameter
CC      := gcc -std=gnu99
LIBS    := -lm -lz -lpthread
INCLUDE := -I ./lib
CFLAGS  := -g -Wall -O3 $(INCLUDE)

.PHONY : everything objs clean veryclean rebuild
  
everything : $(TARGET) $(TARGET_SO)
  
objs : $(OBJS)
  
rebuild: veryclean everything
                
clean :
	rm -fr *.so
	rm -fr *.o
	rm -fr lib/*.so
	rm -fr lib/*.o
    
veryclean : clean
	rm -fr ../bin/$(TARGET)
	rm -fr ../bin/$(TARGET_SO)
  
$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	mkdir -p ../bin
	mv $@ ../bin/

$(TARGET_SO) : $(OBJS)
	$(CC) $(CFLAGS) -shared -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	mkdir -p ../bin
	mv $@ ../bin/

%.o: %.c
	$(CC) $(CFLAGS) -c -fPIC $< -o $@