CC= g++
OPT=-g

BUILD_DATE := "$(shell date)"
LAVA_UPDATE := "17.04.2018"
LAVA_VERSION := "1.0"

CPPFLAGS = -std=c++11 $(OPT) -I /home/ezgi/tools/include -I sonic -DBUILD_DATE=\"$(BUILD_DATE)\" -DLAVA_UPDATE=\"$(LAVA_UPDATE)\" -DLAVA_VERSION=\"$(LAVA_VERSION)\"
LDFLAGS = usr/local/lib/libhts.a usr/local/lib/libsonic.a -lz -lm -lpthread
SOURCES = stringOp.cpp read.cpp sv.cpp cluster.cpp bam.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = lava

all: $(SOURCES) $(EXECUTABLE)
	rm -f *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(OPT) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -c $(CPPFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) *.o *~

libs:
	make -C sonic

