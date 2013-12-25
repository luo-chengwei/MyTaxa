CC=g++
CFLAGS=-c -Wall
LDFLAGS=-g -Wall
SOURCES=src/run.cpp src/algo.cpp src/taxonomy.cpp src/utility.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MyTaxa
all:$(SOURCES) $(EXECUTABLE)
$(EXECUTABLE):$(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@
.cpp.o:
		$(CC) $(CFLAGS) $< -o $@
		
clean:
		rm -rf src/*.o
