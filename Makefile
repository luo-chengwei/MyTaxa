CC=g++
CFLAGS=-c -Wall
LDFLAGS=-g -Wall
SOURCES=run.cpp algo.cpp taxonomy.cpp utility.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=MeTaxa
all:$(SOURCES) $(EXECUTABLE)
$(EXECUTABLE):$(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@
.cpp.o:
		$(CC) $(CFLAGS) $< -o $@
		
clean:
		rm -rf *.o
