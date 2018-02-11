CC=g++
CXXFLAGS=-std=c++11
OBJS=diffevol.o ions.o

ions: $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o ions

all: ions

clean:
	rm $(OBJS)
