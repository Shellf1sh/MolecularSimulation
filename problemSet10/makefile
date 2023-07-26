CC = g++
CFLAGS = -std=c++11 -Wall

SRCS = mainLJ.cpp functionsLJ.cpp classesLJ.cpp
OBJS = $(SRCS:.cpp=.o)

TARGET = myprogramLJ

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

main.o: mainLJ.cpp functionsLJ.h classesLJ.h
	$(CC) $(CFLAGS) -c mainLJ.cpp

functions.o: functionsLJ.cpp functionsLJ.h
	$(CC) $(CFLAGS) -c functionsLJ.cpp

classes.o: classesLJ.cpp classesLJ.h functionsLJ.h
	$(CC) $(CFLAGS) -c classesLJ.cpp

clean:
	rm -f $(OBJS) $(TARGET)
