CXX = g++

CXXFLAGS = -std=c++17 -O2 -Wall $(shell root-config --cflags)
GLIBS = $(shell root-config --glibs)

TARGET = makeTriggerEfficiency.x
SRCDIR = ./src/
OBJDIR = ./obj/
SRC = $(SRCDIR)makeTriggerEfficiency.C
OBJ = $(OBJDIR)makeTriggerEfficiency.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(GLIBS)

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
