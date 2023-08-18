# Directories
SRCDIR = src
BUILDDIR = build

# Compiler and Linker flags
CC = mpicc
CFLAGS = -Wall -Werror -fopenmp 
LDFLAGS = -lm

# Target executable
TARGET = $(BUILDDIR)/program

# Source and Object files
SRC = $(SRCDIR)/main.c
OBJ = $(BUILDDIR)/main.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

$(OBJ): $(SRC)
	mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $(SRC) -o $(OBJ)

clean:
	rm -rf $(BUILDDIR)
