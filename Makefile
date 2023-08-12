CC = mpicc
CFLAGS = -Wall -Wextra -g -std=c11
LDFLAGS = -lm

SRCS = main.c
OBJS = $(SRCS:.c=.o)
EXECUTABLE = program

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJS)
