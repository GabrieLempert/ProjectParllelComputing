CC = gcc
CFLAGS = -Wall -Wextra -std=c11
LDFLAGS =

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
