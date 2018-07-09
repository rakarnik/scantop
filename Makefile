#
# Macros
#

CC = /usr/bin/g++
CC_OPTIONS = -O3 -Wall -Wextra
#CC_OPTIONS = -O0 -ggdb -Wall -Wextra
LNK_OPTIONS =
BIN_DIR = bin

#
# INCLUDE directories
#
INCLUDE = -I. -Isrc


#
# Build scantop
#

all: scantop

scantop: \
		bin/bgmodel.o\
		bin/motif.o\
		bin/scantop.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/bgmodel.o\
		bin/motif.o\
		bin/scantop.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o\
		-o bin/scantop

clean: 
	rm -f bin/*.o bin/scantop

dir_guard=@mkdir -p $(@D)

$(BIN_DIR)/%.o: src/%.cpp
	$(dir_guard)
	$(CC) $(CC_OPTIONS) -c $(INCLUDE) -o $@ $<
