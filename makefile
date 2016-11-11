DEST_DIR = ~/bin

CFLAGS =  -O3 -Wall -Wextra -std=c++11
SPQRFLAGS =  -lOGDF -lCOIN -pthread 

ALL = libcorrect bundler orientcontigs 

all: $(ALL)

libcorrect: 
	g++ $(CFLAGS) -o libcorrect lib_correct_bed.cpp

bundler: 
	g++ $(CFLAGS) -o bundler bundler.cpp

orientcontigs: 
	g++ $(CFLAGS) -o orientcontigs orientcontigs.cpp

clean:
	rm -f $(ALL)

install:
	 cp $(ALL) $(DEST_DIR)
