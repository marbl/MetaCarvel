DEST_DIR = ~/bin

CFLAGS =  -O3 -Wall -Wextra -std=c++11
SPQRFLAGS =  -lOGDF -lCOIN -pthread 

ALL = libcorrect bundler orientcontigs spqr 

all: $(ALL)

libcorrect: 
	g++ $(CFLAGS) -o libcorrect libcorrect.cpp

bundler: 
	g++ $(CFLAGS) -o bundler bundler.cpp

orientcontigs: 
	g++ $(CFLAGS) -o orientcontigs orientcontigs.cpp

spqr:
	g++  -L/cbcb/sw/RedHat-7-x86_64/users/jayg/local/OGDF/2015.05/_release/ -I/cbcb/sw/RedHat-7-x86_64/users/jayg/local/OGDF/2015.05/include -o spqr spqr.cpp $(CFLAGS) $(SPQRFLAGS)

clean:
	rm -f $(ALL)

install:
	 cp $(ALL) $(DEST_DIR)
