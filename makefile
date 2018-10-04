DEST_DIR = ~/bin

CFLAGS =  -O3 -Wall -Wextra -std=c++11
SPQRFLAGS =  -lOGDF -lCOIN -pthread 

#####MODIFY THESE PATHS BASED ON YOUR INSTALLATION LOCATION####
OGDF_INCL = -I OGDF/include/
OGDF_LINK = -L OGDF/_release/
############################


ALL = libcorrect bundler orientcontigs spqr

all: $(ALL)

libcorrect: 
	g++ $(CFLAGS) -o libcorrect libcorrect.cpp

bundler: 
	g++ $(CFLAGS) -o bundler bundler.cpp

orientcontigs: 
	g++ $(CFLAGS) -o orientcontigs orientcontigs.cpp

spqr:
	g++ spqr.cpp $(CFLAGS) $(OGDF_INCL) $(OGDF_LINK) $(SPQRFLAGS) -o spqr

clean:
	rm -f $(ALL)

install:
	 cp $(ALL) $(DEST_DIR)
