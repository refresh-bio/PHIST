all: phist subsystem

phist:  utils/phist.cpp
	g++ -std=c++11 -O3 utils/phist.cpp -o utils/phist
subsystem: 
	$(MAKE) -C kmer-db

clean:
	$(MAKE) clean -C kmer-db
	-rm utils/phist  
