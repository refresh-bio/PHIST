all: phist matcher subsystem

ZLIB_DIR=./kmer-db/libs

phist:  utils/phist.cpp
	g++ -std=c++11 -O3 utils/phist.cpp -o utils/phist

matcher: utils/matcher.cpp utils/input_file.cpp
	g++ -std=c++11 -o utils/matcher -I${ZLIB_DIR} -O3 utils/matcher.cpp utils/input_file.cpp -lz

subsystem: 
	$(MAKE) -C kmer-db

clean:
	$(MAKE) clean -C kmer-db
	-rm utils/phist
	-rm utils/matcher  
