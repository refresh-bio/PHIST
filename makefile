all: phist matcher subsystem ng_zlib

ZLIB_DIR=./3rd_party/zlib-ng

phist: utils/phist.cpp ng_zlib subsystem
	$(CXX) -std=c++11 -O3 utils/phist.cpp -o utils/phist

matcher: utils/matcher.cpp utils/input_file.cpp ng_zlib
	$(CXX) -std=c++11 -o utils/matcher -I${ZLIB_DIR} -O3 utils/matcher.cpp utils/input_file.cpp $(ZLIB_DIR)/libz.a

ng_zlib:
	cd $(ZLIB_DIR) && ./configure --zlib-compat && $(MAKE) libz.a

subsystem: 
	$(MAKE) -C kmer-db

clean:
	$(MAKE) clean -C kmer-db
	cd $(ZLIB_DIR) && $(MAKE) -f Makefile.in clean
	-rm $(ZLIB_DIR)/libz.a
	-rm utils/phist
	-rm utils/matcher  
