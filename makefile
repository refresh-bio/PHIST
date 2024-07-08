all: phist matcher subsystem ng_zlib

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
    uname_M := "x86_64"
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
    uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
endif

CFLAGS=-O3 -std=c++11

ifeq ($(uname_S),Linux)
	CFLAGS+=-fabi-version=6
	CFLAGS+=-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif

ifeq ($(uname_S),Darwin)
	CFLAGS+= -lc -static-libgcc
endif


ZLIB_DIR=./3rd_party/zlib-ng

phist: utils/phist.cpp ng_zlib subsystem
	$(CXX) $(CFLAGS) utils/phist.cpp -o utils/phist

matcher: utils/matcher.cpp utils/input_file.cpp ng_zlib
	$(CXX) $(CFLAGS) -o utils/matcher -I${ZLIB_DIR} utils/matcher.cpp utils/input_file.cpp $(ZLIB_DIR)/libz.a

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
