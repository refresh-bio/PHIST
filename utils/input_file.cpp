/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/
#include "input_file.h"

#include <zlib.h>

#include <memory>
#include <fstream>
#include <string>
#include <cassert>
#include <bitset>
#include <cstdio>



#ifndef WIN32
	#define my_fopen    fopen
	#define my_fseek    fseek
	#define my_ftell    ftell
#else
	#define my_fopen    fopen
	#define my_fseek    _fseeki64
	#define my_ftell    _ftelli64
#endif

// *****************************************************************************************
//
 bool FastaFile::open(const std::string& filename) {

	status = false;

	FILE * in;
	
	// try to open without adding extension
	in = fopen(filename.c_str(), "rb");
	isGzipped = filename.substr(filename.length() - 3) == ".gz";
	
	if (!in) {
		return status;
	}

	my_fseek(in, 0, SEEK_END);
	rawSize = my_ftell(in);
	my_fseek(in, 0, SEEK_SET);

	rawData = reinterpret_cast<char*>(malloc(rawSize + 1));
	size_t blocksRead = fread(rawData, rawSize, 1, in);
	rawData[rawSize] = 0; // add null termination 
	fclose(in);
		
	if (!blocksRead) {
		return status;
	}

	char* data;
	size_t total = 0;

	if (isGzipped) {
		status = unzip(rawData, rawSize, data, totalLen);
		if (!status) {
			return status;
		}
	}
	else {
		data = rawData;
		totalLen = rawSize;
	}

	status = extractSubsequences(data, totalLen, subsequences, lengths, headers);

	return status;
}

 // *****************************************************************************************
 //
 /*
bool GenomeInputFile::load(
	uint32_t kmerLength,
	std::vector<kmer_t>& kmersBuffer,
	std::vector<uint32_t>& positionsBuffer) {
	
	if (!status) {
		return false;
	}

	

	
		// determine max k-mers count
		size_t sum_sizes = 0;
		for (auto e : lengths)
			sum_sizes += e - kmerLength + 1; 

		kmersBuffer.clear();
		kmersBuffer.resize(sum_sizes);

		positionsBuffer.clear();
		positionsBuffer.resize(sum_sizes);

		uint32_t kmersCount = 0;
		kmer_t* currentKmers = kmersBuffer.data();
		uint32_t* currentPositions = positionsBuffer.data();


		for (size_t i = 0; i < chromosomes.size(); ++i) {
			size_t count = extractKmers(chromosomes[i], lengths[i], kmerLength, currentKmers, currentPositions);
			currentKmers += count;
			currentPositions += count;
			kmersCount += count;
		}
	
		//LOG_DEBUG << "Extraction: " << kmersCount << " kmers, " << chromosomes.size() << " chromosomes, " << totalLen << " bases" << endl;
	}
	
	// free memory
	if (data != rawData) {
		free(reinterpret_cast<void*>(data));
	}
	free(reinterpret_cast<void*>(rawData));
	
	return status;
}
*/


// *****************************************************************************************
//
bool FastaFile::unzip(char* compressedData, size_t compressedSize, char*&outData, size_t &outSize) {
	bool ok = true;
	size_t blockSize = 10000000;
	outData = reinterpret_cast<char*>(malloc(blockSize));

	// Init stream structure
	z_stream stream;
	stream.zalloc = Z_NULL;
	stream.zfree = Z_NULL;
	stream.opaque = Z_NULL;
	stream.avail_in = compressedSize;
	stream.next_in = reinterpret_cast<Bytef*>(compressedData);

	if (inflateInit2(&stream, 31) == Z_OK) {

		// read data in portions
		char *ptr = outData;

		size_t allocated = blockSize;

		// decompress file in portions
		for (;;) {
			stream.avail_out = blockSize;
			stream.next_out = reinterpret_cast<Bytef*>(ptr);
			int ret = inflate(&stream, Z_NO_FLUSH);

			switch (ret)
			{
			case Z_NEED_DICT:
			case Z_DATA_ERROR:
			case Z_MEM_ERROR:
				ok = false;
				ret = Z_STREAM_END;
				break;
			}

			if (ret == Z_OK && stream.avail_out == 0) {
				outSize = stream.total_out;
			}

			if (ret == Z_STREAM_END) {
				outSize = stream.total_out;
				//multistream detection
				if (stream.avail_in >= 2 && stream.next_in[0] == 0x1f && stream.next_in[1] == 0x8b) {
					if (inflateReset(&stream) != Z_OK) {
						//LOG_NORMAL << "Error while reading gzip file\n";
						exit(1);
					}
				}
				else
					break;
			}

			// reallocate only when some data left
			allocated += blockSize;
			outData = reinterpret_cast<char*>(realloc(outData, allocated + 1)); // allocate for null termination
			ptr = outData + outSize;
		}

		inflateEnd(&stream);
		outData[outSize] = 0;
	}
	else {
		ok = false;
	}

	return ok;
}

// *****************************************************************************************
//
bool FastaFile::extractSubsequences(
	char* data,
	size_t& totalLen,
	std::vector<char*>& subsequences,
	std::vector<size_t>& lengths,
	std::vector<char*>& headers) {

	// extract contigs
	char * header = nullptr;
	char * ptr = data;
	
	while ((header = strchr(ptr, '>'))) { // find begining of header
		*header = 0; // put 0 as separator (end of previous chromosome)
		if (subsequences.size()) {
			lengths.push_back(header - subsequences.back());
		}

		++header; // omit '<'
		headers.push_back(header);

		ptr = strchr(header, '\n'); // find end of header
		if (*(ptr - 1) == '\r') { // on Windows
			*(ptr - 1) = 0;
		}
		*ptr = 0; // put 0 as separator
		++ptr; // move to next character (begin of chromosome)
		subsequences.push_back(ptr); // store chromosome
	}

	lengths.push_back(data + totalLen - subsequences.back());

	// remove newline characters from chromosome
	totalLen = 0;
	for (size_t i = 0; i < subsequences.size(); ++i) {
		// determine chromosome end
		char* newend = std::remove_if(subsequences[i], subsequences[i] + lengths[i], [](char c) -> bool { return c == '\n' || c == '\r';  });
		*newend = 0;
		lengths[i] = newend - subsequences[i];
		totalLen += lengths[i];
	//	assert(lengths[i] == strlen(subsequences[i]));
	}

	return true;
}






