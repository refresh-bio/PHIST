#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <cassert>

#include "kmer_helper.h"

// *****************************************************************************************
//


class FastaFile {
public:

	const std::vector<char*>& getSubsequences() const { return subsequences; }
	const std::vector<size_t>& getLengths() const { return lengths; }
	const std::vector<char*>& getHeaders() const { return headers; }

	size_t numSubsequences() const { return subsequences.size(); }


	FastaFile() : rawSize(0), rawData(nullptr), totalLen(0), status(true), isGzipped(false) {}
	~FastaFile() { delete [] rawData;  }

	bool open(const std::string& filename);
	
	bool close() { delete[] rawData; return true; }


protected:
	size_t rawSize;
	char* rawData;
	size_t totalLen;
	bool status;
	bool isGzipped;

	std::vector<char*> subsequences;
	std::vector<size_t> lengths;
	std::vector<char*> headers;

	bool unzip(char* compressedData, size_t compressedSize, char*&outData, size_t &outSize);
	
	bool extractSubsequences(
		char* data,
		size_t& totalLen,
		std::vector<char*>& subsequences,
		std::vector<size_t>& lengths,
		std::vector<char*>& headers);

	
};
