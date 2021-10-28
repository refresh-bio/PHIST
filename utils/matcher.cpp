/*******************************************************************************

PHIST
Copyright (C) 2021, A. Zielezinski, S. Deorowicz, and A. Gudys
https://github.com/refresh-bio/PHIST

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program. If not, see https://www.gnu.org/licenses/.

******************************************************************************/
#include "input_file.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <map>
#include <chrono>
#include <iostream>


using namespace std;

union GenomeCoords {
	struct {
		uint32_t pos;
		uint16_t chr;
		uint16_t is_rev;
	};

	uint64_t raw;
};

struct Match {
	uint32_t vir_start;
	uint32_t vir_last;
	GenomeCoords host_start;
	GenomeCoords host_last;

	Match(uint32_t vir_start, uint32_t vir_last, GenomeCoords host_start, GenomeCoords host_last) :
		vir_start(vir_start),
		vir_last(vir_last),
		host_start(host_start),
		host_last(host_last) 
	{}

	void asRanges(std::pair<uint32_t, uint32_t>& vir_range, std::pair<uint32_t, uint32_t>& host_range, int k) {
		
		vir_range.first = vir_start + 1;
		vir_range.second= vir_last + k ; // 1-based indexing
	
		if (host_last.is_rev) {
			host_range.second = host_last.pos + 1; // 1-based indexing
			host_range.first = host_start.pos + k;
		}
		else {
			host_range.first = host_start.pos + 1; // 1-based indexing
			host_range.second = host_last.pos + k;
		}
	}
};


bool findSwitch(std::vector<std::string>& params, const std::string& name) {
	auto it = find(params.begin(), params.end(), name); // verbose mode
	if (it != params.end()) {
		params.erase(it);
		return true;
	}

	return false;
}

template <typename T>
bool findOption(std::vector<std::string>& params, const std::string& name, T& v) {
	auto prevToEnd = std::prev(params.end());
	auto it = find(params.begin(), prevToEnd, name); // verbose mode
	if (it != prevToEnd) {
		std::istringstream iss(*std::next(it));
		if (iss >> v) {
			params.erase(it, it + 2);
			return true;
		}
	}

	return false;
}

int main(int argc, char** argv) {

	cout << "PHIST-Matcher utility 1.0.0" << endl
		<< "A.Zielezinski, S. Deorowicz, A. Gudys (c) 2021" << endl << endl;
	
	std::vector<std::string> params(argc - 1);
	std::transform(argv + 1, argv + argc, params.begin(), [](char* s)->string { return s; });

	int k;
	if (!findOption(params, "-k", k)) {
		k = 25;
	}

	if (params.size() != 3) {
		cout << "USAGE:" << endl
			<< "matcher [-k <length>] <phage> <host> <matches>" << endl << endl
			<< "Parameters:" << endl
			<< "\tlength - minimum match length (25 by default)" << endl
			<< "\tphage - phage FASTA file (gzipped or not)" << endl
			<< "\thost - host FASTA file (gzipped or not)" << endl
			<< "\tmatches - CSV table with all exact matches" << endl;
		return 0;
	}

	auto start = std::chrono::high_resolution_clock::now();

	const std::string& virPath = params[0];
	const std::string& hostPath = params[1];

	cout << "Finding exact matches..." << endl
		<< "minimum length: " << k << endl
		<< "phage FASTA:    " << virPath << endl
		<< "host FASTA:     " << hostPath << endl  << endl;

	FastaFile virFasta;
	FastaFile hostFasta;
	if (!virFasta.open(virPath) || !hostFasta.open(hostPath)) {
		cout << "Unable to open input files" << endl;
		return -1;
	}

	std::vector<std::vector<kmer_t>> virKmerCollections(virFasta.numSubsequences());
	std::unordered_set<kmer_t> uniqueKmers; // this set will be used for filtering host kmers
	AlwaysPassFilter apf;

	// iterate over virus subsequences
	for (int chr_id = 0; chr_id < virKmerCollections.size(); ++chr_id) {
		std::vector<kmer_t>& kmers = virKmerCollections[chr_id];
		
		kmers.resize(virFasta.getLengths()[chr_id] - k + 1);
		
		extract_kmers<KmerMode::Forward, AlwaysPassFilter>(
			virFasta.getSubsequences()[chr_id], 
			virFasta.getLengths()[chr_id], 
			k, 
			apf, 
			kmers.data(), 
			nullptr);

		for (auto kmer : kmers) {
			uniqueKmers.insert(kmer);
		}
	}

	// iterate over host subsequences
	std::multimap<kmer_t, GenomeCoords> hostKmers;


	SetBasedFilter filter(uniqueKmers);

	for (uint16_t chr_id = 0; chr_id < hostFasta.numSubsequences(); ++chr_id) {
		std::vector<kmer_t> kmers(hostFasta.getLengths()[chr_id] - k + 1);
		std::vector<uint32_t> positions(hostFasta.getLengths()[chr_id] - k + 1);

		// forward direction
		size_t count = extract_kmers<KmerMode::Forward, SetBasedFilter>(
			hostFasta.getSubsequences()[chr_id],
			hostFasta.getLengths()[chr_id],
			k,
			filter,
			kmers.data(),
			positions.data());

		for (int i = 0; i < count; ++i) {
			GenomeCoords coords = { positions[i], chr_id, 0 };
			hostKmers.insert(std::make_pair(kmers[i], coords));
		}

		// reverse direction
		count = extract_kmers<KmerMode::Reverse, SetBasedFilter>(
			hostFasta.getSubsequences()[chr_id],
			hostFasta.getLengths()[chr_id],
			k,
			filter,
			kmers.data(),
			positions.data());

		for (int i = 0; i < count; ++i) {
			GenomeCoords coords = { positions[i], chr_id, 1 };
			hostKmers.insert(std::make_pair(kmers[i], coords));
		}
	}

	// perform matching from virus point of view
	ofstream outfile(params[2]);
	outfile << virPath << "," << hostPath << endl;

	std::vector<Match> matches;

	// iterate over virus chromosomes
	for (int vir_cid = 0; vir_cid < virKmerCollections.size(); ++vir_cid) {
		const auto& col = virKmerCollections[vir_cid];
		string vir_header = virFasta.getHeaders()[vir_cid];
		
		// iterate over virus positions
		for (uint64_t vir_pos = 0; vir_pos < col.size(); ++vir_pos) {

			// get host 
			kmer_t kmer = col[vir_pos];
			auto range = hostKmers.equal_range(kmer);
			
			if (range.first != range.second) {
				
				// iterate over hit range of host positions
				for (auto it = range.first; it != range.second; ++it) {
					GenomeCoords host_hit = it->second;
					
					bool consumed = false;

					// iterate over current matches and check if hit extends any of it
					for (auto& match : matches) {
						if (host_hit.raw == match.host_last.raw - 1 && host_hit.is_rev) {
							// continue reverse match
							--match.host_last.pos;
							++match.vir_last;
							consumed = true;
						}
						else if (host_hit.raw == match.host_last.raw + 1 && !host_hit.is_rev) {
							// continue forward match
							++match.host_last.pos;
							++match.vir_last;
							consumed = true;
						}
					}

					// if hit does not extend any existing match - create a new one
					if (!consumed) {
						matches.emplace_back(vir_pos, vir_pos, host_hit, host_hit);
					}

				}
			}

			// save and remove unextended matches
			for (auto it = matches.begin(); it != matches.end(); )
			{
				if (it->vir_last != vir_pos) {
					
					std::pair<uint32_t, uint32_t> vir_range, host_range;
					it->asRanges(vir_range, host_range, k);
					
					outfile
						<< vir_header << ':' << vir_range.first << "-" << vir_range.second << ","  
						<< hostFasta.getHeaders()[it->host_last.chr] << ":" << host_range.first << "-" << host_range.second << endl;

					it = matches.erase(it);
				}
				else {
					++it;
				}
			}	
		}

		// if there are some matches left
		for (auto it = matches.begin(); it != matches.end(); ++it) {
			std::pair<uint32_t, uint32_t> vir_range, host_range;
			it->asRanges(vir_range, host_range, k);

			outfile
				<< vir_header << ':' << vir_range.first << "-" << vir_range.second << ","
				<< hostFasta.getHeaders()[it->host_last.chr] << ":" << host_range.first << "-" << host_range.second << endl;
		}
		matches.clear();
	} 

	outfile.close();

	auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start);
	cout << "Finished in " << time.count() << " seconds" << endl;

	return 0;
}


