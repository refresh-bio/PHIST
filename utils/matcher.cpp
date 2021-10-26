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
#include <unordered_map>

using namespace std;


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
}

int main(int argc, char** argv) {

	std::vector<std::string> params(argc - 1);
	std::transform(argv + 1, argv + argc, params.begin(), [](char* s)->string { return s; });

	int k;
	if (!findOption(params, "-k", k)) {
		return -1;
	}

	const std::string& virPath = params[0];
	const std::string& hostPath = params[1];

	FastaFile virFasta;
	FastaFile hostFasta;
	virFasta.open(virPath);
	hostFasta.open(hostPath);

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
	std::unordered_multimap<kmer_t, uint64_t> hostKmersFov;
	std::unordered_multimap<kmer_t, uint64_t> hostKmersRev;

	SetBasedFilter filter(uniqueKmers);

	for (uint64_t chr_id = 0; chr_id < hostFasta.numSubsequences(); ++chr_id) {
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
			uint64_t global_pos = chr_id << 32 | positions[i];
			hostKmersFov.insert(std::make_pair(kmers[i], global_pos));
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
			uint64_t global_pos = (1LL << 63) | (chr_id << 32) | positions[i];
			hostKmersFov.insert(std::make_pair(kmers[i], global_pos));
		}
	}

	// perform matching from virus point of view
	ofstream outfile(params[2]);
	outfile << "phage,host" << endl;

	for (int vir_cid = 0; vir_cid < virKmerCollections.size(); ++vir_cid) {
		const auto& col = virKmerCollections[vir_cid];
		string vir_header = virFasta.getHeaders()[vir_cid];
		
		uint64_t vir_match_start = 0;
		uint64_t vir_match_last = 0;
		
		for (uint64_t i = 0; i < col.size(); ++i) {

			kmer_t kmer = col[i];
			auto it = hostKmersFov.find(kmer);
			
			if (it != hostKmersFov.end()) {
				uint64_t vir_pos = i + 1; // 1-based indexing

				if (vir_match_start == 0) {
					// start match region
					vir_match_start = vir_match_last = vir_pos;
				}
				else {
					// continue match region
					vir_match_last = vir_pos;
				}

				uint64_t host_pos = it->second;

				bool is_rev = host_pos >> 63;
				host_pos = (host_pos << 1) >> 1; // clear MSB
				uint64_t host_chr_id = host_pos >> 32;
				host_pos = (host_pos << 32) >> 32; // clear upper half
	
			}
			else {
				// check if region has just finished
				if (vir_match_start) {
					outfile
						<< vir_header << ':' << vir_match_start << "-" << vir_match_last + k - 1 << "," << endl;
						//<< hostFasta.getHeaders()[host_chr_id] << ":" << host_pos << "-" << host_pos + k - 1 << endl;
					vir_match_start = vir_match_last = 0;
				}
			}

			
		}
	} 

	outfile.close();

	return 0;
}


