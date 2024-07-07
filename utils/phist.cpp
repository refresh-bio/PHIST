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


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>


using namespace std;


long int strtol(const char* str, char** endptr) 
{
	
	long int val = 0;
	char *p = (char*)str;
	bool is_negative = false;


	if (*p == '-')
	{
		is_negative = true;
		++p;
	}

	while (*p >= '0' && *p <= '9')
	{
		val = val * 10 + (*p++ - '0');
	}

	if (endptr)
		*endptr = p;

	return is_negative ? -val : val;
}



struct Organism {
public:
	string name;
	uint32_t kmer_count;

	template<class Iterator>
	Organism(Iterator name_begin, Iterator name_end) : name(name_begin, name_end), kmer_count(0) {}
};

struct Hit {
	uint32_t host_id;
	uint32_t common_kmers;
	
	Hit(uint32_t host_id, uint32_t common_kmers) : host_id(host_id), common_kmers(common_kmers) {}

};

struct Phage : public Organism {
public:
	

	std::vector<Hit> hits;

	template<class Iterator>
	Phage(Iterator name_begin, Iterator name_end) : Organism(name_begin, name_end) {}
};


int main(int argc, char** argv) {
	
	cout << "PHIST utility 1.0.0" << endl
		<< "A.Zielezinski, S. Deorowicz, A. Gudys (c) 2021" << endl << endl;
	
	vector<string> params;
	
	for (int i = 1; i < argc; ++i) {
		params.push_back(argv[i]);
	}

	if (params.size() != 2) {
		cout << "USAGE:" << endl
			<< "phist <input> <output>" << endl << endl
			<< "Parameters:" << endl
			<< "\tinput - CSV file in a sparse format with a number of common k-mers between phages and bacteria" << endl
			<< "\t        (result of running `kmer-db new2all -sparse phages.db bacteria.list`)," << endl
			<< "\toutput - CSV file with assignments of phages to their most probable hosts" << endl;
	}

	auto start = std::chrono::high_resolution_clock::now();
	
	ifstream input(params[0]);

	if (!input) {
		return 0;
	}

	size_t bufsize = 1 << 30; // 1 GB line buffer
	char* line = new char[bufsize]; 

	vector<Phage> phages;
	vector<Organism> bacteria;

	//
	// Extract phages names
	//
	input.getline(line, bufsize);
	char * end = line + input.gcount() - 1;  // do not count \n

	// get k-mer length
	char * begin = line;
	char * p = std::find(begin, end, ':');
	begin = p+2;
	uint32_t k = strtol(begin, &p);

	begin = line;
	p = std::find(begin, end, ',');
	p = std::find(p + 1, end, ',');

	begin = p + 1;
	do {
		p = std::find(begin, end, ',');
		phages.emplace_back(begin, p);
		begin = p + 1;
	} while (end - begin > 1);

	//
	// Extract phages k-mers count
	// 
	input.getline(line, bufsize);
	end = line + input.gcount() - 1;  // do not count \n
	
	// omit two first cells
	begin = line;
	p = std::find(begin, end, ',');
	p = std::find(p + 1, end, ',');
	
	begin = p + 1;
	p = end;
	int phage_id = 0;
	do {
		phages[phage_id].kmer_count = strtol(begin, &p); // assume no white characters after the number -> p points comma
		++phage_id;
		begin = p + 1;
	} while (end - begin > 1);

	//
	// Process bacteria
	//
	cout << "Processing bacteria from Kmer-db table..." << endl;

	uint32_t bact_id = 0;
	while (input.getline(line, bufsize)) {
		// show progress
		if ((bact_id + 1) % 10 == 0) {
			cout << "\r" << bact_id + 1 << "..." << std::flush;
		}
		
		// extract name
		end = line + input.gcount() - 1; // do not count \n
		begin = line;
		p = std::find(begin, end, ',');
		bacteria.emplace_back(begin, p);
		begin = p + 1;

		// extract kmer count
		Organism & bact = bacteria.back();
		bact.kmer_count = strtol(begin, &p); // assume no white characters after the number -> p points comma
		begin = p + 1;
		
		// extract number of common kmers
		while (end - begin > 1) {
			// each entry is in the form <phage_id>:<common_kmers_count>
			
			uint32_t phage_id = strtol(begin, &p); // assume no white characters after number -> p points colon
			--phage_id; // indexing in file is 1-based

			Phage& phage = phages[phage_id];
			
			begin = p + 1; 
			uint32_t common_kmers = strtol(begin, &p); // assume no white characters after number -> p points comma
			begin = p + 1;

			if (phage.hits.empty() || common_kmers == phage.hits.front().common_kmers) {
				// empty collection or same as current best - add new 
				phage.hits.emplace_back(bact_id, common_kmers);
			}
			else if (common_kmers > phage.hits.front().common_kmers) {
				// better then current best - replace
				phage.hits.clear();
				phage.hits.emplace_back(bact_id, common_kmers);
			}
			

		}

		++bact_id;
	}
	cout << "\r" << bact_id << " [OK]" << endl;
	input.close();
	
	//
	// Store results
	// 
	ofstream output(params[1]);

	output << "phage,host,#common-kmers,pvalue,adj-pvalue" << endl;

	for (Phage& ph : phages) {
		
		// no host
		if (ph.hits.empty()) {
			output << ph.name << endl;
		}
		else {
			// sort increasingly by the host length
			std::stable_sort(ph.hits.begin(), ph.hits.end(), [&bacteria](const Hit& h1, const Hit& h2)->bool {
				return bacteria[h1.host_id].kmer_count < bacteria[h2.host_id].kmer_count;
			});
			
			for (const auto& hit : ph.hits) {
				
				const Organism& host = bacteria[hit.host_id];
				int len_common = hit.common_kmers + k - 1;
				int len_host = host.kmer_count + k - 1;
				int len_phage = ph.kmer_count + k - 1;
				
				long double num_canonical = (len_common % 2)
					? pow(4, len_common) / 2
					: (pow(4, len_common) + pow(4, len_common / 2)) / 2;

				long double lambda = (double)(len_host - len_common + 1) * (len_phage - len_common + 1) / num_canonical;
				long double pval = 1 - std::exp(-lambda);
				
				// adjust by the number of potential hosts
				long double adj_pval = std::min(bact_id * pval, (long double)1.0);

				output << ph.name << ',' << bacteria[hit.host_id].name << ',' << hit.common_kmers << ',' << std::scientific << pval << "," << adj_pval << endl;
			}
		}
	}

	//
	output.close();
	delete[] line;

	auto time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -start);
	cout << "File analyzed in " << time.count() << " seconds" << endl;

}
