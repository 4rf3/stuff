
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include "delaunator.hpp"

using namespace std::literals;

std::vector<double> read_data(std::string fname)
{
	std::ifstream inf{ fname };
	std::string header;
	getline(inf, header);
	std::stringstream ss{ header };
	size_t x_col = 0, y_col = 0;
	size_t cols = 1;
	for (; !ss.eof(); cols++)
	{
		std::string head;
		ss >> head;
		if (head == "XM")
			x_col = cols;
		if (head == "YM")
			y_col = cols;
	}
	std::vector<double> data;
	while (!inf.eof())
	{
		std::string str;
		double x, y;
		for (size_t i = 0; i < cols; i++)
		{
			if (i == x_col)
				inf >> x;
			else if (i == y_col)
				inf >> y;
			else
				inf >> str;
		}
		data.push_back(x);
		data.push_back(y);
	}
	return data;
}

void write_data_nb_counts(std::string fname, const std::vector<double>& coords, const std::vector<size_t>& nb_counts)
{
	std::ofstream ofs{ fname };
	ofs << "X\tY\tneighbors\n";
	for (size_t i = 0; i < coords.size() / 2; i++)
	{
		ofs << coords[2 * i] << "\t" << coords[2 * i + 1] << "\t" << nb_counts[i] << "\n";
	}
}

struct res_data
{
	bool valid{ false };
	std::vector<std::pair<std::string, std::string>> data;
	double thr{ 0.95 };
	bool silent{ false }, benchmark{ false };
};

void process_pair(std::string f0, std::string f1, const res_data& conf)
{
	std::cout << "Processing " << f0 << " + " << f1;
	auto start_timer = std::chrono::high_resolution_clock::now();
	auto data0 = read_data(f0);
	auto data1 = read_data(f1);
	if (!conf.silent) std::cout << " .";
	size_t count0 = (data0.end() - data0.begin()) / 2, count1 = (data1.end() - data1.begin()) / 2;
	auto data{ data0 };
	copy(data1.begin(), data1.end(), back_inserter(data));
	size_t total = data.size() / 2;
	auto triang = delaunator::Delaunator{ data };
	if (!conf.silent) std::cout << " .";
	double avg_sqdist = 0;
	std::vector<double> sqdist;
	for (size_t i = 0; i < triang.triangles.size(); i += 3)
	{
		double dx = data[2 * triang.triangles[i]] - data[2 * triang.triangles[i + 1]];
		double dy = data[2 * triang.triangles[i] + 1] - data[2 * triang.triangles[i + 1] + 1];
		sqdist.push_back(dx * dx + dy * dy);
		dx = data[2 * triang.triangles[i + 1]] - data[2 * triang.triangles[i + 2]];
		dy = data[2 * triang.triangles[i + 1] + 1] - data[2 * triang.triangles[i + 2] + 1];
		sqdist.push_back(dx * dx + dy * dy);
		dx = data[2 * triang.triangles[i]] - data[2 * triang.triangles[i + 2]];
		dy = data[2 * triang.triangles[i] + 1] - data[2 * triang.triangles[i + 2] + 1];
		sqdist.push_back(dx * dx + dy * dy);
	}
	if (!conf.silent) std::cout << " .";
	sort(sqdist.begin(), sqdist.end());
	if (!conf.silent) std::cout << " .";
	double threshold = sqdist[(size_t)sqdist.size() * conf.thr];
	std::vector<size_t> neigh0(data.size() / 2, 0), neigh1(data.size() / 2, 0); //neigh 0 are neighbors with type=0
	for (size_t i = 0; i < triang.triangles.size(); i += 3)
	{
		// [0], [1], [2] are CCW, count INTO the edge
		// [0] -> [1]
		double dx = data[2 * triang.triangles[i]] - data[2 * triang.triangles[i + 1]];
		double dy = data[2 * triang.triangles[i] + 1] - data[2 * triang.triangles[i + 1] + 1];
		double dist = dx * dx + dy * dy;
		if (dist < threshold)
		{
			if (triang.triangles[i] < count0)
				neigh0[triang.triangles[i + 1]] += 1;
			else
				neigh1[triang.triangles[i + 1]] += 1;
		}
		// [1] -> [2]
		dx = data[2 * triang.triangles[i + 1]] - data[2 * triang.triangles[i + 2]];
		dy = data[2 * triang.triangles[i + 1] + 1] - data[2 * triang.triangles[i + 2] + 1];
		dist = dx * dx + dy * dy;
		if (dist < threshold)
		{
			if (triang.triangles[i + 1] < count0)
				neigh0[triang.triangles[i + 2]] += 1;
			else
				neigh1[triang.triangles[i + 2]] += 1;
		}
		// [2] -> [0]
		dx = data[2 * triang.triangles[i]] - data[2 * triang.triangles[i + 2]];
		dy = data[2 * triang.triangles[i] + 1] - data[2 * triang.triangles[i + 2] + 1];
		dist = dx * dx + dy * dy;
		if (dist < threshold)
		{
			if (triang.triangles[i + 2] < count0)
				neigh0[triang.triangles[i]] += 1;
			else
				neigh1[triang.triangles[i]] += 1;
		}
	}
	if (!conf.silent) std::cout << " .\n";
	size_t	cnt0_near0 = 0,
		cnt1_near0 = 0,
		cnt0_near1 = 0,
		cnt1_near1 = 0;
	for (size_t i = 0; i < count0; i++)
	{
		cnt0_near0 += neigh0[i];
		cnt1_near0 += neigh1[i];
	}
	for (size_t i = count0; i < total; i++)
	{
		cnt0_near1 += neigh0[i];
		cnt1_near1 += neigh1[i];
	}
	auto suffix = ".neighbors"s;
	auto suff0 = f0.substr(f0.rfind('.'));
	auto fout0 = f0.replace(f0.rfind(suff0), suff0.length(), suffix + ".txt"s);
	write_data_nb_counts(fout0, data0, neigh0);
	std::vector<size_t> ng1;
	copy(neigh1.begin() + count0, neigh1.end(), back_inserter(ng1));
	auto suff1 = f1.substr(f1.rfind('.'));
	auto fout1 = f1.replace(f1.rfind(suff1), suff1.length(), suffix + ".txt"s);
	write_data_nb_counts(fout1, data1, ng1);
	if (conf.benchmark && !conf.silent)
	{
		auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_timer);
		std::cout << "Processing took " << elapsed.count() / 1000.0 << " ms\n";
	}
	if (!conf.silent)
	{
		std::cout << "(0) is " << f0 << ", written to " << fout0 << ",\n(1) is " << f1 << ", written to " << fout1 << ":\n";
		std::cout << "Average of " << (double)cnt0_near0 / count0 << " #(0) adjacent to (0) in center\n";
		std::cout << "Average of " << (double)cnt1_near1 / count1 << " #(1) adjacent to (1) in center\n";
		std::cout << "Average of " << (double)cnt1_near0 / count0 << " #(1) adjacent to (0) in center\n";
		std::cout << "Average of " << (double)cnt0_near1 / count1 << " #(0) adjacent to (1) in center\n\n";
	}
}

void show_help()
{
	std::cout << R"(Description: Nearest neighbor counter.
Given a set of points, split into two kinds, for every point
count nearest neighbors of the same kind.
Usage:  neighbor_counter [-ish?] [-t threshold] input-file1 input-file2 ...
Options:
  -t  Threshold parameter, a decimal number in (0,1). The defalt value is 0.95.
      Make it less if there are missing data points in the data sets
  -i  Changes how the input files are grouped. E.g., for input files A B C D E F grouping is:
      without -i:  A   B   C   D   E   F
                   \___/   \___/   \___/
      with -i:     A   B   С   D   E   F
                   \___\___\___/  /   /
                        \___\____/   /
                             \______/
  -s  Silent output. Do not write progress and results to console.
  -h, -?
      Show this help and quit)";
}

res_data get_data(size_t argc, char* argv[])
{
	res_data data_sets;
	size_t arg = 1;
	bool interleaved = true;
	while (arg < argc && argv[arg][0] == '-')
	{
		if (strspn(argv[arg] + 1, "tish?b") == 0)
		{
			std::cout << "Invalid parameter: " << argv[arg] << std::endl;
			return {};
		}
		for (size_t i = 1; i < strlen(argv[arg]); i++)
		{
			switch (argv[arg][1])
			{
			case 'i':
				interleaved = false;
				break;
			case 's':
				data_sets.silent = true;
				break;
			case 'h':
			case '?':
				show_help();
				return {};
				break;
			case 't':
			{
				i = strlen(argv[arg]);
				++arg;
				std::istringstream thr_str{ argv[arg] };
				thr_str >> data_sets.thr;
			}
				break;
			case 'b':
				data_sets.benchmark = true;
				break;
			default:
				break;
			}
		}		
		++arg;
	}
	if ((argc - arg) % 2 == 1)
	{
		std::cout << "Odd number of data files: " << argc - arg << std::endl;
		return {};
	}

	auto num_sets = (argc - arg) / 2;
	for (size_t i = 0; i < num_sets; i++)
	{
		int i0, i1;
		if (!interleaved)
		{
			i0 = arg + i;
			i1 = arg + i + num_sets;
		}
		else
		{
			i0 = arg + 2 * i;
			i1 = arg + 2 * i + 1;
		}
		data_sets.data.push_back({argv[i0], argv[i1]});
	}
	data_sets.valid = true;
	return data_sets;
}

int main(int argc, char* argv[])
{
	auto data_sets = get_data(argc, argv);
	if (!data_sets.valid)
		return -1;
	for (auto&& [f0, f1] : data_sets.data)
		process_pair(f0, f1, data_sets);
	return 0;
}