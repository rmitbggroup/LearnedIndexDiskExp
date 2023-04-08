/*
 * This example shows how to use pgm::DynamicPGMIndex, a std::map-like container supporting inserts and deletes.
 * Compile with:
 *   g++ updates.cpp -std=c++17 -I../include -o updates
 * Run with:
 *   ./updates
 */

#include <vector>
#include <cstdlib>
#include <iostream>

#include <random>
#include <fstream>
#include <chrono>
#include <unistd.h>
#include <algorithm>
#include "pgm/pgm_index_dynamic.hpp"
#include "pgm/util.h"
#include "pgm/utils.h"

#define KeyType uint64_t
#define ValueType uint64_t

#define Memory 0

void test_bulk_load(int hybrid_mode, char *index_name, char *key_path, int count, int has_size) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    std::vector<std::pair<KeyType, ValueType>> data(count);
    for (int i = 0; i < count; i++) {
        data[i].first = keys[i];
        data[i].second = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if Memory
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(data.begin(), data.end());
#else
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(true, hybrid_mode==0, data.begin(), data.end());
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;
    std::cout << "inner node size:" << dynamic_pgm.report_main_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << dynamic_pgm.report_disk_file_size() << " bytes" << std::endl;
    return;
}

void test_lookup(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int step, int need_sleep) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    std::vector<std::pair<KeyType, ValueType>> data(count);
    for (int i = 0; i < count; i++) {
        data[i].first = keys[i];
        data[i].second = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if Memory
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(data.begin(), data.end());
#else
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(true, hybrid_mode == 0, data.begin(), data.end());
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;

    int c = 0;
    double tc = 0;
    int nc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
    std::cout << "start to test... " << std::endl;
    if (need_sleep == 1) sleep(7);
    std::cout << "start to record... " << std::endl;
    bool found = false;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
#if Memory
        dynamic_pgm.find(search_keys[j]);
#else
        c = 0;
        int ic = 0;
        int lc = 0;
        found = dynamic_pgm.find_on_disk(search_keys[j], &c, &ic, &lc);
        tc += c;
        if (!found) {
            nc+=1;
            std::cout.precision(17);
            std::cout << "not: " << search_keys[j] << std::endl;
            return;
        }
#endif
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops" << std::endl;
    std::cout << tc/s_count<< " block/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
}

void test_scan(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int r_size , int step) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    std::vector<std::pair<KeyType, ValueType>> data(count);
    for (int i = 0; i < count; i++) {
        data[i].first = keys[i];
        data[i].second = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if Memory
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(data.begin(), data.end());
#else
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(true, hybrid_mode == 0, data.begin(), data.end());
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;

    int c = 0;
    double tc = 0;
    int nc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
//    delete[] keys;
    std::cout << "start to test... " << std::endl;
//    sleep(10);
    std::cout << "start to record... " << std::endl;
    bool found = false;
    KeyType rs[r_size];
    pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk *a = new pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk [102];
    pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk *b = new pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk [102];
    int fsize = 0;
    KeyType hi = keys[count-1];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
#if Memory
#else
        c = 0;
         dynamic_pgm.range_on_disk(search_keys[j], hi, 100, &c,
                                   a, b, &fsize);
        tc += c;
#endif
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops" << std::endl;
    std::cout << tc/s_count<< " block/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    delete []keys;
}

int get_next(std::vector<int>& v, int _seed) {
    int n = v.size();

    srand(_seed);

    // Make sure the number is within
    // the index range
    int index = rand() % n;

    // Get random number from the vector
    int num = v[index];

    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();
    return num;
}

void test_insert(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int insert_count = 10000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    unsigned seed = 0;
    std::shuffle(keys, keys + count,
            std::default_random_engine(seed));
    std::sort(keys, keys+count-insert_count);
    std::vector<std::pair<KeyType, ValueType>> data(count-insert_count);
    for (int i = 0; i < count-insert_count; i++) {
        data[i].first = keys[i];
        data[i].second = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if Memory
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(data.begin(), data.end());
#else
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(true, hybrid_mode == 0, data.begin(), data.end());
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;
    std::cout << "start to insert... " << std::endl;
    std::vector<int> v(insert_count);
    for (int i = 0; i < insert_count; i++)
        v[i] = i;
    std::vector<int> i_o(insert_count);
    for (int i = 0; i < insert_count; i++) {
        i_o[i] = get_next(v, insert_count-i);
    }

#if Profiling
    long long *latency = new long long[insert_count];

    long long *sels = new long long[insert_count];
    long long *inls = new long long[insert_count];
    long long *smls = new long long[insert_count];
#endif
    long long total_upper_latency = 0;
    long long total_leaf_latency = 0;

    int _start = count - insert_count;
    long long sel, inl, sml;
    long long total_updated_tuples = 0;
    double total_smo = 0;
    int smo_c = 0;
    int update_t = 0;
    bulk_s = std::chrono::high_resolution_clock::now();
    for (int i = _start, j = 0; i < count; i++, j++) {
#if Memory
        dynamic_pgm.insert_or_assign(keys[i], keys[i]+1);
#else
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_start_time = std::chrono::high_resolution_clock::now();
#endif
        dynamic_pgm.insert_on_disk(keys[_start + i_o[j]], keys[_start + i_o[j]] + 1, &sel, &inl, &sml, &smo_c, &update_t);
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_end_time = std::chrono::high_resolution_clock::now();
        long long i_batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(i_lookups_end_time - i_lookups_start_time).count();
        latency[j] = i_batch_lookup_time;
        sels[j] = sel;
        smls[j] = sml;
        inls[j] = inl;
        total_updated_tuples += update_t;
        total_smo += smo_c;
#endif
#endif

    }
    bulk_e = std::chrono::high_resolution_clock::now();
    bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<< double(insert_count)*1e9/bulk_t<< " ops"<< std::endl;
    std::cout << "inner node size:" << dynamic_pgm.report_main_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << dynamic_pgm.report_disk_file_size() << " bytes" << std::endl;
#if Profiling
    std::cout << "smos:" << total_smo << std::endl;
    std::cout << "updates:" << total_updated_tuples << std::endl;
    std::ofstream outFile;
    outFile.open("in_latency.txt");
    for (int i = 0; i < insert_count; i++) {
        outFile << latency[i] << std::endl;//","<< sels[i] << ","<< inls[i] << ","<< smls[i] << std::endl;
    }
    outFile.close();
#endif
    return;
}

void test_mixed(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int case_id, int total_op=10000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    unsigned seed = 0;
    std::shuffle(keys, keys + count,
            std::default_random_engine(seed));
    std::sort(keys, keys+count-total_op);
    std::vector<std::pair<KeyType, ValueType>> data(count-total_op);
    for (int i = 0; i < count-total_op; i++) {
        data[i].first = keys[i];
        data[i].second = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if Memory
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(data.begin(), data.end());
#else
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(true, hybrid_mode == 0, data.begin(), data.end());
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;


//    int total_op = count/2;
    int start = count - total_op;
    std::vector<int> v(total_op);
    for (int i = 0; i < total_op; i++)
        v[i] = i;
    std::vector<KeyType> i_o(total_op);
    for (int i = 0; i < total_op; i++) {
        i_o[i] = keys[start + get_next(v, total_op-i)];
    }
    for (int i = 0; i < total_op; i++) {
        keys[start +i] = i_o[i];
    }
    int s_count;
    int step_s;
    int step_w;
    if (case_id == 1) { // read_heavey
        s_count = 18 * total_op / 20;
        step_s = 18;
        step_w = 2;
    } else if (case_id == 2) { // write_heavey
        s_count = 2* total_op / 20;
        step_s = 2;
        step_w = 18;
    } else { // balence
        s_count = total_op / 2;
        step_s = 10;
        step_w = 10;
    }
    KeyType *search_keys = new KeyType[s_count];
    std::mt19937_64 gen(19937);
    for (int i = 0; i < s_count;) {
        std::uniform_int_distribution<int> dis(0, start + (i/step_s + 1)*step_w);
        for (int _i = 0; _i < step_s; _i++) {
            search_keys[i] = keys[dis(gen)];
            i++;
        }
    }

    int w_count = total_op - s_count;
    std::cout << "start to hybrid test... " << std::endl;
    bool found = false;
    int nc = 0;
    int c = 0;
    int _start = count - total_op;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    for (int i = _start, j = 0; i < (_start + w_count) && j < s_count; ) {
        if (i % 100000 == 0){
            std::cout << i/100000 << std::endl;
        }

        for (int _i = 0; _i < step_w; _i++) {
#if Memory
            dynamic_pgm.insert_or_assign(keys[i], keys[i]+1);
#else
            long long sel, ins, smol;
            int smo_c, update_c;
            dynamic_pgm.insert_on_disk(keys[i],keys[i]+1, &sel, &ins, &smol, &smo_c, &update_c);
#endif
            i += 1;
        }
        for (int _j = 0; _j < step_s; _j++) {
#if Memory
            dynamic_pgm.find(search_keys[j]);
#else
            int ic = 0;
            int lc = 0;
            found = dynamic_pgm.find_on_disk(search_keys[j],&c, &ic, &lc);
            if (!found) {
                nc+=1;
            }
            j += 1;
#endif
        }
    }
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    double scount = total_op;
    std::cout << 1e9*scount/bulk_lookup_time<< " ops" << std::endl;
    std::cout << nc << std::endl;
    std::cout << "inner node size:" << dynamic_pgm.report_main_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << dynamic_pgm.report_disk_file_size() << " bytes" << std::endl;
    return;
}


void test_bulk_search(int hybrid_mode, char *key_path, int count, int has_size, int s_count, int case_id, int r_size) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    std::vector<std::pair<KeyType, ValueType>> data(count);
    for (int i = 0; i < count; i++) {
        data[i].first = keys[i];
        data[i].second = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if Memory
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(data.begin(), data.end());
#else
    pgm::DynamicPGMIndex<KeyType, ValueType> dynamic_pgm(true, hybrid_mode==0, data.begin(), data.end());
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;
    std::cout << "inner node size:" << dynamic_pgm.report_main_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << dynamic_pgm.report_disk_file_size() << " bytes" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;

    int c = 0;
    double tc = 0;
    int nc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count - r_size*2, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count - r_size*2, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    int ic = 0;
    int lc = 0;
    KeyType hi = keys[count-1];
    delete[] keys;
    bool found = false;
    double tlc = 0;
    double tic = 0;


    long long *latency = new long long[s_count];
//    int *blocks = new int[s_count];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
#if Memory
        dynamic_pgm.find(search_keys[j]);
#else
        c = 0;
        int lc = 0;
        int ic = 0;
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_start_time = std::chrono::high_resolution_clock::now();
#endif
        found = dynamic_pgm.find_on_disk(search_keys[j], &c, &ic, &lc);
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_end_time = std::chrono::high_resolution_clock::now();
        long long i_batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(i_lookups_end_time - i_lookups_start_time).count();
        latency[j] = i_batch_lookup_time;
#endif
//        std::cout << lc << std::endl;
        tlc += lc;
        tic += ic;
        tc += c;
//        blocks[j] = c;
        if (!found) {
            nc+=1;
        }
#endif
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops - lookup" << std::endl;
    std::cout << tc/s_count<< " block/lookup - lookup" << std::endl;
    std::cout << tic/s_count<< " inner block/lookup - lookup" << std::endl;
    std::cout << tlc/s_count<< " inner level/lookup - lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    std::cout << "\n\n\n\n" << std::endl;

#if Profiling
    std::ofstream outFile;
    outFile.open("s_latency.txt");
    for (int i = 0; i < s_count; i++) {
        outFile << latency[i] << std::endl;
    }
    outFile.close();
#endif
    pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk *a = new pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk [102];
    pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk *b = new pgm::DynamicPGMIndex<KeyType, ValueType>::ItemOnDisk [102];
    int fsize = 0;

    tc = 0;
    lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
#if Memory
#else
        c = 0;
        dynamic_pgm.range_on_disk(search_keys[j], hi, r_size, &c,
                                  a, b, &fsize);
        tc += c;
#endif
    }
    lookups_end_time = std::chrono::high_resolution_clock::now();
    batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops" << std::endl;
    std::cout << tc/s_count<< " block/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    return;
}

int main(int argc, char *argv[]) {
    auto flags = parse_flags(argc, argv);
    std::string key_file_path = get_required(flags, "keys_file");
    std::string op_type = get_required(flags, "op_type");
    std::string index_name = get_required(flags, "index_file");
    std::string utility_file = index_name + "_utility";
    int count = stoi(get_required(flags, "total_count"));
    int search_count = stoi(get_with_default(flags, "search_count", "200000"));
    int has_size = stoi(get_required(flags, "has_size"));
    int case_id = stoi(get_with_default(flags, "case_id", "1"));
    int step = stoi(get_with_default(flags, "step", "999"));
    int r_size = stoi(get_with_default(flags, "r_size", "100"));
    int need_sleep = stoi(get_with_default(flags, "need_sleep", "0"));
    int memory_type = stoi(get_with_default(flags, "memory_type", "0")); // all disk
    if (op_type == "bulk") {
        test_bulk_load(memory_type, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size);
    } else if (op_type == "lookup") {
        test_lookup(memory_type, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, step, need_sleep);
    } else if (op_type == "scan") {
        test_scan(memory_type, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, r_size, step);
    } else if (op_type == "insert") {
        test_insert(memory_type, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size);
    } else if (op_type == "mix_workload") {
        test_mixed(memory_type, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, case_id);
    } else if (op_type == "bulk_search_range") {
        test_bulk_search(memory_type, const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, r_size);
    }
    return 0;
}
