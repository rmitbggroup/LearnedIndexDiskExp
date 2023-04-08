#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <map>
#include "storage_management.h"
#include <unistd.h>
#include <string>
#include "util.h"
#include "utils.h"

using namespace std;

void test_lipp_bulk(char *index_name, char *key_path, int count, int has_size) {
  StorageManager<KeyType, ValueType> sm;
  sm.init(index_name, true);
  std::ifstream fin(key_path, std::ios::binary);
  KeyType *keys = new KeyType[count];
  if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
   }
  fin.read((char*)(keys), sizeof(KeyType)*count);
  fin.close();
  std::sort(keys, keys+count);

  ValueType *values = new ValueType[count];
  int rp = 0;
  for (int i = 0; i < count; i++) {
    values[i] = keys[i] + 1;
  }
  for (int i = rp; i < count-1; i++) {
      keys[i] = keys[i+1];
  }
  count -= 1;
  for (int i = 0; i < count; i++) {
    values[i] = keys[i] + 1;
  }
  cout << "start to build... " << endl;
  std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
  sm.bulk_load_disk(keys, values, count);
  sm.sys_metablock(false);
  std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
  long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
  std::cout << "bulk load time: " << bulk_lookup_time/1e9 << std::endl;
  std::cout << "file size:" << sm.get_file_size() << std::endl;
  return;
}

void test_lipp_lookup(char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int step, int need_sleep) {
    if (case_id > 2) {
        throw std::invalid_argument("not support this query case...");
    }
    StorageManager<KeyType, ValueType> sm;
    sm.init(index_name, false);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);

    bool found;
    int bc;
    double sc = 0;
    ValueType v;
    bc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;

    cout << "start to test... " << endl;
    if (need_sleep == 1) sleep(7);
    cout << "start to record... " << endl;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    int nc = 0;
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        int lc,ic;
        lc = ic = 0;
        found = sm.searchv2(search_keys[j], &v, &bc, &ic, &lc);
        sc += bc;
        if (!found) {
            nc+=1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    cout << 1e9*s_count/batch_lookup_time<< " ops" << endl;
    cout << sc/s_count<< " block/lookup" << endl;
    cout << "not found:" << nc << endl;
}

void test_lipp_scan(char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int step, int r_size) {
    if (case_id > 2) {
        throw std::invalid_argument("not support this query case...");
    }
    StorageManager<KeyType, ValueType> index;
    index.init(index_name, false);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);

    bool found;
    int bc;
    double sc = 0;
    KeyType* search_keys = new KeyType[s_count];
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
    ValueType v;
    KeyType rs[r_size];
    cout << "start to test... " << endl;

    cout << "start to scan... " << endl;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    int nc = 0;
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        int ic = 0;
        int lc = 0;
        index.range_query_len(rs, search_keys[j], r_size, &bc, &ic, &lc);
        sc += bc;
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    cout << 1e9*s_count/batch_lookup_time<< " ops" << endl;
    cout << sc/s_count<< " block/lookup" << endl;
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

void test_lipp_insert(char *index_name, char *key_path, int count, int has_size, int insert_count=10000000) {
    StorageManager<KeyType, ValueType> index;
    index.init(index_name, true);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            default_random_engine(seed));
    std::sort(keys, keys+count-insert_count);
    ValueType *values = new ValueType[count-insert_count];
    for (int i = 0; i < count-insert_count; i++) {
        values[i] = keys[i] + 1;
    }
    int bc = 0;
    int nc = 0;
    double sc = 0;
    bool found;
    ValueType v;
    double t = 0;

    cout << "start to bulk... " << endl;
    index.bulk_load_disk(keys, values, count-insert_count);

    std::vector<int> _v(insert_count);
    for (int i = 0; i < insert_count; i++)
        _v[i] = i;
    std::vector<int> i_o(insert_count);
    for (int i = 0; i < insert_count; i++) {
        i_o[i] = get_next(_v, insert_count-i);
    }
    cout << "start to insert... " << endl;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    int _start = count-insert_count;

    long long *latency = new long long[insert_count];
#if Profiling
    long long *sels = new long long[insert_count];
    long long *inls = new long long[insert_count];
    long long *smls = new long long[insert_count];
    long long *mals = new long long[insert_count];
#endif
    double tsmo = 0;
    long long total_update_tuples = 0;

//    long long *otls = new long long[insert_count];
    for (int i = _start, j = 0; i < count; i++,j++) {
        if (j%100000 == 0) {
            cout << j/100000 << endl;
        }
        //
        // if (i == 801) {
        //     printf("xxx");
        // }
        long long se_l, in_l, sm_l, ma_l;
        int smo = 0;
        int updated_tuples = 0;
#if Profiling
        auto i_start = std::chrono::high_resolution_clock::now();
#endif
        index.insertv2(keys[_start + i_o[j]], keys[_start + i_o[j]] + 1, &se_l,
                       &in_l, &sm_l, &ma_l, &smo, &updated_tuples);

        auto i_end = std::chrono::high_resolution_clock::now();
        auto i_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                i_end - i_start)
                .count();
        latency[j] = i_b;
#if Profiling
        sels[j] = se_l;
        inls[j] = in_l;
        smls[j] = sm_l;
        mals[j] = ma_l;
        tsmo += smo;
        total_update_tuples += updated_tuples;
#endif
//        otls[j] = i_b - se_l - in_l - sm_l - ma_l;
    }
    index.sys_metablock(false);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    cout << double(insert_count)*1e9/bulk_lookup_time<< " ops" << endl;
    cout << tsmo << " total smos" << endl;
    cout << total_update_tuples << " updated_tuples" << endl;
    std::cout << "file size:" << index.get_file_size() << std::endl;
#if Profiling
    std::ofstream outFile;
    outFile.open("i_latency.txt");
    for (int i = 0; i < insert_count; i++) {
        outFile << latency[i] << "," << sels[i] << "," << inls[i] << "," << smls[i] << "," << mals[i] << std::endl;
    }
    outFile.close();
#endif
    return;
}

void test_lipp_mix(char *index_name, char *key_path, int count, int has_size, int case_id, int total_op = 10000000) {
    StorageManager<KeyType, ValueType> index;
    index.init(index_name, true);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            default_random_engine(seed));
    std::sort(keys, keys+count - total_op);
    ValueType *values = new ValueType[count];
    for (int i = 0; i < count; i++) {
        values[i] = keys[i] + 1;
    }
//    int total_op = count/2;
    int start = count - total_op;
    std::vector<int> _v(total_op);
    for (int i = 0; i < total_op; i++)
        _v[i] = i;
    std::vector<KeyType> i_o(total_op);
    for (int i = 0; i < total_op; i++) {
        i_o[i] = keys[start + get_next(_v, total_op-i)];
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
    int bc = 0;
    bool found;
    ValueType v;
    double t = 0;
    int nc = 0;
    cout << "start to bulk... " << endl;
    index.bulk_load_disk(keys, values, count - total_op);
    cout << "start to hybrid test... " << endl;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    long long se_l, in_l, sm_l, ma_l;
    int smo;
    for (int i = start, j = 0; i < (start + w_count) && j < s_count; ) {
        for (int _i = 0; _i < step_w; _i++) {
            int updated = 0;
            index.insertv2(keys[i], values[i], &se_l, &in_l, &sm_l, &ma_l, &smo, &updated);
            i += 1;
        }
        for (int _j = 0; _j < step_s; _j++) {
            int ic = 0; int lc  =0 ;
            found = index.searchv2(search_keys[j], &v, &bc, &ic, &lc);
            j += 1;
            if (!found) nc+=1;
        }
    }
    index.sys_metablock(false);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    double scount = total_op;
    cout << 1e9*scount/bulk_lookup_time<< " ops" << endl;
    std::cout << "file size:" << index.get_file_size() << std::endl;
    cout << nc << endl;
    return;
}

void test_bulk_search(char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int r_size) {
    StorageManager<KeyType, ValueType> sm;
    sm.init(index_name, true);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);

    ValueType *values = new ValueType[count];
    for (int i = 0; i < count; i++) {
        values[i] = keys[i] + 1;
    }
    cout << "start to build... " << endl;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    sm.bulk_load_disk(keys, values, count);
    sm.sys_metablock(false);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    std::cout << "bulk load time: " << bulk_lookup_time/1e9 << std::endl;
    std::cout << "file size:" << sm.get_file_size() << std::endl;
    std::cout << "\n\n\n\n" << std::endl;

    bool found;
    int bc;
    double sc = 0;
    ValueType v;
    bc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count - r_size*2, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count - r_size*2, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
//    delete[] values;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    int nc = 0;
    double tlc = 0;
    long long *latency = new long long[s_count];
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        int lc = 0;
        int ic = 0;
//#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_start_time = std::chrono::high_resolution_clock::now();
//#endif
        found = sm.searchv2(search_keys[j], &v, &bc, &ic, &lc);
        std::chrono::high_resolution_clock::time_point i_lookups_end_time = std::chrono::high_resolution_clock::now();
//#if Profiling
        long long i_batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                i_lookups_end_time - i_lookups_start_time).count();
        latency[j] = i_batch_lookup_time;
//#endif
        sc += bc;
        tlc += lc;
        if (!found) {
            nc+=1;
        }
    }
//#if Profiling
    std::ofstream outFile;
    outFile.open("s_latency.txt");
    for (int i = 0; i < s_count; i++) {
        outFile << latency[i] << std::endl;
    }
    outFile.close();
//#endif
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    cout << 1e9*s_count/batch_lookup_time<< " ops" << endl;
    cout << sc/s_count<< " block/lookup" << endl;
    cout << tlc/s_count<< " level/lookup" << endl;
    cout << "not found:" << nc << endl;
    std::cout << "\n\n\n\n" << std::endl;
    return;
    KeyType rs[r_size];
    sc = 0;
    lookups_start_time = std::chrono::high_resolution_clock::now();
    nc = 0;
    double tic = 0; //item count
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        int ic = 0;
        int lc = 0;
        sm.range_query_len(rs, search_keys[j], r_size, &bc, &ic, &lc);
        tlc += lc;
        tic += ic;
        sc += bc;
    }
    lookups_end_time = std::chrono::high_resolution_clock::now();
    batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    cout << 1e9*s_count/batch_lookup_time<< " ops" << endl;
    cout << sc/s_count<< " block/lookup" << endl;
    cout << tlc/s_count<< " node/lookup" << endl;
    cout << tic/s_count<< " item/lookup" << endl;
    std::cout << "\n\n\n\n" << std::endl;
}


int main(int argc, char *argv[]) {
        auto flags = parse_flags(argc, argv);
        std::string key_file_path = get_required(flags, "keys_file");
        std::string op_type = get_required(flags, "op_type");
        std::string index_name = get_required(flags, "index_file");
        int count = stoi(get_required(flags, "total_count"));
        int search_count = stoi(get_with_default(flags, "search_count", "200000"));
        int has_size = stoi(get_required(flags, "has_size"));
        int case_id = stoi(get_with_default(flags, "case_id", "1"));
        int step = stoi(get_with_default(flags, "step", "999"));
        int r_size = stoi(get_with_default(flags, "r_size", "100"));
        int need_sleep = stoi(get_with_default(flags, "need_sleep", "0"));
        if (op_type == "bulk") {
                test_lipp_bulk(const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size);
        } else if (op_type == "lookup") {
                test_lipp_lookup(const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, step, need_sleep);
        } else if (op_type == "scan") {
                test_lipp_scan(const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, step, r_size);
        } else if (op_type == "insert") {
                test_lipp_insert(const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size);
        } else if (op_type == "mix_workload") {
                test_lipp_mix(const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, case_id);
        } else if (op_type == "bulk_search_range") {
        	test_bulk_search(const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, r_size);
    	}	
        return 0;
}