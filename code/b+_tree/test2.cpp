#include "b_tree.h"
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unistd.h>
#include <algorithm>
#include "util.h"
#include "utils.h"


void test_bulk_load(int hybrid_mode, char *index_name, char *key_path, int count, int has_size) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    BTree bt(hybrid_mode, true, index_name, true);
    LeafNodeIterm *data = new LeafNodeIterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
    bt.bulk_load(data, count);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << "bulk_load time:" << bulk_t / 1e9 << " sec" << std::endl;
    std::cout << "inner node size:" << bt.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << bt.get_file_size() << " bytes" << std::endl;
    return;
}

void test_lookup(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id,
                 int step, int need_sleep) {
    if (case_id > 1) throw std::invalid_argument("not support this query case...");
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    BTree bt(hybrid_mode, true, index_name, true);
    LeafNodeIterm *data = new LeafNodeIterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i] + 1;
    }
    bt.bulk_load(data, count);
    bt.sync_metanode();

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
        found = bt.lookup(search_keys[j], &c);
        tc += c;
        if (!found) {
            nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << 1e9 * s_count / batch_lookup_time << " ops" << std::endl;
    std::cout << tc / s_count << " block/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    delete[]search_keys;
}

void test_scan(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id,
               int r_size, int step) {
    if (case_id > 1) throw std::invalid_argument("not support this query case...");
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    BTree bt(hybrid_mode, true, index_name, true);
    LeafNodeIterm *data = new LeafNodeIterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
    bt.bulk_load(data, count);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << "bulk_load time:" << bulk_t / 1e9 << " sec" << std::endl;
    bt.sync_metanode();

    int c = 0;
    double tc = 0;
    int nc = 0;
    KeyType *search_keys = new KeyType[s_count];
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
    std::cout << "start to test... " << std::endl;

    std::cout << "start to record... " << std::endl;
    bool found = false;
    KeyType rs[r_size];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        found = bt.scan(rs, search_keys[j], r_size, &c);
        tc += c;
        if (!found) {
            nc += 1;
            std::cout << "not: " << j << std::endl;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << 1e9 * s_count / batch_lookup_time << " ops" << std::endl;
    std::cout << tc / s_count << " block/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
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

void test_insert(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int insert_count = 20000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            std::default_random_engine(seed));
    std::sort(keys, keys + (count-insert_count));
    BTree bt(hybrid_mode, true, index_name, true);
    LeafNodeIterm *data = new LeafNodeIterm[count-insert_count];
    for (int i = 0; i < count-insert_count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i] + 1;
    }
    std::cout << "start to bulk... " << std::endl;
    bt.bulk_load(data, count-insert_count);

    // generate a permutation for [0, insert_count-1]
    std::vector<int> v(insert_count);
    for (int i = 0; i < insert_count; i++)
        v[i] = i;
    std::vector<int> i_o(insert_count);
    for (int i = 0; i < insert_count; i++) {
        i_o[i] = get_next(v, insert_count-i);
    }
    std::cout << "start to insert... " << std::endl;
//    sleep(10);

    long long *latency = new long long[insert_count];
#if Profiling
    long long *sels = new long long[insert_count];
    long long *smls = new long long[insert_count];
    long long *inls = new long long[insert_count];
#endif
    long long sel;
    long long sml;
    long long inl;
    double total_smos = 0;
    double total_updated = 0;
    int update_c = 0;
    int smo_c = 0;
    auto bulk_s = std::chrono::high_resolution_clock::now();
    int _start = count-insert_count;


    for (int i = _start, j = 0; i < count; i++, j++) {
        if (j % 100000 == 0)
            std::cout << "---|" <<j/100000<<"|---"<<std::endl;
        update_c = 0, smo_c = 0;
#if Profiling
        auto i_start = std::chrono::high_resolution_clock::now();
#endif
	    bt.insert_key_entry(keys[_start + i_o[j]], keys[_start + i_o[j]] + 1, &sel, &sml, &inl, &smo_c, &update_c);

        auto i_end = std::chrono::high_resolution_clock::now();
        auto i_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                i_end - i_start)
                .count();

        latency[j] = i_b;
#if Profiling
        sels[j] = sel;
        smls[j] = sml;
        inls[j] = inl;
        total_smos += smo_c;
        total_updated += update_c;
#endif
    }
    auto bulk_e = std::chrono::high_resolution_clock::now();
    auto bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << insert_count * 1e9 * (1.0) / bulk_t << " ops" << std::endl;
    std::cout << "inner node size:" << bt.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << bt.get_file_size() << " bytes" << std::endl;
    std::cout << "smos:" << total_smos << std::endl;
    std::cout << "tuples:" << total_updated << std::endl;

#if Profiling
    std::ofstream outFile;
    outFile.open("in_latency.txt");
    for (int i = 0; i < insert_count; i++) {
        outFile << latency[i] << std::endl; //","<< sels[i]<< ","<< inls[i] << "," << smls[i] <<std::endl;
    }
    outFile.close();
#endif
    return;
}

void test_mixed(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int case_id, int total_op = 20000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            std::default_random_engine(seed));
    std::sort(keys, keys + count - total_op);
    LeafNodeIterm *data = new LeafNodeIterm[count - total_op];
    for (int i = 0; i < count - total_op; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i] + 1;
    }
    // generate a permutation for [0, insert_count-1]
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
    BTree bt(hybrid_mode, true, index_name, true);
    std::cout << "start to bulk... " << std::endl;
    bt.bulk_load(data, count - total_op);
    bt.sync_metanode();
    std::cout << "start to hybrid test... " << std::endl;
    bool found = false;
    int nc = 0;
    int c = 0;
    long long sel, sml, inl;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    for (int i =start, j = 0; i < (start + w_count) && j < s_count;) {
        for (int _i = 0; _i < step_w; _i++) {
            int smo_c, update_c;
            bt.insert_key_entry(keys[i], keys[i] + 1, &sel, &sml, &inl, &smo_c, &update_c);
            i += 1;
        }
        for (int _j = 0; _j < step_s; _j++) {
            found = bt.lookup(search_keys[j], &c);
            j += 1;
            if (!found) nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    double scount = total_op;
    std::cout << 1e9 * scount / bulk_lookup_time << " ops" << std::endl;
    std::cout << "not found:"<< nc << std::endl;
    std::cout << "inner node size:" << bt.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << bt.get_file_size() << " bytes" << std::endl;
    return;
}

void test_bulk_search_range(int hybrid_mode, char *index_name, char *key_path, int count, int has_size, int s_count,
                            int case_id, int r_size) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    BTree bt(hybrid_mode, true, index_name, true);
    LeafNodeIterm *data = new LeafNodeIterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i] + 1;
    }
    auto bulk_s = std::chrono::high_resolution_clock::now();
    bt.bulk_load(data, count);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << "bulk_load time:" << bulk_t / 1e9 << " sec" << std::endl;
    std::cout << "inner node size:" << bt.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << bt.get_file_size() << " bytes" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;


    bt.sync_metanode();

    int c = 0; double tc = 0; int nc = 0; KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count - r_size*2, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count - r_size*2, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }

    delete []keys;
    delete []data;

    bool found = false;
    long long *latency = new long long[s_count];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_start_time = std::chrono::high_resolution_clock::now();
#endif
        found = bt.lookup(search_keys[j], &c);
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_end_time = std::chrono::high_resolution_clock::now();
        long long i_batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                i_lookups_end_time - i_lookups_start_time).count();
        latency[j] = i_batch_lookup_time;
#endif
        tc += c;
        if (!found) {
            nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << 1e9 * s_count / batch_lookup_time << " ops - lookup" << std::endl;
    std::cout << tc / s_count << " block/lookup - lookup" << std::endl;
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
    KeyType rs[r_size];
    nc = 0;
    tc = 0;
    lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        found = bt.scan(rs, search_keys[j], r_size, &c);
        tc += c;
        if (!found) {
            nc += 1;
            std::cout << "not: " << j << std::endl;
        }
    }
    lookups_end_time = std::chrono::high_resolution_clock::now();
    batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << 1e9 * s_count / batch_lookup_time << " ops - scan" << std::endl;
    std::cout << tc / s_count << " block/lookup - scan" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    delete[]search_keys;
    std::cout << "\n\n\n\n" << std::endl;
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
        test_bulk_load(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                       count, has_size);
    } else if (op_type == "lookup") {
        test_lookup(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                    count, has_size, search_count, case_id, step, need_sleep);
    } else if (op_type == "scan") {
        test_scan(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()), count,
                  has_size, search_count, case_id, r_size, step);
    } else if (op_type == "insert") {
        test_insert(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                    count, has_size);
    } else if (op_type == "mix_workload") {
        test_mixed(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                   count, has_size, case_id);
    } else if (op_type == "bulk_search_range") {
        test_bulk_search_range(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                               count, has_size, search_count, case_id, r_size);
    }
    return 0;
}
