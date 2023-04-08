#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include "fiting_tree_memory.h"
#include <unistd.h>
#include "util.h"
#include "utils.h"

void test_bulk_load(int memory_type, int eb, char *index_name, char *key_path, int count, int has_size) {
    std::cout << eb << std::endl;
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
   if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
   }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    FITingTree ft(eb, index_name, true, memory_type);
    Iterm *data =new Iterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i];
    }
    #if PGMSegment
    std::vector<KeyType> data2(keys, keys+count);
    #endif
    auto bulk_s = std::chrono::high_resolution_clock::now();
    #if PGMSegment
    ft.bulk_load_pgm(data, count, data2.begin(), data2.end(), eb);
    #else
    ft.bulk_load(data, count);
    #endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            bulk_e - bulk_s)
                            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;

    std::cout << "inner node size:" << ft.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << ft.get_file_size() << " bytes" << std::endl;
    delete []keys;
    delete []data;
    return;
}

void test_lookup(int memory_type, int eb, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int step, int need_sleep) {
    std::cout << eb << std::endl;
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    FITingTree ft(eb, index_name, true, memory_type);
    Iterm *data =new Iterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i];
    }
#if PGMSegment
    std::vector<KeyType> data2(keys, keys+count);
#endif
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if PGMSegment
    ft.bulk_load_pgm(data, count, data2.begin(), data2.end(), eb);
#else
    ft.bulk_load(data, count);
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;


//    std::ifstream fin(key_path, std::ios::binary);
//    KeyType *keys = new KeyType[count];
//    if (has_size == 1) {
//        uint64_t size;
//        fin.read((char*)(&size), sizeof(uint64_t));
//    }
//    fin.read((char*)(keys), sizeof(KeyType)*count);
//    fin.close();
//    std::sort(keys, keys+count);
//
//    FITingTree ft(eb, index_name, false);

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
        c = 0;
        found = ft.lookup(search_keys[j], &c);
        tc += c;
        if (!found) {
            nc+=1;
            // std::cout << j << std::endl;
            // break;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops" << std::endl;
    std::cout << tc/s_count<< " block/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    delete []search_keys;
}

void test_scan(int memory_type, int eb, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id, int step, int r_size) {

    std::cout << eb << std::endl;
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    FITingTree ft(eb, index_name, true, memory_type);
    Iterm *data =new Iterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i];
    }
#if PGMSegment
    std::vector<KeyType> data2(keys, keys+count);
#endif
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if PGMSegment
    ft.bulk_load_pgm(data, count, data2.begin(), data2.end(), eb);
#else
    ft.bulk_load(data, count);
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

    std::cout << "start to record... " << std::endl;
    bool found = false;
    KeyType rs[r_size];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        ft.scan(search_keys[j], rs, &c, r_size);
        tc += c;
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops" << std::endl;
    std::cout << tc/s_count<< " block/lookup" << std::endl;
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

void test_insert(int memory_type, int eb, char *index_name, char *key_path, int count, int has_size, int insert_count = 20000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys + count -insert_count);
    FITingTree ft(eb, index_name, true, memory_type);
    Iterm *data =new Iterm[count -insert_count];
    for (int i = 0; i < count -insert_count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i]+1;
    }
    std::cout << "start to bulk... " << std::endl;

    #if PGMSegment
    std::vector<KeyType> data2(keys, keys+ count -insert_count);
    ft.bulk_load_pgm(data, count -insert_count, data2.begin(), data2.end(), eb);
    #endif
    // ft.bulk_load_pgm(data, count/2);
    int c = 0;
    // generate a permutation for [0, insert_count-1]
    std::vector<int> v(insert_count);
    for (int i = 0; i < insert_count; i++)
        v[i] = i;
    std::vector<int> i_o(insert_count);
    for (int i = 0; i < insert_count; i++) {
        i_o[i] = get_next(v, insert_count-i);
    }

    std::cout << "start to insert... " << std::endl;
    // sleep(10);
    int _start = count-insert_count;
    auto bulk_s = std::chrono::high_resolution_clock::now();
    for (int i = count -insert_count, j=0; i < count; i++, j++) {
        if ((i - (count -insert_count))%100000 == 0) std::cout << (i - (count -insert_count))/100000 << std::endl;
//        ft.insert_key_entry_f(keys[i], keys[i]+1);
        ft.insert_key_entry_f(keys[_start + i_o[j]], keys[_start + i_o[j]] + 1);
    }
    auto bulk_e = std::chrono::high_resolution_clock::now();
    auto bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            bulk_e - bulk_s)
                            .count();
    std::cout<< double(insert_count)*1e9/bulk_t<< " ops"<< std::endl;
    std::cout << "inner node size:" << ft.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << ft.get_file_size() << " bytes" << std::endl;
    delete []data;
    delete []keys;
    return;
}

void test_mixed(int memory_type, int eb, char *index_name, char *key_path, int count, int has_size, int case_id, int total_op = 20000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count - total_op);
    Iterm *data =new Iterm[count - total_op];
    for (int i = 0; i < count - total_op; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i]+1;
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
        std::uniform_int_distribution<int> dis(0, start + (i / step_s + 1) * step_w);
        for (int _i = 0; _i < step_s; _i++) {
            search_keys[i] = keys[dis(gen)];
            i++;
        }
    }
    int w_count = total_op - s_count;
    FITingTree ft(eb, index_name, true, memory_type);
    std::cout << "start to bulk... " << std::endl;
    // ft.bulk_load(data, count/2);
    std::vector<KeyType> data2(keys, keys+count - total_op);
    ft.bulk_load_pgm(data, count - total_op, data2.begin(), data2.end(), eb);
    std::cout << "start to hybrid test... " << std::endl;
    bool found = false;
    int nc = 0;
    int c = 0;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    for (int i = count - total_op, j = 0; i < (count - total_op + w_count) && j < s_count; ) {
        for (int _i = 0; _i < step_w; _i++) {
            ft.insert_key_entry_f(keys[i], keys[i]+1);
            i += 1;
        }
        for (int _j = 0; _j < step_s; _j++) {
            found = ft.lookup(search_keys[j], &c);
            j += 1;
            if (!found){
                nc+=1;
            }
        }
    }
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    double scount = total_op;
    std::cout << 1e9*scount/bulk_lookup_time<< " ops" << std::endl;
    std::cout << nc << std::endl;
    std::cout << "inner node size:" << ft.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << ft.get_file_size() << " bytes" << std::endl;
    delete []keys;
    delete []data;
    delete []search_keys;
    return;
}

void test_bulk_search(int memory_type, int eb, char *index_name, char *key_path, int count, int has_size
        , int s_count, int case_id, int r_size) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char*)(&size), sizeof(uint64_t));
    }
    fin.read((char*)(keys), sizeof(KeyType)*count);
    fin.close();
    std::sort(keys, keys+count);
    FITingTree ft(eb, index_name, true, memory_type);
    Iterm *data =new Iterm[count];
    for (int i = 0; i < count; i++) {
        data[i].key = keys[i];
        data[i].value = keys[i];
    }
#if PGMSegment
    std::vector<KeyType> data2(keys, keys+count);
#endif
    auto bulk_s = std::chrono::high_resolution_clock::now();
#if PGMSegment
    ft.bulk_load_pgm(data, count, data2.begin(), data2.end(), eb);
#else
    ft.bulk_load(data, count);
#endif
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;

    std::cout << "inner node size:" << ft.get_inner_size() << " bytes" << std::endl;
    std::cout << "file size:" << ft.get_file_size() << " bytes" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    delete []data;

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
    bool found = false;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        found = ft.lookup(search_keys[j], &c);
        tc += c;
        if (!found) {
            nc+=1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops - lookup" << std::endl;
    std::cout << tc/s_count<< " block/lookup - scan" << std::endl;
    std::cout << "not found:" << nc << std::endl;
    std::cout << "\n\n\n\n" << std::endl;

    found = false;
    KeyType rs[r_size];
    tc = 0;
    lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        ft.scan(search_keys[j], rs, &c, r_size);
        tc += c;
    }
    lookups_end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration_cast<std::chrono::nanoseconds>(lookups_end_time - lookups_start_time).count();
    std::cout << 1e9*s_count/batch_lookup_time<< " ops - scan" << std::endl;
    std::cout << tc/s_count<< " block/lookup - scan" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    delete []search_keys;

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
        int eb = stoi(get_with_default(flags, "error", "64"));
        int memory_type = stoi(get_with_default(flags, "memory_type", "0")); // all disk
        if (op_type == "bulk") {
                test_bulk_load(memory_type, eb, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size);
        } else if (op_type == "lookup") {
                test_lookup(memory_type, eb, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, step, need_sleep);
        } else if (op_type == "scan") {
                test_scan(memory_type, eb, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, step, r_size);
        } else if (op_type == "insert") {
                test_insert(memory_type, eb, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size);
        } else if (op_type == "mix_workload") {
                test_mixed(memory_type, eb, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, case_id);
        } else if (op_type == "bulk_search_range") {
                test_bulk_search(memory_type, eb, const_cast<char*>(index_name.c_str()), const_cast<char*>(key_file_path.c_str()), count, has_size, search_count, case_id, r_size);
        }

        return 0;
}