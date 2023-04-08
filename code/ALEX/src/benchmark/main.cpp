#include "../core/alex.h"

#include <iomanip>

#include "flags.h"
#include "utils.h"
#include <unistd.h>

// Modify these if running your own workload
#define KeyType uint64_t
#define ValueType uint64_t

void test_bulk_load(int hybrid_mode, char *index_name, char *key_path, char *inner_file, int count, int has_size) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    alex::Alex<KeyType, ValueType> index(hybrid_mode, true, index_name, inner_file);
    auto values = new std::pair<KeyType, ValueType>[count];
    for (int i = 0; i < count; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<ValueType>(keys[i] + 1);
    }
    std::cout << "start bulk_load..." << std::endl;
    auto bulk_s = std::chrono::high_resolution_clock::now();
    index.bulk_load(values, count);
    index.sync_metanode(true);
    index.sync_metanode(false);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << "bulk_load time:" << bulk_t / 1e9 << " sec" << std::endl;
    std::cout << "inner node size:" << index.get_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << index.get_file_size() << " bytes" << std::endl;
    delete[]keys;
    delete[]values;
    return;
}

void
test_lookup(int hybrid_mode, char *index_name, char *key_path, char *inner_file, int count, int has_size, int s_count,
            int case_id, int step, int need_sleep) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    alex::Alex<KeyType, ValueType> index(hybrid_mode, true, index_name, inner_file);
    auto values = new std::pair<KeyType, ValueType>[count];
    for (int i = 0; i < count; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<ValueType>(keys[i]+1);
    }
    std::cout<<"start bulk_load..." << std::endl;
    auto bulk_s = std::chrono::high_resolution_clock::now();
    index.bulk_load(values, count);
    index.sync_metanode(true);
    index.sync_metanode(false);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout<<"bulk_load time:"<< bulk_t/1e9 << " sec"<< std::endl;


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
    ValueType v;
    int c = 0;
    double tc = 0;
    int nc = 0;
    double tl = 0;
    int l = 0;
    double tci = 0;
    int ci = 0;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        l = 0;
        ci = 0;
        found = index.get_leaf_disk(search_keys[j], &v, &c, &l, &ci);
        tc += c;
        tl += l;
        tci += ci;
        if (!found) {
            nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << 1e9 * s_count / batch_lookup_time << " ops" << std::endl;
    std::cout << tc / s_count << " block/lookup" << std::endl;
    std::cout << tl / s_count << " levels/lookup" << std::endl;
    std::cout << tci / s_count << " inners/lookup" << std::endl;
    std::cout << "not found:" << nc << std::endl;

}

void
test_scan(int hybrid_mode, char *index_name, char *key_path, char *inner_file, int count, int has_size, int s_count,
          int case_id, int r_size, int step) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    alex::Alex<KeyType, ValueType> index(hybrid_mode, true, index_name, inner_file);
     auto values = new std::pair<KeyType, ValueType>[count];
     for (int i = 0; i < count; i++) {
         values[i].first = keys[i];
         values[i].second = static_cast<ValueType>(keys[i]+1);
     }
     std::cout<<"start bulk_load..." << std::endl;
     auto bulk_s = std::chrono::high_resolution_clock::now();
     index.bulk_load(values, count);
     index.sync_metanode(true);
     index.sync_metanode(false);
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
        search_keys = get_search_keys(keys, count - s_count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count - s_count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }

    std::cout << "start to test... " << std::endl;

    std::cout << "start to record... " << std::endl;
    std::cout << r_size << std::endl;
    bool found = false;
    KeyType rs[r_size];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        index.scan(rs, search_keys[j], &c, r_size);
        tc += c;
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << 1e9 * s_count / batch_lookup_time << " ops" << std::endl;
    std::cout << tc / s_count << " block/lookup" << std::endl;
    delete[] keys;
}

int get_next(std::vector<int> &v, int _seed) {
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

void test_insert(int hybrid_mode, char *index_name, char *key_path, char *inner_file, int count, int has_size,
                 int insert_count = 10000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    // To obtain a time-based seed
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            std::default_random_engine(seed));
    std::sort(keys, keys + count - insert_count);
    alex::Alex<KeyType, ValueType> index(hybrid_mode, true, index_name, inner_file);
    auto values = new std::pair<KeyType, ValueType>[count - insert_count];
    for (int i = 0; i < count - insert_count; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<ValueType>(keys[i] + 1);
    }
    std::cout << "start bulk_load..." << std::endl;
    index.bulk_load(values, count - insert_count);
    index.sync_metanode(true);
    index.sync_metanode(false);

    std::vector<int> v(insert_count);
    for (int i = 0; i < insert_count; i++)
        v[i] = i;
    std::vector<int> i_o(insert_count);
    for (int i = 0; i < insert_count; i++) {
        i_o[i] = get_next(v, insert_count - i);
    }
    std::cout << "start to insert... " << std::endl;

    int _start = count - insert_count;

    long long *latency = new long long[insert_count];
#if Profiling
    long long *sels = new long long[insert_count];
    long long *inls = new long long[insert_count];
    long long *smls = new long long[insert_count];
    long long *mals = new long long[insert_count];
#endif
    long long sel;
    long long sml;
    long long inl;
    long long mal;
    double tsmo = 0;
    int smo = 0;
    auto bulk_s = std::chrono::high_resolution_clock::now();
    for (int i = _start, j = 0; i < count; i++, j++) {
//        if ((i - _start) % 100000 == 0) std::cout << (i - _start) / 100000 << std::endl;
#if Profiling
        auto i_start = std::chrono::high_resolution_clock::now();
#endif
        index.insert_disk(keys[_start + i_o[j]], keys[_start + i_o[j]] + 1,
                          &sel, &inl, &sml, &mal, &smo);
#if Profiling
        auto i_end = std::chrono::high_resolution_clock::now();
        auto i_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                i_end - i_start)
                .count();
        latency[j] = i_b;

        sels[j] = sel;
        inls[j] = inl;
        smls[j] = sml;
        mals[j] = mal;
#endif
        tsmo += smo;
        smo = 0;

    }
    index.sync_metanode(true);
    index.sync_metanode(false);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    auto bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << double(count - insert_count) * 1e9 / bulk_t << " ops" << std::endl;
    std::cout << "inner node size:" << index.get_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << index.get_file_size() << " bytes" << std::endl;
    std::cout << "total smo:" << tsmo << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    delete[]keys;
    delete[]values;
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

void
test_mixed(int hybrid_mode, char *index_name, char *key_path, char *inner_file, int count, int has_size, int case_id,
           int total_op = 20000000) {
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    // To obtain a time-based seed
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            std::default_random_engine(seed));
    std::sort(keys, keys + count - total_op);
    auto values = new std::pair<KeyType, ValueType>[count - total_op];
    for (int i = 0; i < count - total_op; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<ValueType>(keys[i] + 1);
    }

    int start = count - total_op;
    std::vector<int> v(total_op);
    for (int i = 0; i < total_op; i++)
        v[i] = i;
    std::vector<KeyType> i_o(total_op);
    for (int i = 0; i < total_op; i++) {
        i_o[i] = keys[start + get_next(v, total_op - i)];
    }
    for (int i = 0; i < total_op; i++) {
        keys[start + i] = i_o[i];
    }
    int s_count;
    int step_s;
    int step_w;
    if (case_id == 1) { // read_heavey
        s_count = 18 * total_op / 20;
        step_s = 18;
        step_w = 2;
    } else if (case_id == 2) { // write_heavey
        s_count = 2 * total_op / 20;
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

#if DiskSetting
    alex::Alex<KeyType, ValueType> index(hybrid_mode, true, index_name, inner_file);
#else
    alex::Alex<KeyType, ValueType> index;
#endif

    std::cout << "start to bulk... " << std::endl;
    index.bulk_load(values, count - total_op);

#if DiskSetting
    index.sync_metanode(true);
    index.sync_metanode(false);
#endif

    std::cout << "start to hybrid test... " << std::endl;
    bool found = false;
    int nc = 0;
    int c = 0;
    int l = 0;
    int ci = 0;
    ValueType _v;
    long long sel;
    long long sml;
    long long inl;
    long long mal;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    int _start = count - total_op;
    for (int i = _start, j = 0; i < (_start + w_count) && j < s_count;) {
        for (int _i = 0; _i < step_w; _i++) {
#if DiskSetting
            int smo;
            index.insert_disk(keys[i], keys[i] + 1, &sel, &inl, &sml, &mal, &smo);
#else
            index.insert(keys[i], keys[i]+1);
#endif
            i += 1;
        }
        for (int _j = 0; _j < step_s; _j++) {
#if DiskSetting
            found = index.get_leaf_disk(search_keys[j], &_v, &c, &l, &ci);
            if (!found) {
                nc += 1;
            }
#else
            index.find(search_keys[j]);
#endif
            j += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    double scount = total_op;
    std::cout << 1e9 * scount / bulk_lookup_time << " ops" << std::endl;
    std::cout << nc << std::endl;
    std::cout << "inner node size:" << index.get_memory_size() << " bytes" << std::endl;
    std::cout << "file size:" << index.get_file_size() << " bytes" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    return;
}

void test_disk_bulk_search(int hybrid_mode, char *index_name, char *key_path, char *inner_file, int count, int has_size,
                           int s_count, int case_id, int step, int r_size) {
    ////////bulk///////////
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
    alex::Alex<KeyType, ValueType> index(hybrid_mode, true, index_name, inner_file);
    auto values = new std::pair<KeyType, ValueType>[count];
    for (int i = 0; i < count; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<ValueType>(keys[i] + 1);
    }
    std::cout << "start bulk_load..." << std::endl;
    auto bulk_s = std::chrono::high_resolution_clock::now();
    index.bulk_load(values, count);
    index.sync_metanode(true);
    index.sync_metanode(false);
    auto bulk_e = std::chrono::high_resolution_clock::now();
    long long bulk_t = std::chrono::duration_cast<std::chrono::nanoseconds>(
            bulk_e - bulk_s)
            .count();
    std::cout << "bulk_load_time:" << bulk_t / 1e9 << " sec" << std::endl;
    std::cout << "inner_node_size:" << index.get_memory_size() << " bytes" << std::endl;
    std::cout << "file_size:" << index.get_file_size() << " bytes" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    delete[]values;
    ///////////lookup////////////
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count - r_size * 2, s_count); //trick for scan
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count - r_size * 2, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
    bool found = false;
    ValueType v;
    int c = 0;
    double tc = 0;
    int nc = 0;
    double tl = 0;
    int l = 0;
    double tci = 0;
    int ci = 0;
#if Profiling
    int *_level_ = new int[s_count];
    int *_i_bs = new int[s_count];
#endif
    long long *latency = new long long[s_count];

    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        l = 0;
        ci = 0;
#if Profiling
        std::chrono::high_resolution_clock::time_point i_lookups_start_time = std::chrono::high_resolution_clock::now();
#endif
        found = index.get_leaf_disk(search_keys[j], &v, &c, &l, &ci);

        std::chrono::high_resolution_clock::time_point i_lookups_end_time = std::chrono::high_resolution_clock::now();
        long long i_batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                i_lookups_end_time - i_lookups_start_time).count();
        latency[j] = i_batch_lookup_time;
#if Profiling
        _level_[j] = l;
        _i_bs[j] = ci;
#endif
        tc += c;
        tl += l;
        tci += ci;

        if (!found) {
            nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << "lookup_ops:" << 1e9 * s_count / batch_lookup_time << std::endl;
    std::cout << "block/lookup:" << tc / s_count << std::endl;
    std::cout << "levels/lookup:" << tl / s_count << std::endl;
    std::cout << "inners/lookup:" << tci / s_count << std::endl;
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
    ///////////range//////////
    KeyType rs[r_size];
    tc = 0;
    lookups_start_time = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < s_count; j++) {
        c = 0;
        index.scan(rs, search_keys[j], &c, r_size);
        tc += c;
    }
    lookups_end_time = std::chrono::high_resolution_clock::now();
    batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    std::cout << "scan_ops:" << 1e9 * s_count / batch_lookup_time << std::endl;
    std::cout << "block/scan:" << tc / s_count << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
}

int main(int argc, char *argv[]) {
    auto flags = parse_flags(argc, argv);
    std::string key_file_path = get_required(flags, "keys_file");
    std::string op_type = get_required(flags, "op_type");
    std::string index_name = get_required(flags, "index_file");
    std::string data_file = index_name + "_data";
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
                       const_cast<char *>(data_file.c_str()), count, has_size);
    } else if (op_type == "lookup") {
        test_lookup(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                    const_cast<char *>(data_file.c_str()), count, has_size, search_count, case_id, step, need_sleep);
    } else if (op_type == "scan") {
        test_scan(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                  const_cast<char *>(data_file.c_str()), count, has_size, search_count, case_id, r_size, step);
    } else if (op_type == "insert") {
        test_insert(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                    const_cast<char *>(data_file.c_str()), count, has_size);
    } else if (op_type == "mix_workload") {
        test_mixed(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                   const_cast<char *>(data_file.c_str()), count, has_size, case_id);
    } else if (op_type == "bulk_search_range") {
        test_disk_bulk_search(memory_type, const_cast<char *>(index_name.c_str()),
                              const_cast<char *>(key_file_path.c_str()), const_cast<char *>(data_file.c_str()), count,
                              has_size, search_count, case_id, step, r_size);
    }
    return 0;
}
