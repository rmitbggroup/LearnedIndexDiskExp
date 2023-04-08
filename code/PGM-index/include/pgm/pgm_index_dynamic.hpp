// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2019 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "pgm_index.hpp"
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <sstream>

#define Profiling 0

namespace pgm {
/**
 * A sorted associative container that contains key-value pairs with unique keys.
 * @tparam K the type of a key
 * @tparam V the type of a value
 * @tparam PGMType the type of @ref PGMIndex to use in the container
 */
template<typename K, typename V, typename PGMType = PGMIndex<K, 64>>
class DynamicPGMIndex {
    class ItemA;
    class ItemB;
    class Iterator;

    using Item = std::conditional_t<std::is_pointer_v<V> || std::is_arithmetic_v<V>, ItemA, ItemB>;
    using Level = std::vector<Item>;



    /// just for re-load remove 'const' (disk version)
    uint8_t base;            ///< base^i is the maximum size of the ith level.
    uint8_t min_level;       ///< Levels 0..min_level are combined into one large level.
    uint8_t min_index_level; ///< Minimum level on which an index is constructed.
    size_t buffer_max_size;        ///< Size of the combined upper levels, i.e. max_size(0) + ... + max_size(min_level).
    uint8_t used_levels;           ///< Equal to 1 + last level whose size is greater than 0, or = min_level if no data.
    std::vector<Level> levels;     ///< (i-min_level)th element is the data array at the ith level.
    std::vector<PGMType> pgms;     ///< (i-min_index_level)th element is the index at the ith level.

    /// ***** for disk version *****///
    char *suffix = "test";
    typedef struct {
        uint8_t base;
        uint8_t min_level;
        uint8_t min_index_level;
        size_t buffer_max_size;
        uint8_t used_levels;
        int32_t item_count[12];
    } Metadata;
    const int MetadataSize = sizeof(Metadata);
public:
    typedef struct {
        K key;
        V value;
    } ItemOnDisk;
    bool inner_disk = false;
private:
    const int32_t ItemOnDiskSize = sizeof(ItemOnDisk);
    const int32_t ItemCountPerBlock = int32_t(BLOCK_SIZE/ItemOnDiskSize);

    Metadata metadata;


    std::vector<FILE *> data_file_handlers;
    std::vector<FILE *> index_file_handlers;
    FILE *metadata_handler;



    // supported internal operation
    int32_t ACCESSED_BLOCK_COUNT_Data = 0;
    void ADD_BLOCK_COUNT_Data() {ACCESSED_BLOCK_COUNT_Data += 1;}
    #define BLOCK_RECORDING 1

    void initial_metadata() {
        metadata.min_index_level = min_index_level;
        metadata.min_level = min_level;
        metadata.base = base;
        metadata.buffer_max_size = buffer_max_size;
        metadata.used_levels = used_levels;
        for (auto i = 0; i < 12; i++) {
            metadata.item_count[i] = 0;
        }
    }

    void load_metadata() {
        char *block = new char[BLOCK_SIZE];
        read_block(metadata_handler, block, 0);
        memcpy(&metadata, block, MetadataSize);
        min_index_level = metadata.min_index_level;
        min_level = metadata.min_level;
        base = metadata.base;
        buffer_max_size = metadata.buffer_max_size;
        used_levels = metadata.used_levels;
        delete []block;
    }

    void flush_metadata() {
        char *block = new char[BLOCK_SIZE];
        memcpy(block,&metadata,MetadataSize);
        write_block(metadata_handler, block, 0);
        delete []block;
        return;
    }

    int32_t get_item_count(int level) {
        if (level > 11) return -1;
        return metadata.item_count[level-min_level];
    }

    void set_item_count(int level, int new_count) {
        metadata.item_count[level-min_level] = new_count;
        return;
    }

    void create_related_empty_files() {
        for (auto i = 0; i < 12 - min_level; i++) {
            std::ostringstream oss;
            oss << i << "." << suffix << "_data";
            FILE *fp = fopen(oss.str().c_str(), "wb");
            fclose(fp);
        }
        for (auto i = 0; i < 12 - min_index_level; i++) {
            std::ostringstream oss;
            oss << i << "." << suffix << "_index";
            FILE *fp = fopen(oss.str().c_str(), "wb");
            fclose(fp);
        }

        {
            std::ostringstream oss;
            oss << suffix << ".metadata";
            FILE *fp = fopen(oss.str().c_str(), "wb");
            fclose(fp);
        }
        return;
    }

    void close_files() {
        for (auto i = 0; i < 12 - min_level; i++) {
            fclose(data_file_handlers[i]);
        }

        for (auto i = 0; i < 12 - min_index_level; i++) {
            fclose(index_file_handlers[i]);
        }

        fclose(metadata_handler);
    }

    void obtain_handlers() {
        for (auto i = 0; i < 12 - min_level; i++) {
            std::ostringstream oss;
            oss << i << "." << suffix << "_data";
            data_file_handlers[i] = fopen(oss.str().c_str(), "r+b");

        }
        for (auto i = 0; i < 12 - min_index_level; i++) {
            std::ostringstream oss;
            oss << i << "." << suffix << "_index";
            index_file_handlers[i] = fopen(oss.str().c_str(), "r+b");

        }

        {
            std::ostringstream oss;
            oss << suffix << ".metadata";
            metadata_handler = fopen(oss.str().c_str(), "r+b");

        }
        return;
    }

    FILE *get_data_file_handler(uint8_t level) {
        return data_file_handlers[level-min_level];
    }

    FILE *get_index_file_handler(uint8_t level) {
        return index_file_handlers[level-min_index_level];
    }

    void read_block_or_not(FILE * file_handler, int32_t *last_block, int32_t block_to_fetch, char *block_data) {
            if (block_to_fetch != *last_block) {
                read_block(file_handler, block_data, block_to_fetch);
                *last_block = block_to_fetch;
#if BLOCK_RECORDING
                ADD_BLOCK_COUNT_Data();
#endif
            }
        }

    void item_block_offset(size_t it, int *block_to_fetch, int *item_offset) {
            *block_to_fetch = int(it / ItemCountPerBlock);
            *item_offset = it % ItemCountPerBlock;
        }

    ItemOnDisk obtain_item_on_disk(size_t pos, FILE *file_handler, int *last_block, char *block_data) {
            int block_to_fetch = -1;
            int item_offset = -1;
            item_block_offset(pos, &block_to_fetch, &item_offset);
            read_block_or_not(file_handler, last_block, block_to_fetch, block_data);
            ItemOnDisk *_it_seg = (ItemOnDisk *)(block_data + sizeof(ItemOnDisk) * item_offset);
            return {_it_seg->key, _it_seg->value};
        }

    size_t lower_bound_bl_disk(FILE *file_handler, size_t first, size_t last, const K &x, size_t N, ItemOnDisk *item) {
            int last_block = -1;
            char *block_data = new char[BLOCK_SIZE];
            if (first == last) return first;
            auto n = last - first;
            while (n > 1) {
                auto half = n / 2;
                first = obtain_item_on_disk(first+half, file_handler, &last_block, block_data).key < x
                        ? first + half : first;
                n -= half;
            }
            // the second part should have a lower overhead
            size_t pos = first + (obtain_item_on_disk(first, file_handler, &last_block, block_data).key < x);
            if (pos != N) {
                *item = obtain_item_on_disk(pos, file_handler, &last_block, block_data);
            }
            delete []block_data;
            return pos;
        }

    void pairwise_merge_on_disk(ItemOnDisk new_item, uint8_t target,
                                size_t size_hint, size_t insertion_point) {
        ItemOnDisk *tmp_a = new ItemOnDisk [size_hint + get_item_count(target)];
        ItemOnDisk *tmp_b = new ItemOnDisk [size_hint + get_item_count(target)];
        char *block_data = new char[BLOCK_SIZE];

        // read the data in min_level
        {
            size_t item_count = get_item_count(min_level);
            int block_count = item_count / ItemCountPerBlock;
            int item_count_last = item_count % ItemCountPerBlock;
            if (item_count_last > 0) block_count += 1;
            size_t cp_offset = 0;
            for (auto i = 0; i < block_count; i++) {
                read_block(get_data_file_handler(min_level), block_data, i);
                int cp_size = ItemCountPerBlock;
                if (i == block_count - 1 && item_count_last > 0) cp_size = item_count_last;
                memcpy(tmp_b + cp_offset, block_data, ItemOnDiskSize * cp_size);
                cp_offset += cp_size;
            }

            memcpy(tmp_a, tmp_b, ItemOnDiskSize * insertion_point);
            tmp_a[insertion_point] = new_item;
            // check whether the insertion point is the end of b
            // todo: check the `index'
            if (insertion_point < item_count)
                memcpy(tmp_a + insertion_point + 1, tmp_b + insertion_point, ItemOnDiskSize * (item_count - insertion_point));
        }

        size_t tmp_size = get_item_count(min_level) + 1;
        // check order
//        for (int i = 0; i < tmp_size-1; i++) {
//            if (tmp_a[i].key > tmp_a[i+1].key) {
//                std::cout << "wrong" << std::endl;
//            }
//        }
        uint8_t merge_limit = get_item_count(target) == 0 ? target - 1 : target;
        auto alternate = true;
        for (uint8_t i = 1 + min_level; i <= merge_limit; ++i, alternate = !alternate) {
            auto tmp = (alternate ? tmp_a : tmp_b);
            auto out = (alternate ? tmp_b : tmp_a);
            auto item_count_i = get_item_count(i);

            {
                ItemOnDisk *tmp_i = new ItemOnDisk [item_count_i];
                int block_count = item_count_i / ItemCountPerBlock;
                int item_count_last = item_count_i % ItemCountPerBlock;
                if (item_count_last > 0) block_count += 1;
                size_t cp_offset = 0;
                for (auto bi = 0; bi < block_count; bi++) {
                    read_block(get_data_file_handler(i), block_data, bi);
                    int cp_size = ItemCountPerBlock;
                    if (bi == block_count - 1 && item_count_last > 0) cp_size = item_count_last;
                    memcpy(tmp_i + cp_offset, block_data, ItemOnDiskSize * cp_size);
                    cp_offset += cp_size;
                }
//                for (int _i = 0; _i < item_count_i-1; _i++) {
//                    if (tmp_i[_i].key > tmp_i[_i+1].key) {
//                        std::cout << "wrong" << std::endl;
//                    }
//                }
                // merge process
                int l_i = 0;
                int r_i = 0;
                int o_i = 0;
                while (l_i < tmp_size && r_i < item_count_i) {
                    if (tmp[l_i].key < tmp_i[r_i].key) {
                        out[o_i] = tmp[l_i];
                        l_i += 1;
                    } else if (tmp[l_i].key > tmp_i[r_i].key) {
                        out[o_i] = tmp_i[r_i];
                        r_i += 1;
                    } else {
                        printf("LOG - repeated key found");
                        out[o_i] = tmp[l_i];
                        l_i += 1;
                        r_i += 1;
                    }
                    o_i += 1;
                }
                // todo: copy the rest to the results
                if (l_i < tmp_size) {
                    memcpy(out+o_i, tmp + l_i, (tmp_size - l_i) * ItemOnDiskSize);
                    o_i += (tmp_size - l_i);
                } else if (r_i<item_count_i) {
                    memcpy(out+o_i, tmp_i + r_i, (item_count_i - r_i) * ItemOnDiskSize);
                    o_i += (item_count_i - r_i);
                }
                tmp_size = o_i;
                delete []tmp_i;
            }

            // clear item count
            set_item_count(i, 0);
            // todo: do we need to truncate the file
        }
        set_item_count(min_level, 0);

        // write the data into target-file & set the size of target file
        ItemOnDisk *items = (ItemOnDisk *)block_data;
        size_t block = 0;
        size_t start = 0;
        ItemOnDisk *final = (alternate ? tmp_a : tmp_b);
        for (;start < tmp_size;) {
            size_t end = std::min(start + ItemCountPerBlock, tmp_size);
            // todo: optimize this
            for (auto i = start; i < end; i++) {
                items[i-start].key = final[i].key;
                items[i-start].value = final[i].value;
            }
            write_block(get_data_file_handler(target), block_data, block);
            start += ItemCountPerBlock;
            block += 1;
        }
        delete []block_data;
        set_item_count(target, tmp_size);

        // build index if needed
        if (has_pgm(target)) {
            // todo: check the item type
            std::vector<K> item_vec(tmp_size);
//            item_vec.reserve(tmp_size);
            ItemOnDisk *in = alternate ? tmp_a : tmp_b;
            for (int i = 0; i < tmp_size; i++) {
                item_vec[i] = in[i].key;
            }
//            auto item_vec = alternate ? std::vector<ItemOnDisk> (tmp_a, tmp_a + sizeof(tmp_a))
//                                      :  std::vector<ItemOnDisk> (tmp_b, tmp_b + sizeof(tmp_b));
            // flush index as needed
            pgm(target) = PGMType(item_vec.begin(), item_vec.end(), inner_disk, get_index_file_handler(target));
        }

        delete []tmp_a;
        delete []tmp_b;
    }

    void insert_item_on_file(int level, size_t position, ItemOnDisk item) {
        FILE *file_handler = get_data_file_handler(level);
        size_t item_count = get_item_count(level);

        char *block_data = new char[BLOCK_SIZE];
        int block_id = position / ItemCountPerBlock;
        int offset = position % ItemCountPerBlock;
        // read all data >= position into main memory, and then rewrite them on disk
        if (position == item_count) {
            if (offset == 0) { // we just write a block at block_id
                memcpy(block_data, &item, ItemOnDiskSize);
                write_block(file_handler, block_data, block_id);
            } else { // read the original block and add the new item
                read_block(file_handler, block_data, block_id);
                ItemOnDisk *temp = (ItemOnDisk *)block_data;
                temp[offset] = item;
                write_block(file_handler, block_data, block_id);
            }
            delete []block_data;
            return;
        }
        size_t change_count = item_count - position;
        ItemOnDisk *items = new ItemOnDisk[offset + change_count + 1];

        int block_end = item_count / ItemCountPerBlock;
        int item_in_last_block = item_count % ItemCountPerBlock;
        if (item_in_last_block > 0) block_end += 1; else item_in_last_block = ItemCountPerBlock;

        // read all the tuples to be moved
        size_t cp_count = 0;
        for (auto i = block_id; i < block_end; i++) {
            read_block(file_handler, block_data, i);
            if (i == block_id) {
                memcpy(items, block_data, ItemOnDiskSize * offset);
                items[offset] = item;
                cp_count = offset + 1;
                int remain = i == block_end - 1 ? item_in_last_block - offset : ItemCountPerBlock - offset;
                memcpy(items + cp_count, block_data + ItemOnDiskSize * offset, ItemOnDiskSize * remain);
                cp_count += remain;
                continue;
            }
            memcpy(items + cp_count, block_data, ItemOnDiskSize * ( i == block_end - 1 ? item_in_last_block : ItemCountPerBlock));
            cp_count += (i == block_end - 1 ? item_in_last_block : ItemCountPerBlock);
        }

        assert(cp_count == offset + change_count + 1);

        // write all the tuples in
        ItemOnDisk *to_write_items = (ItemOnDisk *)block_data;
        size_t block = block_id;
        size_t start = 0;
        for (;start < cp_count;) {
            size_t end = std::min(start + ItemCountPerBlock, cp_count);
            for (auto i = start; i < end; i++) {
                to_write_items[i-start].key = items[i].key;
                to_write_items[i-start].value = items[i].value;
            }
            write_block(file_handler, block_data, block);
            start += ItemCountPerBlock;
            block += 1;
        }

        delete []block_data;
        delete []items;
    }

    void insert_on_disk_internal(ItemOnDisk new_item, long long *search_latency,
                                 long long *insert_latency, long long *smo_latency,
                                 int *smo_count, int *updated_tuples) {
#if Profiling
        *smo_count = 0;
        *updated_tuples = 0;
        std::chrono::high_resolution_clock::time_point search_start = std::chrono::high_resolution_clock::now();
        std::chrono::high_resolution_clock::time_point search_end = search_start;
        std::chrono::high_resolution_clock::time_point insert_start = std::chrono::high_resolution_clock::now();
        std::chrono::high_resolution_clock::time_point insert_end = insert_start;
        std::chrono::high_resolution_clock::time_point smo_start = std::chrono::high_resolution_clock::now();
        std::chrono::high_resolution_clock::time_point smo_end = smo_start;
#endif
        ItemOnDisk item;
        auto item_count = get_item_count(min_level);
#if Profiling
        search_start = std::chrono::high_resolution_clock::now();
#endif
        auto insertion_point = lower_bound_bl_disk(get_data_file_handler(min_level),
                                                   0, item_count,
                                                   new_item.key, item_count, &item);
#if Profiling
        search_end = std::chrono::high_resolution_clock::now();
#endif
        if (insertion_point != item_count && item.key == new_item.key) {
            // todo: update the value on disk; we hope we do this thing in `lower_bound_bl_disk`
            //  to reduce the i/o cost
//            return;
        } else if (get_item_count(min_level) < buffer_max_size) {
#if Profiling
            insert_start = std::chrono::high_resolution_clock::now();
#endif
//            std::cout<<"insert into the data file" << std::endl;
            insert_item_on_file(min_level, insertion_point, new_item);
            used_levels = used_levels == min_level ? min_level + 1 : used_levels;
            set_item_count(min_level, item_count + 1);
            // todo: update the metadata
#if Profiling
            insert_end = std::chrono::high_resolution_clock::now();
#endif
//            return;
        } else {
#if Profiling
            *smo_count += 1;
            smo_start = std::chrono::high_resolution_clock::now();
#endif
            size_t slots_required = buffer_max_size + 1;
            uint8_t i;
            for (i = min_level + 1; i < used_levels; ++i) {
                auto item_count_i = get_item_count(i);
                auto slots_left_in_level = max_size(i) - item_count_i;
                if (slots_required <= slots_left_in_level)
                    break;
                slots_required += item_count_i;
            }
#if Profiling
            *updated_tuples += slots_required;
#endif
            auto need_new_level = i == used_levels;
            if (need_new_level) {
                ++used_levels;
                if (i - min_index_level >= int(pgms.size()))
                    pgms.emplace_back();
                // we suppose we have allocated these files.
            }
//        std::cout<<"merge..." << std::endl;
            pairwise_merge_on_disk(new_item, i, slots_required, insertion_point);
#if Profiling
            smo_end = std::chrono::high_resolution_clock::now();
#endif
        }
#if Profiling
        *search_latency = std::chrono::duration_cast<std::chrono::nanoseconds>(search_end - search_start).count();
        *insert_latency = std::chrono::duration_cast<std::chrono::nanoseconds>(insert_end - insert_start).count();
        *smo_latency = std::chrono::duration_cast<std::chrono::nanoseconds>(smo_end - smo_start).count();
#endif
        return;
    }
    // support interface
public:
    size_t report_disk_file_size() {
        size_t total_size = 0;
        for (auto i = min_level; i < used_levels; ++i) {
            auto item_count = get_item_count(i);
            if (item_count == 0) continue;
            FILE *_temp = get_data_file_handler(i);
            fseek(_temp,0,SEEK_END);
            total_size += ftell(_temp);
            if (inner_disk && has_pgm(i)) {
                _temp = get_index_file_handler(i);
                fseek(_temp,0,SEEK_END);
                total_size += ftell(_temp);
            }
        }
        return total_size;
    }

    size_t report_main_memory_size() {
        size_t total_size = 0;
        if (inner_disk) return total_size;
        for (auto i = min_level; i < used_levels; ++i) {
            auto item_count = get_item_count(i);
            if (item_count == 0) continue;
            if (has_pgm(i)) {
                total_size += pgm(i).size_in_bytes();
            }
        }
        return total_size;
    }
    int find_on_disk(const K &key, int *c, int *ic, int *lc) {
        ACCESSED_BLOCK_COUNT_Data = 0;
        ItemOnDisk item;
            for (auto i = min_level; i < used_levels; ++i) {
                auto item_count = get_item_count(i);
                if (item_count == 0) continue;

                size_t first = 0;
                size_t last = item_count;
                if (has_pgm(i)) {
                    // pgm search on disk
                    auto range = inner_disk ? pgm(i).search_disk(get_index_file_handler(i), key, c, lc) : pgm(i).search(key);
                    first = range.lo;
                    last = range.hi;
                }
                // binary search on data file
                *ic = *c;
                auto pos = lower_bound_bl_disk(get_data_file_handler(i), first, last, key, item_count, &item);
                *c += ACCESSED_BLOCK_COUNT_Data;
                if (pos != item_count && item.key == key) {
                    return 1;
                }
            }
            return 0;
        }

    int range_on_disk(const K &lo, const K &hi, int len, int *c, ItemOnDisk *tmp_a, ItemOnDisk *tmp_b, int *fsize, ItemOnDisk *results = nullptr) {
        ACCESSED_BLOCK_COUNT_Data = 0;
            if (lo > hi)
                throw std::invalid_argument("lo > hi");

            int used_a = 0;
            int used_b = 0;
            auto alternate = true;
            char *block_data = new char[BLOCK_SIZE];
            for (auto i = min_level; i < used_levels; i++) {
                int item_count_i = get_item_count(i);
                if (item_count_i == 0) continue;

                size_t lo_first = 0;
                size_t lo_last = item_count_i;
                if (has_pgm(i)) {
                    int lc = 0;
                    auto range = inner_disk ? pgm(i).search_disk(get_index_file_handler(i), lo, c, &lc) : pgm(i).search(lo);
                    lo_first = range.lo;
                    lo_last = range.hi;
                    // todo: Currently, I just remove the search for the upper bound
                } else {
                    ItemOnDisk  _item;
                    lo_first = lower_bound_bl_disk(get_data_file_handler(i), lo_first, lo_last, lo, lo_last, &_item);
                }
                // todo: check we use `begin' or `clear'
                auto tmp_size = (alternate ? used_a : used_b);
                ItemOnDisk *tmp_it = (alternate ? tmp_a : tmp_b);
                ItemOnDisk *out_it = (alternate ? tmp_b : tmp_a);

                int l_i = 0;
                int o_i = 0;
                int block_id = lo_first / ItemCountPerBlock;
                int offset = lo_first % ItemCountPerBlock;
                //ItemOnDisk *item_r = new ItemOnDisk[BLOCK_SIZE];
                int r_i = lo_first;
                bool is_out_of_range = false;
                // we set l_i `<=' is to handle the case that the remaining in
                while (l_i <= tmp_size && r_i < item_count_i && !is_out_of_range && offset < ItemCountPerBlock) {
                    int last_block = -1;
                    read_block_or_not(get_data_file_handler(i), &last_block, block_id, block_data);
                    ItemOnDisk *item_r = (ItemOnDisk *) (block_data);

                    while (l_i < tmp_size && r_i < item_count_i && offset < ItemCountPerBlock) {
                        // skip the useless ones
                        if (o_i >= len) {is_out_of_range = true; break;}
                        if (item_r[offset].key < lo) {
                            offset += 1;
                            r_i += 1;
                            continue;
                        } else if (item_r[offset].key > hi) {
                            is_out_of_range = true;
                            break;
                        }


                        if (tmp_it[l_i].key < item_r[offset].key) {


                            out_it[o_i].key = tmp_it[l_i].key;
                            out_it[o_i].value = tmp_it[l_i].value;
                            l_i += 1;
                        } else if (tmp_it[l_i].key > item_r[offset].key) {


                            out_it[o_i].key = item_r[offset].key;
                            out_it[o_i].value = item_r[offset].value;
                            offset += 1;
                            r_i += 1;
                        } else {


                            out_it[o_i].key = tmp_it[l_i].key;
                            out_it[o_i].value = tmp_it[l_i].value;
                            l_i += 1;
                            r_i += 1;
                            offset += 1;
                        }
                        o_i += 1;
                    }

                    // if l_i == tmp_size and we still can fetch from right table
                    while (offset < ItemCountPerBlock && l_i == tmp_size && r_i < item_count_i && !is_out_of_range) {
                        if (o_i >= len) {is_out_of_range = true; break;}
                        if (item_r[offset].key < lo) {
                            offset += 1;
                            r_i += 1;
                            continue;
                        } else if (item_r[offset].key > hi) {
                            is_out_of_range = true;
                            break;
                        }

                        out_it[o_i].key = item_r[offset].key;
                        out_it[o_i].value = item_r[offset].value;
                        offset += 1;
                        r_i += 1;
                        o_i += 1;
                    }

                    offset = 0;
                    block_id += 1;
                }
                // here, when only to handle the case that the remaining tuples in left table
                // append the remaining from one side
                for (; l_i < tmp_size && (is_out_of_range || r_i == item_count_i); l_i++) {

                    out_it[o_i].key = tmp_it[l_i].key;
                    out_it[o_i].value = tmp_it[l_i].value;
                    o_i += 1;
                }
                alternate ? used_b = o_i : used_a = o_i;
//                tmp_it.clear();
                alternate = !alternate;
            }
            *c += ACCESSED_BLOCK_COUNT_Data;
//            ItemOnDisk *out = (alternate ? tmp_a : tmp_b);
            *fsize = std::min(alternate ? used_a: used_b, len);
        return (alternate ? 0 : 1);

        }

    void insert_on_disk(K key, V value, long long *search_latency,
                        long long *insert_latency, long long *smo_latency,int *smo_count, int *updated_tuples) {
        ItemOnDisk x;
        x.key = key;
        x.value = value;
        insert_on_disk_internal(x, search_latency, insert_latency, smo_latency, smo_count, updated_tuples);
    }
    /// ***** End for disk version *****///

private:
    const Level &level(uint8_t level) const { return levels[level - min_level]; }
    const PGMType &pgm(uint8_t level) const { return pgms[level - min_index_level]; }
    Level &level(uint8_t level) { return levels[level - min_level]; }
    PGMType &pgm(uint8_t level) { return pgms[level - min_index_level]; }
    bool has_pgm(uint8_t level) const { return level >= min_index_level; }
    size_t max_size(uint8_t level) const { return size_t(1) << (level * ceil_log2(base)); }
    uint8_t max_fully_allocated_level() const { return min_level + 2; }
    uint8_t ceil_log_base(size_t n) const { return (ceil_log2(n) + ceil_log2(base) - 1) / ceil_log2(base); }
    constexpr static uint8_t ceil_log2(size_t n) { return n <= 1 ? 0 : sizeof(long long) * 8 - __builtin_clzll(n - 1); }

    void pairwise_merge(const Item &new_item,
                        uint8_t target,
                        size_t size_hint,
                        typename Level::iterator insertion_point) {
        Level tmp_a(size_hint + level(target).size());
        Level tmp_b(size_hint + level(target).size());

        // Insert new_item in sorted order in the first level
        auto alternate = true;
        auto it = std::move(level(min_level).begin(), insertion_point, tmp_a.begin());
        *it++ = new_item;
        it = std::move(insertion_point, level(min_level).end(), it);
        auto tmp_size = std::distance(tmp_a.begin(), it);

        // Merge subsequent levels
        uint8_t merge_limit = level(target).empty() ? target - 1 : target;
        // todo: we may do a multiple-source merge
        for (uint8_t i = 1 + min_level; i <= merge_limit; ++i, alternate = !alternate) {
            auto tmp_begin = (alternate ? tmp_a : tmp_b).begin();
            auto tmp_end = tmp_begin + tmp_size;
            auto out_begin = (alternate ? tmp_b : tmp_a).begin();
            decltype(out_begin) out_end;

            auto can_delete_permanently = i == used_levels - 1;
            if (can_delete_permanently)
                out_end = merge<true, true>(tmp_begin, tmp_end, level(i).begin(), level(i).end(), out_begin);
            else
                out_end = merge<false, true>(tmp_begin, tmp_end, level(i).begin(), level(i).end(), out_begin);
            tmp_size = std::distance(out_begin, out_end);

            // Empty this level and the corresponding index
            // todo: truncate file or set item count with 0
            level(i).clear();
            if (i >= max_fully_allocated_level())
                level(i).shrink_to_fit();
            if (has_pgm(i))
                pgm(i) = PGMType();
        }

        // todo: truncate file or set item count with 0
        level(min_level).clear();
        level(target) = std::move(alternate ? tmp_a : tmp_b);
        level(target).resize(tmp_size);

        // Rebuild index, if needed
        // todo: build index file if needed
        if (has_pgm(target))
            pgm(target) = PGMType(level(target).begin(), level(target).end());
    }

    void insert(const Item &new_item) {
        auto insertion_point = lower_bound_bl(level(min_level).begin(), level(min_level).end(), new_item);
        if (insertion_point != level(min_level).end() && *insertion_point == new_item) {
            *insertion_point = new_item;
            return;
        }

        if (level(min_level).size() < buffer_max_size) {
            level(min_level).insert(insertion_point, new_item);
            used_levels = used_levels == min_level ? min_level + 1 : used_levels;
            return;
        }

        size_t slots_required = buffer_max_size + 1;
        uint8_t i;
        for (i = min_level + 1; i < used_levels; ++i) {
            auto slots_left_in_level = max_size(i) - level(i).size();
            if (slots_required <= slots_left_in_level)
                break;
            slots_required += level(i).size();
        }

        auto need_new_level = i == used_levels;
        if (need_new_level) {
            ++used_levels;
            levels.emplace_back();
            if (i - min_index_level >= int(pgms.size()))
                pgms.emplace_back();
        }

        pairwise_merge(new_item, i, slots_required, insertion_point);
    }

public:

    using key_type = K;
    using mapped_type = V;
    using value_type = Item;
    using size_type = size_t;
    using iterator = Iterator;
    /**
     * Constructs an empty container.
     * @param base determines the size of the ith level as base^i
     * @param buffer_level determines the size of level 0, equal to the sum of base^i for i = 0, ..., buffer_level
     * @param index_level the minimum level at which an index is constructed to speed up searches
     */
    DynamicPGMIndex(uint8_t base = 8, uint8_t buffer_level = 0, uint8_t index_level = 0)
        : base(base),
          min_level(buffer_level ? buffer_level : ceil_log_base(128) - (base == 2)),
          min_index_level(std::max<size_t>(min_level + 1, index_level ? index_level : ceil_log_base(size_t(1) << 21))),
          buffer_max_size(),
          used_levels(min_level),
          levels(),
          pgms() {
        if (base < 2 || (base & (base - 1u)) != 0)
            throw std::invalid_argument("base must be a power of two");

        for (auto j = 0; j <= min_level; ++j)
            buffer_max_size += max_size(j);

        levels.resize(32 - used_levels);
        level(min_level).reserve(buffer_max_size);
        data_file_handlers.resize(12 - min_level);
        index_file_handlers.resize(12 - min_index_level);
        for (uint8_t i = min_level + 1; i < max_fully_allocated_level(); ++i)
            level(i).reserve(max_size(i));
    }


    /// in memory code
    template<typename Iterator>
    DynamicPGMIndex(Iterator first, Iterator last, uint8_t base = 8, uint8_t buffer_level = 0, uint8_t index_level = 0)
            : DynamicPGMIndex(base, buffer_level, index_level) {
        size_t n = std::distance(first, last);
        used_levels = std::max<uint8_t>(ceil_log_base(n), min_level) + 1;
        levels.resize(std::max<uint8_t>(used_levels, 32) - min_level + 1);
        level(min_level).reserve(buffer_max_size);
        for (uint8_t i = min_level + 1; i < max_fully_allocated_level(); ++i)
            level(i).reserve(max_size(i));

        if (n == 0) {
            used_levels = min_level;
            return;
        }

        // Copy only the first of each group of pairs with same key value
        auto &target = level(used_levels - 1);
        target.resize(n);
        auto out = target.begin();
        *out++ = Item(first->first, first->second);
        while (++first != last) {
            if (first->first < std::prev(out)->first)
                throw std::invalid_argument("Range is not sorted");
            if (first->first != std::prev(out)->first)
                *out++ = Item(first->first, first->second);
        }
        target.resize(std::distance(target.begin(), out));

        if (has_pgm(used_levels - 1)) {
            pgms = decltype(pgms)(used_levels - min_index_level);
            pgm(used_levels - 1) = PGMType(target.begin(), target.end());
        }
    }

    /**
     * Constructs the container on the sorted data in the range [first, last).
     * @tparam Iterator
     * @param first, last the range containing the sorted elements to be indexed
     * @param base determines the size of the ith level as base^i
     * @param buffer_level determines the size of level 0, equal to the sum of base^i for i = 0, ..., buffer_level
     * @param index_level the minimum level at which an index is constructed to speed up searches
     */
    template<typename Iterator>
    DynamicPGMIndex(bool is_first_time, bool _inner_disk, Iterator first, Iterator last, uint8_t base = 8, uint8_t buffer_level = 0, uint8_t index_level = 0)
        : DynamicPGMIndex(base, buffer_level, index_level) {

        size_t n = std::distance(first, last);
        used_levels = std::max<uint8_t>(ceil_log_base(n), min_level) + 1;
        levels.resize(std::max<uint8_t>(used_levels, 32) - min_level + 1);
        level(min_level).reserve(buffer_max_size);

        inner_disk = _inner_disk;
        std::cout << "is inner:" << inner_disk << std::endl;

        /// create/load related file; initial metadata
        if (is_first_time) create_related_empty_files();
        obtain_handlers();
        // we can flush the metadata at last instead of each operation
        if (is_first_time) initial_metadata();
        else load_metadata();
        // if is_first_time is false, we suppose `first' and `last' are null, we
        // do not need them
        if (!is_first_time) {
            // load related pgm tree into main memory
            pgms = decltype(pgms)(12);
            for (auto i  = min_level; i < 12; i++) {
                if (get_item_count(i) > 0 && has_pgm(i)){
                    pgm(i - min_index_level) = PGMType(get_index_file_handler(i), inner_disk);
                }
            }
            return;
        }
        for (uint8_t i = min_level + 1; i < max_fully_allocated_level(); ++i)
            level(i).reserve(max_size(i));

        // todo: set metadata
        if (n == 0) {
            used_levels = min_level;
            return;
        }

        // Copy only the first of each group of pairs with same key value
        auto &target = level(used_levels - 1);
        target.resize(n);
        auto out = target.begin();
        *out++ = Item(first->first, first->second);
        while (++first != last) {
            if (first->first < std::prev(out)->first)
                throw std::invalid_argument("Range is not sorted");
            if (first->first != std::prev(out)->first)
                *out++ = Item(first->first, first->second);
        }
        target.resize(std::distance(target.begin(), out));
        // flush the data into corresponding data file
        assert(get_item_count(used_levels - 1) == 0);
        char *block_data = new char[BLOCK_SIZE];
        ItemOnDisk *items = (ItemOnDisk *)block_data;
        int32_t block = 0;
        size_t start = 0;
        for (;start < target.size();) {
            size_t end = std::min(start + ItemCountPerBlock, target.size());
            for (auto i = start; i < end; i++) {
                items[i-start].key = target[i].first;
                items[i-start].value = target[i].second;
            }
            write_block(get_data_file_handler(used_levels-1), block_data, block);
            start += ItemCountPerBlock;
            block += 1;
        }
        delete []block_data;

        // update the metadata
        // todo: delete the data in level
        set_item_count(used_levels - 1, target.size());

        if (has_pgm(used_levels - 1)) {
            pgms = decltype(pgms)(used_levels - min_index_level);
            // flush the index into index file as needed
            pgm(used_levels - 1) = PGMType(target.begin(), target.end(), inner_disk, get_index_file_handler(used_levels - 1));
        }
    }

    /**
     * Inserts an element into the container if @p key does not exists in the container. If @p key already exists, the
     * corresponding value is updated with @p value.
     * @param key element key to insert or update
     * @param value element value to insert
     */
    void insert_or_assign(const K &key, const V &value) { insert(Item(key, value)); }

    /**
     * Removes the specified element from the container.
     * @param key key value of the element to remove
     */
    void erase(const K &key) { insert(Item(key)); }

    /**
     * Finds an element with key equivalent to @p key.
     * @param key key value of the element to search for
     * @return an iterator to an element with key equivalent to @p key. If no such element is found, end() is returned
     */
    iterator find(const K &key) const {
        for (auto i = min_level; i < used_levels; ++i) {
            if (level(i).empty())
                continue;

            auto first = level(i).begin();
            auto last = level(i).end();
            if (has_pgm(i)) {
                auto range = pgm(i).search(key);
                first = level(i).begin() + range.lo;
                last = level(i).begin() + range.hi;
            }
            auto it = lower_bound_bl(first, last, key);
            if (it != level(i).end() && it->first == key)
                return it->deleted() ? end() : iterator(this, i, it);
        }

        return end();
    }

    /**
     * Returns a copy of the elements with key between and including @p lo and @p hi.
     * @param lo lower endpoint of the range query
     * @param hi upper endpoint of the range query, must be greater than or equal to @p lo
     * @return a vector of key-value pairs satisfying the range query
     */
    std::vector<std::pair<K, V>> range(const K &lo, const K &hi) const {
        if (lo > hi)
            throw std::invalid_argument("lo > hi");

        Level tmp_a;
        Level tmp_b;
        auto alternate = true;

        for (auto i = min_level; i < used_levels; ++i) {
            if (level(i).empty())
                continue;

            auto lo_first = level(i).begin();
            auto lo_last = level(i).end();
            auto hi_first = level(i).begin();
            auto hi_last = level(i).end();
            if (has_pgm(i)) {
                auto range = pgm(i).search(lo);
                lo_first = level(i).begin() + range.lo;
                lo_last = level(i).begin() + range.hi;
                range = pgm(i).search(hi);
                hi_first = level(i).begin() + range.lo;
                hi_last = level(i).begin() + range.hi;
            }

            auto it_lo = lower_bound_bl(lo_first, lo_last, lo);
            auto it_hi = std::upper_bound(std::max(it_lo, hi_first), hi_last, hi);
            auto range_size = std::distance(it_lo, it_hi);
            if (range_size == 0)
                continue;

            auto tmp_size = (alternate ? tmp_a : tmp_b).size();
            (alternate ? tmp_b : tmp_a).resize(tmp_size + range_size);
            auto tmp_it = (alternate ? tmp_a : tmp_b).begin();
            auto out_it = (alternate ? tmp_b : tmp_a).begin();
            tmp_size = std::distance(out_it, merge<false, false>(tmp_it, tmp_it + tmp_size, it_lo, it_hi, out_it));
            (alternate ? tmp_b : tmp_a).resize(tmp_size);
            alternate = !alternate;
        }

        std::vector<std::pair<K, V>> result;
        result.reserve((alternate ? tmp_a : tmp_b).size());
        auto first = (alternate ? tmp_a : tmp_b).begin();
        auto last = (alternate ? tmp_a : tmp_b).end();
        for (auto it = first; it != last; ++it)
            if (!it->deleted())
                result.emplace_back(it->first, it->second);
        return result;
    }

    /**
     * Returns an iterator pointing to the first element that is not less than (i.e. greater or equal to) @p key.
     * @param key key value to compare the elements to
     * @return an iterator to an element with key not less than @p key. If no such element is found, end() is returned
     */
    iterator lower_bound(const K &key) const {
        typename Level::const_iterator lb;
        auto lb_set = false;
        uint8_t lb_level;
        std::set<K> deleted;

        for (auto i = min_level; i < used_levels; ++i) {
            if (level(i).empty())
                continue;

            auto first = level(i).begin();
            auto last = level(i).end();
            if (has_pgm(i)) {
                auto range = pgm(i).search(key);
                first = level(i).begin() + range.lo;
                last = level(i).begin() + range.hi;
            }

            for (auto it = lower_bound_bl(first, last, key);
                 it != level(i).end() && (!lb_set || it->first < lb->first); ++it) {
                if (it->deleted())
                    deleted.emplace(it->first);
                else if (deleted.find(it->first) == deleted.end()) {
                    if (it->first == key)
                        return iterator(this, i, it);
                    lb = it;
                    lb_level = i;
                    lb_set = true;
                    break;
                }
            }
        }

        if (lb_set)
            return iterator(this, lb_level, lb);
        return end();
    }

    /**
     * Checks if the container has no elements, i.e. whether begin() == end().
     * @return true if the container is empty, false otherwise
     */
    bool empty() const { return begin() == end(); }

    /**
     * Returns an iterator to the beginning.
     * @return an iterator to the beginning
     */
    iterator begin() const { return lower_bound(std::numeric_limits<K>::min()); }

    /**
     * Returns an iterator to the end.
     * @return an iterator to the end
     */
    iterator end() const { return iterator(this, levels.size() - 1, levels.back().end()); }

    /**
     * Returns the number of elements with key that compares equal to the specified argument key, which is either 1
     * or 0 since this container does not allow duplicates.
     * @param key key value of the elements to count
     * @return number of elements with the given key, which is either 1 or 0.
     */
    size_t count(const K &key) const { return find(key) == end() ? 0 : 1; }

    /**
     * Returns the number of elements in the container.
     * @return the number of elements in the container
     */
    size_t size() const {
        // TODO: scanning the levels and using a hash table for the encountered keys could be more time efficient
        return std::distance(begin(), end());
    }

    /**
     * Returns the size of the container (data + index structure) in bytes.
     * @return the size of the container in bytes
     */
    size_t size_in_bytes() const {
        size_t bytes = levels.size() * sizeof(Level);
        for (auto &l: levels)
            bytes += l.size() * sizeof(Item);
        return index_size_in_bytes() + bytes;
    }

    /**
     * Returns the size of the index used in this container in bytes.
     * @return the size of the index used in this container in bytes
     */
    size_t index_size_in_bytes() const {
        size_t bytes = 0;
        for (auto &p: pgms)
            bytes += p.size_in_bytes();
        return bytes;
    }

private:

    template<bool SkipDeleted, bool Move, typename In1, typename In2, typename OutIterator>
    static OutIterator merge(In1 first1, In1 last1, In2 first2, In2 last2, OutIterator result) {
        while (first1 != last1 && first2 != last2) {
            if (*first2 < *first1) {
                if constexpr (Move) *result = std::move(*first2);
                else *result = *first2;
                ++first2;
                ++result;
            } else if (*first1 < *first2) {
                if constexpr (Move) *result = std::move(*first1);
                else *result = *first1;
                ++first1;
                ++result;
            } else if (SkipDeleted && first1->deleted()) {
                ++first1;
                ++first2;
            } else {
                if constexpr (Move) *result = std::move(*first1);
                else *result = *first1;
                ++first1;
                ++first2;
                ++result;
            }
        }
        if constexpr (Move)
            return std::move(first2, last2, std::move(first1, last1, result));
        return std::copy(first2, last2, std::copy(first1, last1, result));
    }

    template<class RandomIt>
    static RandomIt lower_bound_bl(RandomIt first, RandomIt last, const K &x) {
        if (first == last)
            return first;
        auto n = std::distance(first, last);
        while (n > 1) {
            auto half = n / 2;
            __builtin_prefetch(&*(first + half / 2), 0, 0);
            __builtin_prefetch(&*(first + half + half / 2), 0, 0);
            first = first[half] < x ? first + half : first;
            n -= half;
        }
        return first + (*first < x);
    }
};

namespace internal {

/* LoserTree implementation adapted from Timo Bingmann's https://tlx.github.io and http://stxxl.org, and from
 * Johannes Singler's http://algo2.iti.uni-karlsruhe.de/singler/mcstl. These three libraries are distributed under the
 * Boost Software License 1.0. */
template<typename T>
class LoserTree {
    using Source = uint8_t;

    struct Loser {
        T key;         ///< Copy of the current key in the sequence.
        Source source; ///< Index of the sequence.
    };

    Source k;                  ///< Smallest power of 2 greater than the number of nodes.
    std::vector<Loser> losers; ///< Vector of size 2k containing loser tree nodes.

    static uint64_t next_pow2(uint64_t x) {
        return x == 1 ? 1 : uint64_t(1) << (sizeof(unsigned long long) * 8 - __builtin_clzll(x - 1));
    }

    /** Called recursively to build the initial tree. */
    Source init_winner(const Source &root) {
        if (root >= k)
            return root;

        auto left = init_winner(2 * root);
        auto right = init_winner(2 * root + 1);
        if (losers[right].key >= losers[left].key) {
            losers[root] = losers[right];
            return left;
        } else {
            losers[root] = losers[left];
            return right;
        }
    }

public:

    LoserTree() = default;

    explicit LoserTree(const Source &ik) : k(next_pow2(ik)), losers(2 * k) {
        for (auto i = ik - 1u; i < k; ++i) {
            losers[i + k].key = std::numeric_limits<T>::max();
            losers[i + k].source = std::numeric_limits<Source>::max();
        }
    }

    /** Returns the index of the sequence with the smallest element. */
    Source min_source() const {
        assert(losers[0].source != std::numeric_limits<Source>::max());
        return losers[0].source;
    }

    /** Inserts the initial element of the sequence source. */
    void insert_start(const T *key_ptr, const Source &source) {
        Source pos = k + source;
        assert(pos < losers.size());
        losers[pos].source = source;
        losers[pos].key = *key_ptr;
    }

    /** Deletes the smallest element and insert a new element in its place. */
    void delete_min_insert(const T *key_ptr) {
        auto source = losers[0].source;
        auto key = key_ptr ? *key_ptr : std::numeric_limits<T>::max();

        for (auto pos = (k + source) / 2; pos > 0; pos /= 2) {
            if (losers[pos].key < key || (key >= losers[pos].key && losers[pos].source < source)) {
                std::swap(losers[pos].source, source);
                std::swap(losers[pos].key, key);
            }
        }

        losers[0].source = source;
        losers[0].key = key;
    }

    /** Initializes the tree. */
    void init() { losers[0] = losers[init_winner(1)]; }
};

} // namespace internal

template<typename K, typename V, typename PGMType>
class DynamicPGMIndex<K, V, PGMType>::Iterator {
    friend class DynamicPGMIndex;

    using level_iterator = typename Level::const_iterator;
    using dynamic_pgm_type = DynamicPGMIndex<K, V, PGMType>;

    struct Cursor {
        uint8_t level_number;
        level_iterator iterator;
        Cursor() = default;
        Cursor(uint8_t level_number, const level_iterator iterator) : level_number(level_number), iterator(iterator) {}
    };

    const dynamic_pgm_type *super;  ///< Pointer to the container that is being iterated.
    Cursor current;                 ///< Pair (level number, iterator to the current element).
    bool initialized;               ///< true iff the members tree and iterators have been initialized.
    uint8_t unconsumed_count;       ///< Number of iterators that have not yet reached the end.
    internal::LoserTree<K> tree;    ///< Tournament tree with one leaf for each iterator.
    std::vector<Cursor> iterators;  ///< Vector with pairs (level number, iterator).

    void lazy_initialize() {
        if (initialized)
            return;

        // For each level create and position an iterator to the first key > current
        iterators.reserve(super->used_levels - super->min_level);
        for (uint8_t i = super->min_level; i < super->used_levels; ++i) {
            auto &level = super->level(i);
            if (level.empty())
                continue;

            size_t lo = 0;
            size_t hi = level.size();
            if (super->has_pgm(i)) {
                auto range = super->pgm(i).search(current.iterator->first);
                lo = range.lo;
                hi = range.hi;
            }

            auto pos = std::upper_bound(level.begin() + lo, level.begin() + hi, current.iterator->first);
            if (pos != level.end())
                iterators.emplace_back(i, pos);
        }

        tree = decltype(tree)(iterators.size());
        for (size_t i = 0; i < iterators.size(); ++i)
            tree.insert_start(&iterators[i].iterator->first, i);
        tree.init();

        initialized = true;
        unconsumed_count = iterators.size();
    }

    void advance() {
        if (unconsumed_count == 0) {
            *this = super->end();
            return;
        }

        auto step = [&] {
            auto &it_min = iterators[tree.min_source()];
            auto level_number = it_min.level_number;
            auto result = it_min.iterator;
            ++it_min.iterator;
            if (it_min.iterator == super->level(level_number).end()) {
                tree.delete_min_insert(nullptr);
                --unconsumed_count;
            } else
                tree.delete_min_insert(&it_min.iterator->first);
            return Cursor(level_number, result);
        };

        Cursor tmp;
        do {
            tmp = step();
            while (unconsumed_count > 0 && iterators[tree.min_source()].iterator->first == tmp.iterator->first)
                step();
        } while (unconsumed_count > 0 && tmp.iterator->deleted());

        if (tmp.iterator->deleted())
            *this = super->end();
        else
            current = tmp;
    }

    Iterator(const dynamic_pgm_type *p, uint8_t level_number, const level_iterator it)
        : super(p), current(level_number, it), initialized(), unconsumed_count(), tree(), iterators() {};

public:

    using difference_type = typename decltype(levels)::difference_type;
    using value_type = const Item;
    using pointer = const Item *;
    using reference = const Item &;
    using iterator_category = std::forward_iterator_tag;

    Iterator &operator++() {
        lazy_initialize();
        advance();
        return *this;
    }

    Iterator operator++(int) {
        Iterator i(current);
        ++i;
        return i;
    }

    reference operator*() const { return *current.iterator; }
    pointer operator->() const { return &*current.iterator; };

    bool operator==(const Iterator &rhs) const {
        return current.level_number == rhs.current.level_number && current.iterator == rhs.current.iterator;
    }

    bool operator!=(const Iterator &rhs) const { return !(*this == rhs); }
};

#pragma pack(push, 1)

template<typename K, typename V, typename PGMType>
class DynamicPGMIndex<K, V, PGMType>::ItemA {
    const static V tombstone;

    template<typename T = V, std::enable_if_t<std::is_pointer_v<T>, int> = 0>
    static V get_tombstone() { return new std::remove_pointer_t<V>(); }

    template<typename T = V, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
    static V get_tombstone() { return std::numeric_limits<V>::max(); }

public:
    K first;
    V second;

    ItemA() { /* do not (default-)initialize for a more efficient std::vector<ItemA>::resize */ }
    explicit ItemA(const K &key) : first(key), second(tombstone) {}
    explicit ItemA(const K &key, const V &value) : first(key), second(value) {
        if (second == tombstone)
            throw std::invalid_argument("The specified value is reserved and cannot be used.");
    }

    operator K() const { return first; }
    bool deleted() const { return this->second == tombstone; }
};

template<typename K, typename V, typename PGMType>
const V DynamicPGMIndex<K, V, PGMType>::ItemA::tombstone = get_tombstone<V>();

template<typename K, typename V, typename PGMType>
class DynamicPGMIndex<K, V, PGMType>::ItemB {
    bool flag;

public:
    K first;
    V second;

    ItemB() = default;
    explicit ItemB(const K &key) : flag(true), first(key), second() {}
    explicit ItemB(const K &key, const V &value) : flag(false), first(key), second(value) {}

    operator K() const { return first; }
    bool deleted() const { return flag; }
};

#pragma pack(pop)

}