// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
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

#include "piecewise_linear_model.hpp"
#include "storage_utility.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cstring>
namespace pgm {

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

/**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */
struct ApproxPos {
    size_t pos; ///< The approximate position of the key.
    size_t lo;  ///< The lower bound of the range.
    size_t hi;  ///< The upper bound of the range.
};

/**
 * A space-efficient index that enables fast search operations on a sorted sequence of @c n numbers.
 *
 * A search returns a struct @ref ApproxPos containing an approximate position of the sought key in the sequence and
 * the bounds of a range of size 2*Epsilon+1 where the sought key is guaranteed to be found if present.
 * If the key is not present, the range is guaranteed to contain a key that is not less than (i.e. greater or equal to)
 * the sought key, or @c n if no such key is found.
 * In the case of repeated keys, the index finds the position of the first occurrence of a key.
 *
 * The @p Epsilon template parameter should be set according to the desired space-time trade-off. A smaller value
 * makes the estimation more precise and the range smaller but at the cost of increased space usage.
 *
 * Internally the index uses a succinct piecewise linear mapping from keys to their position in the sorted order.
 * This mapping is represented as a sequence of linear models (segments) which, if @p EpsilonRecursive is not zero, are
 * themselves recursively indexed by other piecewise linear mappings.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam EpsilonRecursive controls the size of the search range in the internal structure
 * @tparam Floating the floating-point type to use for slopes
 */
template<typename K, size_t Epsilon = 64, size_t EpsilonRecursive = 4, typename Floating = float>
class PGMIndex {
protected:
    template<typename, size_t, size_t, uint8_t, typename>
    friend class BucketingPGMIndex;

    template<typename, size_t, typename>
    friend class EliasFanoPGMIndex;

    static_assert(Epsilon > 0);
    struct Segment;

    size_t n;                           ///< The number of elements this index was built on.
    K first_key;                        ///< The smallest element.
    std::vector<Segment> segments;      ///< The segments composing the index.
    std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[], in reverse order.

    long SegmentSize = sizeof(Segment);
    long SegmentCountPerBlock = long(BLOCK_SIZE/SegmentSize);
    // todo: implement a disk version
    template<typename RandomIt>
    static void build(RandomIt first, RandomIt last,
                      size_t epsilon, size_t epsilon_recursive,
                      std::vector<Segment> &segments,
                      std::vector<size_t> &levels_offsets) {
        auto n = (size_t) std::distance(first, last);
        if (n == 0)
            return;

        levels_offsets.push_back(0);
        segments.reserve(n / (epsilon * epsilon));

        auto ignore_last = *std::prev(last) == std::numeric_limits<K>::max(); // max() is the sentinel value
        auto last_n = n - ignore_last;
        last -= ignore_last;

        auto build_level = [&](auto epsilon, auto in_fun, auto out_fun) {
            auto n_segments = internal::make_segmentation_par(last_n, epsilon, in_fun, out_fun);
            std::cout << "n-" << n_segments << std::endl;
            if (segments.back().slope == 0 && last_n > 1) {
                // Here we need to ensure that keys > *(last-1) are approximated to a position == prev_level_size
                segments.emplace_back(*std::prev(last) + 1, 0, last_n);
                ++n_segments;
            }
            segments.emplace_back(last_n); // Add the sentinel segment
            return n_segments;
        };

        // Build first level
        auto in_fun = [&](auto i) {
            auto x = first[i];
            // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
            // to x=first[i] such that x+1!=first[i+1], we map the values x+1,...,first[i+1]-1 to their correct rank i
            auto flag = i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1];
            return std::pair<K, size_t>(x + flag, i);
        };
        auto out_fun = [&](auto cs) { segments.emplace_back(cs); };
        last_n = build_level(epsilon, in_fun, out_fun);
        levels_offsets.push_back(levels_offsets.back() + last_n + 1);

        // Build upper levels
        while (epsilon_recursive && last_n > 1) {
            auto offset = levels_offsets[levels_offsets.size() - 2];
            auto in_fun_rec = [&](auto i) { return std::pair<K, size_t>(segments[offset + i].key, i); };
            last_n = build_level(epsilon_recursive, in_fun_rec, out_fun);
            levels_offsets.push_back(levels_offsets.back() + last_n + 1);
        }
    }

    /**
     * Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key.
     * @param key the value of the element to search for
     * @return an iterator to the segment responsible for the given key
     */
    auto segment_for_key(const K &key) const {
        if constexpr (EpsilonRecursive == 0) {
            return std::prev(std::upper_bound(segments.begin(), segments.begin() + segments_count(), key));
        }

        auto it = segments.begin() + *(levels_offsets.end() - 2);
        for (auto l = int(height()) - 2; l >= 0; --l) {
            auto level_begin = segments.begin() + levels_offsets[l];
            auto pos = std::min<size_t>((*it)(key), std::next(it)->intercept);
            auto lo = level_begin + PGM_SUB_EPS(pos, EpsilonRecursive + 1);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Segment);
            if constexpr (EpsilonRecursive <= linear_search_threshold) {
                for (; std::next(lo)->key <= key; ++lo)
                    continue;
                it = lo;
            } else {
                auto level_size = levels_offsets[l + 1] - levels_offsets[l] - 1;
                auto hi = level_begin + PGM_ADD_EPS(pos, EpsilonRecursive, level_size);
                it = std::prev(std::upper_bound(lo, hi, key));
            }
        }
        return it;
    }

    int32_t ACCESSED_BLOCK_COUNT_Index = 0;
    void ADD_BLOCK_COUNT_Index() {ACCESSED_BLOCK_COUNT_Index += 1;}
    #define BLOCK_RECORDING 1

    void read_block_or_not(FILE * file_handler, int32_t *last_block, int32_t block_to_fetch, char *block_data) {
        if (block_to_fetch != *last_block) {
            read_block(file_handler, block_data, block_to_fetch);
            *last_block = block_to_fetch;
            #if BLOCK_RECORDING
            ADD_BLOCK_COUNT_Index();
            #endif
        }
    }

    void item_block_offset(size_t it, int *block_to_fetch, int *item_offset) {
        *block_to_fetch = int(it / SegmentCountPerBlock) + 1;
        *item_offset = it % SegmentCountPerBlock;
    }

    size_t predict_pos(size_t it, FILE* file_handler, int *last_block, char *block_data, K key) {
        int block_to_fetch = -1;
        int item_offset = -1;
        item_block_offset(it, &block_to_fetch, &item_offset);
        read_block_or_not(file_handler, last_block, block_to_fetch, block_data);
//        if (ACCESSED_BLOCK_COUNT_Index == 4) {
//            int x = 0;
//        }
        Segment *it_seg = (Segment *)(block_data + sizeof(Segment) * item_offset);
        size_t est_pos = (*it_seg)(key);

        int next_offset = item_offset + 1;
        if (item_offset == SegmentCountPerBlock - 1) {
            next_offset = 0;
            block_to_fetch = *last_block + 1;
            read_block_or_not(file_handler, last_block, block_to_fetch, block_data);
        }
        Segment *it_next_seg = (Segment *)(block_data + sizeof(Segment) * next_offset);

        size_t pos = std::min<size_t>(est_pos, it_next_seg->intercept);
        return pos;
    }

    K obtain_key_from_disk(size_t it, FILE *file_handler, int *last_block, char *block_data) {
        int block_to_fetch = -1;
        int item_offset = -1;
        item_block_offset(it, &block_to_fetch, &item_offset);
        read_block_or_not(file_handler, last_block, block_to_fetch, block_data);
        Segment *_it_seg = (Segment *)(block_data + sizeof(Segment) * item_offset);
        return _it_seg->key;
    }

    size_t segment_for_key_disk(FILE *file_handler, const K &key, int *lc) {
        // do not consider the case `EpsilonRecursive == 0'
        int last_block = -1;
        char *block_data = new char[BLOCK_SIZE];
        size_t it = *(levels_offsets.end() - 2);
        for (auto l = int(height() -2); l >= 0; l--) {
            *lc += 1;
            auto pos = predict_pos(it, file_handler, &last_block, block_data, key);
            size_t level_begin = levels_offsets[l];
            size_t lo = level_begin + PGM_SUB_EPS(pos, EpsilonRecursive + 1);

            static constexpr size_t linear_search_threshold = 8 * 64 / sizeof(Segment);
            if constexpr (EpsilonRecursive <= linear_search_threshold) {
                for (;;++lo) {
                    auto lo_next_key = obtain_key_from_disk(lo + 1, file_handler, &last_block, block_data);
                    if (lo_next_key> key) break;
                }
                it = lo;
            } else {
                size_t level_size = levels_offsets[l+1] - levels_offsets[l] - 1;
                size_t hi = level_begin + PGM_ADD_EPS(pos, EpsilonRecursive, level_size);
                size_t count = hi - lo;
                // https://cplusplus.com/reference/algorithm/upper_bound/
                while (count > 0) {
                    auto step = count / 2;
                    auto _it_ = lo + step;
                    auto _it_key = obtain_key_from_disk(_it_, file_handler, &last_block, block_data);
                    if (!(key < _it_key)) {
                        lo = _it_ + 1;
                        count -= step + 1;
                    } else count = step;
                }
                it = lo - 1; //prev()
            }
        }
        auto pos = predict_pos(it, file_handler, &last_block, block_data, key);
        delete []block_data;
        return pos;
    }

public:

    static constexpr size_t epsilon_value = Epsilon;

    /**
     * Constructs an empty index.
     */
    PGMIndex() = default;

    /**
     * Constructs the index on the given sorted vector.
     * @param data the vector of keys to be indexed, must be sorted
     */
    explicit PGMIndex(const std::vector<K> &data) : PGMIndex(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted keys in the range [first, last).
     * @param first, last the range containing the sorted keys to be indexed
     */
    template<typename RandomIt>
    PGMIndex(RandomIt first, RandomIt last)
        : n(std::distance(first, last)),
          first_key(n ? *first : K(0)),
          segments(),
          levels_offsets() {
        build(first, last, Epsilon, EpsilonRecursive, segments, levels_offsets);
    }

    template<typename RandomIt>
    PGMIndex(RandomIt first, RandomIt last, bool _inner_disk, FILE *file_handler)
        : n(std::distance(first, last)),
        first_key(n ? *first : K(0)),
        segments(),
        levels_offsets() {
        build(first, last, Epsilon, EpsilonRecursive, segments, levels_offsets);
        if (_inner_disk) {
            write_file(file_handler);
        }
    }

    PGMIndex(FILE *file_handler, bool inner_disk) {
        load_file(file_handler, inner_disk);
    }

    void load_file(FILE *file_handler, bool inner_disk) {
        char *block_data = new char[BLOCK_SIZE];
        read_block(file_handler, block_data, 0);
        int offset = 0;
        memcpy(&first_key, block_data + offset, sizeof(K));
        offset += sizeof(K);
        memcpy(&n, block_data + offset, sizeof(size_t));
        offset += sizeof(size_t);
        size_t level_size;
        memcpy(&level_size, block_data + offset, sizeof(size_t));
        offset += sizeof(size_t);
        size_t seg_size;
        memcpy(&seg_size, block_data + offset, sizeof(size_t));
        offset += sizeof(size_t);
        levels_offsets = decltype(levels_offsets)(level_size);
        size_t *p = (size_t *)(block_data + offset);
        for (int i = 0; i < level_size; i++) {
            levels_offsets[i] = p[i];
        }

        if (!inner_disk) {
            segments = decltype(segments)(seg_size);
            int block_count = seg_size / SegmentCountPerBlock;
            int item_count_last = seg_size % SegmentCountPerBlock;
            int cp_offset = 0;
            for (auto i = 0; i < block_count; i++) {
                read_block(file_handler, block_data, i+1);
                int cp_size = SegmentCountPerBlock;
                if (i == block_count - 1) cp_size = item_count_last;
                Segment *p = (Segment *)(block_data);
                for (auto j = 0; j < cp_size; j++, cp_size++) {
                    segments[cp_offset] = p[j];
                }
            }
        }
    }
    /**
     * Returns the approximate position and the range where @p key can be found.
     * @param key the value of the element to search for
     * @return a struct with the approximate position and bounds of the range
     */

    ApproxPos search(const K &key) const {
        auto k = std::max(first_key, key);
        auto it = segment_for_key(k);
        auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
        auto lo = PGM_SUB_EPS(pos, Epsilon);
        auto hi = PGM_ADD_EPS(pos, Epsilon, n);
        return {pos, lo, hi};
    }

    // todo: implement a disk version
    ApproxPos search_disk(FILE *file_handler, const K&key, int *c, int *lc) {
        auto k = std::max(first_key, key);
        // on disk version, we calculate the pos
        ACCESSED_BLOCK_COUNT_Index = 0;
        auto pos = segment_for_key_disk(file_handler, key, lc);
        *c = ACCESSED_BLOCK_COUNT_Index;
        return {pos, PGM_SUB_EPS(pos, Epsilon), PGM_ADD_EPS(pos, Epsilon, n)};
    }

    void write_file(FILE *file) {
        int32_t block = 0;
        char *block_data = new char[BLOCK_SIZE];
        memcpy(block_data, &first_key, sizeof(K));
        // write offsets info
        size_t *pointer = (size_t *)(block_data+sizeof(K));
        size_t offset_size = levels_offsets.size();
        pointer[0] = n;
        pointer[1] = offset_size;
        pointer[2] = (size_t)(segments.size());
        for (auto i = 0; i < offset_size; i++) {
            pointer[i+3] = levels_offsets[i];
        }
        write_block(file, block_data, block);
        block += 1;

        // write segments info
        int32_t start = 0;
        Segment *items = (Segment *)block_data;
        for (;start < segments.size();) {
            int32_t end = std::min<size_t>(start + SegmentCountPerBlock, segments.size());
            for (auto i = start; i < end; i++) {
                items[i-start] = segments[i];
            }
            write_block(file, block_data, block);
            start += SegmentCountPerBlock;
            block += 1;
        }
        delete []block_data;
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t segments_count() const { return segments.empty() ? 0 : levels_offsets[1] - 1; }

    /**
     * Returns the number of levels of the index.
     * @return the number of levels of the index
     */
    size_t height() const { return levels_offsets.size() - 1; }

    /**
     * Returns the size of the index in bytes.
     * @return the size of the index in bytes
     */
    size_t size_in_bytes() const { return segments.size() * sizeof(Segment) + levels_offsets.size() * sizeof(size_t); }
};

#pragma pack(push, 1)

template<typename K, size_t Epsilon, size_t EpsilonRecursive, typename Floating>
struct PGMIndex<K, Epsilon, EpsilonRecursive, Floating>::Segment {
    K key;             ///< The first key that the segment indexes.
    Floating slope;    ///< The slope of the segment.
    int32_t intercept; ///< The intercept of the segment.

    Segment() = default;

    Segment(K key, Floating slope, int32_t intercept) : key(key), slope(slope), intercept(intercept) {};

    explicit Segment(size_t n) : key(std::numeric_limits<K>::max()), slope(), intercept(n) {};

    explicit Segment(const typename internal::OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment &cs)
        : key(cs.get_first_x()) {
        auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
        if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
            throw std::overflow_error("Change the type of Segment::intercept to int64");
        slope = cs_slope;
        intercept = cs_intercept;
    }

    friend inline bool operator<(const Segment &s, const K &k) { return s.key < k; }
    friend inline bool operator<(const K &k, const Segment &s) { return k < s.key; }
    friend inline bool operator<(const Segment &s, const Segment &t) { return s.key < t.key; }

    operator K() { return key; };

    /**
     * Returns the approximate position of the specified key.
     * @param k the key whose position must be approximated
     * @return the approximate position of the specified key
     */
    inline size_t operator()(const K &k) const {
        auto pos = int64_t(slope * (k - key)) + intercept;
        return pos > 0 ? size_t(pos) : 0ull;
    }
};

#pragma pack(pop)

}
