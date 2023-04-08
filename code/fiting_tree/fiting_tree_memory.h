#include "limits"
#include <vector>
#include <iostream>
#include "storage_management.h"
#include <algorithm>
#include <cstring>
#include "pgm_seg.h"
#include "stx/btree_multimap.h"

#define FInnerNodeType 1
#define FLeafNodeType 0
#define Disk 1
#define Binary 0
#define BufferSize 256
const long MaxNodeSize = BlockSize;

#define PGMSegment 1

typedef struct {
    KeyType key;
    ValueType value;
} Iterm;
bool ItermCompare(Iterm a, Iterm b){
	if(a.key < b.key)
		return 1;
	else 
		return 0;
}
#define ItermSize sizeof(Iterm)
#define ItermCountOneBlock (MaxNodeSize/ItermSize)

class Segment {
    Iterm *data;
    int count;
    public:
    Segment(Iterm *_data,int _count ) {
        data = _data;
        count = _count;
    }
    ~Segment() {
    }
    bool search(int start, int end, KeyType key, ValueType *v) {
        bool found = false;
        int l = start;
        int r = end;
        while (l <= r) {
            int mid = l + (r - l) / 2;
            if (data[mid].key == key) {
                found = true;
                *v = data[mid].value;
                return found;
            } else if (data[mid].key < key) l = mid + 1;
            else r = mid - 1;
        }
        return found;
    }
    int get_size() {
        return count;
    }

};

typedef struct FLeafNodeIterm {
    KeyType key;
    double slope;
    Segment *segment;
} FLeafNodeIterm;
#define FLeafNodeItermSize sizeof(LeafNodeIterm)
#define FMaxLeafNodeItemCount (MaxNodeSize/FLeafNodeItermSize) 

#if Disk
    typedef struct {
        KeyType key;
        double slope;
        int block_id;
        int iterm_count;
        int added_in_buffer;
        #if PGMSegment
        int32_t intercept;
        #endif
    } LeafNodeItermOnDisk;
    #define LeafNodeItermSizeOnDisk sizeof(LeafNodeItermOnDisk)


    // used for when we store all inner nodes in main memory
    typedef struct {
        double slope;
        int block_id;
        int iterm_count;
        int added_in_buffer;
    #if PGMSegment
        int32_t intercept;
    #endif
    } LeafNodeValueOnDisk;
    #define LeafNodeValueOnDiskSize sizeof(LeafNodeValueOnDisk)
    #define MaxItemLeafValueCount (BlockSize/LeafNodeValueOnDiskSize)
    #define MaxItemInnerNodeLeafValue (BlockSize/(sizeof(KeyType) + sizeof(void *)))
    struct traits_inner : stx::btree_default_map_traits<KeyType, LeafNodeValueOnDisk>
    {
        static const bool       selfverify = false;
        static const bool       debug = false;

        static const int        leafslots = MaxItemLeafValueCount;
        static const int        innerslots = MaxItemInnerNodeLeafValue;
    };
    typedef stx::btree_multimap<KeyType, LeafNodeValueOnDisk, std::less<KeyType>, traits_inner> stx_btree;


    typedef struct {
        char tag;
        int item_count;
    } LeaftNodeHeaderOnDisk;
    #define LeaftNodeHeaderSizeOnDisk sizeof(LeaftNodeHeaderOnDisk)
    #define MaxLeafNodeItemCountOnDisk ((MaxNodeSize - LeaftNodeHeaderSizeOnDisk)/LeafNodeItermSizeOnDisk)

    typedef struct {
        KeyType key;
        int block_id;
    } InnerNodeItermOnDisk;
    #define InnerNodeItermSizeOnDisk sizeof(InnerNodeItermOnDisk)

    typedef struct {
        char tag;
        int item_count;
        int left_child;
    } InnerNodeHeaderOnDisk;
    #define InnerNodeHeaderSizeOnDisk sizeof(InnerNodeHeaderOnDisk)
    #define MaxInnerNodeItemCountOnDisk ((MaxNodeSize - InnerNodeHeaderSizeOnDisk)/InnerNodeItermSize)

    typedef struct {
        KeyType key;
        int status;
        int added_block;
        int level;
        InnerNodeItermOnDisk inid;
    } FBuildStatus;  
#endif

typedef struct {
    char tag;
    int item_count;
    FLeafNodeIterm items[FMaxLeafNodeItemCount];
} FLeafNode;
#define FLeafNodeSize sizeof(FLeafNode)

typedef struct {
    KeyType key;
    void *address;
} FInnerNodeIterm;
#define FInnerNodeItermSize sizeof(FInnerNodeIterm)
#define FMaxInnerNodeItemCount (MaxNodeSize/FInnerNodeItermSize)

typedef struct {
    char tag;
    int item_count;
    FInnerNodeIterm items[FMaxInnerNodeItemCount];
} FInnerNode;
#define FInnerNodeSize sizeof(FInnerNode)

class FITingTree {
    private:
        #if Disk
        StorageManager *sm;
        MetaNode metanode;
        #define ALL_DISK 0
        #define LEAF_DISK 1
        int hybrid_mode = ALL_DISK;
        void load_metanode() {
            Block block = sm->get_block(0);
            memcpy(&metanode, block.data, MetaNodeSize);
        }

        void sync_metanode() {
            Block block;
            memcpy(block.data, &metanode, MetaNodeSize);
            sm->write_block(0, block);
        }
        stx_btree inner_btree;
        #endif
        int error_bound = 32;
        std::vector<Segment*> segments;
        void *root_node = nullptr;

    public:
    FITingTree() = default;
    FITingTree(int error) {
        error_bound = error;
    }
    #if Disk
    FITingTree(int error, char *file_name, bool first) {
        sm = new StorageManager(first, file_name);
        error_bound = error;
        load_metanode();
    }

    FITingTree(int error, char *file_name, bool first, int _hybrid_mode) {
        sm = new StorageManager(first, file_name);
        error_bound = error;
        load_metanode();
        hybrid_mode = _hybrid_mode;
    }
    #endif
    ~FITingTree() {
        free(root_node);
        delete sm;
    }
    // Suppose data is sorted

    class Seg {
        public:
            KeyType key;             ///< The first key that the segment indexes.
            float slope;    ///< The slope of the segment.
            int32_t intercept; ///< The intercept of the segment.
            int32_t start;
            int32_t end;
        
        Seg(KeyType key, float slope, int32_t intercept) : key(key), slope(slope), intercept(intercept) {};

        explicit Seg(size_t n) : key(std::numeric_limits<KeyType>::max()), slope(), intercept(n) {};

        explicit Seg(const typename OptimalPiecewiseLinearModel<KeyType, size_t>::CanonicalSegment &cs)
            : key(cs.get_first_x()) {
            // auto[cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
            auto m = cs.get_floating_point_segment(key);
            auto cs_slope = m.first;
            auto cs_intercept = m.second - cs.get_start_index();
            start = cs.get_start_index();
            end = cs.get_end_index();
            if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
                throw std::overflow_error("Change the type of Segment::intercept to int64");
            slope = cs_slope;
            intercept = cs_intercept;
        }
    };

    template<typename RandomIt>
    void _ts(RandomIt first, RandomIt last, size_t epsilon, std::vector<Seg> &segments) {
        auto n = (size_t) std::distance(first, last);
        auto ignore_last = *std::prev(last) == std::numeric_limits<KeyType>::max(); // max() is the sentinel value
            auto last_n = n - ignore_last;
            last -= ignore_last;
            // std::cout << last_n << std::endl;
            // std::cout << *last << std::endl;
            auto build_level = [&](auto epsilon, auto in_fun, auto out_fun) {
                auto n_segments = make_segmentation_par(last_n, epsilon, in_fun, out_fun);
                // if (segments.back().slope == 0 && last_n > 1) {
                //     // Here we need to ensure that keys > *(last-1) are approximated to a position == prev_level_size
                //     segments.emplace_back(*std::prev(last) + 1, 0, last_n);
                //     ++n_segments;
                // }
                // segments.emplace_back(last_n); // Add the sentinel segment
                return n_segments;
            };

            // Build first level
            auto in_fun = [&](auto i) {
                auto x = first[i];
                // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
                // to x=first[i] such that x+1!=first[i+1], we map the values x+1,...,first[i+1]-1 to their correct rank i
                auto flag = i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1];
                return std::pair<KeyType, size_t>(x + flag, i);
            };
            auto out_fun = [&](auto cs) { segments.emplace_back(cs); };
            last_n = build_level(epsilon, in_fun, out_fun);
            // std::cout << last_n << std::endl;
            // std::cout << "segment size:" << segments.size() << std::endl;
            // std::cout << segments.at(segments.size()-1).start<< "," << segments.at(segments.size()-1).end << std::endl;
    }

    #if Disk
    std::vector<LeafNodeItermOnDisk> write_segment_pgm(Iterm *data, long item_count, std::vector<Seg> segments, bool is_first, int last_block, int next_block) {
        std::vector<LeafNodeItermOnDisk> leaf_iterms;
        if (is_first) {
            Iterm _data; 
            _data.key = KeyType(0);
            _data.value = ValueType(metanode.block_count + 1);
            
            char _data_[BlockSize];
            memcpy(_data_, &_data, ItermSize);
            sm->write_with_size(metanode.block_count, _data_, BlockSize);
            metanode.block_id = metanode.block_count;
            last_block = metanode.block_count;
            metanode.block_count += 1;
            metanode.original_key = segments.at(0).key;
            metanode.added_in_buffer = 0;
            // is_first =false;
        }
        long current_pos = 0;
//        long total_write_size = 0;
//        long total_offset = 0;
//        long total_blocks = 0;
        for (int i = 0; i < segments.size(); i++) {
            Seg _seg = segments.at(i);
            int copy_size = _seg.end - _seg.start + 1;
            Iterm *_data = new Iterm[2 +  copy_size + BufferSize]; //(Iterm *)malloc(copy_size * ItermSize);
                // add one more block for each segment for the new insertion.
            memcpy(_data + 2, data + current_pos, copy_size * ItermSize);
            current_pos += copy_size; 

            long write_size = (2 + copy_size + BufferSize) * ItermSize;
//            total_write_size += write_size;
            int wroten = write_size / MaxNodeSize;
            int offset = write_size % MaxNodeSize;
            if (offset > 0) {
                wroten += 1;

//                total_blocks += wroten;
            }
            _data[0].key = KeyType(last_block);
            _data[0].value = i == (segments.size() - 1)? ValueType(next_block) : ValueType(metanode.block_count + wroten);
            _data[1].key = KeyType(copy_size);
            sm->write_with_size(metanode.block_count, _data, write_size);

            delete[] _data;
            if (offset > 0) {
                char empty[MaxNodeSize - offset];
                long start = metanode.block_count * MaxNodeSize + write_size;
                sm->write_arbitrary(start, empty, MaxNodeSize - offset);

//                total_offset += (MaxNodeSize - offset);
            }

            //update with new links
            if (!is_first && i == 0 && last_block != 0) {
                char _temp[MaxNodeSize];
                sm->get_block(last_block, _temp);
                Iterm *_data = (Iterm *)_temp;
                _data[0].value = ValueType(metanode.block_count);
                sm->write_with_size(last_block, _temp, MaxNodeSize);
            }
            // must move it here ...
            last_block = metanode.block_count;

            if (!is_first && i == (segments.size() - 1) && next_block != 0) {
                char _temp[MaxNodeSize];
                sm->get_block(next_block, _temp);
                Iterm *_data = (Iterm *)_temp;
                _data[0].key = KeyType(metanode.block_count);
                sm->write_with_size(next_block, _temp, MaxNodeSize);
            }

            LeafNodeItermOnDisk lnid;
            lnid.key = _seg.key;
            lnid.slope = _seg.slope;
            lnid.iterm_count = copy_size;
            lnid.block_id = metanode.block_count;
            lnid.added_in_buffer = 0;
            lnid.intercept = _seg.intercept;
            leaf_iterms.push_back(lnid);
            metanode.block_count += wroten;
            sync_metanode();
        }
//        std::cout << total_write_size << "," << total_blocks << "," << total_offset << std::endl;
        return leaf_iterms;
    }


    std::vector<LeafNodeItermOnDisk> write_segment(Iterm *data, long item_count, bool is_first) {
        std::vector<LeafNodeItermOnDisk> leaf_iterms;
    #else
    std::vector<FLeafNodeIterm> write_segment(Iterm *data, long item_count) {
         std::vector<FLeafNodeIterm> leaf_iterms;
    #endif
        double high_slope = std::numeric_limits<double>::max();
        double low_slope = 0;
        KeyType original_key = data[0].key;
        double original_location = 0;
        // bool is_first = true;
        long i = 1;
        for (; i < item_count; i++) {
            KeyType temp_key = data[i].key;
            double temp_slope = (i - original_location) / (temp_key - original_key);
            if (temp_slope <= high_slope && temp_slope >= low_slope) {
                double temp_high_slope = ((i + error_bound) - original_location) / (temp_key - original_key);
                double temp_low_slope = ((i - error_bound) - original_location) / (temp_key - original_key);
                high_slope = std::min<double>(temp_high_slope, high_slope);
                low_slope = std::max<double>(temp_low_slope, low_slope);
            } else {
                double slope = (high_slope + low_slope) / 2;
                long copy_size = i - original_location;
                if (copy_size == 1) slope = 1;
                Iterm *_data = new Iterm[1+ copy_size + BufferSize]; //(Iterm *)malloc(copy_size * ItermSize);
                // add one more block for each segment for the new insertion.
                memcpy(_data+1, data + int(original_location), copy_size * ItermSize);

                #if Disk
                    //metanode.block_count
                    
                    long write_size = (1 + copy_size + BufferSize) * ItermSize;
                    int wroten = write_size / MaxNodeSize;
                    int offset = write_size % MaxNodeSize;
                    if (offset > 0) {
                        wroten += 1;
                    }
                    _data[0].key = (KeyType) (metanode.block_count + wroten);
                    _data[0].value = (ValueType) (copy_size);
                    sm->write_with_size(metanode.block_count, _data, write_size);
                    delete[] _data;
                    if (offset > 0) {
                        char empty[MaxNodeSize - offset];
                        long start = metanode.block_count * MaxNodeSize + write_size;
                        sm->write_arbitrary(start, empty, MaxNodeSize - offset);
                    }
                    LeafNodeItermOnDisk lnid;
                    lnid.key = original_key;
                    lnid.slope = slope;
                    lnid.iterm_count = copy_size;
                    lnid.block_id = metanode.block_count;
                    lnid.added_in_buffer = 0;
                    leaf_iterms.push_back(lnid);
                    if (is_first) {
                        // int _block = (copy_size+1) / ItermCountOneBlock;
                        // int _offset = (copy_size+1) % ItermCountOneBlock;
                        // long offset = (metanode.block_count + _block) * MaxNodeSize + _offset * ItermSize;
                        metanode.block_id = metanode.block_count + 1;
                        metanode.original_key = original_key;
                        metanode.added_in_buffer = 0;
                        is_first =false;
                    }
                    metanode.block_count += wroten;
                    sync_metanode();
                #else
                    Segment *s = new Segment(_data, copy_size);
                    segments.push_back(s); //针票了
                    FLeafNodeIterm lni;
                    lni.key = original_key;
                    lni.slope = slope;
                    // TODO: check here whether we just store index...
                    lni.segment = s;
                    leaf_iterms.push_back(lni);
                #endif
                original_location = i;
                original_key = data[i].key;
                high_slope = std::numeric_limits<double>::max();
                low_slope = 0;
            }
        }
        double slope = (high_slope + low_slope) / 2;
        long copy_size = i - original_location;
        if (copy_size == 1) slope = 1;
        Iterm *_data = new Iterm[copy_size + BufferSize]; //(Iterm *)malloc(copy_size * ItermSize);
        memcpy(_data, data + int(original_location), copy_size * ItermSize);
        #if Disk
            long write_size = (copy_size + BufferSize) * ItermSize;
            sm->write_with_size(metanode.block_count, _data, write_size);
            delete[] _data;
            int wroten = write_size / MaxNodeSize;
            int offset = write_size % MaxNodeSize;
            if (offset > 0) {
                wroten += 1;
                char empty[MaxNodeSize - offset];
                long start = metanode.block_count * MaxNodeSize + write_size;
                sm->write_arbitrary(start, empty, MaxNodeSize - offset);
            }
            LeafNodeItermOnDisk lnid;
            lnid.key = original_key;
            lnid.slope = slope;
            lnid.iterm_count = copy_size;
            lnid.block_id = metanode.block_count;
            lnid.added_in_buffer = 0;
            leaf_iterms.push_back(lnid);
            metanode.block_count += wroten;
            sync_metanode();
        #else
            Segment *s = new Segment(_data, copy_size);
            segments.push_back(s);
            FLeafNodeIterm lni;
            lni.key = original_key;
            lni.slope = slope;
            // TODO: check here whether we just store index...
            lni.segment = s;
            leaf_iterms.push_back(lni);
        #endif
       
        #if Disk
            metanode.level += 1;
        #endif

        return leaf_iterms;
    }

    template<typename RandomIt>
    void obtain_segments_pgm(RandomIt first, RandomIt last, size_t epsilon, std::vector<Seg> &segments) {
        auto n = (size_t) std::distance(first, last);
        auto ignore_last = *std::prev(last) == std::numeric_limits<KeyType>::max(); // max() is the sentinel value
            auto last_n = n - ignore_last;
            last -= ignore_last;
            // std::cout << last_n << std::endl;
            // std::cout << *last << std::endl;
            auto build_level = [&](auto epsilon, auto in_fun, auto out_fun) {
                auto n_segments = make_segmentation_par(last_n, epsilon, in_fun, out_fun);
                // if (segments.back().slope == 0 && last_n > 1) {
                //     // Here we need to ensure that keys > *(last-1) are approximated to a position == prev_level_size
                //     segments.emplace_back(*std::prev(last) + 1, 0, last_n);
                //     ++n_segments;
                // }
                // segments.emplace_back(last_n); // Add the sentinel segment
                return n_segments;
            };

            // Build first level
            auto in_fun = [&](auto i) {
                auto x = first[i];
                // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
                // to x=first[i] such that x+1!=first[i+1], we map the values x+1,...,first[i+1]-1 to their correct rank i
                auto flag = i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1];
                return std::pair<KeyType, size_t>(x + flag, i);
            };
            auto out_fun = [&](auto cs) { segments.emplace_back(cs); };
            last_n = build_level(epsilon, in_fun, out_fun);
            // std::cout << last_n << std::endl;
            // std::cout << segments.size() << std::endl;
            // std::cout << segments.at(segments.size()-1).start<< "," << segments.at(segments.size()-1).end << std::endl;
    }

    template<typename RandomIt>
    void bulk_load_pgm(Iterm *data, long item_count, RandomIt first, RandomIt last, size_t epsilon) {
        std::vector<Seg> segments; 
        _ts(first, last, epsilon, segments);
//        std::cout << segments.size() << std::endl;
        std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment_pgm(data, item_count, segments, true, 0, 0);
//        std::cout << sm->get_file_size() << std::endl;
        if (hybrid_mode == ALL_DISK)
            _build_tree(data, item_count, leaf_iterms);
        else if (hybrid_mode == LEAF_DISK) {
            std::vector<std::pair<KeyType, LeafNodeValueOnDisk> > pairs(leaf_iterms.size());
            for (int i = 0; i < leaf_iterms.size(); i++) {
                pairs[i].first = leaf_iterms[i].key;
                pairs[i].second = LeafNodeValueOnDisk {leaf_iterms[i].slope, leaf_iterms[i].block_id,
                                                       leaf_iterms[i].iterm_count,
                                                       leaf_iterms[i].added_in_buffer, leaf_iterms[i].intercept};
            }
            inner_btree.bulk_load(pairs.begin(), pairs.end());
        }
    }

    void bulk_load(Iterm *data, long item_count) {
        std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment(data, item_count, true);
        _build_tree(data, item_count, leaf_iterms);
    }

    void _build_tree(Iterm *data, long item_count, std::vector<LeafNodeItermOnDisk> leaf_iterms) {
        int MAX_COUNT_LEAF = MaxLeafNodeItemCountOnDisk;
        int MAX_COUNT_INNER = MaxInnerNodeItemCountOnDisk;
        int LHEADER = LeaftNodeHeaderSizeOnDisk;
        int IHEADER = InnerNodeHeaderSizeOnDisk;
        std::cout << "all segments!" << std::endl;
        std::cout << leaf_iterms.size() << std::endl;
        int leaf_node_count = leaf_iterms.size() / MAX_COUNT_LEAF;
        int m = leaf_iterms.size() % MAX_COUNT_LEAF;

        if (m > 0) leaf_node_count += 1;
        int address[leaf_node_count];
        
        KeyType _keys[leaf_node_count];

        int _index = 0;
        for (int i = 0; i < leaf_iterms.size(); ) {
            int end_index = std::min<int>(leaf_iterms.size(), i + MAX_COUNT_LEAF) - 1;
            LeaftNodeHeaderOnDisk lnhd;
            lnhd.item_count = end_index - i + 1;
            lnhd.tag = FLeafNodeType;
            LeafNodeItermOnDisk _lnid[lnhd.item_count];
            for (int j = 0; j < lnhd.item_count; j++) {
                _lnid[j] = leaf_iterms.at(j+i);
            }
            Block block;
            memcpy(block.data, &lnhd, LeaftNodeHeaderSizeOnDisk);
            long _isize = lnhd.item_count*LeafNodeItermSizeOnDisk;
            memcpy(block.data + LeaftNodeHeaderSizeOnDisk, _lnid, _isize);
            sm->write_block(metanode.block_count, block);
            address[_index] = metanode.block_count;
            _keys[_index] = _lnid[0].key;
            metanode.block_count += 1;
            i += MAX_COUNT_LEAF;
            sync_metanode();
            _index += 1;
        }

        metanode.level += 1;

        int c = 0;
        while (_index > 1) {
            metanode.level += 1;
            c = (_index - 1) / MAX_COUNT_INNER;
            m = (_index - 1) % MAX_COUNT_INNER;
            int node_count = m == 0 ? c : c + 1;
            int wroten_item = 1;
            for (int i = 0; i < node_count; i++) {
                int end_index = std::min<int>(_index, wroten_item + MAX_COUNT_INNER) - 1;
                InnerNodeHeaderOnDisk inhd;
                inhd.left_child = i == 0 ? address[0] : -1; // we do not need to record the left sibling here currently
                inhd.item_count = end_index - wroten_item + 1;
                inhd.tag = FInnerNodeType;
                InnerNodeItermOnDisk inid[inhd.item_count];
                for (int j = 0; j < inhd.item_count; j++) {
                    inid[j].block_id = address[wroten_item+j];
                    inid[j].key = _keys[wroten_item+j];
                }
                Block block;
                memcpy(block.data, &inhd, InnerNodeHeaderSizeOnDisk);
                long _isize = inhd.item_count * InnerNodeItermSizeOnDisk;
                memcpy(block.data + InnerNodeHeaderSizeOnDisk, inid, _isize);
                sm->write_block(metanode.block_count, block);
                address[i] = metanode.block_count;
                _keys[i] = inid[0].key;
                metanode.block_count += 1;
                wroten_item += inhd.item_count;
                sync_metanode();
            }
            _index = node_count;
        }
        metanode.root_block_id = address[0];
        sync_metanode();
    }

    bool _lookup(KeyType key, void *address, ValueType *v) {
        if (((FInnerNode*)address)->tag == FInnerNodeType) {
            FInnerNode *temp = (FInnerNode*)address;
            int i = 0;
            // todo binary
            for (; i < temp->item_count-1; i++) {
                if (temp->items[i].key <= key && temp->items[i+1].key > key) break;
            }
            return _lookup(key, temp->items[i].address, v);
        } else {
            FLeafNode *temp = (FLeafNode*)address;
            int i = 0;
             // todo binary
            for (; i < temp->item_count-1; i++) {
                if (temp->items[i].key <= key && temp->items[i+1].key > key) break;
            }
            Segment *s = temp->items[i].segment;
            int pos = int(temp->items[i].slope*(key - temp->items[i].key));
            // binary search
            int start = std::max<int>(pos - error_bound, 0);
            int end = std::min<int>(pos + error_bound, s->get_size() - 1);
            return s->search(start, end, key, v); 
        }
    }

    #if Disk
    bool _get_from_segment_buffer(KeyType key, int block_id, ValueType *v, int *block_count) {
        char *data = new char[MaxNodeSize];
        sm->get_block(block_id, data);
        *block_count += 1;
        Iterm * _temp = (Iterm * )data;
        bool found = false;
        // the first slot records the next block id
        for (int _i = 1; _i < metanode.added_in_buffer+1; _i++) {
            if (_temp[_i].key == key) {
                found = true;
                *v = _temp[_i].value;
                break;
            }
        }
        return found;
    }

    void read_block_or_not(int32_t *last_block, int32_t block_to_fetch, char *block_data, int *block_count) {
        if (block_to_fetch != *last_block) {
            sm->get_block(block_to_fetch, block_data);
            *last_block = block_to_fetch;
            *block_count += 1;
        }
    }

    void scan_in_leaf_node(KeyType key, int start_block, int pos, ValueType *v, int *block_count, int item_count, int item_buffer) {
        int start = std::max<int>(pos - error_bound, 0) + 2;
        int end = std::min<int>(pos + error_bound + 2, item_count) + 2;

        // where to start to scan
//        bool found = false;
        int last_block = -1;
        char *_data2 = new char[MaxNodeSize];
        while (start <= end) {
            int mid = start + (end - start) / 2;
            int _block = mid / ItermCountOneBlock;
            int _offset = mid % ItermCountOneBlock;
            read_block_or_not(&last_block, start_block + _block, _data2, block_count);
            Iterm * _temp = (Iterm * )_data2;
            if (_temp[_offset].key == key) {
//                found = true;
//                *v =_temp[_offset].value;
                break;
            } else if (_temp[_offset].key < key) start = mid + 1;
            else end = mid - 1;
        }
        // start search from mid


    }
    bool search_in_leaf_node(KeyType key, int start_block, int pos, ValueType *v, int *block_count, int item_count, int item_buffer) {
        int start = std::max<int>(pos - error_bound, 0) + 2;
        int end = std::min<int>(pos + error_bound + 2, item_count) + 2;

        bool found = false;
        int last_block = -1;
        char *_data2 = new char[MaxNodeSize];
        while (start <= end) {
            int mid = start + (end - start) / 2;
            int _block = mid / ItermCountOneBlock;
            int _offset = mid % ItermCountOneBlock;
            read_block_or_not(&last_block, start_block + _block, _data2, block_count);
            Iterm * _temp = (Iterm * )_data2;
            if (_temp[_offset].key == key) {
                found = true;
                *v =_temp[_offset].value;
                break;
            } else if (_temp[_offset].key < key) start = mid + 1;
            else end = mid - 1;
        }

        if (!found && item_buffer > 0) {
            int _block = (item_count+2) / ItermCountOneBlock;
            int _offset = (item_count+2) % ItermCountOneBlock;
            long offset = (start_block + _block) * MaxNodeSize + _offset * ItermSize;
            sm->read_block_arbitrary(_data2, offset);
            *block_count += 1;
            Iterm * _temp = (Iterm * )_data2;
            for (int _i = 0; _i < item_buffer; _i++) {
                if (_temp[_i].key == key) {
                    found = true;
                    *v = _temp[_i].value;
                    break;
                }
            }
        }
        delete[] _data2;
        return found;
    }

    bool _lookup_disk(KeyType key, int block_id, ValueType *v, int *block_count) {
        if (key < metanode.original_key) { // must in the first segment buffer if have
            if (metanode.added_in_buffer > 0)
                return _get_from_segment_buffer(key, metanode.block_id, v, block_count);
            else return false;
        }
        char *data = new char[MaxNodeSize];

        sm->get_block(block_id, data);
        *block_count += 1;

        bool is_inner = data[0] == FInnerNodeType;
        while (is_inner) {
            InnerNodeHeaderOnDisk *inhd;
            inhd = (InnerNodeHeaderOnDisk *)data;
            InnerNodeItermOnDisk* inids = (InnerNodeItermOnDisk*) (data + InnerNodeHeaderSizeOnDisk);
            int next_block = 0;
            if (inids[0].key > key) {
                next_block = inhd->left_child;
            } else {
                int i = 0;
                #if Binary
                    int _start = 0;
                    int _end = inhd->item_count - 1;
                    while (_start <= _end) {
                        int mid = _start + (_end - _start) / 2;
                        if (inids[mid].key <= key) _start = mid + 1;
                        else _end = mid - 1;
                    }
                    i = _start  - 1;
                #else
                    for (; i < inhd->item_count - 1; i++) {
                        if (inids[i].key <= key && inids[i+1].key > key) break;
                    }
                #endif
                next_block = inids[i].block_id;
            }
            
            sm->get_block(next_block, data);
            *block_count += 1;

            is_inner = data[0] == FInnerNodeType;
        }
        LeaftNodeHeaderOnDisk *lnhd = (LeaftNodeHeaderOnDisk *)data;
        LeafNodeItermOnDisk *lnids = (LeafNodeItermOnDisk *)(data + LeaftNodeHeaderSizeOnDisk);
        int i = 0;
        #if Binary
            int _start = 0;
            int _end = lnhd->item_count - 1;
            while (_start <= _end) {
                int mid = _start + (_end - _start) / 2;
                if (lnids[mid].key <= key) _start = mid + 1;
                else _end = mid - 1;
            }
            i = _start  - 1;
        #else
            for (; i < lnhd->item_count - 1; i++) {
                if (lnids[i].key <= key && lnids[i+1].key > key) break;
            }
        #endif
        #if PGMSegment
        int pos = int64_t(lnids[i].slope * (key- lnids[i].key)) + lnids[i].intercept;
        #else
        int pos = int(lnids[i].slope * (key - lnids[i].key));
        #endif

        bool found = search_in_leaf_node(key, lnids[i].block_id, pos, v, block_count, lnids[i].iterm_count, lnids[i].added_in_buffer);
        delete []data;
        return found;
    }
    #endif

    bool lookup_disk(KeyType key, int *block_count, ValueType *v) {
        bool found = false;
        found = _lookup_disk(key, metanode.root_block_id, v, block_count);
        return found;
    }

    bool lookup_leaf_disk(KeyType key, int *block_count, ValueType *v) {
        // search in inner nodes
        stx_btree::iterator it = inner_btree.find_for_disk(key);
        // find the right leaf segment to search
        if ((it != inner_btree.begin() && it.key() >= key) || (it == inner_btree.begin() && it.key() == key) || it == inner_btree.end()) {
            if (it.key() > key || it == inner_btree.end())
                it--;
            LeafNodeValueOnDisk lfv = it.data();

            #if PGMSegment
            int pos = int64_t(lfv.slope * (key- it.key())) + lfv.intercept;
            #else
            int pos = int(lfv.slope * (key - it.key()));
            #endif
            bool found = search_in_leaf_node(key, lfv.block_id, pos, v, block_count, lfv.iterm_count, lfv.added_in_buffer);
            return found;
        } else if (it == inner_btree.begin() && it.key() > key) {
            return _get_from_segment_buffer(key, metanode.block_id, v, block_count);
        } else {
            // should not happen
            throw std::invalid_argument("cannot find a key that larger than `it.key()'...");
//            std::cout << "should not happen" << std::endl;
        }
    }

    size_t get_inner_size() {
        if (hybrid_mode == LEAF_DISK)
            return inner_btree.get_inner_size();
        else return 0;
    }

    size_t get_file_size() {
        return sm->get_file_size();
    }
    int lookup(KeyType key, int *block_count) {
        ValueType v;
        if (hybrid_mode == ALL_DISK) return lookup_disk(key, block_count, &v);
        else if (hybrid_mode == LEAF_DISK) return lookup_leaf_disk(key, block_count, &v);
        return 0;
    }

    std::vector<LeafNodeItermOnDisk> add_left_buffer(bool *no_op, KeyType key, ValueType value, int next_block_id) {
        char *_data2 = new char[MaxNodeSize];
        sm->get_block(metanode.block_id, _data2);
        Iterm *_temp = (Iterm *)_data2;
        int item_i = metanode.added_in_buffer+1;
        for (; item_i > 1 && _temp[item_i - 1].key > key; item_i--) { // the first one is meta
            _temp[item_i] = _temp[item_i - 1];
        }
        _temp[item_i].key = key;
        _temp[item_i].value = value;
        if (metanode.added_in_buffer == ItermCountOneBlock - 2) { // the first one is meta for next and last block
            // resegment
            #if PGMSegment
            KeyType *_xdata20 = new KeyType[ItermCountOneBlock-1];
            for (int _xdatai = 1; _xdatai < ItermCountOneBlock; _xdatai++) {
                _xdata20[_xdatai-1] = _temp[_xdatai].key;
            }
            std::vector<KeyType> _xdata21(_xdata20, _xdata20+(ItermCountOneBlock-1));
            std::vector<Seg> segments; 
            _ts(_xdata21.begin(), _xdata21.end(), error_bound, segments);
            std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment_pgm(_temp, ItermCountOneBlock-1, segments, false, metanode.block_id, next_block_id);
            #else
            std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment(buffer, buffer_size, false);
            #endif
            metanode.added_in_buffer = 0;
            *no_op = false;
            return leaf_iterms;
        }

        metanode.added_in_buffer += 1;
        sm->write_with_size(metanode.block_id, _data2, MaxNodeSize);
        std::vector<LeafNodeItermOnDisk> leaf_iterms;
        *no_op = true;
        return leaf_iterms;
    }

    void merge(Iterm *buffer, Iterm *old_buffer, Iterm *old_segment, int tl, int tr) {
        int m = 0, ml = 0, mr = 0;
        while (ml < tl && mr < tr) {
            buffer[m++] = old_buffer[ml].key <= old_segment[mr].key ? old_buffer[ml++] : old_segment[mr++];
        }
        while (ml < tl) {
            buffer[m++] = old_buffer[ml++];
        }
        while (mr < tr) {
            buffer[m++] = old_segment[mr++];
        }
        return;
    }

    std::vector<LeafNodeItermOnDisk> add_segment_buffer_for_leaf_disk(bool *no_op, KeyType true_key, ValueType value, stx_btree::iterator it) {
//        LeafNodeItermOnDisk *lnids = (LeafNodeItermOnDisk *)(data + LeaftNodeHeaderSizeOnDisk);
        char *_data2 = new char[MaxNodeSize];
        int iterm_count = it.data().iterm_count;
        int add_buffer = it.data().added_in_buffer;
        int start_block = it.data().block_id;
        int _block = (iterm_count+2) / ItermCountOneBlock;
        int _offset = (iterm_count+2)% ItermCountOneBlock;
        long offset = (start_block + _block) * MaxNodeSize + _offset * ItermSize;
        sm->read_block_arbitrary(_data2, offset);
        Iterm * _temp = (Iterm * )_data2;
        int _i = add_buffer - 1;
        for (; _i >= 0; _i--) {
            if (_temp[_i].key > true_key) {
                _temp[_i+1].key = _temp[_i].key;
                _temp[_i+1].value = _temp[_i].value;
                continue;
            }
            break;
        }
        _temp[_i+1].key = true_key;
        _temp[_i+1].value = value;

        add_buffer += 1;
        if (add_buffer == BufferSize) {
            int buffer_size = add_buffer + iterm_count;
            Iterm *buffer = new Iterm[buffer_size];
            Iterm *old_buffer = new Iterm[BufferSize];
            memcpy(old_buffer, _temp, BufferSize * ItermSize);
            Iterm *old_segment = new Iterm[iterm_count];
            int cp_c = 0;
            // memcpy(buffer + cp_c, _temp, BufferSize * ItermSize);
            // cp_c += BufferSize;
            int last_block = 0;
            int next_block = 0;
            for (int block_i = 0; block_i < _block; block_i++) {
                sm->get_block(start_block + block_i, _data2);
                if (block_i == 0) {
                    _temp = (Iterm *)_data2;
                    last_block = int(_temp[0].key);
                    next_block = int(_temp[0].value);
                    // BUG!!!!!
                    memcpy(old_segment + cp_c,  _temp + 2, (ItermCountOneBlock-2) * ItermSize);
                    cp_c += (ItermCountOneBlock-2);
                } else {
                    memcpy(old_segment + cp_c, _temp, ItermCountOneBlock * ItermSize);
                    cp_c += ItermCountOneBlock;
                }

            }
            if(_offset > 0 && _block > 0) {
                sm->get_block(start_block + _block, _data2);
                memcpy(old_segment + cp_c, _temp, _offset * ItermSize);
            } else if (_offset > 0 && _block == 0){
                sm->get_block(start_block + _block, _data2);
                _temp = (Iterm *)_data2;
                last_block = int(_temp[0].key);
                next_block = int(_temp[0].value);
                memcpy(old_segment + cp_c, _temp+2, (_offset-2) * ItermSize);
            }
            // change it to merge sort...
            // std::sort(buffer, buffer + buffer_size, ItermCompare);
            merge(buffer, old_buffer, old_segment, BufferSize, iterm_count);
            delete []old_buffer;
            delete []old_segment;
            // resegment and write
#if Disk
#if PGMSegment
            KeyType *_xdata20 = new KeyType[buffer_size];
            for (int _xdatai = 0; _xdatai < buffer_size; _xdatai++) {
                _xdata20[_xdatai] = buffer[_xdatai].key;
            }
            std::vector<KeyType> _xdata21(_xdata20, _xdata20+buffer_size);
            std::vector<Seg> segments;
            _ts(_xdata21.begin(), _xdata21.end(), error_bound, segments);
            std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment_pgm(buffer, buffer_size, segments, false, last_block, next_block);
            delete[] _data2;
            // leak
            delete[] buffer;
            delete[] _xdata20;
            *no_op = false;
            return leaf_iterms;
#else
            std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment(buffer, buffer_size, false);
#endif
#else
            std::vector<FLeafNodeIterm> leaf_iterms = write_segment(buffer, lnids[i].added_in_buffer + lnids[i].iterm_count);
#endif
        } else {
            // flush data
            LeafNodeValueOnDisk lfv = it.data();
            lfv.added_in_buffer += 1;
            it.set_value(lfv);
            sm->write_arbitrary(offset, _data2, BufferSize * ItermSize);
            *no_op = true;
            delete[] _data2;
            std::vector<LeafNodeItermOnDisk> leaf_iterms;
            return leaf_iterms;
        }
    }

    std::vector<LeafNodeItermOnDisk> add_segment_buffer(bool *no_op, int start_block, KeyType true_key, ValueType value, char *data, int i, int current_block_id) {
            LeafNodeItermOnDisk *lnids = (LeafNodeItermOnDisk *)(data + LeaftNodeHeaderSizeOnDisk);
            char *_data2 = new char[MaxNodeSize];
            int _block = (lnids[i].iterm_count+2) / ItermCountOneBlock;
            int _offset = (lnids[i].iterm_count+2)% ItermCountOneBlock;
            long offset = (start_block + _block) * MaxNodeSize + _offset * ItermSize;
            sm->read_block_arbitrary(_data2, offset);
            Iterm * _temp = (Iterm * )_data2;
            int _i = lnids[i].added_in_buffer - 1;
            for (; _i >= 0; _i--) {
                if (_temp[_i].key > true_key) {
                    _temp[_i+1].key = _temp[_i].key;
                    _temp[_i+1].value = _temp[_i].value;
                    continue;
                }
                break;
            }
            _temp[_i+1].key = true_key;
            _temp[_i+1].value = value;

            lnids[i].added_in_buffer+=1;
            if (lnids[i].added_in_buffer == BufferSize) {
                int buffer_size = lnids[i].added_in_buffer + lnids[i].iterm_count;
                Iterm *buffer = new Iterm[buffer_size];
                Iterm *old_buffer = new Iterm[BufferSize];
                memcpy(old_buffer, _temp, BufferSize * ItermSize);
                Iterm *old_segment = new Iterm[lnids[i].iterm_count];
                int cp_c = 0;
                // memcpy(buffer + cp_c, _temp, BufferSize * ItermSize);
                // cp_c += BufferSize;
                int last_block = 0;
                int next_block = 0;
                for (int block_i = 0; block_i < _block; block_i++) {
                    sm->get_block(start_block + block_i, _data2);
                    if (block_i == 0) {
                        _temp = (Iterm *)_data2;
                        last_block = int(_temp[0].key);
                        next_block = int(_temp[0].value);
                        // BUG!!!!!
                        memcpy(old_segment + cp_c,  _temp + 2, (ItermCountOneBlock-2) * ItermSize);
                        cp_c += (ItermCountOneBlock-2);
                    } else {
                        memcpy(old_segment + cp_c, _temp, ItermCountOneBlock * ItermSize);
                        cp_c += ItermCountOneBlock;
                    }
                    
                }
                if(_offset > 0 && _block > 0) {
                    sm->get_block(start_block + _block, _data2);
                    memcpy(old_segment + cp_c, _temp, _offset * ItermSize);
                } else if (_offset > 0 && _block == 0){
                    sm->get_block(start_block + _block, _data2);
                    _temp = (Iterm *)_data2;
                    last_block = int(_temp[0].key);
                    next_block = int(_temp[0].value);
                    memcpy(old_segment + cp_c, _temp+2, (_offset-2) * ItermSize);
                }
                // change it to merge sort...
                // std::sort(buffer, buffer + buffer_size, ItermCompare);
                merge(buffer, old_buffer, old_segment, BufferSize, lnids[i].iterm_count);
                delete []old_buffer;
                delete []old_segment;
                // resegment and write
                #if Disk
                    #if PGMSegment
                    KeyType *_xdata20 = new KeyType[buffer_size];
                    for (int _xdatai = 0; _xdatai < buffer_size; _xdatai++) {
                        _xdata20[_xdatai] = buffer[_xdatai].key;
                    }
                    std::vector<KeyType> _xdata21(_xdata20, _xdata20+buffer_size);
                    std::vector<Seg> segments; 
                    _ts(_xdata21.begin(), _xdata21.end(), error_bound, segments);
                    std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment_pgm(buffer, buffer_size, segments, false, last_block, next_block);
                    delete[] _data2;
                    // leak
                    delete[] buffer;
                    delete[] _xdata20;
                    *no_op = false;
                    return leaf_iterms;
                    #else
                    std::vector<LeafNodeItermOnDisk> leaf_iterms = write_segment(buffer, buffer_size, false);
                    #endif
                #else
                    std::vector<FLeafNodeIterm> leaf_iterms = write_segment(buffer, lnids[i].added_in_buffer + lnids[i].iterm_count);
                #endif
            } else {
                // flush header // must in this order
                sm->write_with_size(current_block_id, data, MaxNodeSize);
                // flush data
                sm->write_arbitrary(offset, _data2, BufferSize * ItermSize);
                *no_op = true;
                delete[] _data2;
                std::vector<LeafNodeItermOnDisk> leaf_iterms;
                return leaf_iterms;
            }
    }

    // here we suppose we have added the smallest key into the index when bulk loading.
    FBuildStatus insert_key(KeyType key, ValueType value, int *block_count, int current_block_id, KeyType true_key) {
        // load the leaf node
        FBuildStatus bs;
        char *data = new char[MaxNodeSize];

        sm->get_block(current_block_id, data);
        *block_count += 1;

        bool is_inner = data[0] == FInnerNodeType;
        if (is_inner) {
            InnerNodeHeaderOnDisk *inhd;
            inhd = (InnerNodeHeaderOnDisk *)data;
            InnerNodeItermOnDisk* inids = (InnerNodeItermOnDisk*) (data + InnerNodeHeaderSizeOnDisk);
            int i = 0;
            int current_block = 0;
            bool is_left_most =false;
            if (key < inids[0].key) {
                current_block = inhd->left_child;
                if (current_block == 0) std::cout << "not right block id" << std::endl;
                is_left_most = true;
            }
            else {
                #if Binary
                int _start = 0;
                int _end = inhd->item_count - 1;
                while (_start <= _end) {
                    int mid = _start + (_end - _start) / 2;
                    if (inids[mid].key <= key) _start = mid + 1;
                    else _end = mid - 1;
                }
                i = _start  - 1;
                #else
                    for (; i < inhd->item_count - 1; i++) {
                        if (inids[i].key <= key && inids[i+1].key > key) break;
                    }
                #endif
                current_block = inids[i].block_id;
            }
            
            // int next_block_id = inids[i].block_id;
            // delete[] data;

            bs = insert_key(key, value, block_count, current_block, true_key);
            if (bs.status == AddNewEntry) {
                // must behind i before i+1
                
                int _i = inhd->item_count - 1;
                // if it is in special case, i = 0 is still right
                int _step_ = is_left_most ? 1 : 0;
                for (; _i > (i - _step_); _i--) {
                    inids[_i + 1] = inids[_i];
                }
                inids[_i + 1] = bs.inid;
                inhd->item_count += 1;
                if (inhd->item_count == MaxInnerNodeItemCountOnDisk) {
                    InnerNodeHeaderOnDisk _inhd;
                    _inhd.tag = FInnerNodeType;
                    _inhd.left_child = 0;
                    int mid = inhd->item_count / 2;
                    int n_size = inhd->item_count - mid;
                    InnerNodeItermOnDisk* _inids = new InnerNodeItermOnDisk[n_size];
                    _inhd.item_count = n_size;
                    memcpy(_inids, inids + mid, n_size * InnerNodeItermSizeOnDisk);
                    inhd->item_count = mid;
                    sm->write_with_size(current_block_id, data, MaxNodeSize);
                    bs.status = AddNewEntry;

                    Block block;
                    memcpy(block.data, &_inhd, InnerNodeHeaderSizeOnDisk);
                    long _isize = _inhd.item_count*InnerNodeItermSizeOnDisk;
                    memcpy(block.data + InnerNodeHeaderSizeOnDisk, _inids, _isize);
                    sm->write_block(metanode.block_count, block);
                    bs.inid.block_id = metanode.block_count;
                    bs.inid.key = _inids[0].key;
                    bs.key = inids[0].key;
                    metanode.block_count += 1;
                    sync_metanode();
                    delete []_inids;
                } else {
                    sm->write_with_size(current_block_id, data, MaxNodeSize);
                    bs.status = NoOperation;
                }
            }
        } else {
            LeaftNodeHeaderOnDisk *lnhd = (LeaftNodeHeaderOnDisk *)data;
            LeafNodeItermOnDisk *lnids = (LeafNodeItermOnDisk *)(data + LeaftNodeHeaderSizeOnDisk);
            int i = 0;
            #if Binary
                int _start = 0;
                int _end = lnhd->item_count - 1;
                while (_start <= _end) {
                    int mid = _start + (_end - _start) / 2;
                    if (lnids[mid].key <= key) _start = mid + 1;
                    else _end = mid - 1;
                }
                i = _start  - 1;
            #else
                for (; i < lnhd->item_count - 1; i++) {
                    if (lnids[i].key <= key && lnids[i+1].key > key) break;
                }
            #endif
            // handle special case
            bool no_op = false;
            std::vector<LeafNodeItermOnDisk> leaf_iterms;
            bool is_specal = false;
            if (i == 0 && lnids[i].key > key) {
                leaf_iterms = add_left_buffer(&no_op, key, value, lnids[i].block_id);
                is_specal = true;
            } else {
                leaf_iterms = add_segment_buffer(&no_op, lnids[i].block_id, true_key, value, data, i, current_block_id);
            }
            if (no_op) {
                bs.status = NoOperation;
                // leak
                delete []data;
                return bs;
            }
            
            int _step_ = is_specal ? 0 : 1;
            int _lower_ = is_specal ? i - 1 : i;
            // update the leaf node
            // right for both left most or other case
            if (lnhd->item_count + leaf_iterms.size() - _step_ < MaxLeafNodeItemCountOnDisk) { // just add it
                // i-th
                for (int _j = lnhd->item_count  - 1; _j > _lower_; _j--) {
                    lnids[_j + leaf_iterms.size() - _step_] = lnids[_j];
                }
                for (int _j = 0; _j < leaf_iterms.size(); _j++) {
                    lnids[i + _j] = leaf_iterms.at(_j);
                }
                lnhd->item_count = lnhd->item_count + leaf_iterms.size() - _step_;
                sm->write_with_size(current_block_id, data, MaxNodeSize);
                bs.status = NoOperation;
                // leak
                delete []data;
                return bs;
            }
            //copy and split
            LeafNodeItermOnDisk all_lnids[lnhd->item_count + leaf_iterms.size() - _step_];
            int all_i = 0;
            for (int _j = 0; _j < _lower_; _j++, all_i++) {
                all_lnids[all_i] = lnids[_j];
            }
            for (int _j = 0; _j < leaf_iterms.size(); _j++, all_i++) {
                all_lnids[all_i] = leaf_iterms.at(_j);
            }
            for (int _j = _lower_ + 1; _j < lnhd->item_count; _j++, all_i++) {
                all_lnids[all_i] = lnids[_j];
            }
            int mid = (leaf_iterms.size() + lnhd->item_count  - _step_) / 2;
            LeaftNodeHeaderOnDisk r_lnhd;
            r_lnhd.item_count = leaf_iterms.size() + lnhd->item_count  - _step_ - mid;
            r_lnhd.tag = FLeafNodeType;
            LeafNodeItermOnDisk r_lnids[r_lnhd.item_count];
            for (int _i = 0; _i < r_lnhd.item_count; _i++) {
                r_lnids[_i] = all_lnids[_i + mid];
            }
            for (int _i = 0; _i < mid; _i++) {
                lnids[_i] = all_lnids[_i];
            }
            lnhd->item_count = mid;
            Block block;
            memcpy(block.data, &r_lnhd, LeaftNodeHeaderSizeOnDisk);
            long _isize = r_lnhd.item_count*LeafNodeItermSizeOnDisk;
            memcpy(block.data + LeaftNodeHeaderSizeOnDisk, r_lnids, _isize);
            sm->write_block(metanode.block_count, block);
            bs.status = AddNewEntry;
            bs.inid.key = r_lnids[0].key;
            bs.inid.block_id = metanode.block_count;
            
            metanode.block_count += 1;
            sync_metanode();
            sm->write_with_size(current_block_id, data, MaxNodeSize);
        }
        delete[] data;
        return bs;

    }

    void insert_disk(KeyType key, ValueType value) {
        //int *block_count, int current_block_id
        int block_count = 0;
        KeyType search_key = key;
        // if (key < metanode.original_key) search_key = metanode.original_key;
        FBuildStatus bs = insert_key(search_key, value, &block_count, metanode.root_block_id, key);
        if (bs.status == AddNewEntry) {
            Block block;
            InnerNodeHeaderOnDisk inh;
            inh.item_count = 1;
            inh.tag = FInnerNodeType;
            inh.left_child = metanode.root_block_id;
            InnerNodeItermOnDisk inid;
            inid.block_id =bs.inid.block_id;
            inid.key = bs.inid.key;
            memcpy(block.data, &inh, InnerNodeHeaderSizeOnDisk);
            memcpy(block.data + InnerNodeHeaderSizeOnDisk, &inid, InnerNodeItermSizeOnDisk);
            sm->write_block(metanode.block_count, block);
            metanode.root_block_id = metanode.block_count;
            metanode.level += 1;
            metanode.block_count += 1;
            sync_metanode();
        }
    }

    void insert_key_leaf_disk(KeyType key, ValueType value) {
        // search in inner nodes
        stx_btree::iterator it = inner_btree.find_for_disk(key);
        // find the right leaf segment to search
        if (it != inner_btree.begin() && it.key() >= key || (it == inner_btree.begin() && it.key() == key ) || it == inner_btree.end()) {// insert into the related segment's buffer
            if (it.key() > key || it == inner_btree.end())
                it--;
            LeafNodeValueOnDisk lfv = it.data();
            bool no_op = false;
            std::vector<LeafNodeItermOnDisk> added = add_segment_buffer_for_leaf_disk(&no_op, key, value, it);
            if (no_op) return;
            for (int i = 0; i < added.size(); ++i) {
                LeafNodeValueOnDisk lfv;
                lfv.block_id = added[i].block_id; lfv.iterm_count = added[i].iterm_count;
                lfv.added_in_buffer = added[i].added_in_buffer, lfv.intercept = added[i].intercept;
                lfv.slope = added[i].slope;
                if (i == 0) it.set_value(lfv);
                else inner_btree.insert(added[i].key, lfv);
            }
        } else if (it == inner_btree.begin() && it.key() > key) { //insert into the left-most buffer
            bool no_op = false;
            int next_block_id = metanode.block_id;
//            metanode.block_id;
            std::vector<LeafNodeItermOnDisk> added = add_left_buffer(&no_op, key, value, next_block_id);
            if (no_op) return;
            for (int i = 0; i < added.size(); i++) {
                inner_btree.insert(added[i].key, {added[i].slope, added[i].block_id,
                                                  added[i].iterm_count,
                                                  added[i].added_in_buffer, added[i].intercept});
            }
        } else {
            // should not happen
            std::invalid_argument("cannot find a key that larger than `it.key()'...");
        }
    }
    void insert_key_entry_f(KeyType key, ValueType value) {
        if (hybrid_mode == ALL_DISK) return insert_disk(key, value);
        else if (hybrid_mode == LEAF_DISK) return insert_key_leaf_disk(key, value);
    }

    void scan(KeyType key, KeyType *v, int *block_count, int len) {
        *block_count = 0;
        if (hybrid_mode == ALL_DISK)
            simple_scan_disk(key, metanode.root_block_id, v, block_count, len);
        else if (hybrid_mode == LEAF_DISK)
            scan_leaf_disk(key, v, block_count, len);
        return;
    }

    void scan_leaf_disk(KeyType key, KeyType *v, int *block_count, int len) {
        // search in inner nodes
        stx_btree::iterator it = inner_btree.find_for_disk(key);
        // find the right leaf segment to search
        if (it != inner_btree.begin() && it.key() >= key) { // search from a segment
            if (it.key() > key)
                it--;
            LeafNodeValueOnDisk lfv = it.data();

            #if PGMSegment
            int pos = int64_t(lfv.slope * (key- it.key())) + lfv.intercept;
            #else
            int pos = int(lfv.slope * (key - it.key()));
            #endif
            int start_block = lfv.block_id;
            int start = std::max<int>(pos - error_bound, 0) + 2;
            int end = std::min<int>(pos + error_bound+2, lfv.iterm_count) + 2;
            return scan_in_a_segment(start, end, key, start_block, block_count,
                                     lfv.iterm_count, len, v, false);

        } else if (it == inner_btree.begin() && it.key() > key) { // search from left-most buffer
            return scan_in_a_segment(1, metanode.added_in_buffer, key, metanode.block_id, block_count,
                                     metanode.added_in_buffer, len, v, true);
        } else {
            // should not happen
            std::invalid_argument("cannot find a key that larger than `it.key()'...");
        }
    }

    // todo: we do not test scan after insert; add code later
    void scan_in_a_segment(int start, int end, int key, int start_block, int *block_count,
                           int iterm_count, int len, ValueType *v, bool is_left_most) {
//        int start = std::max<int>(pos - error_bound, 0) + 2;
//        int end = std::min<int>(pos + error_bound+2, iterm_count) + 2;
        int last_block = -1;
        char *_data2 = new char[MaxNodeSize];
        Iterm *_temp = (Iterm*)_data2;

        // handle the case that, the left-most buffer is empty.
        // we will not go into the `while loop' start = 1, end = 0
        int mid = iterm_count + (is_left_most ? 1 : 2);
        int _offset;

        while (start <= end) {
            mid = start + (end - start) / 2;
            int _block = mid / ItermCountOneBlock;
            _offset = mid % ItermCountOneBlock;
            read_block_or_not(&last_block, start_block+_block, _data2, block_count);
            _temp = (Iterm * )_data2;
            if (_temp[_offset].key == key) {
                break;
            } else if (_temp[_offset].key < key) start = mid + 1;
            else end = mid - 1;
        }

        int _i = 0;
        int tc = iterm_count + (is_left_most ? 1 : 2);
        _temp = (Iterm * )_data2;
        int next_segement_block = 0;

        for (;_i < len;) {
            if (_offset < ItermCountOneBlock && mid < tc) {
                v[_i] = _temp[_offset].key;
                _offset += 1;
                mid += 1;
                _i += 1;
            } else if (_offset == ItermCountOneBlock && mid < tc) {
                // read next block
                sm->get_block(last_block + 1, _data2);
                *block_count += 1;
                last_block += 1;
                _temp = (Iterm * )_data2;
                _offset = 0;
            } else if (mid == tc) {
                // reade next segment;
                sm->get_block(start_block, _data2);
                *block_count += 1;
                // last_block = 0;
                _temp = (Iterm*)_data2;
                if (next_segement_block == 0) //lazy fetch next block id
                    next_segement_block = int(_temp[0].value);
                if (next_segement_block == 0) {
                    break;
                }
                sm->get_block(next_segement_block, _data2);
                *block_count += 1;
                _temp = (Iterm * )_data2;
                next_segement_block = int(_temp[0].value);
                tc = int(_temp[1].key) + 2;
                mid = 2;
                _offset = 2;

            }
        }
        delete[] _data2;
        return;
    }

    void simple_scan_disk(KeyType key, int block_id, KeyType *v, int *block_count, int len) {
        if (key < metanode.original_key) { // must in the first segment buffer if have
            return scan_in_a_segment(1, metanode.added_in_buffer, key, metanode.block_id, block_count,
                              metanode.added_in_buffer, len, v, true);
        }
        char *data = new char[MaxNodeSize];

        sm->get_block(block_id, data);
        *block_count += 1;

        bool is_inner = data[0] == FInnerNodeType;
        while (is_inner) {
            InnerNodeHeaderOnDisk *inhd;
            inhd = (InnerNodeHeaderOnDisk *)data;
            InnerNodeItermOnDisk* inids = (InnerNodeItermOnDisk*) (data + InnerNodeHeaderSizeOnDisk);
            int i = 0;
            int next_id = 0;
            if (inids[0].key > key) {
                next_id = inhd->left_child;
            } else {
                #if Binary
                int _start = 0;
                int _end = inhd->item_count - 1;
                while (_start <= _end) {
                    int mid = _start + (_end - _start) / 2;
                    if (inids[mid].key <= key) _start = mid + 1;
                    else _end = mid - 1;
                }
                i = _start  - 1;
                #else
                    for (; i < inhd->item_count - 1; i++) {
                        if (inids[i].key <= key && inids[i+1].key > key) break;
                    }
                #endif
                next_id = inids[i].block_id;
            }
            
            sm->get_block(next_id, data);
            *block_count += 1;
            is_inner = data[0] == FInnerNodeType;
        }
        LeaftNodeHeaderOnDisk *lnhd = (LeaftNodeHeaderOnDisk *)data;
        LeafNodeItermOnDisk *lnids = (LeafNodeItermOnDisk *)(data + LeaftNodeHeaderSizeOnDisk);
        int i = 0;
        #if Binary
            int _start = 0;
            int _end = lnhd->item_count - 1;
            while (_start <= _end) {
                int mid = _start + (_end - _start) / 2;
                if (lnids[mid].key <= key) _start = mid + 1;
                else _end = mid - 1;
            }
            i = _start  - 1;
        #else
            for (; i < lnhd->item_count - 1; i++) {
                if (lnids[i].key <= key && lnids[i+1].key > key) break;
            }
        #endif

        #if PGMSegment
        int pos = int64_t(lnids[i].slope * (key- lnids[i].key)) + lnids[i].intercept;
        #else
        int pos = int(lnids[i].slope * (key - lnids[i].key));
        #endif
        int start_block = lnids[i].block_id;
        int start = std::max<int>(pos - error_bound, 0) + 2;
        int end = std::min<int>(pos + error_bound+2, lnids[i].iterm_count) + 2;
        return scan_in_a_segment(start, end, key, start_block, block_count,
                          lnids[i].iterm_count, len, v, false);


    }

    #if Disk
        int get_level() {
            return metanode.level;
        }
    #endif
};