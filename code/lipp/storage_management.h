#define KeyType uint64_t
#define ValueType uint64_t

#include<map>
#include<cstring>
#include<stack>
long BLOCK_SIZE = 8192/2;

typedef uint8_t bitmap_t;
#define BITMAP_WIDTH (sizeof(bitmap_t) * 8)
#define BITMAP_SIZE(num_items) (((num_items) + BITMAP_WIDTH - 1) / BITMAP_WIDTH)
#define BITMAP_GET(bitmap, pos) (((bitmap)[(pos) / BITMAP_WIDTH] >> ((pos) % BITMAP_WIDTH)) & 1)
#define BITMAP_SET(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] |= 1 << ((pos) % BITMAP_WIDTH))
#define BITMAP_CLEAR(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] &= ~bitmap_t(1 << ((pos) % BITMAP_WIDTH)))
#define BITMAP_NEXT_1(bitmap_item) __builtin_ctz((bitmap_item))

#define Profiling 1

template<class T, class P>
class StorageManager {

    const double BUILD_LR_REMAIN = 0;

    inline int compute_gap_count(int size) {
        if (size >= 1000000) return 1;
        if (size >= 100000) return 2;
        return 5;
    }

    FILE *fp;
    typedef struct {
        int is_two;
        int size;
        int build_size;
        int fixed;
        int number_inserts;
        int num_insert_to_data;
        int num_items;
        // int empty_size1; // for header
        int empty_size2; // for items
        double slope;
        long double intercept;
    } NodeHeaderD;
    int NodeHeaderDSize = sizeof(NodeHeaderD);

    typedef struct {
        // 1 none
        // 2 child
        // 3 data
        char tag; // we know the type from the bitmap in header
        union {
            struct {
                T key;
                P value;
            } data;
            struct {
                int block;
                int offset;
            } addr;
        } comp;
    } ItemD;
    int ItemDSize = sizeof(ItemD);
    int MaxItemCount = BLOCK_SIZE / ItemDSize;

    typedef struct {
        int root_block;
        int root_offset;
        int NEXT_BLOCK;
        int NEXT_OFFSET;
    } MetaNode;
    int MetaNodeSize = sizeof(MetaNode);

    int NEXT_BLOCK;
    int NEXT_OFFSET;
    MetaNode mb;

    private:
        void write_data(void *data, long offset, int len) {
            fseek(fp, offset, SEEK_SET);
            fwrite(data, len, 1, fp);
            return;
        }

        void read_block(void *data, int block_id) {
            fseek(fp, block_id * BLOCK_SIZE, SEEK_SET);
            fread(data, BLOCK_SIZE, 1, fp);
            return;
        }

        void write_block(void *data, int block_id) {
            fseek(fp, block_id * BLOCK_SIZE, SEEK_SET);
            fwrite(data, BLOCK_SIZE, 1, fp);
            return;
        }

        void read_data(void *data, long offset, int len) {
            fseek(fp, offset, SEEK_SET);
            fread(data, len, 1, fp);
            return;
        }


    public:
        //StorageManager (
        void sys_metablock(bool update_root, int block = 0, int offset = 0) {
            char empty_block[BLOCK_SIZE];
            mb.NEXT_BLOCK = NEXT_BLOCK;
            mb.NEXT_OFFSET = NEXT_OFFSET; 
            if (update_root) {
                mb.root_block = block;
                mb.root_offset = offset;
            }
            memcpy(empty_block, &mb, MetaNodeSize);
            write_data(empty_block, 0, BLOCK_SIZE);
        }

        size_t get_file_size() {
            fseek(fp, 0, SEEK_END);
            return ftell(fp);
        }

        void load_metablock() {
            char data[BLOCK_SIZE];
            read_block(data, 0);
            MetaNode *_mb = (MetaNode *)(data);
            mb.NEXT_BLOCK = _mb->NEXT_BLOCK;
            mb.NEXT_OFFSET = _mb->NEXT_OFFSET;
            mb.root_block = _mb->root_block;
            mb.root_offset = _mb->root_offset;
        }

        void init(char *fn, bool is_first) {
            if (is_first) {
                fp = fopen(fn,"wb");
                // char empty_block[BLOCK_SIZE];
                // MetaNode mb;
                NEXT_BLOCK = 1;
                NEXT_OFFSET = 0;
                // mb.root_block = 1;
                // mb.root_offset = 0; 
                // memcpy(empty_block, &mb, MetaNodeSize);
                // write_data(empty_block, 0, BLOCK_SIZE);
                sys_metablock(true, 1, 0);
                fclose(fp);
            }
            fp = fopen(fn,"r+b");
            load_metablock();
        } 

        void _write_node_internal(NodeHeaderD nhd, ItemD *its, int *block, int *offset) {
            if (BLOCK_SIZE - NEXT_OFFSET < NodeHeaderDSize) {
                char empty[BLOCK_SIZE - NEXT_OFFSET];
                write_data(empty, NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET, BLOCK_SIZE - NEXT_OFFSET);
                NEXT_BLOCK += 1;
                NEXT_OFFSET = 0;
            }
            *block = NEXT_BLOCK;
            *offset = NEXT_OFFSET;
            // Write Header
            long start_offset = NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET;
            long _start_header = start_offset;
            write_data(&nhd, start_offset, NodeHeaderDSize);
            start_offset += NodeHeaderDSize;

            nhd.empty_size2 = (BLOCK_SIZE - (start_offset % BLOCK_SIZE)) % ItemDSize;
            if (start_offset % BLOCK_SIZE == 0) {
                nhd.empty_size2 = 0;
            }
            write_data(&nhd, _start_header, NodeHeaderDSize);
            char empty[nhd.empty_size2];
            write_data(empty, start_offset, nhd.empty_size2);
            start_offset += nhd.empty_size2;
            
            
            int c1 = (BLOCK_SIZE - (start_offset % BLOCK_SIZE)) / ItemDSize;
            if (c1 < MaxItemCount) {
                if (c1 > nhd.num_items) c1 = nhd.num_items;
                write_data(its, start_offset, ItemDSize * c1);
                start_offset += ItemDSize * c1;
            } else {
                c1 = 0;   
            }
            
            // NEXT_BLOCK = start_offset / BLOCK_SIZE;
            char data[BLOCK_SIZE];
            for (int i = c1; i < nhd.num_items; ) {
                int _c = nhd.num_items - i;
                if (_c > MaxItemCount) _c = MaxItemCount;
                if (_c == MaxItemCount) {
                    memcpy(data, its + i,  _c * ItemDSize);
                    write_data(data, start_offset, BLOCK_SIZE);
                    start_offset += BLOCK_SIZE;
                } else {
                    write_data(its + i, start_offset, _c * ItemDSize);
                    start_offset += _c * ItemDSize;
                }
                i += _c;
            }
            
            NEXT_BLOCK = start_offset / BLOCK_SIZE;
            NEXT_OFFSET = start_offset % BLOCK_SIZE;
            // todo update the metanode
            // sys_metablock(false);
        }

        void update_child_address(int block, int offset, int i, int child_block, int child_offset) {
            char data[BLOCK_SIZE];
            read_block(data, block);
            NodeHeaderD *nhd = (NodeHeaderD *) (data + offset);
            //long i_offset = block * BLOCK_SIZE + offset + NodeHeaderDSize + nhd->empty_size2 + i * ItemDSize;

            long _offset_ = (block * BLOCK_SIZE + offset + NodeHeaderDSize + nhd->empty_size2) % BLOCK_SIZE;
            int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
            if (_offset_ == 0) {
                c1 = 0;
            }
            ItemD id;
            id.tag = 2;
            id.comp.addr.block = child_block;
            id.comp.addr.offset = child_offset;
            long i_offset = 0;
            if (i < c1) {
                // itd = (ItemD *) (data + _offset_ + pos * ItemDSize);
                i_offset = block * BLOCK_SIZE + offset + NodeHeaderDSize + nhd->empty_size2 + i * ItemDSize;
            } else {
                block += 1;
                int _b = (i - c1) / MaxItemCount;
                block += _b;
                int _offset = (i - c1) % MaxItemCount;
                i_offset = block * BLOCK_SIZE + _offset * ItemDSize;
                // read_block(data, block2);
                // itd = (ItemD *) (data + offset2*ItemDSize);
            }
            
            write_data(&id, i_offset, ItemDSize);
        }

        bool searchv1(T key, P *value, int *bc){
            NodeHeaderD nhd;
            ItemD itd;
            long _offset1 = 0;
            long _offset2 = 0;
            int block = mb.root_block;
            int offset = mb.root_offset;
            *bc=0;
            while (true) {
                // read node header
                _offset1 = block * BLOCK_SIZE + offset;
                read_data(&nhd, _offset1, NodeHeaderDSize);
                double v = nhd.slope * static_cast<long double>(key) + nhd.intercept;
                int pos = 0;
                if (v > std::numeric_limits<int>::max() / 2) {
                    pos = nhd.num_items - 1;
                } else if (v < 0) {
                    pos = 0;
                } else {
                    pos = std::min(nhd.num_items - 1, static_cast<int>(v));
                }
                // read iterm
                long _offset_ = (block * BLOCK_SIZE + offset + NodeHeaderDSize + nhd.empty_size2) % BLOCK_SIZE;
                int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                if (_offset_ == 0) {
                    c1 = 0;
                }
                
                if (pos < c1) {
                    //itd = (ItemD *) (data + _offset_ + pos * ItemDSize);
                    _offset2 = block*BLOCK_SIZE + _offset_ + pos * ItemDSize;
                } else {
                    // block2 = block + 1;
                    int _b = (pos - c1) / MaxItemCount;
                    // block2 += _b;
                    (pos - c1) % MaxItemCount;
                    // read_block(data, block2);
                    // *bc += 1;
                    // itd = (ItemD *) (data + offset2*ItemDSize);
                    _offset2 = (block + 1 + _b) * BLOCK_SIZE + ((pos - c1) % MaxItemCount) * ItemDSize;
                }


                // _offset2 = _offset1 + NodeHeaderDSize + nhd.empty_size2 + ItemDSize * pos;
                read_data(&itd, _offset2, ItemDSize);
                if (itd.tag == 1) {
                    return false;
                } else if (itd.tag == 2) {
                    block = itd.comp.addr.block;
                    offset = itd.comp.addr.offset;
                } else {
                    if (itd.comp.data.key == key) {
                        *value = itd.comp.data.value;
                        return true;
                    }
                    return false;
                }
            }
        }

        bool searchv2(T key, P *value, int *bc, int *ic, int *lc){
            NodeHeaderD *nhd;
            char data[BLOCK_SIZE];
            ItemD *itd;
            int block1 = mb.root_block;
            int offset1 = mb.root_offset;
            int block2 = 0;
            int offset2 = 0;
            *bc = 0;
            // load_metablock();
            // *bc += 1;
            while (true) {
                // read node header
                *lc += 1;
                read_block(data, block1);
                *bc += 1;
                nhd = (NodeHeaderD *) (data + offset1);
                double v = nhd->slope * static_cast<long double>(key) + nhd->intercept;
                int pos = 0;
                if (v > std::numeric_limits<int>::max() / 2) {
                    pos = nhd->num_items - 1;
                } else if (v < 0) {
                    pos = 0;
                } else {
                    pos = std::min(nhd->num_items - 1, static_cast<int>(v));
                }
                // read iterm
                long _offset_ = (block1 * BLOCK_SIZE + offset1 + NodeHeaderDSize + nhd->empty_size2) % BLOCK_SIZE;
                int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                if (_offset_ == 0) {
                    c1 = 0;
                }
                if (pos < c1) {
                    itd = (ItemD *) (data + _offset_ + pos * ItemDSize);
                } else {
                    block2 = block1 + 1;
                    int _b = (pos - c1) / MaxItemCount;
                    block2 += _b;
                    offset2 = (pos - c1) % MaxItemCount;
                    read_block(data, block2);
                    *bc += 1;
                    itd = (ItemD *) (data + offset2*ItemDSize);
                }
                
                if (itd->tag == 1) {
                    return false;
                } else if (itd->tag == 2) {
                    block1 = itd->comp.addr.block;
                    offset1 = itd->comp.addr.offset;
                } else {
                    if (itd->comp.data.key == key) {
                        *value = itd->comp.data.value;
                        return true;
                    }
                    return false;
                }
            }
        }

        int range_query_len(T *results, T lower, int len, int *BC, int *ic, int *lc) {
            char data[BLOCK_SIZE];
            return range_core_len(data, false, false, results, 0, mb.root_block, mb.root_offset, lower, len, BC, ic, lc);
        }

        int range_core_len(char *o_data, bool same, bool SAT, T *results, int pos, int block, int offset, T lower, int len, int *BC, int *ic, int *lc) {
            NodeHeaderD node;
            char data[BLOCK_SIZE];
            ItemD *itd;
            if (!same) {
                read_block(data, block);
            } else {
                memcpy(data, o_data, BLOCK_SIZE);
            }
            
            *BC += 1;
            memcpy(&node, data + offset, NodeHeaderDSize);
            int _init_block = block;
            // Profile Node Count
            *lc += 1;
            if (SAT) {
                long _offset_ = (block * BLOCK_SIZE + offset + NodeHeaderDSize + node.empty_size2) % BLOCK_SIZE;
                int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                if (_offset_ == 0) {
                    block += 1;
                }
                for (int i = 0; i < node.num_items;) {
                    if (block != _init_block) {
                        read_block(data, block);
                        *BC += 1;
                    }
                    ItemD *_items = (ItemD *) (data + _offset_);
                    if (c1 > (node.num_items - i)) {
                        c1 = node.num_items - i;
                    }
                    for (int j = 0; j < c1; j++) {
                        // Profile item count
                        *ic += 1;
                        if (_items[j].tag == 3) {
                            results[pos] = _items[j].comp.data.key;
                            pos ++;
                        } else if (_items[j].tag == 2) {
                            int _b = _items[j].comp.addr.block;
                            int _o = _items[j].comp.addr.offset;
                            pos = range_core_len(data, _items[j].comp.addr.block == block, true, results, pos,_items[j].comp.addr.block, _items[j].comp.addr.offset, lower, len, BC, ic, lc);
                        }
                        if (pos >= len) {
                            return pos;
                        }
                    }
                
                    i += c1;
                    block += 1;
                    _offset_ = 0;
                    c1 = MaxItemCount;
                }
                return pos;
            }
            else {
                int lower_pos = predict_pos(node, lower);
                long _offset_ = (block * BLOCK_SIZE + offset + NodeHeaderDSize + node.empty_size2) % BLOCK_SIZE;
                int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                if (_offset_ == 0) {
                    c1 = 0;
                }
                int remains = 0;
                if (lower_pos < c1) {
                    itd = (ItemD *) (data + _offset_ + lower_pos * ItemDSize);
                    remains = (c1 > node.num_items ? node.num_items : c1) - (lower_pos + 1);
                } else {
                    // int block2;
                    int offset2;
                    block = block + 1;
                    int _b = (lower_pos - c1) / MaxItemCount;
                    block += _b;
                    offset2 = (lower_pos - c1) % MaxItemCount;
                    read_block(data, block);
                    *BC += 1;
                    itd = (ItemD *) (data + offset2*ItemDSize);
                    remains = MaxItemCount - (offset2 + 1);
                    if (remains + pos + 1> node.num_items) {
                        remains = node.num_items - pos - 1;
                    }
                }
                // Profile Item Count
                *ic += 1;
                if (itd->tag == 2) {
                    pos = range_core_len(data, itd->comp.addr.block == block, false, results, pos,
                                         itd->comp.addr.block, itd->comp.addr.offset, lower, len, BC, ic ,lc);
                } else if (itd->tag == 3) {
                    if (itd->comp.data.key >= lower) {
                        results[pos] = itd->comp.data.key;
                        pos ++;
                    }
                }
                if (pos >= len) return pos;
                if (lower_pos + 1 >= node.num_items) return pos;

                int i = lower_pos + 1;
                itd += 1;
                for (;;) {
                    if (remains > (node.num_items - i)) {
                        remains = node.num_items  - i;
                    }
                    for (int j = 0; j < remains; j++) {
                        // Profile Item Count
                        *ic += 1;
                        if (itd[j].tag == 3) {
                            results[pos] = itd[j].comp.data.key;
                            pos ++;
                        } else if (itd[j].tag == 2) {
                            pos = range_core_len(data, itd[j].comp.addr.block == block, true, results, pos,
                                                 itd[j].comp.addr.block, itd[j].comp.addr.offset, lower, len, BC, ic, lc);
                        }
                        if (pos >= len) return pos;
                    }
                    block += 1;
                    i += remains;
                    remains = MaxItemCount;
                    if (i >= node.num_items) break;
                    else {
                        read_block(data, block);
                        *BC += 1;
                        itd = (ItemD *) data;
                    }
                }
                return pos;
            }
        }

        int predict_pos(NodeHeaderD nhd, T key) {
            double v = nhd.slope * static_cast<long double>(key) + nhd.intercept;
            int pos = 0;
            if (v > std::numeric_limits<int>::max() / 2) {
                pos = nhd.num_items - 1;
            } else if (v < 0) {
                pos = 0;
            } else {
                pos = std::min(nhd.num_items - 1, static_cast<int>(v));
            } 
            return pos;
        }

        void build_node_two(T key1, P value1, T key2, P value2, int *block, int *offset) {
            if (key1 > key2) {
                std::swap(key1, key2);
                std::swap(value1, value2);
            }
            NodeHeaderD node;
            node.is_two = 1;
            node.build_size = 2;
            node.size = 2;
            node.fixed = 0;
            node.number_inserts = node.num_insert_to_data = 0;
            node.num_items = 8;
            ItemD items[8];
            for (int i = 0; i < node.num_items; i++) {
                items[i].tag = 1;
            }
            const long double mid1_key = key1;
            const long double mid2_key = key2;

            const double mid1_target = node.num_items / 3;
            const double mid2_target = node.num_items * 2 / 3;
            node.slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
            node.intercept = mid1_target - node.slope * mid1_key;
            int pos = predict_pos(node, key1);
            items[pos].tag = 3;
            items[pos].comp.data.key = key1;
            items[pos].comp.data.value = value1;

            pos = predict_pos(node, key2);
            items[pos].tag = 3;
            items[pos].comp.data.key = key2;
            items[pos].comp.data.value = value2;
            _write_node_internal(node, items, block, offset);
            return;
        }

        void bulk_load_disk(T *_keys, P *_values, int _size, int *r_block = nullptr, int *r_offset = nullptr) {

            typedef struct {
                int segment_id;
                int begin;
                int end;
                int level; // top level = 1
                NodeHeaderD node;
            } Segment;
            std::stack<Segment> s;

            
            typedef struct {
                int parent_id;
                int i; //i-th child of its parent;
            } _map;

            typedef struct {
                int block;
                int offset;
            } _addr;
           

            std::vector<_map> mappings;
            //std::vector<_addr> addresses;
            std::map<int, _addr> addresses;
        
            int SEGMENTID = 0;

            NodeHeaderD ret;
            s.push((Segment){SEGMENTID, 0, _size, 1, ret});
            mappings.push_back(_map{-1, 0});
            SEGMENTID += 1;
            bool is_first = true;
            while (!s.empty()) {
                int _block = 0;
                int _offset = 0;
                const int begin = s.top().begin;
                const int end = s.top().end;
                const int level = s.top().level;
                NodeHeaderD node = s.top().node;
                int current_segid = s.top().segment_id;
                s.pop();
                if (end - begin == 2) {
                    build_node_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1], &_block, &_offset);
                    /////// TODO
                } else {
                    T* keys = _keys + begin;
                    P* values = _values + begin;
                    const int size = end - begin;
                    const int BUILD_GAP_CNT = compute_gap_count(size);

                    node.is_two = 0;
                    node.build_size = size;
                    node.size = size;
                    node.fixed = 0;
                    node.number_inserts = node.num_insert_to_data = 0;


                    // FMCD method
                    // Here the implementation is a little different with Algorithm 1 in our paper.
                    // In Algorithm 1, U_T should be (keys[size-1-D] - keys[D]) / (L - 2).
                    // But according to the derivation described in our paper, M.A should be less than 1 / U_T.
                    // So we added a small number (1e-6) to U_T.
                    // In fact, it has only a negligible impact of the performance.
                    {
                        const int L = size * static_cast<int>(BUILD_GAP_CNT + 1);
                        int i = 0;
                        int D = 1;
                        double Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                                    (static_cast<double>(L - 2)) + 1e-6;
                        while (i < size - 1 - D) {
                            while (i + D < size && keys[i + D] - keys[i] >= Ut) {
                                i ++;
                            }
                            if (i + D >= size) {
                                break;
                            }
                            D = D + 1;
                            if (D * 3 > size) break;
                            Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                                (static_cast<double>(L - 2)) + 1e-6;
                        }
                        if (D * 3 <= size) {

                            node.slope = 1.0 / Ut;
                            node.intercept = (L - node.slope * (static_cast<long double>(keys[size - 1 - D]) +
                                                                static_cast<long double>(keys[D]))) / 2;

                            node.num_items = L;
                        } else {
                            int mid1_pos = (size - 1) / 3;
                            int mid2_pos = (size - 1) * 2 / 3;

                            const long double mid1_key = (static_cast<long double>(keys[mid1_pos]) +
                                                        static_cast<long double>(keys[mid1_pos + 1])) / 2;
                            const long double mid2_key = (static_cast<long double>(keys[mid2_pos]) +
                                                        static_cast<long double>(keys[mid2_pos + 1])) / 2;

                            node.num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                            const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                            const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

                            node.slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                            node.intercept = mid1_target - node.slope * mid1_key;
                        }
                    }
                    const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
                    node.intercept+= lr_remains;
                    node.num_items += lr_remains * 2;

                    if (size > 1e6) {
                        node.fixed = 1;
                    }
                    ItemD * items = new ItemD[node.num_items];
                    for (int i = 0; i < node.num_items; i++) {
                        items[i].tag = 1;
                    }

                    for (int item_i = predict_pos(node, keys[0]), offset = 0; offset < size; ) {
                        int next = offset + 1, next_i = -1;
                        while (next < size) {
                            next_i = predict_pos(node, keys[next]);
                            if (next_i == item_i) {
                                next ++;
                            } else {
                                break;
                            }
                        }
                        if (next == offset + 1) {
                            items[item_i].tag = 3;
                            items[item_i].comp.data.key = keys[offset];
                            items[item_i].comp.data.value = values[offset];
                        } else {
                            items[item_i].tag = 2;
                            NodeHeaderD _nhd;
                            s.push((Segment){SEGMENTID, begin + offset, begin + next, level + 1, _nhd});
                            _map m;
                            m.parent_id = current_segid;
                            m.i = item_i;
                            mappings.push_back(m);
                            SEGMENTID += 1;
                        }
                        if (next >= size) {
                            break;
                        } else {
                            item_i = next_i;
                            offset = next;
                        }
                    }


                    _write_node_internal(node, items, &_block, &_offset);
                    

                }
                _addr add;
                add.block = _block;
                add.offset = _offset;
                addresses[current_segid] = add;
                if (r_block == nullptr && is_first) {
                    sys_metablock(true, _block, _offset);
                    is_first = false;
                } else if (r_block != nullptr && is_first){
                    *r_block = _block;
                    *r_offset = _offset;
                    is_first = false;
                }
                _map m = mappings.at(current_segid);
                if (m.parent_id != -1) {
                    _addr fadd = addresses[m.parent_id];
                    update_child_address(fadd.block, fadd.offset, m.i, _block, _offset);
                }
            }
        }

        void scan_and_destory_tree(NodeHeaderD _root, int block, int offset, T* keys, P* values) {
            typedef std::pair<int, NodeHeaderD> Segment;
            typedef std::pair<int, int> NodeAddress;
            std::stack<Segment> s;
            std::stack<NodeAddress> addr;

            s.push(Segment(0, _root));
            addr.push(NodeAddress(block, offset));
            char data[BLOCK_SIZE];
            char data2[BLOCK_SIZE];
            while (!s.empty()) {
                int begin = s.top().first;
                NodeHeaderD node = s.top().second;
                const int SHOULD_END_POS = begin + node.size;
                s.pop();
                int _block = addr.top().first;
                int _offset = addr.top().second;
                addr.pop();
                long _offset_ = (_block * BLOCK_SIZE + _offset + NodeHeaderDSize + node.empty_size2) % BLOCK_SIZE;
                int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                if (_offset_ == 0) {
                    _block += 1;
                }
                for (int i = 0; i < node.num_items;) {
                    read_block(data, _block);
                    ItemD *_items = (ItemD *) (data + _offset_);
                    if (c1 > (node.num_items - i)) {
                        c1 = node.num_items - i;
                    }
                    for (int j = 0; j < c1; j++) {
                        if (_items[j].tag == 3) {
                            keys[begin] = _items[j].comp.data.key;
                            values[begin] = _items[j].comp.data.value;
                            begin ++;
                        } else if (_items[j].tag == 2) {
                            int _b = _items[j].comp.addr.block;
                            int _o = _items[j].comp.addr.offset;
                            read_block(data2, _b);
                            NodeHeaderD _n;
                            memcpy(&_n, data2 + _o, NodeHeaderDSize);
                            s.push(Segment(begin, _n));
                            addr.push(NodeAddress(_b, _o));
                            begin += _n.size;
                        }
                    }
                
                    i += c1;
                    _block += 1;
                    _offset_ = 0;
                    c1 = MaxItemCount;
                }
                // TODO: reclaim the storage space
            }
        }
        typedef struct {
            int block;
            int offset;
        } Address;

        void count_root() {
            int num_items = 0;
            int num_nodes = 0;
            int insert_count = 0;
            count_one_node(mb.root_block, mb.root_offset, &num_items, &num_nodes, &insert_count);
            std::cout << num_items << std::endl;
            std::cout << num_nodes << std::endl;
            std::cout << insert_count << std::endl;
            // std::cout << count_one_node(mb.root_block, mb.root_offset) << std::endl;
        }

        void count_one_node(int block, int offset, int *NumItems, int *NumNodes, int *InsertCount) {
            char data[BLOCK_SIZE];
            char data2[BLOCK_SIZE];
            read_block(data2, block);
                //read one node header
            int _block = block;
            
            // int count  = 0;
            // int node_count = 1;
            NodeHeaderD *nhd = (NodeHeaderD *) (data2 + offset); 
            long _offset_ = (_block * BLOCK_SIZE + offset + NodeHeaderDSize + nhd->empty_size2) % BLOCK_SIZE;
            int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
            if (_offset_ == 0) {
                _block += 1;
            }
            
            *NumItems += nhd->num_items;
            *NumNodes += 1;
            for (int i = 0; i < nhd->num_items;) {
                read_block(data, _block);
                ItemD *_items = (ItemD *) (data + _offset_);
                if (c1 > (nhd->num_items - i)) {
                    c1 = nhd->num_items - i;
                }
                for (int j = 0; j < c1; j++) {
                    if (_items[j].tag == 2) {
                        count_one_node(_items[j].comp.addr.block, _items[j].comp.addr.offset, NumItems, NumNodes, InsertCount);
                    } else if (_items[j].tag == 3) {
                        *InsertCount += 1;
                    }
                }
                i += c1;
                _block += 1;
                _offset_ = 0;
                c1 = MaxItemCount;
            }
        }

        void insertv2(const T& key, const P& value, long long* search_latency, long long* insert_latency,
                      long long *smo_latency, long long *maintain_latency, int *smo, int *updated_tuples) {
#if Profiling
            *smo = 0;
            *updated_tuples = 0;
            *smo_latency = 0;
            std::chrono::high_resolution_clock::time_point search_start = std::chrono::high_resolution_clock::now();
            std::chrono::high_resolution_clock::time_point search_end = search_start;
            std::chrono::high_resolution_clock::time_point insert_start = std::chrono::high_resolution_clock::now();
            std::chrono::high_resolution_clock::time_point insert_end = insert_start;
            std::chrono::high_resolution_clock::time_point smo_start = std::chrono::high_resolution_clock::now();
            std::chrono::high_resolution_clock::time_point smo_end = smo_start;
            std::chrono::high_resolution_clock::time_point maintain_start = std::chrono::high_resolution_clock::now();
            std::chrono::high_resolution_clock::time_point maintain_end = maintain_start;
#endif

            int MAX_DEPTH = 128;
            NodeHeaderD path1[MAX_DEPTH];
            Address path2[MAX_DEPTH];
            int path_size = 0;
            int insert_to_data = 0;

            NodeHeaderD *nhd;
            ItemD *itd;
            char data[BLOCK_SIZE];
            int block1 = mb.root_block;
            int offset1 = mb.root_offset;
            int block2 = 0;
            int offset2 = 0;
#if Profiling
            search_start = std::chrono::high_resolution_clock::now();
#endif
            for (;;) {
                read_block(data, block1);
                //read one node header
                nhd = (NodeHeaderD *) (data + offset1);   
                path2[path_size].block = block1;
                path2[path_size].offset = offset1;
                memcpy(path1 + path_size, nhd, NodeHeaderDSize);
                path1[path_size].size ++;
                path1[path_size].number_inserts ++;
                path_size += 1;
                int pos = predict_pos(*nhd, key);

                // read item
                long _offset_ = (block1 * BLOCK_SIZE + offset1 + NodeHeaderDSize + nhd->empty_size2) % BLOCK_SIZE;
                int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                if (_offset_ == 0) {
                    c1 = 0;
                }
                if (pos < c1) {
                    offset2 = _offset_ + pos * ItemDSize;
                    itd = (ItemD *) (data + offset2);
                    block2 = block1;
                } else {
                    block2 = block1 + 1;
                    int _b = (pos - c1) / MaxItemCount;
                    block2 += _b;
                    offset2 = ((pos - c1) % MaxItemCount) * ItemDSize;
                    read_block(data, block2);
                    itd = (ItemD *) (data + offset2);
                }
                // handle insert
                if (itd->tag == 1) {
#if Profiling
                    search_end = std::chrono::high_resolution_clock::now();
                    insert_start = std::chrono::high_resolution_clock::now();
#endif
                    itd->tag = 3;
                    itd->comp.data.key = key;
                    itd->comp.data.value = value;
                    write_block(data, block2);
#if Profiling
                    insert_end = std::chrono::high_resolution_clock::now();
#endif
                    // printf("i1\n");
                    break;
                } else if (itd->tag == 3) {
#if Profiling
                    search_end = std::chrono::high_resolution_clock::now();
                    // create new node
                    *smo += 1;
                    smo_start = std::chrono::high_resolution_clock::now();
#endif
                    int _block = 0;
                    int _offset = 0;
                    build_node_two(itd->comp.data.key, itd->comp.data.value, key, value, &_block, &_offset);
                    // printf("i2\n");
                    itd->tag = 2;
                    itd->comp.addr.block = _block;
                    itd->comp.addr.offset = _offset;
                    // The new node may be in the same block here
                    // we just rewrite data back to disk, which will overwrite the new node info.
                    // write_block(data, block2); // we need to re-read the block and re-write the data
                    write_data(itd, block2 * BLOCK_SIZE + offset2, ItemDSize);
                    insert_to_data = 1;
#if Profiling
                    smo_end = std::chrono::high_resolution_clock::now();
                    *smo_latency += std::chrono::duration_cast<std::chrono::nanoseconds>(smo_end - smo_start).count();
#endif
                    break;
                } else if (itd->tag == 2){
                    block1 = itd->comp.addr.block;
                    offset1 = itd->comp.addr.offset;
                } else {
                    printf("wroing\n");
                }
            }

#if Profiling
            maintain_start = std::chrono::high_resolution_clock::now();
#endif
            int last_block = 0;
            for (int i = 0; i < path_size; i++) {
                path1[i].num_insert_to_data += insert_to_data;
                if (path2[i].block != last_block) {
                    if (i > 1)
                        write_block(data, last_block);
                    read_block(data, path2[i].block);
                    last_block =  path2[i].block;
                }
                memcpy(data + path2[i].offset, path1+i, NodeHeaderDSize);
            }
            if (last_block > 0)
                write_block(data, last_block);
#if Profiling
            maintain_end = std::chrono::high_resolution_clock::now();

            smo_start = std::chrono::high_resolution_clock::now();
#endif
            for (int i = 0; i < path_size; i++) {
                NodeHeaderD _nhd = path1[i];
                const int num_inserts = _nhd.number_inserts;
                const int num_insert_to_data = _nhd.num_insert_to_data;
                const bool need_rebuild = _nhd.fixed == 0 && _nhd.size >= _nhd.build_size * 4 && _nhd.size >= 64 && num_insert_to_data * 10 >= num_inserts;
                if (need_rebuild) {
                    *smo += 1;
                    // printf("re:%d;%d;%d\n",i,_nhd.size,path_size);
                    const int ESIZE = _nhd.size;
                    *updated_tuples = _nhd.size;
                    T* keys = new T[ESIZE];
                    P* values = new P[ESIZE];
                    scan_and_destory_tree(_nhd, path2[i].block, path2[i].offset, keys, values);
                    int _b = 0;
                    int _o = 0;
                    bulk_load_disk(keys, values, ESIZE, &_b, &_o);
                    delete[] keys;
                    delete[] values;

                    if (i > 0) {
                        int pos = predict_pos(path1[i-1], key);
                        int block = path2[i-1].block;
                        int offset = path2[i-1].offset;

                        long _offset_ = (block * BLOCK_SIZE + offset + NodeHeaderDSize + path1[i-1].empty_size2) % BLOCK_SIZE;
                        int c1 = (BLOCK_SIZE - _offset_)/ItemDSize;
                        if (_offset_ == 0) {
                            c1 = 0;
                        }
                        if (pos < c1) {
                            //itd = (ItemD *) (data + _offset_ + pos * ItemDSize);
                            _offset_ = _offset_ + pos * ItemDSize;
                            read_block(data, block);
                        } else {
                            int _b = (pos - c1) / MaxItemCount;
                            _offset_ = ((pos - c1) % MaxItemCount) * ItemDSize;
                            block = block + 1 + _b;
                            read_block(data, block);
                        }
                        ItemD* _it = (ItemD *)(data + _offset_);
                        _it->comp.addr.block = _b;
                        _it->comp.addr.offset = _o;
                        write_block(data, block);
                    } else {
                        sys_metablock(true, _b, _o);
                    }
                    break;
                }
            }
#if Profiling
            smo_end = std::chrono::high_resolution_clock::now();
            *search_latency = std::chrono::duration_cast<std::chrono::nanoseconds>(search_end - search_start).count();
            *insert_latency = std::chrono::duration_cast<std::chrono::nanoseconds>(insert_end - insert_start).count();
            *smo_latency += std::chrono::duration_cast<std::chrono::nanoseconds>(smo_end - smo_start).count();
            *maintain_latency = std::chrono::duration_cast<std::chrono::nanoseconds>(maintain_end - maintain_start).count();
#endif
            return;
        }
};