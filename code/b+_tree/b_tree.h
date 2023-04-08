/**
 * @file b_tree.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-02-12
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "storage_management.h"
#include <limits> 
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include "stx/btree_multimap.h"
#include "vector"
#include <chrono>

#define Profiling 1

class BTree {
    /// Start - In-main memory B+-tree for the case, all inners nodes are in main memory
    struct traits_inner : stx::btree_default_map_traits<KeyType, int>
    {
        static const bool       selfverify = false;
        static const bool       debug = false;

        static const int        leafslots = MaxItemInLeafNode;
        static const int        innerslots = MaxItemInInnerNode;
    };
    typedef stx::btree_multimap<KeyType, int, std::less<KeyType>, traits_inner> stx_btree;
    stx_btree inner_btree;
    /// END
    StorageManager *sm;
    MetaNode metanode;
    char *utility_file;

    #define ALL_DISK 0
    #define LEAF_DISK 1
    #define ROOT_MEMORY 2
    int hybrid_mode = ALL_DISK;
    private:
        void load_metanode() {
            Block block = sm->get_block(0);
            memcpy(&metanode, block.data, MetaNodeSize);
        }



        // return the first one less or equal than it...
        template<typename NodeIterm>
        int _search_in_node(NodeIterm *data, int item_count, KeyType key) {
            if (item_count < SearchThreshold) {
                for (int i = 0; i < item_count; i++) {
                    if (data[i].key >= key) {
                        return i-1;
                    }
                }
            }

            // binary search
            int l = 0;
            int r = item_count - 1;
            while (l <= r) {
                int mid = l + (r - l) / 2;
                if (data[mid].key < key) l = mid + 1;
                else r = mid - 1;
            }
            return l - 1;
        }

public:
    size_t get_inner_size() {
            if (hybrid_mode == LEAF_DISK) {
                return inner_btree.get_tree_size();
            } else
                return 0;
    }

    size_t get_file_size() {
            return sm->get_file_size();
        }
    void sync_metanode() {
        Block block;
        memcpy(block.data, &metanode, MetaNodeSize);
        sm->write_block(0, block);
    }
        BTree(bool first, char *file_name, bool bulk_load =false) {
            sm = new StorageManager(file_name, first, bulk_load);
            load_metanode();
        }

        // utility_file is used for all inner nodes are in main memory case to quickly
        // reconstruct the main memory tree, which is `main_file'+`.utility'
        BTree(int _hybrid_mode, bool first, char *main_file, bool bulk_load = false, char *_utility_file = nullptr) {
            hybrid_mode = _hybrid_mode;
            sm = new StorageManager(main_file, first, bulk_load);
            load_metanode();
            utility_file = _utility_file;
        }

        BTree() = default;
        ~BTree() {
            delete sm;
        }
        void print_next_block() {
            std::cout<<metanode.block_count<<std::endl;
        }
        template<typename NodeHeader, typename NodeIterm, int NodeHeaderSize = InnerNodeHeaderSize, 
                int ItermSize = InnerNodeItermSize, int MaxItermCount = MaxItemInInnerNode>
        InnerNodeIterm split_node(int POS, int block_id, Block block_to_split, NodeIterm *new_entry, bool is_right_most_inner = false) {
            //std::cout << "split" << std::endl;
            InnerNodeIterm ini;
            NodeHeader inh_l;
            NodeHeader inh_r;
            memcpy(&inh_l, block_to_split.data, NodeHeaderSize);
            NodeIterm inis[inh_l.item_count + 1];
            memcpy(&inis, block_to_split.data + NodeHeaderSize, inh_l.item_count * ItermSize);
            KeyType old_key;
            int old_right_most = inh_l.next_block_id;
            inh_r.next_block_id = old_right_most;
            if (POS > -1 && !is_right_most_inner) { // call by inner node case
                old_key = inis[POS].key;
               inis[POS].key = ((InnerNodeIterm *)new_entry)->key;
                ((InnerNodeIterm *)new_entry)->key = old_key; 
            } else if (is_right_most_inner) {
                inh_r.next_block_id =  ((InnerNodeIterm *)new_entry)->block_id;
                ((InnerNodeIterm *)new_entry)->block_id = old_right_most;
            }
            int ins_pos = _search_in_node<NodeIterm>(inis, inh_l.item_count, new_entry->key) + 1;

            int pivot = MaxItermCount / 2;
            inh_l.item_count = pivot + 1;

            
            inh_l.next_block_id = metanode.block_count;
            for (int i = MaxItermCount-1; i >= ins_pos; i--) {
                inis[i+1] = inis[i];

            }
            inis[ins_pos] = *new_entry;
            KeyType pivot_key = inis[pivot].key;
            
            memcpy(block_to_split.data, &inh_l, NodeHeaderSize);
            memcpy(block_to_split.data + NodeHeaderSize, inis, (pivot + 1) * ItermSize);
            sm->write_block(block_id, block_to_split);
            
            Block block_r;
            inh_r.node_type = inh_l.node_type;
            inh_r.item_count = MaxItermCount - pivot;
            inh_r.level = inh_l.level; 
            memcpy(block_r.data + NodeHeaderSize, inis + pivot+1, inh_r.item_count * ItermSize);
            memcpy(block_r.data, &inh_r, NodeHeaderSize);
            sm->write_block(metanode.block_count, block_r);
            ini.block_id = metanode.block_count;
            ini.key = pivot_key;
            metanode.block_count += 1;
//            sync_metanode();
            return ini;
        }

        template<typename NodeHeader, typename NodeIterm, int NodeHeaderSize = InnerNodeHeaderSize, int NodeItermSize = InnerNodeItermSize>
        void add_entry_in_node(int POS, int block_id, Block block_to_append, NodeIterm *ini, bool is_right_most_inner = false) {
            NodeHeader inh;
            KeyType old_key;
            
            memcpy(&inh, block_to_append.data, NodeHeaderSize);
            NodeIterm inis[inh.item_count+1];
            memcpy(&inis, block_to_append.data + NodeHeaderSize, NodeItermSize * inh.item_count);
            if (POS > -1 && !is_right_most_inner) {
                old_key = inis[POS].key;
                inis[POS].key = ((InnerNodeIterm *)ini)->key;
                ((InnerNodeIterm *)ini)->key = old_key; 
            } else if (is_right_most_inner) {
                int temp = inh.next_block_id;
                inh.next_block_id =  ((InnerNodeIterm *)ini)->block_id;
                ((InnerNodeIterm *)ini)->block_id = temp;
                memcpy(block_to_append.data + NodeHeaderSize + NodeItermSize * inh.item_count, ini, NodeItermSize);
                inh.item_count += 1;
                memcpy(block_to_append.data, &inh, NodeHeaderSize);
                sm->write_block(block_id, block_to_append);
                return;
            }
            int ins_pos = _search_in_node<NodeIterm>(inis, inh.item_count, ini->key) + 1;
            for (int i = inh.item_count-1; i >= ins_pos; i--) {
                inis[i + 1] = inis[i];
            }
            inis[ins_pos] = *ini;
            inh.item_count += 1;
            memcpy(block_to_append.data + NodeHeaderSize, inis, NodeItermSize * inh.item_count); 
            memcpy(block_to_append.data, &inh, NodeHeaderSize);
            sm->write_block(block_id, block_to_append);
            return;
        }


#if Profiling
    long long SEARCH_USED_TIME = 0;
    long long SMO_USED_TIME = 0;
    long long INSERT_USED_TIME = 0;
    int SMO_COUNT = 0;
    int SMO_UPDATE = 0;
#endif
        BuildStatus insert_key(int current_block_id, LeafNodeIterm lni_item) {
#if Profiling
            std::chrono::high_resolution_clock::time_point i_start = std::chrono::high_resolution_clock::now();
#endif
            BuildStatus bs;
            Block block;
            block = sm->get_block(current_block_id);
            bool is_inner = block.data[0] == InnerNodeType;
            if (is_inner) {
                InnerNodeHeader *inh = (InnerNodeHeader *)(block.data);
                // memcpy(&inh, block.data, InnerNodeHeaderSize);

                // InnerNodeIterm inis[inh.item_count];
                InnerNodeIterm *inis = (InnerNodeIterm *)(block.data + InnerNodeHeaderSize);
                // memcpy(&inis, block.data + InnerNodeHeaderSize, InnerNodeItermSize * inh.item_count);
                // int ins_pos = _search_in_node<InnerNodeIterm>(inis, inh.item_count, lni_item.key);
                int ins_pos = _search_in_node<InnerNodeIterm>(inis, inh->item_count, lni_item.key);
#if Profiling
                std::chrono::high_resolution_clock::time_point i_end = std::chrono::high_resolution_clock::now();
                SEARCH_USED_TIME += std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start).count();
#endif
                if (ins_pos == inh->item_count - 1) {
                    bs = insert_key(inh->next_block_id, lni_item);
                }  else {
                    bs = insert_key(inis[ins_pos + 1].block_id, lni_item);
                }
#if Profiling
                std::chrono::high_resolution_clock::time_point ismo_start = std::chrono::high_resolution_clock::now();
#endif
                bool is_right_most_inner = ins_pos == inh->item_count - 1;
                int level = inh->level; // in case of block update
                if (bs.status == AddNewEntry) {
                    if (inh->item_count == MaxItemInInnerNode) { // we need to split
#if Profiling
                        SMO_COUNT += 1;
                        SMO_UPDATE += inh->item_count;
#endif
                        InnerNodeIterm to_add_iterm = split_node<InnerNodeHeader, InnerNodeIterm>(ins_pos + 1, current_block_id, block, &bs.ini, is_right_most_inner);
                        bs.ini = to_add_iterm;
                        bs.added_block += 1;
                    } else { // we insert directly
                        add_entry_in_node<InnerNodeHeader, InnerNodeIterm>(ins_pos + 1, current_block_id, block, &(bs.ini), is_right_most_inner);
                        bs.status = NoOperation;
                        //bs.added_block = 0;
                    }
                }
#if Profiling
                std::chrono::high_resolution_clock::time_point ismo_end = std::chrono::high_resolution_clock::now();
                SMO_USED_TIME += std::chrono::duration_cast<std::chrono::nanoseconds>(ismo_end - ismo_start).count();
#endif
                bs.level = level;


            } else {
#if Profiling
                std::chrono::high_resolution_clock::time_point i_end = std::chrono::high_resolution_clock::now();
                SEARCH_USED_TIME += std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start).count();
#endif

                LeaftNodeHeader *lfn_l = (LeaftNodeHeader *)(block.data);
                int level = lfn_l->level;
                if (lfn_l->item_count == MaxItemInLeafNode) {
#if Profiling
                    SMO_COUNT += 1;
                    SMO_UPDATE += lfn_l->item_count;
                    std::chrono::high_resolution_clock::time_point ii_start = std::chrono::high_resolution_clock::now();
#endif
                    InnerNodeIterm to_add_iterm = split_node<LeaftNodeHeader, LeafNodeIterm, LeaftNodeHeaderSize, LeafNodeItemSize, MaxItemInLeafNode>(-1, current_block_id, block, &lni_item);
#if Profiling
                    std::chrono::high_resolution_clock::time_point ii_end = std::chrono::high_resolution_clock::now();
                    SMO_USED_TIME += std::chrono::duration_cast<std::chrono::nanoseconds>(ii_end - ii_start).count();
#endif
                    bs.ini = to_add_iterm;
                    bs.status = AddNewEntry;
                    bs.added_block = 1;
                } else { // insert directly
#if Profiling
                    std::chrono::high_resolution_clock::time_point ii_start = std::chrono::high_resolution_clock::now();
#endif
                    add_entry_in_node<LeaftNodeHeader, LeafNodeIterm, LeaftNodeHeaderSize, LeafNodeItemSize>(-1, current_block_id, block, &(lni_item));
                    bs.status = NoOperation;
                    bs.added_block = 0;
#if Profiling
                    std::chrono::high_resolution_clock::time_point ii_end = std::chrono::high_resolution_clock::now();
                    INSERT_USED_TIME += std::chrono::duration_cast<std::chrono::nanoseconds>(ii_end - ii_start).count();
#endif
                }
                bs.level = level;

            }
            return bs;
        }

        int insert_key_leaf_disk(LeafNodeIterm lni) {
            // find the leaf node to insert the key
            stx_btree::iterator it = inner_btree.find_for_disk(lni.key);
            int block_id = metanode.last_block;
            bool special_case = false;
            KeyType indexed_key;
            if (it != inner_btree.end()) {
                block_id = it.data();
                indexed_key = it.key();
            } else
                special_case = true;
            // do insert
            Block block;
            sm->get_block(block_id, block.data);
            LeaftNodeHeader *lnh = (LeaftNodeHeader *)(block.data);

            if (lnh->item_count == MaxItemInLeafNode) {
                InnerNodeIterm to_add_iterm = split_node<LeaftNodeHeader, LeafNodeIterm, LeaftNodeHeaderSize, LeafNodeItemSize, MaxItemInLeafNode>(-1, block_id, block, &lni);
                if (special_case) {
                    inner_btree.insert(to_add_iterm.key, metanode.last_block);
                    metanode.last_block = to_add_iterm.block_id;
                } else {
//                    inner_btree.erase_one(indexed_key);
                    it.set_value(to_add_iterm.block_id);
//                    inner_btree.insert(indexed_key, to_add_iterm.block_id);
                    inner_btree.insert(to_add_iterm.key, block_id);
                }
                return 1;
            } else { // insert directly
                add_entry_in_node<LeaftNodeHeader, LeafNodeIterm, LeaftNodeHeaderSize, LeafNodeItemSize>(-1, block_id, block, &(lni));
                return 0;
            }
        }

        int insert_key_disk(LeafNodeIterm lni) {
            BuildStatus bs = insert_key(metanode.root_block_id, lni);
#if Profiling
            std::chrono::high_resolution_clock::time_point i_start = std::chrono::high_resolution_clock::now();
#endif
            if (bs.status == AddNewEntry) { // we need to create a new root node
#if Profiling
                SMO_COUNT += 1;
                SMO_UPDATE += 1;
#endif
                Block block;
                InnerNodeHeader inh;
                inh.item_count = 1;
                // The last one of one level points to child node while others point to the sibling.
                inh.next_block_id = bs.ini.block_id;
                inh.node_type = InnerNodeType;
                inh.level = bs.level + 1;
                memcpy(block.data, &inh, InnerNodeHeaderSize);
                bs.ini.block_id = metanode.root_block_id;
                memcpy(block.data + InnerNodeHeaderSize, &bs.ini, InnerNodeItermSize);
                sm->write_block(metanode.block_count, block);
                metanode.root_block_id = metanode.block_count;
                metanode.block_count += 1;
                metanode.level = inh.level;
                bs.added_block += 1;

//                sync_metanode();
            }
#if Profiling
            std::chrono::high_resolution_clock::time_point i_end = std::chrono::high_resolution_clock::now();
            SMO_USED_TIME += std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start).count();
#endif
            return bs.added_block;
        }


        int insert_key_entry(KeyType key, ValueType value, long long *search_time, long long *smo_time,
                             long long *insert_time, int *smo_c, int *update_c) {
#if Profiling
            SEARCH_USED_TIME = 0;
            SMO_USED_TIME = 0;
            INSERT_USED_TIME = 0;
            SMO_COUNT = 0;
            SMO_UPDATE = 0;
#endif
            LeafNodeIterm lni;
            lni.key = key; lni.value = value;
            if (hybrid_mode == ALL_DISK) insert_key_disk(lni);
            else if (hybrid_mode == LEAF_DISK) insert_key_leaf_disk(lni);
#if Profiling
            *search_time = SEARCH_USED_TIME;
            *smo_time = SMO_USED_TIME;
            *insert_time = INSERT_USED_TIME;
            *update_c = SMO_UPDATE;
            *smo_c = SMO_COUNT;
#endif
            return 0;
        }

        typedef struct Condition
        {
            KeyType min;
            bool include_min;
            KeyType max;
            bool include_max;
        } Condition;
        

        class IndexIterator {
            //char *data = nullptr;
            LeaftNodeHeader *current_leaf_node_header;
            LeafNodeIterm *current_leaf_node_items;
            int next_point;
            Condition cond;
            StorageManager *sm;

            public:
                char *data = nullptr;
                IndexIterator() = default;
                IndexIterator(LeaftNodeHeader *start_header, LeafNodeIterm *start_iterms, int start_point, Condition _cond, StorageManager *_sm) {
                    current_leaf_node_header = start_header;
                    current_leaf_node_items = start_iterms;
                    next_point = start_point;
                    cond = _cond;
                    sm = _sm;
                }

                IndexIterator(Condition _cond, StorageManager *_sm) {
                    data = (char *) malloc(BlockSize);
                    cond = _cond;
                    sm = _sm;
                    current_leaf_node_header = (LeaftNodeHeader *)(data);
                    current_leaf_node_items = (LeafNodeIterm *)(data + LeaftNodeHeaderSize);
                }

                ~IndexIterator() {
                    if (data != nullptr) free(data);
                }

                void set_start(int start) {
                    next_point = start;
                }

                bool has_next() {
                    if (next_point == current_leaf_node_header->item_count) {
                        sm->get_block(current_leaf_node_header->next_block_id, data);
                        next_point = 0;
                    }
                    bool min_c = cond.include_min ? current_leaf_node_items[next_point].key >= cond.min : current_leaf_node_items[next_point].key > cond.min;
                    bool max_c = cond.include_max ? current_leaf_node_items[next_point].key <= cond.max : current_leaf_node_items[next_point].key < cond.max;
                    return min_c && max_c;
                }

                bool has_next_range(int *c) {
                    if (next_point == current_leaf_node_header->item_count) {
                        sm->get_block(current_leaf_node_header->next_block_id, data);
                        *c += 1;
                        next_point = 0;
                    }
                    return current_leaf_node_items[next_point].key >= cond.min;
                    
                }
                bool has_next(int *c) {
                    if (next_point == current_leaf_node_header->item_count) {
                        sm->get_block(current_leaf_node_header->next_block_id, data);
                        *c += 1;
                        next_point = 0;
                    }
                    bool min_c = cond.include_min ? current_leaf_node_items[next_point].key >= cond.min : current_leaf_node_items[next_point].key > cond.min;
                    bool max_c = cond.include_max ? current_leaf_node_items[next_point].key <= cond.max : current_leaf_node_items[next_point].key < cond.max;
                    return min_c && max_c;
                }

                ValueType next() {
                    ValueType v = current_leaf_node_items[next_point].value;
                    next_point += 1;
                    return v;
                }
        };

        IndexIterator get_index_iterator(Condition cond, int *c) {
            IndexIterator ii(cond, sm);
            // find the first block
            int current_block_id = metanode.root_block_id;
            // Block block = sm->get_block(current_block_id);
            
            sm->get_block(current_block_id, ii.data);
            *c += 1;

            bool is_inner = ii.data[0] == InnerNodeType;
            // traversal to leaf node
            while (is_inner) {
                InnerNodeHeader *inh;
                inh = (InnerNodeHeader *)(ii.data);
                InnerNodeIterm* inis = (InnerNodeIterm*) (ii.data + InnerNodeHeaderSize); //[inh->item_count];

                int pos = _search_in_node<InnerNodeIterm>(inis, inh->item_count, cond.min);
                if (pos == inh->item_count - 1) current_block_id = inh->next_block_id;
                else if (pos > -1 && inis[pos].key == cond.min) current_block_id = inis[pos].block_id; // return is the last one <=
                else current_block_id = inis[pos + 1].block_id;
                sm->get_block(current_block_id, ii.data);

                *c += 1;

                is_inner = ii.data[0] == InnerNodeType;
            }

            LeaftNodeHeader *lnh = (LeaftNodeHeader *)(ii.data);
            LeafNodeIterm *lnis = (LeafNodeIterm *) (ii.data + LeaftNodeHeaderSize);
            int pos = _search_in_node<LeafNodeIterm>(lnis, lnh->item_count, cond.min);
            int start = lnis[pos].key == cond.min ? pos : pos + 1;
            ii.set_start(start);
            return ii;
        }

        bool lookup_disk(Condition cond, int *c) {
            IndexIterator idx_it = get_index_iterator(cond, c);
            bool found = false;
            while (idx_it.has_next(c)) {
                idx_it.next();
                found = true;
            }
            return found;
        }

        bool bookup_leaf_disk(Condition cond, int *c) {
            IndexIterator ii = obtain_for_leaf_disk(cond, c);
            // scan forward
            bool found = false;
            while (ii.has_next(c)) {
                ii.next();
                found = true;
            }
            return found;
        }

        bool lookup(KeyType key, int *c) {
            Condition cond;
            cond.include_max = cond.include_min = true;
            cond.max = cond.min = key;
            *c = 0;
            if (hybrid_mode == ALL_DISK) return lookup_disk(cond, c);
            else if (hybrid_mode == LEAF_DISK) return bookup_leaf_disk(cond, c);
        }

        bool scan_disk(ValueType *results, Condition cond, int len, int *c) {
            bool found = false;
            IndexIterator idx_it = get_index_iterator(cond, c);
            int pos = 0;
            while (idx_it.has_next_range(c) && pos < len) {
                results[pos] = idx_it.next();
                pos ++;
                found = true;
            }
            return found;
        }

        IndexIterator obtain_for_leaf_disk(Condition cond, int *c) {
            // search in the inner part
            stx_btree::iterator it = inner_btree.find_for_disk(cond.min);
            int block_id = metanode.last_block;
            if (it != inner_btree.end())
                block_id = it.data();

            // build iterator
            IndexIterator ii(cond, sm);
            sm->get_block(block_id, ii.data);
            *c += 1;

            // find the start position
            LeaftNodeHeader *lnh = (LeaftNodeHeader *)(ii.data);
            LeafNodeIterm *lnis = (LeafNodeIterm *) (ii.data + LeaftNodeHeaderSize);
            int pos = _search_in_node<LeafNodeIterm>(lnis, lnh->item_count, cond.min);
            int start = lnis[pos].key == cond.min ? pos : pos + 1;
            ii.set_start(start);
            return ii;
        }

        bool scan_leaf_disk(ValueType *results, Condition cond, int len, int *c) {
            IndexIterator ii = obtain_for_leaf_disk(cond, c);
            int pos = 0;
            bool found = false;
            while (ii.has_next_range(c) && pos < len) {
                results[pos] = ii.next();
                pos ++;
                found = true;
            }
            return found;
        }

        bool scan(ValueType *results, KeyType lower_key, int len, int *c) {
            Condition cond;
            cond.include_min = true;
            cond.min = lower_key;
            *c = 0;
            if (hybrid_mode == ALL_DISK) return scan_disk(results, cond, len, c);
            else if (hybrid_mode == LEAF_DISK)  return scan_leaf_disk(results, cond, len, c);
        }

        void bulk_load(LeafNodeIterm *data, long item_count, double per=0.8) {
            // sort the data
            // build from bottom to up
            int INSERT_SIZE_LEAF = int(MaxItemInLeafNode*per);
            int _c = item_count / INSERT_SIZE_LEAF;
            int _m = item_count % INSERT_SIZE_LEAF;
            int leaf_node_count = _m == 0 ? _c : _c + 1;
            int valid_item_count = leaf_node_count;
            InnerNodeIterm *inis = new InnerNodeIterm [leaf_node_count];
            std::cout << "start leaf" << std::endl;
            std::vector< std::pair<KeyType, int> > pairs (leaf_node_count-1);
            for (int i = 0; i < leaf_node_count; i++) {
                LeaftNodeHeader lnh;
                Block block;
                lnh.item_count = _m > 0 && i == leaf_node_count  - 1 ? _m : INSERT_SIZE_LEAF;
                lnh.level = 1;
                lnh.next_block_id = metanode.block_count + 1;
                lnh.node_type = LeafNodeType;
                memcpy(block.data, &lnh, LeaftNodeHeaderSize);
                memcpy(block.data + LeaftNodeHeaderSize, data + i * INSERT_SIZE_LEAF, LeafNodeItemSize * lnh.item_count);
                sm->write_block(metanode.block_count, block);
                inis[i].block_id = metanode.block_count;
                inis[i].key = data[i * INSERT_SIZE_LEAF + lnh.item_count - 1].key;

                // preprocessing for leaf_disk mode
                if (i == leaf_node_count - 1) {
                    metanode.last_block = metanode.block_count;
                }
                else {// we do not store the last one in the inner nodes
                    pairs[i].first = data[i * INSERT_SIZE_LEAF + lnh.item_count - 1].key;
                    pairs[i].second = metanode.block_count;
                }

                metanode.block_count += 1;
//                sync_metanode();
            }
            std::cout << "end leaf" << std::endl;
            // inner nodes in main memory case
            if (hybrid_mode == LEAF_DISK) {
                std::cout << "build inner" << std::endl;
                inner_btree.bulk_load(pairs.begin(), pairs.end());
                return;
            }
            int next_level = 2;
            // start to build upper levels
            int INSERT_SIZE_INNER = int(MaxItemInInnerNode*per);
            while (valid_item_count > 1) {
                // the last child block will be put as the next_block_id of the last one
                _c = (valid_item_count - 1) / INSERT_SIZE_INNER;
                _m = (valid_item_count - 1) % INSERT_SIZE_INNER;
                int node_count = _m == 0 ? _c : _c + 1;
                for (int i = 0; i < node_count; i++) {
                    InnerNodeHeader inh;
                    inh.item_count = _m > 0 && i == node_count  - 1 ? _m : INSERT_SIZE_INNER;
                    inh.level = next_level;
                    inh.next_block_id = ((i == node_count - 1) ? inis[valid_item_count - 1].block_id : metanode.block_count + 1);
                    inh.node_type = InnerNodeType;

                    Block block;
                    memcpy(block.data, &inh, InnerNodeHeaderSize);
                    memcpy(block.data + InnerNodeHeaderSize, inis + i * INSERT_SIZE_INNER, InnerNodeItermSize * inh.item_count);
                    // if (metanode.block_count == 196374) {
                    //     std::cout << "next_block_id:" << inh.next_block_id << std::endl;
                    // }
                    sm->write_block(metanode.block_count, block);
                    inis[i].block_id = metanode.block_count;
                    inis[i].key = inis[i * INSERT_SIZE_INNER + inh.item_count - 1].key;
                    metanode.block_count += 1;
//                    sync_metanode();
                }
                valid_item_count = node_count;
                next_level += 1;
            }
            metanode.level = (next_level - 1);
            metanode.root_block_id = inis[0].block_id;
            delete []inis;
        }

        bool delete_key(KeyType key) {
            return true;
        }

        void print_tree_level() {
            printf("level: %d\n", metanode.level);
        }
};
