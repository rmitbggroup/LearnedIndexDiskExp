// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

/*
 * ALEX with key type T and payload type P, combined type V=std::pair<T, P>.
 * Iterating through keys is done using an "Iterator".
 * Iterating through tree nodes is done using a "NodeIterator".
 *
 * Core user-facing API of Alex:
 * - Alex()
 * - void bulk_load(V values[], int num_keys)
 * - void insert(T key, P payload)
 * - int erase_one(T key)
 * - int erase(T key)
 * - Iterator find(T key)  // for exact match
 * - Iterator begin()
 * - Iterator end()
 * - Iterator lower_bound(T key)
 * - Iterator upper_bound(T key)
 *
 * User-facing API of Iterator:
 * - void operator ++ ()  // post increment
 * - V operator * ()  // does not return reference to V by default
 * - const T& key ()
 * - P& payload ()
 * - bool is_end()
 * - bool operator == (const Iterator & rhs)
 * - bool operator != (const Iterator & rhs)
 */

#pragma once

#include <fstream>
#include <iostream>
#include <stack>
#include <type_traits>

#include "alex_base.h"
#include "alex_fanout_tree.h"
#include "alex_nodes.h"
#include "storage_management.h"

#define DiskSetting 1
#define _Debug 0
#define Profiling 1
// Whether we account for floating-point precision issues when traversing down
// ALEX.
// These issues rarely occur in practice but can cause incorrect behavior.
// Turning this on will cause slight performance overhead due to extra
// computation and possibly accessing two data nodes to perform a lookup.
#define ALEX_SAFE_LOOKUP 1
#define MaxModeCount (BlockSize / PhyscialAddrSize)
namespace alex {

template <class T, class P, class Compare = AlexCompare,
          class Alloc = std::allocator<std::pair<T, P>>,
          bool allow_duplicates = true>
class Alex {
  static_assert(std::is_arithmetic<T>::value, "ALEX key type must be numeric.");
  static_assert(std::is_same<Compare, AlexCompare>::value,
                "Must use AlexCompare.");

 public:
  // Value type, returned by dereferencing an iterator
  typedef std::pair<T, P> V;

  // ALEX class aliases
  typedef Alex<T, P, Compare, Alloc, allow_duplicates> self_type;
  typedef AlexModelNode<T, P, Alloc> model_node_type;
  typedef AlexDataNode<T, P, Compare, Alloc, allow_duplicates> data_node_type;

  // Forward declaration for iterators
  class Iterator;
  class ConstIterator;
  class ReverseIterator;
  class ConstReverseIterator;
  class NodeIterator;  // Iterates through all nodes with pre-order traversal

  AlexNode<T, P>* root_node_ = nullptr;
  model_node_type* superroot_ =
      nullptr;  // phantom node that is the root's parent

  /* User-changeable parameters */
  struct Params {
    // When bulk loading, Alex can use provided knowledge of the expected
    // fraction of operations that will be inserts
    // For simplicity, operations are either point lookups ("reads") or inserts
    // ("writes)
    // i.e., 0 means we expect a read-only workload, 1 means write-only
    double expected_insert_frac = 1;
    // Maximum node size, in bytes. By default, 16MB.
    // Higher values result in better average throughput, but worse tail/max
    // insert latency
    int max_node_size = 1 << 24;
    // Approximate model computation: bulk load faster by using sampling to
    // train models
    bool approximate_model_computation = true;
    // Approximate cost computation: bulk load faster by using sampling to
    // compute cost
    bool approximate_cost_computation = false;
  };
  Params params_;

  /* Setting max node size automatically changes these parameters */
  struct DerivedParams {
    // The defaults here assume the default max node size of 16MB
    int max_fanout = 1 << 21;  // assumes 8-byte pointers
    int max_data_node_slots = (1 << 24) / sizeof(V);
//      int max_data_node_slots = (1 << 15) / sizeof(V);

  };
  DerivedParams derived_params_;

  /* Counters, useful for benchmarking and profiling */
  struct Stats {
    int num_keys = 0;
    int num_model_nodes = 0;  // num model nodes
    int num_data_nodes = 0;   // num data nodes
    int num_expand_and_scales = 0;
    int num_expand_and_retrains = 0;
    int num_downward_splits = 0;
    int num_sideways_splits = 0;
    int num_model_node_expansions = 0;
    int num_model_node_splits = 0;
    long long num_downward_split_keys = 0;
    long long num_sideways_split_keys = 0;
    long long num_model_node_expansion_pointers = 0;
    long long num_model_node_split_pointers = 0;
    mutable long long num_node_lookups = 0;
    mutable long long num_lookups = 0;
    long long num_inserts = 0;
    double splitting_time = 0;
    double cost_computation_time = 0;
  };
  Stats stats_;

  /* These are for research purposes, a user should not change these */
  struct ExperimentalParams {
    // Fanout selection method used during bulk loading: 0 means use bottom-up
    // fanout tree, 1 means top-down
    int fanout_selection_method = 0;
    // Policy when a data node experiences significant cost deviation.
    // 0 means always split node in 2
    // 1 means decide between no splitting or splitting in 2
    // 2 means use a full fanout tree to decide the splitting strategy
    int splitting_policy_method = 1;
    // Splitting upwards means that a split can propagate all the way up to the
    // root, like a B+ tree
    // Splitting upwards can result in a better RMI, but has much more overhead
    // than splitting sideways
    bool allow_splitting_upwards = false;
  };
  ExperimentalParams experimental_params_;

  /* Structs used internally */

 private:
  /* Statistics related to the key domain.
   * The index can hold keys outside the domain, but lookups/inserts on those
   * keys will be inefficient.
   * If enough keys fall outside the key domain, then we expand the key domain.
   */
  struct InternalStats {
    T key_domain_min_ = std::numeric_limits<T>::max();
    T key_domain_max_ = std::numeric_limits<T>::lowest();
    int num_keys_above_key_domain = 0;
    int num_keys_below_key_domain = 0;
    int num_keys_at_last_right_domain_resize = 0;
    int num_keys_at_last_left_domain_resize = 0;
    T min_key_in_left = std::numeric_limits<T>::max();
    T max_key_in_right = std::numeric_limits<T>::lowest();
  };
  InternalStats istats_;

  /* Save the traversal path down the RMI by having a linked list of these
   * structs. */
  struct TraversalNode {
    model_node_type* node = nullptr;
    int bucketID = -1;
  };

  /* Used when finding the best way to propagate up the RMI when splitting
   * upwards.
   * Cost is in terms of additional model size created through splitting
   * upwards, measured in units of pointers.
   * One instance of this struct is created for each node on the traversal path.
   * User should take into account the cost of metadata for new model nodes
   * (base_cost). */
  struct SplitDecisionCosts {
    static constexpr double base_cost =
        static_cast<double>(sizeof(model_node_type)) / sizeof(void*);
    // Additional cost due to this node if propagation stops at this node.
    // Equal to 0 if redundant slot exists, otherwise number of new pointers due
    // to node expansion.
    double stop_cost = 0;
    // Additional cost due to this node if propagation continues past this node.
    // Equal to number of new pointers due to node splitting, plus size of
    // metadata of new model node.
    double split_cost = 0;
  };

  // At least this many keys must be outside the domain before a domain
  // expansion is triggered.
  static const int kMinOutOfDomainKeys = 5;
  // After this many keys are outside the domain, a domain expansion must be
  // triggered.
  static const int kMaxOutOfDomainKeys = 1000;
  // When the number of max out-of-domain (OOD) keys is between the min and
  // max, expand the domain if the number of OOD keys is greater than the
  // expected number of OOD due to randomness by greater than the tolereance
  // factor.
  static const int kOutOfDomainToleranceFactor = 2;

  Compare key_less_ = Compare();
  Alloc allocator_ = Alloc();

  /*** Constructors and setters ***/
 public:
  StorageManager *sm;
  StorageManager *data_node_sm;
  MetaNode metanode;
  MetaNode metanode_data_node;
  void load_metanode() {
    Block block = sm->get_block(0);
    memcpy(&metanode, block.data, MetaNodeSize);
    memcpy(&stats_, block.data + MetaNodeSize, sizeof(Stats));
    block = data_node_sm->get_block(0);
    memcpy(&metanode_data_node, block.data, MetaNodeSize);
  }

  void sync_metanode(bool is_inner = true) {
      Block block;
      if (is_inner) {
        memcpy(block.data, &metanode, MetaNodeSize);
        memcpy(block.data + MetaNodeSize, &stats_, sizeof(Stats));
        sm->write_block(0, block);
      } else {
        memcpy(block.data, &metanode_data_node, MetaNodeSize);
        data_node_sm->write_block(0, block);
      }
  }

  Alex() {
    // Set up root as empty data node
    auto empty_data_node = new (data_node_allocator().allocate(1))
        data_node_type(key_less_, allocator_);
    empty_data_node->bulk_load(nullptr, 0);
    root_node_ = empty_data_node;
    stats_.num_data_nodes++;
    create_superroot();
  }

  Alex(bool is_first, char * index_name, char *data_file) {
    auto empty_data_node = new (data_node_allocator().allocate(1))
        data_node_type(key_less_, allocator_);
    empty_data_node->bulk_load(nullptr, 0);
    root_node_ = empty_data_node;
    stats_.num_data_nodes++;
    create_superroot();
    sm = new StorageManager(index_name, is_first);
    data_node_sm = new StorageManager(data_file, is_first);
    load_metanode();
  }

  #define ALL_DISK 0
  #define LEAF_DISK 1
  int hybrid_mode = ALL_DISK;
  Alex(int _hybrid_mode, bool is_first, char * index_name, char *data_file) {
      auto empty_data_node = new (data_node_allocator().allocate(1))
                data_node_type(key_less_, allocator_);
      empty_data_node->bulk_load(nullptr, 0);
      root_node_ = empty_data_node;
      stats_.num_data_nodes++;
      create_superroot();
        sm = new StorageManager(index_name, is_first);
      // when we test hybrid mode, every time we will start again.
      data_node_sm = new StorageManager(data_file, is_first);
      load_metanode();
      hybrid_mode = _hybrid_mode;
  }

    size_t get_file_size() {
      if (hybrid_mode == ALL_DISK) {
          return sm->get_file_size() + data_node_sm->get_file_size();
      } else {
          return data_node_sm->get_file_size();
      }
    }

    size_t get_memory_size() {
      if (hybrid_mode == ALL_DISK) {
          return 0;
      } else {
          size_t total_size = 0;
          // NOTE: next() function has been updated
          for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
               node_it.next()) {
              if (!node_it.current()->is_leaf_ && hybrid_mode == LEAF_DISK)
                 total_size += node_it.current()->node_size();
          }
          return total_size;
      }
  }

  Alex(const Compare& comp, const Alloc& alloc = Alloc())
      : key_less_(comp), allocator_(alloc) {
    // Set up root as empty data node
    auto empty_data_node = new (data_node_allocator().allocate(1))
        data_node_type(key_less_, allocator_);
    empty_data_node->bulk_load(nullptr, 0);
    root_node_ = empty_data_node;
    stats_.num_data_nodes++;
    create_superroot();
  }

  Alex(const Alloc& alloc) : allocator_(alloc) {
    // Set up root as empty data node
    auto empty_data_node = new (data_node_allocator().allocate(1))
        data_node_type(key_less_, allocator_);
    empty_data_node->bulk_load(nullptr, 0);
    root_node_ = empty_data_node;
    stats_.num_data_nodes++;
    create_superroot();
  }

  ~Alex() {
#if DiskSetting
      delete sm;
      delete data_node_sm;
#else
    for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
         node_it.next()) {
      delete_node(node_it.current());
    }
    delete_node(superroot_);
#endif
  }

  // Initializes with range [first, last). The range does not need to be
  // sorted. This creates a temporary copy of the data. If possible, we
  // recommend directly using bulk_load() instead.
  template <class InputIterator>
  explicit Alex(InputIterator first, InputIterator last, const Compare& comp,
                const Alloc& alloc = Alloc())
      : key_less_(comp), allocator_(alloc) {
    std::vector<V> values;
    for (auto it = first; it != last; ++it) {
      values.push_back(*it);
    }
    std::sort(values.begin(), values.end(),
              [this](auto const& a, auto const& b) {
                return key_less_(a.first, b.first);
              });
    bulk_load(values.data(), static_cast<int>(values.size()));
  }

  // Initializes with range [first, last). The range does not need to be
  // sorted. This creates a temporary copy of the data. If possible, we
  // recommend directly using bulk_load() instead.
  template <class InputIterator>
  explicit Alex(InputIterator first, InputIterator last,
                const Alloc& alloc = Alloc())
      : allocator_(alloc) {
    std::vector<V> values;
    for (auto it = first; it != last; ++it) {
      values.push_back(*it);
    }
    std::sort(values.begin(), values.end(),
              [this](auto const& a, auto const& b) {
                return key_less_(a.first, b.first);
              });
    bulk_load(values.data(), static_cast<int>(values.size()));
  }

  explicit Alex(const self_type& other)
      : params_(other.params_),
        derived_params_(other.derived_params_),
        stats_(other.stats_),
        experimental_params_(other.experimental_params_),
        istats_(other.istats_),
        key_less_(other.key_less_),
        allocator_(other.allocator_) {
    superroot_ =
        static_cast<model_node_type*>(copy_tree_recursive(other.superroot_));
    root_node_ = superroot_->children_[0];
  }

  Alex& operator=(const self_type& other) {
    if (this != &other) {
      for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
           node_it.next()) {
        delete_node(node_it.current());
      }
      delete_node(superroot_);
      params_ = other.params_;
      derived_params_ = other.derived_params_;
      experimental_params_ = other.experimental_params_;
      istats_ = other.istats_;
      stats_ = other.stats_;
      key_less_ = other.key_less_;
      allocator_ = other.allocator_;
      superroot_ =
          static_cast<model_node_type*>(copy_tree_recursive(other.superroot_));
      root_node_ = superroot_->children_[0];
    }
    return *this;
  }

  void swap(const self_type& other) {
    std::swap(params_, other.params_);
    std::swap(derived_params_, other.derived_params_);
    std::swap(experimental_params_, other.experimental_params_);
    std::swap(istats_, other.istats_);
    std::swap(stats_, other.stats_);
    std::swap(key_less_, other.key_less_);
    std::swap(allocator_, other.allocator_);
    std::swap(superroot_, other.superroot_);
    std::swap(root_node_, other.root_node_);
  }

  ////////////// Code for Disk Setting /////////////
  #define ALEXModelType 0
  #define ALEXDataType 1
  typedef struct {
    char node_type;
    char empty_byte_count;
    int children_number;
    double slope;
    double intercept;
    // Need to add other stats later
    double cost;
    int level;
    uint8_t duplication_factor_;
  } ModelHeaderOnDisk;
  #define ModelHeaderOnDiskSize sizeof(ModelHeaderOnDisk)

  typedef struct {
    char node_type;
	  char empty_byte_count_1;
    char empty_byte_count_2;
    int pre_leaf_block;
    int pre_leaf_offset;
    int next_leaf_block;
    int next_leaf_offset;
    int data_capacity;
    int num_keys;
    int bitmap_size;
    double slope;
    double intercept;
    T first_key;
    T last_key;

     // stats
    double cost;
    int level;
    uint8_t duplication_factor_;
    long long num_shifts_;
    long long num_exp_search_iterations_;
    int num_lookups_;
    int num_inserts_;
    int num_resizes_;
    int num_right_out_of_bounds_inserts_;
    int num_left_out_of_bounds_inserts_;
    double expected_avg_exp_search_iterations_;
    double expected_avg_shifts_;
    T max_key_;
    T min_key_;
    double expansion_threshold_;
    double contraction_threshold_;
    int  max_slots_;
  } DataHeaderOnDisk;
  #define DataHeaderOnDiskSize sizeof(DataHeaderOnDisk)

  typedef struct {
    T key;
    P value;
  } ItemOnDisk;
  #define ItemOnDiskSize sizeof(ItemOnDisk)

void init_model_node_from_disk(model_node_type *model_node, ModelHeaderOnDisk mh, PhyscialAddr *children) {
  model_node->num_children_ = mh.children_number;
  model_node->level_ = mh.level;
  model_node->cost_ = mh.cost;
  model_node->model_.a_ = mh.slope;
  model_node->model_.b_ = mh.intercept;
  model_node->childrenAddrs = new PhyscialAddr[model_node->num_children_];
  for (int i = 0; i < model_node->num_children_; i++) {
    model_node->childrenAddrs[i] = children[i];
  }
  model_node->duplication_factor_ = mh.duplication_factor_;
}

  void init_from_disk(data_node_type *data_node ,DataHeaderOnDisk dh, uint64_t *bitmaps_, ItemOnDisk *items) {
    data_node->data_capacity_ = dh.data_capacity;
    data_node->num_keys_ = dh.num_keys;
    data_node->bitmap_size_ = dh.bitmap_size;
    data_node->num_shifts_ = dh.num_shifts_;
    data_node->num_exp_search_iterations_ = dh.num_exp_search_iterations_;
    data_node->num_lookups_ = dh.num_lookups_;
    data_node->num_inserts_ = dh.num_inserts_;
    data_node->num_resizes_ = dh.num_resizes_;
    data_node->max_key_ = dh.max_key_;
    data_node->min_key_ = dh.min_key_;
    data_node->num_right_out_of_bounds_inserts_ = dh.num_right_out_of_bounds_inserts_;
    data_node->num_left_out_of_bounds_inserts_ = dh.num_left_out_of_bounds_inserts_;
    data_node->expected_avg_exp_search_iterations_ = dh.expected_avg_exp_search_iterations_;
    data_node->expected_avg_shifts_ = dh.expected_avg_shifts_;
    data_node->duplication_factor_ = dh.duplication_factor_;
    data_node->level_ = dh.level;
    data_node->cost_ = dh.cost;
    data_node->model_.a_ = dh.slope;
    data_node->model_.b_ = dh.intercept;
    data_node->expansion_threshold_ = dh.expansion_threshold_;
    data_node->contraction_threshold_ = dh.contraction_threshold_;
    data_node->max_slots_ = dh.max_slots_;
    #if ALEX_DATA_NODE_SEP_ARRAYS
    data_node->key_slots_ = new T[data_node->data_capacity_];
    for (int i = 0; i < data_node->data_capacity_; i++) {
      data_node->key_slots_[i] = items[i].key;
    }
    // std::copy(other.key_slots_, other.key_slots_ + other.data_capacity_,
    //           key_slots_);
    data_node->payload_slots_ = new P[data_node->data_capacity_];
    // std::copy(other.payload_slots_, other.payload_slots_ + other.data_capacity_,
    //           payload_slots_);
    for (int i = 0; i < data_node->data_capacity_; i++) {
      data_node->payload_slots_[i] = items[i].value;
    }
    #else
        data_node->data_slots_ = new (value_allocator().allocate(data_node->data_capacity_))
            V[data_node->data_capacity_];
        std::copy(other.data_slots_, other.data_slots_ + other.data_capacity_,
                  data_slots_);
    #endif
        data_node->bitmap_ = new uint64_t[data_node->bitmap_size_];
        std::copy(bitmaps_, bitmaps_ + data_node->bitmap_size_, data_node->bitmap_);
  }

  // ==== Meta Data ====
  // node type@char
  // model: a, b@double
  // children number@int
  // = alignment char =
  // ====== Data ======
  // pointer_array@2*int
  // ==================
  ModelHeaderOnDisk write_model_node(model_node_type *model_node, PhyscialAddr *r_addr, bool write_on_r_ddr = false) {
    ModelHeaderOnDisk mhd;
    mhd.node_type = ALEXModelType;
    mhd.slope = (model_node->model_).a_;
    mhd.intercept = (model_node->model_).b_;
    mhd.children_number = model_node->num_children_;
    mhd.duplication_factor_ = model_node->duplication_factor_;
    mhd.cost = model_node->cost_;
    mhd.level = model_node->level_;
    return write_model_node(mhd, model_node->childrenAddrs, r_addr, write_on_r_ddr);
  }

  //==== Meta Data =====
  // node type@char
  // model: a, b@double
  // pre_leaf@int*2
  // next_leaf@int*2
  // data_capacity@int
  // bitmap_size@int
  // num_keys_@int
  // first_key@T
  // last_key@T
  // Empty Bytes
  //==== Data ==========
  // bitmap@uint_64 array
  // Empty Bytes
  // key_payload@T+P arrayslots@T
  //====================
  DataHeaderOnDisk write_data_node(data_node_type *data_node, PhyscialAddr *r_addr, bool write_on_r_ddr, T first_key, T last_key, int next_node_block = -1, int next_node_offset = -1) {
    // int last_data_node_block = metanode.last_data_node_block;
    // int last_data_node_offset = metanode.last_data_node_offset;
    DataHeaderOnDisk dnh;
    dnh.node_type = ALEXDataType;

    dnh.slope = (data_node->model_).a_;
    dnh.intercept = (data_node->model_).b_;
    dnh.first_key = first_key;
    dnh.last_key = last_key;
    dnh.data_capacity = data_node->data_capacity_;
    dnh.num_keys = data_node->num_keys_;
    dnh.bitmap_size = data_node->bitmap_size_;
    dnh.cost = data_node->cost_;
    dnh.level = data_node->level_;
    dnh.duplication_factor_ = data_node->duplication_factor_;

    // dnh.expansion_threshold_ = data_node.expansion_threshold_;
    // dnh.contraction_threshold_ = data_node.contraction_threshold_;
    dnh.num_shifts_ = data_node->num_shifts_;
    dnh.num_exp_search_iterations_ = data_node->num_exp_search_iterations_;
    dnh.num_lookups_ = data_node->num_lookups_;
    dnh.num_inserts_ = data_node->num_inserts_;
    dnh.num_resizes_ = data_node->num_resizes_;
    dnh.max_key_ = data_node->max_key_;
    dnh.min_key_ = data_node->min_key_;
    dnh.num_right_out_of_bounds_inserts_ = data_node->num_right_out_of_bounds_inserts_;
    dnh.num_left_out_of_bounds_inserts_ = data_node->num_left_out_of_bounds_inserts_;
    dnh.expected_avg_exp_search_iterations_ = data_node->expected_avg_exp_search_iterations_;
    dnh.expected_avg_shifts_ = data_node->expected_avg_shifts_;
    dnh.contraction_threshold_ = data_node->contraction_threshold_;
    dnh.max_slots_ = data_node->max_slots_;

    dnh.expansion_threshold_ = data_node->expansion_threshold_;

    uint64_t *bitmap_arr = data_node->bitmap_;

    ItemOnDisk *all_data = new ItemOnDisk[dnh.data_capacity];
    for (int i = 0; i < dnh.data_capacity; i++) {
      all_data[i].key = data_node->key_slots_[i];
      all_data[i].value = data_node->payload_slots_[i];
    }


    if (!write_on_r_ddr) { // re-write the original node at the same space
      if (BlockSize - metanode_data_node.next_offset < DataHeaderOnDiskSize) {
        metanode_data_node.next_block += 1;
        metanode_data_node.next_offset = 0;
      }
      r_addr->block = metanode_data_node.next_block;
      r_addr->offset = metanode_data_node.next_offset;
      r_addr->flag = ALEXDataType;
    }
    r_addr->duplication_factor_ =  data_node->duplication_factor_;
//    if (r_addr->block == 5618 && r_addr->offset == 1280) {
//        std::cout << 0 << std::endl;
//    }
    dnh.pre_leaf_block = metanode_data_node.last_data_node_block;
    dnh.pre_leaf_offset = metanode_data_node.last_data_node_offset;
    dnh.next_leaf_block = next_node_block;
    dnh.next_leaf_offset = next_node_offset;
    metanode_data_node.last_data_node_block = r_addr->block;
    metanode_data_node.last_data_node_offset = r_addr->offset;
    if (dnh.pre_leaf_block != -1) {
      char data[BlockSize];
      data_node_sm->_read_block(data, dnh.pre_leaf_block);
      DataHeaderOnDisk *_dhd = (DataHeaderOnDisk *) (data + dnh.pre_leaf_offset);
      _dhd->next_leaf_block = r_addr->block;
      _dhd->next_leaf_offset = r_addr->offset;
      data_node_sm->write_with_size(dnh.pre_leaf_block, data, BlockSize);
    } else {
      // Note we put them into the metanode
      metanode.left_most_block = r_addr->block;
      metanode.left_most_offset = r_addr->offset;
      metanode.dup_left = data_node->duplication_factor_;
      istats_.min_key_in_left = first_key;
    }
    if (dnh.next_leaf_block > 0) {
      char data[BlockSize];
      data_node_sm->_read_block(data, dnh.next_leaf_block);
      DataHeaderOnDisk *_dhd = (DataHeaderOnDisk *) (data + dnh.next_leaf_offset);
      _dhd->pre_leaf_block = r_addr->block;
      _dhd->pre_leaf_offset = r_addr->offset;
      data_node_sm->write_with_size(dnh.next_leaf_block, data, BlockSize);
    } else if (dnh.next_leaf_block == -1){
      metanode.right_most_block = r_addr->block;
      metanode.right_most_offset = r_addr->offset;
      metanode.dup_right = data_node->duplication_factor_;
      istats_.max_key_in_right = last_key;
    }

    // simu write process
    int next_offset_simu = r_addr->offset;
    next_offset_simu += DataHeaderOnDiskSize;
    dnh.empty_byte_count_1 = (BlockSize - next_offset_simu) % sizeof(uint64_t);
    next_offset_simu += dnh.empty_byte_count_1;
    next_offset_simu += dnh.bitmap_size * sizeof(uint64_t);
    next_offset_simu = next_offset_simu % BlockSize;

    dnh.empty_byte_count_2 = (BlockSize - next_offset_simu) % ItemOnDiskSize;

    // do write header
    long start_offset = r_addr->block * BlockSize + r_addr->offset;
    data_node_sm->write_arbitrary(start_offset, &dnh, DataHeaderOnDiskSize);
    start_offset += DataHeaderOnDiskSize;

    if (dnh.empty_byte_count_1 > 0) {
      char _data[dnh.empty_byte_count_1];
      data_node_sm->write_arbitrary(start_offset, _data, dnh.empty_byte_count_1);
      start_offset += dnh.empty_byte_count_1;
    }

    data_node_sm->write_arbitrary(start_offset, bitmap_arr, dnh.bitmap_size * sizeof(uint64_t));
    start_offset += dnh.bitmap_size * sizeof(uint64_t);

    if (dnh.empty_byte_count_2 > 0) {
      char _data[dnh.empty_byte_count_2];
      data_node_sm->write_arbitrary(start_offset, _data, dnh.empty_byte_count_2);
      start_offset += dnh.empty_byte_count_2;
    }

    data_node_sm->write_arbitrary(start_offset, all_data, dnh.data_capacity * ItemOnDiskSize);
    start_offset += dnh.data_capacity * ItemOnDiskSize;

    if (!write_on_r_ddr) {
      metanode_data_node.next_block = start_offset / BlockSize;
      metanode_data_node.next_offset = start_offset % BlockSize;
    }
//      if (r_addr->block == 5618 && r_addr->offset >= 1280) {
//          int b_block = 46024112 / BlockSize;
//          int b_offset = 46024112 % BlockSize;
//          char data[BlockSize];
//          data_node_sm->_read_block(data, b_block);
//          uint64_t x = *(uint64_t *)((data + b_offset));
//          std::cout << x << std::endl;
//      }

    delete []all_data;
    return dnh;
  }

  // can be mergerd with write_model_node above
  ModelHeaderOnDisk write_model_node(ModelHeaderOnDisk mhd, PhyscialAddr *children, PhyscialAddr *r_addr, bool write_on_r_ddr = false) {
    if (!write_on_r_ddr) { // re-write the original node at the same space
      if (BlockSize - metanode.next_offset < DataHeaderOnDiskSize) {
        char _data[BlockSize - metanode.next_offset];
        sm->write_arbitrary(metanode.next_block * BlockSize + metanode.next_offset, _data, BlockSize - metanode.next_offset);
        metanode.next_block += 1;
        metanode.next_offset = 0;
      }
      r_addr->block = metanode.next_block;
      r_addr->offset = metanode.next_offset;
      r_addr->flag = ALEXModelType;
    }
    r_addr->duplication_factor_ = mhd.duplication_factor_;
    // if (BlockSize - metanode.next_offset < ModelHeaderOnDiskSize) {
    //   char _data[BlockSize - metanode.next_offset];
    //   sm->write_arbitrary(metanode.next_block * BlockSize + metanode.next_offset, _data, BlockSize - metanode.next_offset);
    //   metanode.next_block += 1;
    //   metanode.next_offset = 0;
    // }
    // r_addr->block = metanode.next_block;
    // r_addr->offset = metanode.next_offset;
    // r_addr->flag = ALEXModelType;
    mhd.empty_byte_count = (BlockSize - r_addr->offset - ModelHeaderOnDiskSize) % PhyscialAddrSize;

    // write header
    long start_offset = r_addr->block * BlockSize + r_addr->offset;
    sm->write_arbitrary(start_offset, &mhd, ModelHeaderOnDiskSize);
    start_offset += ModelHeaderOnDiskSize;

    // write empty bytes
    if (mhd.empty_byte_count > 0) {
      char _data[mhd.empty_byte_count];
      sm->write_arbitrary(start_offset, _data, mhd.empty_byte_count);
      start_offset += mhd.empty_byte_count;
    }

    // write children pointers
    // old we can add `test' attribute to make sure it is right!
    // sm->write_arbitrary(start_offset, model_node->childrenAddrs, model_node->num_children_ * PhyscialAddrSize);
    // start_offset += model_node->num_children_ * PhyscialAddrSize;
    //
    int remained = BlockSize - (start_offset % BlockSize);
    int c1 = remained / PhyscialAddrSize;

    if (c1 >= mhd.children_number ) {
      sm->write_arbitrary(start_offset, children, mhd.children_number * PhyscialAddrSize);
      start_offset += mhd.children_number * PhyscialAddrSize;
    } else {
      sm->write_arbitrary(start_offset, children, c1 * PhyscialAddrSize);
      start_offset += remained;
      char data[BlockSize];
      for (int i = c1; i <  mhd.children_number; ) {
        int _c = mhd.children_number - i;
        if (_c > MaxModeCount) _c = MaxModeCount;
        if (_c == MaxModeCount) {
          memcpy(data, children + i, _c * PhyscialAddrSize);
          sm->write_arbitrary(start_offset, data, BlockSize);
          start_offset += BlockSize;
        } else {
          sm->write_arbitrary(start_offset, children + i, _c * PhyscialAddrSize);
          start_offset += _c * PhyscialAddrSize;
        }
        i += _c;
      }
    }

    // update the metanode info
    if (!write_on_r_ddr) {
      metanode.next_block = start_offset / BlockSize;
      metanode.next_offset = start_offset % BlockSize;
    }
    return mhd;
  }

  T get_key(DataHeaderOnDisk *data_node, int pred_pos, int _block_, int _offset_, int *last_block, char *data, int *r_block_count) {
      long key_offset = _block_ * BlockSize + _offset_ + DataHeaderOnDiskSize
                      + data_node->empty_byte_count_1 + data_node->bitmap_size * sizeof(uint64_t)
                      + data_node->empty_byte_count_2 + pred_pos * ItemOnDiskSize;
      int key_block = key_offset / BlockSize;
      int key_offset_in_block = key_offset % BlockSize;
      if (key_block != *last_block) {
          //*b_data = BC->get_block(key_block);
          // BC->get_block(key_block, data);
          data_node_sm->get_block(key_block, data);
          *last_block = key_block;
          *r_block_count += 1;
      }
      T pred_key = ((ItemOnDisk *)(data + key_offset_in_block))->key;
      return pred_key;
  }

  P get_value(DataHeaderOnDisk *data_node, int pred_pos, int _block_, int _offset_, char *data) {
      long key_offset = _block_ * BlockSize + _offset_ + DataHeaderOnDiskSize
                      + data_node->empty_byte_count_1 + data_node->bitmap_size * sizeof(uint64_t)
                      + data_node->empty_byte_count_2 + pred_pos * ItemOnDiskSize;
      int key_offset_in_block = key_offset % BlockSize;
      P value = ((ItemOnDisk *)(data + key_offset_in_block))->value;
      return value;
  }

  bool _key_greater(T key1, T key2) {
    return key1 > key2;
  }

  bool _key_lessequal(T key1, T key2) {
    return key1 <= key2;
  }

  bool _key_equal(T key1, T key2) {
    return key1 == key2;
  }

  P* exp_search(DataHeaderOnDisk *data_node, int pred_pos, T key, int _block_, int _offset_, char *data, int *POS, int *r_block_count, int *last_block, bool is_insert = false) {
      int bound = 1;
      int l, r;
//      int last_block = _block_;
      if (_key_greater(get_key(data_node, pred_pos, _block_, _offset_, last_block, data, r_block_count), key)) {
          int size = pred_pos;
          while (bound < size && _key_greater(get_key(data_node, pred_pos - bound, _block_, _offset_, last_block, data, r_block_count), key)) {
              bound *= 2;
              data_node->num_exp_search_iterations_ ++;
          }
          l = pred_pos - std::min<int>(bound, size);
          r = pred_pos - bound / 2;
      } else {
          int size = data_node->data_capacity - pred_pos;
          while (bound < size && _key_lessequal(get_key(data_node, pred_pos + bound, _block_, _offset_, last_block, data, r_block_count), key)) {
              bound *= 2;
              data_node->num_exp_search_iterations_ ++;
          }
          l = pred_pos + bound / 2;
          r = pred_pos + std::min<int>(bound, size);
      }
      // binary search
      while (l < r) {
          int mid = l + (r -l) / 2;
          if (_key_lessequal(get_key(data_node, mid, _block_, _offset_, last_block, data, r_block_count), key)) {
              l = mid + 1;
          } else {
              r = mid;
          }
      }
      if (is_insert) {
        *POS = l;
        return nullptr;
      }
      l -= 1;
      if (l < 0 || !_key_equal(get_key(data_node, l, _block_, _offset_, last_block, data, r_block_count), key)) {
        if (POS != nullptr)
          *POS = *POS = l > 0 ? l : 0;
        return nullptr;
      }
      P _v = get_value(data_node, l, _block_, _offset_, data);

      if (POS != nullptr) {
          *POS = l;
          return nullptr;
      }
      P *v = (P *)malloc(sizeof(P));
      *v = _v;
      return v;
  }

  P *get_payload_on_disk(T key, int block, int offset, char *data, int *r_block_count, int *POS = nullptr) {
    // DataHeaderOnDisk
    DataHeaderOnDisk _data_node;
    memcpy(&_data_node, data + offset, DataHeaderOnDiskSize);
    double pred_pos_d = _data_node.slope * static_cast<double>(key) +  _data_node.intercept;
    int pred_pos_i = static_cast<int>(pred_pos_d);
    pred_pos_i = std::min<int>(std::max<int>(pred_pos_i, 0), _data_node.data_capacity - 1);
    int last_block = block;
    return exp_search(&_data_node, pred_pos_i, key, block, offset, data, POS, r_block_count, &last_block);
  }

  bool get_leaf_disk_hybrid(T key, P *v, int *r_block_count, int *r_level) {
      char *data = new char[BlockSize];
      AlexNode<T, P>* cur = root_node_;
      int leaf_node_block;
      int leaf_node_offset;
      bool go_while = true;
      if (metanode.is_leaf == ALEXDataType) {
        leaf_node_block = metanode.root_block_id;
        leaf_node_offset = metanode.root_offset;
        P *_v = get_payload_on_disk(key, leaf_node_block, leaf_node_offset, data, r_block_count);
        if (_v == nullptr) {
            return false;
        }
        *v = *_v;
        free(_v);
        return true;
      }
      while (true) {
          auto node = static_cast<model_node_type*>(cur);
          double bucketID_prediction = node->model_.predict_double(key);
          int bucketID = static_cast<int>(bucketID_prediction);
          bucketID = std::min<int>(std::max<int>(bucketID, 0), node->num_children_ - 1);
          if (node->is_leaf[bucketID] == 0) {
              cur = node->children_[bucketID];
              continue;
          }
          leaf_node_block = node->childrenAddrs[bucketID].block;
          leaf_node_offset = node->childrenAddrs[bucketID].offset;
//          std::cout << leaf_node_block << "," << leaf_node_offset << std::endl;
          bool found = check_leaf_node(false, key, v, leaf_node_block, leaf_node_offset, data, bucketID_prediction, r_block_count, r_level);
          delete []data;
          return found;
      }

  }

  bool check_leaf_node(bool same, T key, P *v, int block,int offset, char *data, double pred_pos_d, int *r_block_count, int *r_level) {
//      int block = child_addr.block;
//      int offset = child_addr.offset;
//      if (!same) {
          data_node_sm->get_block(block, data);
          *r_block_count += 1;
//      }
      *r_level += 1;
      DataHeaderOnDisk *dh = (DataHeaderOnDisk *) (data + offset);
      int pred_pos_i_rounded = static_cast<int>(pred_pos_d + 0.5);
      double tolerance = 10 * std::numeric_limits<double>::epsilon() * pred_pos_d;
      if (std::abs(pred_pos_d - pred_pos_i_rounded) <= tolerance) {
          if (pred_pos_i_rounded <= pred_pos_d) {
              if (dh->pre_leaf_block != -1) {
                  char pre_node_data[BlockSize];
                  data_node_sm->get_block(dh->pre_leaf_block, pre_node_data);
                  *r_block_count += 1;
                  DataHeaderOnDisk *pre_dh = (DataHeaderOnDisk *) (pre_node_data + dh->pre_leaf_offset);
                  if (pre_dh->last_key >= key) {
                      P *_v = get_payload_on_disk(key, dh->pre_leaf_block, dh->pre_leaf_offset, pre_node_data, r_block_count);
                      if (_v == nullptr) {
                          return false;
                      }
                      *v = *_v;
                      free(_v);
                      //delete []data;
                      return true;
                  }
              }
          } else {
              if (dh->next_leaf_block != -1) {
                  char next_node_data[BlockSize];
                  data_node_sm->get_block(dh->next_leaf_block, next_node_data);
                  DataHeaderOnDisk *next_dh = (DataHeaderOnDisk *) (next_node_data + dh->next_leaf_offset);
                  if (next_dh->first_key <= key) {
                      P *_v = get_payload_on_disk(key, dh->next_leaf_block, dh->next_leaf_offset, next_node_data, r_block_count);
                      if (_v == nullptr) {
                          return false;
                      }
                      *v = *_v;
                      free(_v);
                      //delete []data;
                      return true;
                  }
              }
          }
      }
      P *_v = get_payload_on_disk(key, block, offset, data, r_block_count);
      if (_v == nullptr) {
          return false;
      }
      *v = *_v;
      free(_v);
      return true;
  }

  bool get_leaf_disk(T key, P *v, int *r_block_count, int *r_level, int *r_inner_block) {
    if (hybrid_mode == ALL_DISK) return get_leaf_disk_all(key, v, r_block_count, r_level, r_inner_block);
    else if (hybrid_mode == LEAF_DISK) return get_leaf_disk_hybrid(key, v, r_block_count, r_level);
    else {std::invalid_argument("unclear model... in get_leaf_disk()");}
  }

  bool get_leaf_disk_all(T key, P *v, int *r_block_count, int *r_level, int *r_inner_block) {
    char *data = new char[BlockSize];
    int block = metanode.root_block_id;
    int offset = metanode.root_offset;
    sm->get_block(block, data);
    *r_block_count += 1;
    *r_inner_block += 1;
    *r_level += 1;

    ModelHeaderOnDisk *mhd;
    mhd = (ModelHeaderOnDisk *)(data + offset);
    PhyscialAddr child_addr;
    while (true) {
      double pred_pos_d = mhd->slope * static_cast<double>(key) + mhd->intercept;
      int pred_pos_i = static_cast<int>(pred_pos_d);
      pred_pos_i = std::min<int>(std::max<int>(pred_pos_i, 0), mhd->children_number - 1);

      long _offset = (block * BlockSize + offset + ModelHeaderOnDiskSize + mhd->empty_byte_count) % BlockSize;
      int c1 = (BlockSize - _offset)/PhyscialAddrSize;
      int block_in_inode  = 0;
      int offset_in_inode = 0;
      if (_offset == 0) c1 = 0;
      if (pred_pos_i < c1) {
        block_in_inode = block;
        offset_in_inode = _offset + pred_pos_i * PhyscialAddrSize;
      } else {
        block_in_inode = block;
        block_in_inode += 1;
        block_in_inode += (pred_pos_i - c1) / MaxModeCount;
        offset_in_inode = ((pred_pos_i - c1) % MaxModeCount) * PhyscialAddrSize;
      }

      if (block_in_inode != block) {
        sm->get_block(block_in_inode, data);
        *r_block_count += 1;
        *r_inner_block += 1;
        block = block_in_inode;
      }
      //PhyscialAddr child_addr = *((PhyscialAddr *) (data + offset_in_inode));
      memcpy(&child_addr, data + offset_in_inode, PhyscialAddrSize);
      if (child_addr.flag == ALEXModelType) {
          if (child_addr.block != block) {
            sm->get_block(child_addr.block, data);
            *r_block_count += 1;
            *r_inner_block += 1;
            block = child_addr.block;
          }
          offset = child_addr.offset;
          mhd = (ModelHeaderOnDisk *)(data + child_addr.offset);
          *r_level += 1;
          continue;
      }
//      std::cout << "child_addr:"<< child_addr.block << ";"<<child_addr.offset<<std::endl;
      bool found = check_leaf_node(block== child_addr.block, key, v, child_addr.block, child_addr.offset, data, pred_pos_d, r_block_count, r_level);
      delete []data;
      return found;
    }

  }

    class IndexIterator {
        uint64_t *_current_bitmap_;
        char bitmap_data[BlockSize];
        int current_access_block;
        int _current_block; // point to block_data
        StorageManager *ii_sm;
        T min_key;
        int cur_bitmap_idx_ = 0;
        uint64_t cur_bitmap_data_ = 0;
        // int data_meta_size = sizeof(_meta_data_node);
        int _current_bitmap_start_idx = 0;
        int _current_read_len_bitmap = 0;
        int _current_bitmap_end_idx = 0;

        public:
            char block_data[BlockSize];
            DataHeaderOnDisk current_leaf_node_header;
            int current_start_block;
            int current_start_offset;
            int next_point;
            IndexIterator() = default;
            IndexIterator(T _min_key, StorageManager *_sm) {
                // block_data = new char;
                min_key = _min_key;
                ii_sm = _sm;
            }
            void init(int *r_block_count) {
                obtain_bitmap(r_block_count, next_point >> 6);
                obtain_current_block();
            }

            void obtain_bitmap(int *r_block_count, int start_index) {
                _current_bitmap_start_idx = start_index;
                _current_read_len_bitmap = (current_leaf_node_header.bitmap_size - _current_bitmap_start_idx) * sizeof(uint64_t);
                if (_current_read_len_bitmap > BlockSize) {
                  _current_read_len_bitmap = BlockSize;
                }
                _current_bitmap_end_idx = _current_bitmap_start_idx + _current_read_len_bitmap / sizeof(uint64_t);

                long offset = current_start_block * BlockSize + current_start_offset + DataHeaderOnDiskSize
                + current_leaf_node_header.empty_byte_count_1 + sizeof(uint64_t)*_current_bitmap_start_idx;

                // ii_sm->get_arbitrary(offset, bitmap_data, len);
                ii_sm->read_arbitrary(bitmap_data, offset, _current_read_len_bitmap);
                *r_block_count += 1;
                // printf("1-1\n");
                _current_bitmap_ = (uint64_t *) bitmap_data;
                cur_bitmap_idx_ = start_index;
                cur_bitmap_data_ = _current_bitmap_[cur_bitmap_idx_ - _current_bitmap_start_idx];
            }

            void obtain_current_block() {
                long key_offset = current_start_block * BlockSize + current_start_offset + DataHeaderOnDiskSize
                            + current_leaf_node_header.empty_byte_count_1 + current_leaf_node_header.bitmap_size * sizeof(uint64_t)
                            + current_leaf_node_header.empty_byte_count_2 + next_point * ItemOnDiskSize;
                _current_block = key_offset / BlockSize;
            }

            uint64_t extract_rightmost_one(uint64_t value) {
                return value & -static_cast<int64_t>(value);
            }
            uint64_t remove_rightmost_one(uint64_t value) {
                return value & (value - 1);
            }

            int get_offset(int word_id, uint64_t bit) {
              return (word_id << 6) + count_ones(bit - 1);
            }

            int count_ones(uint64_t value) {
               int count  = 0;
               while (value)
               {
                 count += value & 1;
                 value >>= 1;
               }
               return static_cast<int>(count);
//            int x = __builtin_popcount(value);
//            return __builtin_popcount(value);
            }

            bool next(int *r_block_count, bool is_first = false) {
		//std::cout<<"enter next:"<<std::endl;
                while (cur_bitmap_data_ == 0) {
                    cur_bitmap_idx_ ++;
                    if (cur_bitmap_idx_ >= current_leaf_node_header.bitmap_size) {
                        // obtain the next node
                        current_start_block = current_leaf_node_header.next_leaf_block;
                        current_start_offset = current_leaf_node_header.next_leaf_offset;

//                        if (current_start_block == 5618 && current_start_offset == 1280) {
//                            std::cout << 0 << std::endl;
////                            if (r_addr->block == 5618 && r_addr->offset >= 1280) {
//                                int b_block = 46024112 / BlockSize;
//                                int b_offset = 46024112 % BlockSize;
//                                char data[BlockSize];
//                                ii_sm->_read_block(data, b_block);
//                                uint64_t x = *(uint64_t *)((data + b_offset));
//                                std::cout << x << std::endl;
////                            }
//                        }


                        //std::cout << "read block:" << current_start_block << std::endl;
                        if (_current_block != current_start_block) {
                            ii_sm->get_block(current_start_block, block_data);
                            *r_block_count += 1;
                            _current_block = current_start_block;
                        }

                        memcpy(&current_leaf_node_header, block_data + current_start_offset, DataHeaderOnDiskSize);
                        int len = current_leaf_node_header.bitmap_size * sizeof(uint64_t);

                        long offset = current_start_block * BlockSize + current_start_offset + DataHeaderOnDiskSize
                        + current_leaf_node_header.empty_byte_count_1;
                        int b_block = offset / BlockSize;
                        int b_offset = offset % BlockSize;
                        next_point = 0;
                        if (b_block == current_start_block && (b_offset + len) <= BlockSize) {
//                            uint64_t x = *(uint64_t *)((block_data + b_offset));
                            memcpy(bitmap_data, block_data+ b_offset, len);
                            cur_bitmap_idx_ = next_point >> 6;
                            _current_bitmap_ = (uint64_t *) bitmap_data;
                            cur_bitmap_data_ = _current_bitmap_[cur_bitmap_idx_];
                            _current_bitmap_start_idx = next_point >> 6;
                        } else obtain_bitmap(r_block_count, next_point >> 6);
                        _current_block = current_start_block;
                    } else if (cur_bitmap_idx_ >= _current_bitmap_end_idx){
                       obtain_bitmap(r_block_count, cur_bitmap_idx_);
                    }
                    cur_bitmap_data_ = _current_bitmap_[cur_bitmap_idx_ - _current_bitmap_start_idx];
                }
                uint64_t bit = extract_rightmost_one(cur_bitmap_data_);
                int temp_next_point = get_offset(cur_bitmap_idx_, bit);
                if (is_first) {
                   while (temp_next_point < next_point && cur_bitmap_data_ != 0) {
                        cur_bitmap_data_ = remove_rightmost_one(cur_bitmap_data_);
                        bit = extract_rightmost_one(cur_bitmap_data_);
                        temp_next_point = get_offset(cur_bitmap_idx_, bit);
                    }
                    next_point = temp_next_point;
                    cur_bitmap_data_ = remove_rightmost_one(cur_bitmap_data_);
                } else {
                    next_point = temp_next_point;//get_offset(cur_bitmap_idx_, bit);
                    cur_bitmap_data_ = remove_rightmost_one(cur_bitmap_data_);
                }
		return true;
            }

            T get_key(int *r_block_count) {
                long key_offset = current_start_block * BlockSize + current_start_offset + DataHeaderOnDiskSize
                            + current_leaf_node_header.empty_byte_count_1 + current_leaf_node_header.bitmap_size * sizeof(uint64_t)
                            + current_leaf_node_header.empty_byte_count_2 + next_point * (sizeof(T) + sizeof(P));
                int key_block = key_offset / BlockSize;
                int key_offset_in_block = key_offset % BlockSize;
//                std::cout<<"block_id:"<<key_block<<std::endl;
                if (key_block != _current_block) {
                //*b_data = BC->get_block(key_block);
                    ii_sm->get_block(key_block, block_data);
                    // printf("1-3\n");
                    _current_block = key_block;
                    *r_block_count += 1;
                }
                T pred_key = ((ItemOnDisk *)(block_data + key_offset_in_block))->key;
                return pred_key;
            }
  };

  IndexIterator get_it(T key, int *r_block_count) {
        IndexIterator ii(key, data_node_sm);
	get_leaf_node_scan(key, &ii.current_leaf_node_header, ii.block_data, &ii.current_start_block, &ii.current_start_offset, &ii.next_point, r_block_count);
        ii.init(r_block_count);
        return ii;
  }

  void get_leaf_node_scan_common(T key, int block, int offset, char *block_data,
                                 char *block_data2, int *r_block_count, double pred_pos_d,
                                 int *l_block, int *l_offset, int *POS, DataHeaderOnDisk *LEAF_NODE) {
//      int block = child_addr.block;
//      int offset = child_addr.offset;
      data_node_sm->get_block(block, block_data2);
      *r_block_count += 1;
      // data_node = (_meta_data_node*)(block_data + _block_if.offset);
      DataHeaderOnDisk *data_node = (DataHeaderOnDisk *)( block_data2 + offset);
      // memcpy(&data_node, block_data + _block_if.offset, data_meta_size);
      int pred_pos_i_rounded = static_cast<int>(pred_pos_d + 0.5);
      double tolerance = 10 * std::numeric_limits<double>::epsilon() * pred_pos_d;
      if (std::abs(pred_pos_d - pred_pos_i_rounded) <= tolerance) {
          if (pred_pos_i_rounded <= pred_pos_d) {
              if (data_node->pre_leaf_block != -1) {
                  // char pre_block_data[BLOCK_SIZE];
                  data_node_sm->get_block(data_node->pre_leaf_block, block_data);
                  *r_block_count += 1;
                  DataHeaderOnDisk *pre_data_node = (DataHeaderOnDisk *) (block_data + data_node->pre_leaf_offset);

                  if (pre_data_node->last_key >= key) {
                      memcpy(LEAF_NODE, pre_data_node, DataHeaderOnDiskSize);
                      get_payload_on_disk(key, data_node->pre_leaf_block,data_node->pre_leaf_offset, block_data, r_block_count, POS);
                      *l_block = data_node->pre_leaf_block;
                      *l_offset = data_node->pre_leaf_offset;
                      memcpy(block_data2, block_data, BlockSize);
                      return;
                  }
              }
          } else {
              if (data_node->next_leaf_block != -1) {
                  data_node_sm->get_block(data_node->next_leaf_block, block_data);
                  *r_block_count += 1;
                  DataHeaderOnDisk *next_data_node = (DataHeaderOnDisk *)(block_data + data_node->next_leaf_offset);

                  // memcpy(&next_data_node, block_data + data_node.next_leaf_offset, data_meta_size);
                  if (next_data_node->first_key <= key) {
                      memcpy(LEAF_NODE, next_data_node, DataHeaderOnDiskSize);
                      get_payload_on_disk(key, data_node->next_leaf_block, data_node->next_leaf_offset, block_data, r_block_count, POS);
                      *l_block = data_node->next_leaf_block;
                      *l_offset = data_node->next_leaf_offset;
                      memcpy(block_data2, block_data, BlockSize);
                      return;
                  }
              }
          }
      }
      memcpy(LEAF_NODE, data_node, DataHeaderOnDiskSize);
      get_payload_on_disk(key, block, offset, block_data2, r_block_count ,POS);
      *l_block = block;
//      std::cout << block << std::endl;
      *l_offset = offset;
      return;
  }

    // we can merge it with get_leaf_disk later...
    void get_leaf_node_scan(T key, DataHeaderOnDisk *LEAF_NODE, char *block_data2, int *l_block, int *l_offset, int *POS, int *r_block_count) {
    // void init_it() {
        char block_data[BlockSize];
        int block = metanode.root_block_id;
        int offset = metanode.root_offset;
        sm->get_block(block, block_data2);
        *r_block_count += 1;
        // if (block_data2[offset] == ALEXDataType) {
        //     memcpy(LEAF_NODE, block_data2 + offset, DataHeaderOnDiskSize);
        //     get_payload_on_disk(key, block, offset, block_data2, r_block_count, POS);
        //     *l_block = block;
        //     *l_offset = offset;
        //     return;
        // }
        ModelHeaderOnDisk *cur;
        cur = (ModelHeaderOnDisk *)(block_data2 + offset);
        while (true) {
            double pred_pos_d = cur->slope * static_cast<double>(key) + cur->intercept;
            int pred_pos_i = static_cast<int>(pred_pos_d);
            pred_pos_i = std::min<int>(std::max<int>(pred_pos_i, 0), cur->children_number - 1);
            long _offset = (block * BlockSize + offset + ModelHeaderOnDiskSize + cur->empty_byte_count) % BlockSize;
            int c1 = (BlockSize - _offset)/PhyscialAddrSize;
            int block_in_inode  = 0;
            int offset_in_inode = 0;
            if (_offset == 0) c1 = 0;
            if (pred_pos_i < c1) {
              block_in_inode = block;
              offset_in_inode = _offset + pred_pos_i * PhyscialAddrSize;
            } else {
              block_in_inode = block;
              block_in_inode += 1;
              block_in_inode += (pred_pos_i - c1) / MaxModeCount;
              offset_in_inode = ((pred_pos_i - c1) % MaxModeCount) * PhyscialAddrSize;
            }

            if (block_in_inode != block) {
                sm->get_block(block_in_inode, block_data2);
                *r_block_count += 1;
                block = block_in_inode;
            }

            PhyscialAddr child_addr = *(PhyscialAddr *)(block_data2 + offset_in_inode);
            if (child_addr.flag == ALEXModelType) {
                if (child_addr.block != block) {
                  sm->get_block(child_addr.block, block_data2);
                  *r_block_count += 1;
                  block = child_addr.block;
                }
                offset = child_addr.offset;
                cur = (ModelHeaderOnDisk *)(block_data2 + offset);
                continue;
            }
            return get_leaf_node_scan_common(key, child_addr.block, child_addr.offset, block_data, block_data2, r_block_count,
                                             pred_pos_d, l_block, l_offset, POS, LEAF_NODE);
//            void get_leaf_node_scan_common(int key, PhyscialAddr child_addr, char *block_data,
//                                           char *block_data2, int *r_block_count, double pred_pos_d,
//                                           int *l_block, int *l_offset, int *POS, DataHeaderOnDisk *LEAF_NODE)
        }
  }

    IndexIterator get_ii_hybrid(T key, int *r_block_count) {

      IndexIterator ii(key, data_node_sm);
      char *data = new char[BlockSize];
      AlexNode<T, P>* cur = root_node_;

      int leaf_node_block;
      int leaf_node_offset;

      if (metanode.is_leaf == ALEXDataType) {
          leaf_node_block = metanode.root_block_id;
          leaf_node_offset = metanode.root_offset;
//          get_leaf_node_scan_common(key, {leaf_node_block, leaf_node_offset}, data, ii.block_data, r_block_count, );
          data_node_sm->get_block(leaf_node_block, ii.block_data);
          DataHeaderOnDisk *data_node = (DataHeaderOnDisk *)( ii.block_data + leaf_node_offset);
          memcpy(&ii.current_leaf_node_header, data_node, DataHeaderOnDiskSize);
          get_payload_on_disk(key, leaf_node_block, leaf_node_offset, ii.block_data, r_block_count ,&ii.next_point);
          ii.current_start_block = leaf_node_block;
          ii.current_start_offset = leaf_node_offset;
          *r_block_count += 1;
          ii.init(r_block_count);
          return ii;
      }
      while (true) {
          auto node = static_cast<model_node_type*>(cur);
          double bucketID_prediction = node->model_.predict_double(key);
          int bucketID = static_cast<int>(bucketID_prediction);
          bucketID = std::min<int>(std::max<int>(bucketID, 0), node->num_children_ - 1);
          if (node->is_leaf[bucketID] == 0) {
              cur = node->children_[bucketID];
              continue;
          }
          leaf_node_block = node->childrenAddrs[bucketID].block;
          leaf_node_offset = node->childrenAddrs[bucketID].offset;
          get_leaf_node_scan_common(key, leaf_node_block, leaf_node_offset, data, ii.block_data, r_block_count,
                                           bucketID_prediction, &ii.current_start_block, &ii.current_start_offset, &ii.next_point, &ii.current_leaf_node_header);
          delete []data;
          ii.init(r_block_count);
          return ii;
      }

  }

  void scan(T *keys, T key, int *r_block_count, int len) {
      IndexIterator ii;
      if (hybrid_mode == ALL_DISK)
        ii = get_it(key, r_block_count);
      else if (hybrid_mode == LEAF_DISK)
        ii = get_ii_hybrid(key, r_block_count);
      else
          std::invalid_argument("not support the mode... in scan");
      bool is_first = true;
      for (int i = 0; i < len;) {
          ii.next(r_block_count, is_first);
          is_first = false;
          T _key = ii.get_key(r_block_count);

          if (_key >= key) {
              keys[i] = _key;
              i++;
          } else {
//		    std::cout << "small" << std::endl;
	    }
      }
      return;
  }


 private:
  // Deep copy of tree starting at given node
  AlexNode<T, P>* copy_tree_recursive(const AlexNode<T, P>* node) {
    if (!node) return nullptr;
    if (node->is_leaf_) {
      return new (data_node_allocator().allocate(1))
          data_node_type(*static_cast<const data_node_type*>(node));
    } else {
      auto node_copy = new (model_node_allocator().allocate(1))
          model_node_type(*static_cast<const model_node_type*>(node));
      int cur = 0;
      while (cur < node_copy->num_children_) {
        AlexNode<T, P>* child_node = node_copy->children_[cur];
        AlexNode<T, P>* child_node_copy = copy_tree_recursive(child_node);
        int repeats = 1 << child_node_copy->duplication_factor_;
        for (int i = cur; i < cur + repeats; i++) {
          node_copy->children_[i] = child_node_copy;
        }
        cur += repeats;
      }
      return node_copy;
    }
  }

 public:
  // When bulk loading, Alex can use provided knowledge of the expected fraction
  // of operations that will be inserts
  // For simplicity, operations are either point lookups ("reads") or inserts
  // ("writes)
  // i.e., 0 means we expect a read-only workload, 1 means write-only
  // This is only useful if you set it before bulk loading
  void set_expected_insert_frac(double expected_insert_frac) {
    assert(expected_insert_frac >= 0 && expected_insert_frac <= 1);
    params_.expected_insert_frac = expected_insert_frac;
  }

  // Maximum node size, in bytes.
  // Higher values result in better average throughput, but worse tail/max
  // insert latency.
  void set_max_node_size(int max_node_size) {
    assert(max_node_size >= sizeof(V));
    params_.max_node_size = max_node_size;
    derived_params_.max_fanout = params_.max_node_size / sizeof(void*);
    derived_params_.max_data_node_slots = params_.max_node_size / sizeof(V);
  }

  // Bulk load faster by using sampling to train models.
  // This is only useful if you set it before bulk loading.
  void set_approximate_model_computation(bool approximate_model_computation) {
    params_.approximate_model_computation = approximate_model_computation;
  }

  // Bulk load faster by using sampling to compute cost.
  // This is only useful if you set it before bulk loading.
  void set_approximate_cost_computation(bool approximate_cost_computation) {
    params_.approximate_cost_computation = approximate_cost_computation;
  }

  /*** General helpers ***/

 public:
// Return the data node that contains the key (if it exists).
// Also optionally return the traversal path to the data node.
// traversal_path should be empty when calling this function.
// The returned traversal path begins with superroot and ends with the data
// node's parent.
#if ALEX_SAFE_LOOKUP
  forceinline data_node_type* get_leaf(
      T key, std::vector<TraversalNode>* traversal_path = nullptr) const {
    if (traversal_path) {
      traversal_path->push_back({superroot_, 0});
    }
    AlexNode<T, P>* cur = root_node_;
    if (cur->is_leaf_) {
      return static_cast<data_node_type*>(cur);
    }

    while (true) {
      auto node = static_cast<model_node_type*>(cur);
      double bucketID_prediction = node->model_.predict_double(key);
      int bucketID = static_cast<int>(bucketID_prediction);
      bucketID =
          std::min<int>(std::max<int>(bucketID, 0), node->num_children_ - 1);
      if (traversal_path) {
        traversal_path->push_back({node, bucketID});
      }
      cur = node->children_[bucketID];
      if (cur->is_leaf_) {
        stats_.num_node_lookups += cur->level_;
        auto leaf = static_cast<data_node_type*>(cur);
        // Doesn't really matter if rounding is incorrect, we just want it to be
        // fast.
        // So we don't need to use std::round or std::lround.
        int bucketID_prediction_rounded =
            static_cast<int>(bucketID_prediction + 0.5);
        double tolerance =
            10 * std::numeric_limits<double>::epsilon() * bucketID_prediction;
        // https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
        if (std::abs(bucketID_prediction - bucketID_prediction_rounded) <=
            tolerance) {
          if (bucketID_prediction_rounded <= bucketID_prediction) {
            if (leaf->prev_leaf_ && leaf->prev_leaf_->last_key() >= key) {
              if (traversal_path) {
                // Correct the traversal path
                correct_traversal_path(leaf, *traversal_path, true);
              }
              return leaf->prev_leaf_;
            }
          } else {
            if (leaf->next_leaf_ && leaf->next_leaf_->first_key() <= key) {
              if (traversal_path) {
                // Correct the traversal path
                correct_traversal_path(leaf, *traversal_path, false);
              }
              return leaf->next_leaf_;
            }
          }
        }
        return leaf;
      }
    }
  }
#else
  data_node_type* get_leaf(
      T key, std::vector<TraversalNode>* traversal_path = nullptr) const {
    if (traversal_path) {
      traversal_path->push_back({superroot_, 0});
    }
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      auto node = static_cast<model_node_type*>(cur);
      int bucketID = node->model_.predict(key);
      bucketID =
          std::min<int>(std::max<int>(bucketID, 0), node->num_children_ - 1);
      if (traversal_path) {
        traversal_path->push_back({node, bucketID});
      }
      cur = node->children_[bucketID];
    }

    stats_.num_node_lookups += cur->level_;
    return static_cast<data_node_type*>(cur);
  }
#endif

 private:
  // Make a correction to the traversal path to instead point to the leaf node
  // that is to the left or right of the current leaf node.
  inline void correct_traversal_path(data_node_type* leaf,
                                     std::vector<TraversalNode>& traversal_path,
                                     bool left, int _dup = 0) const {
    if (left) {
      int repeats = hybrid_mode == LEAF_DISK ? 1 << _dup : 1 << leaf->duplication_factor_;
      TraversalNode& tn = traversal_path.back();
      model_node_type* parent = tn.node;
      // First bucket whose pointer is to leaf
      int start_bucketID = tn.bucketID - (tn.bucketID % repeats);
      if (start_bucketID == 0) {
        // Traverse back up the traversal path to make correction
        while (start_bucketID == 0) {
          traversal_path.pop_back();
          repeats = 1 << parent->duplication_factor_;
          tn = traversal_path.back();
          parent = tn.node;
          start_bucketID = tn.bucketID - (tn.bucketID % repeats);
        }
        int correct_bucketID = start_bucketID - 1;
        tn.bucketID = correct_bucketID;
        AlexNode<T, P>* cur = parent->children_[correct_bucketID];
        bool _is_leaf = hybrid_mode != LEAF_DISK ?
                cur->is_leaf_ : parent->is_leaf[correct_bucketID] == 1;
        while (!_is_leaf) {
          auto node = static_cast<model_node_type*>(cur);
          traversal_path.push_back({node, node->num_children_ - 1});
          _is_leaf = hybrid_mode != LEAF_DISK ?
                  (node->children_[node->num_children_ - 1])->is_leaf_ : node->is_leaf[node->num_children_ - 1] == 1;
          if (!_is_leaf) cur = node->children_[node->num_children_ - 1];
        }
        //assert(cur == leaf->prev_leaf_);
      } else {
        tn.bucketID = start_bucketID - 1;
      }
    } else {
      int repeats = hybrid_mode == LEAF_DISK ? 1 << _dup : 1 << leaf->duplication_factor_;
      TraversalNode& tn = traversal_path.back();
      model_node_type* parent = tn.node;
      // First bucket whose pointer is not to leaf
      int end_bucketID = tn.bucketID - (tn.bucketID % repeats) + repeats;
      if (end_bucketID == parent->num_children_) {
        // Traverse back up the traversal path to make correction
        while (end_bucketID == parent->num_children_) {
          traversal_path.pop_back();
          repeats = 1 << parent->duplication_factor_;
          tn = traversal_path.back();
          parent = tn.node;
          end_bucketID = tn.bucketID - (tn.bucketID % repeats) + repeats;
        }
        int correct_bucketID = end_bucketID;
        tn.bucketID = correct_bucketID;
        AlexNode<T, P>* cur = parent->children_[correct_bucketID];
          bool _is_leaf = hybrid_mode != LEAF_DISK ?
                          cur->is_leaf_ : parent->is_leaf[correct_bucketID] == 1;
        while (!_is_leaf) {
          auto node = static_cast<model_node_type*>(cur);
          traversal_path.push_back({node, 0});
          _is_leaf = hybrid_mode != LEAF_DISK ?
                  (node->children_[0])->is_leaf_ : node->is_leaf[0] == 1;
            if (!_is_leaf) cur = node->children_[0];
        }
        //assert(cur == leaf->next_leaf_);
      } else {
        tn.bucketID = end_bucketID;
      }
    }
  }

  // Return left-most data node
  data_node_type* first_data_node() const {
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      cur = static_cast<model_node_type*>(cur)->children_[0];
    }
    return static_cast<data_node_type*>(cur);
  }

  // Return right-most data node
  data_node_type* last_data_node() const {
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      auto node = static_cast<model_node_type*>(cur);
      cur = node->children_[node->num_children_ - 1];
    }
    return static_cast<data_node_type*>(cur);
  }

  // Returns minimum key in the index
  T get_min_key() const { return first_data_node()->first_key(); }

  // Returns maximum key in the index
  T get_max_key() const { return last_data_node()->last_key(); }

  // Link all data nodes together. Used after bulk loading.
  void link_all_data_nodes() {
    data_node_type* prev_leaf = nullptr;
    for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
         node_it.next()) {
      AlexNode<T, P>* cur = node_it.current();
      if (cur->is_leaf_) {
        auto node = static_cast<data_node_type*>(cur);
        if (prev_leaf != nullptr) {
          prev_leaf->next_leaf_ = node;
          node->prev_leaf_ = prev_leaf;
        }
        prev_leaf = node;
      }
    }
  }

  // Link the new data nodes together when old data node is replaced by two new
  // data nodes.
  void link_data_nodes(const data_node_type* old_leaf,
                       data_node_type* left_leaf, data_node_type* right_leaf) {
    if (old_leaf->prev_leaf_ != nullptr) {
      old_leaf->prev_leaf_->next_leaf_ = left_leaf;
    }
    left_leaf->prev_leaf_ = old_leaf->prev_leaf_;
    left_leaf->next_leaf_ = right_leaf;
    right_leaf->prev_leaf_ = left_leaf;
    right_leaf->next_leaf_ = old_leaf->next_leaf_;
    if (old_leaf->next_leaf_ != nullptr) {
      old_leaf->next_leaf_->prev_leaf_ = right_leaf;
    }
  }

  /*** Allocators and comparators ***/

 public:
  Alloc get_allocator() const { return allocator_; }

  Compare key_comp() const { return key_less_; }

 private:
  typename model_node_type::alloc_type model_node_allocator() {
    return typename model_node_type::alloc_type(allocator_);
  }

  typename data_node_type::alloc_type data_node_allocator() {
    return typename data_node_type::alloc_type(allocator_);
  }

  typename model_node_type::pointer_alloc_type pointer_allocator() {
    return typename model_node_type::pointer_alloc_type(allocator_);
  }

  void delete_node(AlexNode<T, P>* node) {
    if (node == nullptr) {
      return;
    } else if (node->is_leaf_) {
      data_node_allocator().destroy(static_cast<data_node_type*>(node));
      data_node_allocator().deallocate(static_cast<data_node_type*>(node), 1);
    } else {
      model_node_allocator().destroy(static_cast<model_node_type*>(node));
      model_node_allocator().deallocate(static_cast<model_node_type*>(node), 1);
    }
  }

  // True if a == b
  template <class K>
  forceinline bool key_equal(const T& a, const K& b) const {
    return !key_less_(a, b) && !key_less_(b, a);
  }

  /*** Bulk loading ***/

 public:
  // values should be the sorted array of key-payload pairs.
  // The number of elements should be num_keys.
  // The index must be empty when calling this method.
  void bulk_load(const V values[], int num_keys) {
    if (stats_.num_keys > 0 || num_keys <= 0) {
      return;
    }
    delete_node(root_node_);  // delete the empty root node from constructor

    stats_.num_keys = num_keys;

    // Build temporary root model, which outputs a CDF in the range [0, 1]
    root_node_ =
        new (model_node_allocator().allocate(1)) model_node_type(0, allocator_);
    T min_key = values[0].first;
    T max_key = values[num_keys - 1].first;
    root_node_->model_.a_ = 1.0 / (max_key - min_key);
    root_node_->model_.b_ = -1.0 * min_key * root_node_->model_.a_;

    // Compute cost of root node
    LinearModel<T> root_data_node_model;
    data_node_type::build_model(values, num_keys, &root_data_node_model,
                                params_.approximate_model_computation);
    DataNodeStats stats;
    root_node_->cost_ = data_node_type::compute_expected_cost(
        values, num_keys, data_node_type::kInitDensity_,
        params_.expected_insert_frac, &root_data_node_model,
        params_.approximate_cost_computation, &stats);

    // Recursively bulk load
    PhyscialAddr RootNodeDisk;
    bulk_load_node(values, num_keys, root_node_, num_keys, &RootNodeDisk,
                   0, 0, 0,
                   &root_data_node_model);

    if (root_node_->is_leaf_) {
      static_cast<data_node_type*>(root_node_)
          ->expected_avg_exp_search_iterations_ = stats.num_search_iterations;
      static_cast<data_node_type*>(root_node_)->expected_avg_shifts_ =
          stats.num_shifts;
    }
    create_superroot();
    update_superroot_key_domain();
//    link_all_data_nodes();
#if DiskSetting
    std::cout << RootNodeDisk.flag << std::endl;
    metanode.root_block_id = RootNodeDisk.block;
    metanode.root_offset = RootNodeDisk.offset;
    metanode.is_leaf = RootNodeDisk.flag;
    metanode.dup_root = RootNodeDisk.duplication_factor_;
    sync_metanode();
    for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
        node_it.next()) {
        if (node_it.current()->is_leaf_ && hybrid_mode == LEAF_DISK)
            delete_node(node_it.current());
        else if (hybrid_mode == ALL_DISK) delete_node(node_it.current());
    }
#endif

  }

 private:
  // Only call this after creating a root node
  void create_superroot() {
    if (!root_node_) return;
    delete_node(superroot_);
    superroot_ = new (model_node_allocator().allocate(1))
        model_node_type(static_cast<short>(root_node_->level_ - 1), allocator_);
    superroot_->num_children_ = 1;
    superroot_->children_ =
        new (pointer_allocator().allocate(1)) AlexNode<T, P>*[1];
    update_superroot_pointer();
  }

  // Updates the key domain based on the min/max keys and retrains the model.
  // Should only be called immediately after bulk loading or when the root node
  // is a data node.
  void update_superroot_key_domain() {
    assert(stats_.num_inserts == 0 || root_node_->is_leaf_);
    istats_.key_domain_min_ = get_min_key();
    istats_.key_domain_max_ = get_max_key();
    istats_.num_keys_at_last_right_domain_resize = stats_.num_keys;
    istats_.num_keys_at_last_left_domain_resize = stats_.num_keys;
    istats_.num_keys_above_key_domain = 0;
    istats_.num_keys_below_key_domain = 0;
    superroot_->model_.a_ =
        1.0 / (istats_.key_domain_max_ - istats_.key_domain_min_);
    superroot_->model_.b_ =
        -1.0 * istats_.key_domain_min_ * superroot_->model_.a_;
  }

  void update_superroot_pointer() {
    superroot_->children_[0] = root_node_;
    superroot_->level_ = static_cast<short>(root_node_->level_ - 1);
  }

  // Recursively bulk load a single node.
  // Assumes node has already been trained to output [0, 1), has cost.
  // Figures out the optimal partitioning of children.
  // node is trained as if it's a model node.
  // data_node_model is what the node's model would be if it were a data node of
  // dense keys.

  // the return flag shows the built node is model or leaf
  bool bulk_load_node(const V values[], int num_keys, AlexNode<T, P>*& node,
                      int total_keys, PhyscialAddr *r_addr, uint8_t duplication_factor_,
                      double expected_avg_exp_search_iterations_, double expected_avg_shifts_,
                      const LinearModel<T>* data_node_model = nullptr) {
    // Automatically convert to data node when it is impossible to be better
    // than current cost
    if (num_keys <= derived_params_.max_data_node_slots *
                        data_node_type::kInitDensity_ &&
        (node->cost_ < kNodeLookupsWeight || node->model_.a_ == 0)) {
      stats_.num_data_nodes++;
      auto data_node = new (data_node_allocator().allocate(1))
          data_node_type(node->level_, derived_params_.max_data_node_slots,
                         key_less_, allocator_);
      data_node->bulk_load(values, num_keys, data_node_model,
                           params_.approximate_model_computation);
      data_node->cost_ = node->cost_;
      data_node->duplication_factor_ = duplication_factor_;
      data_node->expected_avg_exp_search_iterations_ = expected_avg_exp_search_iterations_;
      data_node->expected_avg_shifts_ = expected_avg_shifts_;
      delete_node(node);
      node = data_node;

      #if DiskSetting
      write_data_node(data_node, r_addr, false, values[0].first, values[num_keys-1].first);
      #endif
      return true;
    }

    // Use a fanout tree to determine the best way to divide the key space into
    // child nodes
    std::vector<fanout_tree::FTNode> used_fanout_tree_nodes;
    std::pair<int, double> best_fanout_stats;
    if (experimental_params_.fanout_selection_method == 0) {
      int max_data_node_keys = static_cast<int>(
          derived_params_.max_data_node_slots * data_node_type::kInitDensity_);
      best_fanout_stats = fanout_tree::find_best_fanout_bottom_up<T, P>(
          values, num_keys, node, total_keys, used_fanout_tree_nodes,
          derived_params_.max_fanout, max_data_node_keys,
          params_.expected_insert_frac, params_.approximate_model_computation,
          params_.approximate_cost_computation, key_less_);
    } else if (experimental_params_.fanout_selection_method == 1) {
      best_fanout_stats = fanout_tree::find_best_fanout_top_down<T, P>(
          values, num_keys, node, total_keys, used_fanout_tree_nodes,
          derived_params_.max_fanout, params_.expected_insert_frac,
          params_.approximate_model_computation,
          params_.approximate_cost_computation, key_less_);
    }
    int best_fanout_tree_depth = best_fanout_stats.first;
    double best_fanout_tree_cost = best_fanout_stats.second;

    // Decide whether this node should be a model node or data node
    if (best_fanout_tree_cost < node->cost_ ||
        num_keys > derived_params_.max_data_node_slots *
                       data_node_type::kInitDensity_) {
      // Convert to model node based on the output of the fanout tree
      stats_.num_model_nodes++;
      auto model_node = new (model_node_allocator().allocate(1))
          model_node_type(node->level_, allocator_);
      if (best_fanout_tree_depth == 0) {
        // slightly hacky: we assume this means that the node is relatively
        // uniform but we need to split in
        // order to satisfy the max node size, so we compute the fanout that
        // would satisfy that condition
        // in expectation
        best_fanout_tree_depth =
            static_cast<int>(std::log2(static_cast<double>(num_keys) /
                                       derived_params_.max_data_node_slots)) +
            1;
        used_fanout_tree_nodes.clear();
        int max_data_node_keys = static_cast<int>(
            derived_params_.max_data_node_slots * data_node_type::kInitDensity_);
        fanout_tree::compute_level<T, P>(
            values, num_keys, node, total_keys, used_fanout_tree_nodes,
            best_fanout_tree_depth, max_data_node_keys,
            params_.expected_insert_frac, params_.approximate_model_computation,
            params_.approximate_cost_computation);
      }
      int fanout = 1 << best_fanout_tree_depth;
      model_node->model_.a_ = node->model_.a_ * fanout;
      model_node->model_.b_ = node->model_.b_ * fanout;
      model_node->num_children_ = fanout;
      model_node->children_ =
          new (pointer_allocator().allocate(fanout)) AlexNode<T, P>*[fanout];

      //#if DiskSetting
      model_node->childrenAddrs = new PhyscialAddr[fanout];
      //#endif
      if (hybrid_mode == LEAF_DISK) {
          model_node->is_leaf = new int[fanout];;
          for (int i = 0; i < fanout; i++) {
              // 0 is model 1 is leaf
              model_node->is_leaf[i] = 0;
          }
      }
      // Instantiate all the child nodes and recurse
      int cur = 0;
      for (fanout_tree::FTNode& tree_node : used_fanout_tree_nodes) {
        auto child_node = new (model_node_allocator().allocate(1))
            model_node_type(static_cast<short>(node->level_ + 1), allocator_);
        child_node->cost_ = tree_node.cost;
        child_node->duplication_factor_ =
            static_cast<uint8_t>(best_fanout_tree_depth - tree_node.level);
        int repeats = 1 << child_node->duplication_factor_;
        double left_value = static_cast<double>(cur) / fanout;
        double right_value = static_cast<double>(cur + repeats) / fanout;
        double left_boundary = (left_value - node->model_.b_) / node->model_.a_;
        double right_boundary =
            (right_value - node->model_.b_) / node->model_.a_;
        child_node->model_.a_ = 1.0 / (right_boundary - left_boundary);
        child_node->model_.b_ = -child_node->model_.a_ * left_boundary;
        model_node->children_[cur] = child_node;
        LinearModel<T> child_data_node_model(tree_node.a, tree_node.b);
        bool is_leaf = bulk_load_node(values + tree_node.left_boundary,
                       tree_node.right_boundary - tree_node.left_boundary,
                       model_node->children_[cur], total_keys, // here we write the pointer for model node
                       &(model_node->childrenAddrs[cur]), // here we write the location for leaf node
                       static_cast<uint8_t>(best_fanout_tree_depth - tree_node.level),
                       tree_node.expected_avg_search_iterations,
                       tree_node.expected_avg_shifts,
                       &child_data_node_model);
        model_node->children_[cur]->duplication_factor_ =
            static_cast<uint8_t>(best_fanout_tree_depth - tree_node.level);
        if (model_node->children_[cur]->is_leaf_) {
          static_cast<data_node_type*>(model_node->children_[cur])
              ->expected_avg_exp_search_iterations_ =
              tree_node.expected_avg_search_iterations;
          static_cast<data_node_type*>(model_node->children_[cur])
              ->expected_avg_shifts_ = tree_node.expected_avg_shifts;
        }
        if (is_leaf && hybrid_mode == LEAF_DISK) {
            model_node->is_leaf[cur] = 1;
        }
        // if the returned one is leaf,
        for (int i = cur + 1; i < cur + repeats; i++) {
            // set flag
            if (is_leaf && hybrid_mode == LEAF_DISK) model_node->is_leaf[i] = 1;
            model_node->children_[i] = model_node->children_[cur];
          #if DiskSetting
            model_node->childrenAddrs[i] = model_node->childrenAddrs[cur];
          #endif
        }
        cur += repeats;
      }

      delete_node(node);
      node = model_node;

      #if DiskSetting
      if (hybrid_mode == ALL_DISK)
        write_model_node(model_node, r_addr);
      #endif

      return false;
    } else {
      // Convert to data node
      stats_.num_data_nodes++;
      auto data_node = new (data_node_allocator().allocate(1))
          data_node_type(node->level_, derived_params_.max_data_node_slots,
                         key_less_, allocator_);
      data_node->bulk_load(values, num_keys, data_node_model,
                           params_.approximate_model_computation);
      data_node->cost_ = node->cost_;
      data_node->duplication_factor_ = duplication_factor_;
      data_node->expected_avg_exp_search_iterations_ = expected_avg_exp_search_iterations_;
      data_node->expected_avg_shifts_ = expected_avg_shifts_;
      delete_node(node);
      node = data_node;

      #if DiskSetting
      write_data_node(data_node, r_addr, false, values[0].first, values[num_keys-1].first);
      #endif

      return true;
    }
  }

  // Caller needs to set the level, duplication factor, and neighbor pointers of
  // the returned data node
  data_node_type* bulk_load_leaf_node_from_existing(
      const data_node_type* existing_node, int left, int right,
      bool compute_cost = true, const fanout_tree::FTNode* tree_node = nullptr,
      bool reuse_model = false, bool keep_left = false,
      bool keep_right = false) {
    auto node = new (data_node_allocator().allocate(1))
        data_node_type(key_less_, allocator_);
    stats_.num_data_nodes++;
    if (tree_node) {
      // Use the model and num_keys saved in the tree node so we don't have to
      // recompute it
      LinearModel<T> precomputed_model(tree_node->a, tree_node->b);
      node->bulk_load_from_existing(existing_node, left, right, keep_left,
                                    keep_right, &precomputed_model,
                                    tree_node->num_keys);
    } else if (reuse_model) {
      // Use the model from the existing node
      // Assumes the model is accurate
      int num_actual_keys = existing_node->num_keys_in_range(left, right);
      LinearModel<T> precomputed_model(existing_node->model_);
      precomputed_model.b_ -= left;
      precomputed_model.expand(static_cast<double>(num_actual_keys) /
                               (right - left));
      node->bulk_load_from_existing(existing_node, left, right, keep_left,
                                    keep_right, &precomputed_model,
                                    num_actual_keys);
    } else {
      node->bulk_load_from_existing(existing_node, left, right, keep_left,
                                    keep_right);
    }
    node->max_slots_ = derived_params_.max_data_node_slots;
    if (compute_cost) {
      node->cost_ = node->compute_expected_cost(existing_node->frac_inserts());
    }
    return node;
  }

  /*** Lookup ***/

 public:
  // Looks for an exact match of the key
  // If the key does not exist, returns an end iterator
  // If there are multiple keys with the same value, returns an iterator to the
  // right-most key
  // If you instead want an iterator to the left-most key with the input value,
  // use lower_bound()
  typename self_type::Iterator find(const T& key) {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_key(key);
    if (idx < 0) {
      return end();
    } else {
      return Iterator(leaf, idx);
    }
  }

  typename self_type::ConstIterator find(const T& key) const {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_key(key);
    if (idx < 0) {
      return cend();
    } else {
      return ConstIterator(leaf, idx);
    }
  }

  size_t count(const T& key) {
    ConstIterator it = lower_bound(key);
    size_t num_equal = 0;
    while (!it.is_end() && key_equal(it.key(), key)) {
      num_equal++;
      ++it;
    }
    return num_equal;
  }

  // Returns an iterator to the first key no less than the input value
  typename self_type::Iterator lower_bound(const T& key) {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_lower(key);
    return Iterator(leaf, idx);  // automatically handles the case where idx ==
                                 // leaf->data_capacity
  }

  typename self_type::ConstIterator lower_bound(const T& key) const {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_lower(key);
    return ConstIterator(leaf, idx);  // automatically handles the case where
                                      // idx == leaf->data_capacity
  }

  // Returns an iterator to the first key greater than the input value
  typename self_type::Iterator upper_bound(const T& key) {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_upper(key);
    return Iterator(leaf, idx);  // automatically handles the case where idx ==
                                 // leaf->data_capacity
  }

  typename self_type::ConstIterator upper_bound(const T& key) const {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_upper(key);
    return ConstIterator(leaf, idx);  // automatically handles the case where
                                      // idx == leaf->data_capacity
  }

  std::pair<Iterator, Iterator> equal_range(const T& key) {
    return std::pair<Iterator, Iterator>(lower_bound(key), upper_bound(key));
  }

  std::pair<ConstIterator, ConstIterator> equal_range(const T& key) const {
    return std::pair<ConstIterator, ConstIterator>(lower_bound(key),
                                                   upper_bound(key));
  }

  // Directly returns a pointer to the payload found through find(key)
  // This avoids the overhead of creating an iterator
  // Returns null pointer if there is no exact match of the key
  P* get_payload(const T& key) const {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    int idx = leaf->find_key(key);
    if (idx < 0) {
      return nullptr;
    } else {
      return &(leaf->get_payload(idx));
    }
  }

  // Looks for the last key no greater than the input value
  // Conceptually, this is equal to the last key before upper_bound()
  typename self_type::Iterator find_last_no_greater_than(const T& key) {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    const int idx = leaf->upper_bound(key) - 1;
    if (idx >= 0) {
      return Iterator(leaf, idx);
    }

    // Edge case: need to check previous data node(s)
    while (true) {
      if (leaf->prev_leaf_ == nullptr) {
        return Iterator(leaf, 0);
      }
      leaf = leaf->prev_leaf_;
      if (leaf->num_keys_ > 0) {
        return Iterator(leaf, leaf->last_pos());
      }
    }
  }

  // Directly returns a pointer to the payload found through
  // find_last_no_greater_than(key)
  // This avoids the overhead of creating an iterator
  P* get_payload_last_no_greater_than(const T& key) {
    stats_.num_lookups++;
    data_node_type* leaf = get_leaf(key);
    const int idx = leaf->upper_bound(key) - 1;
    if (idx >= 0) {
      return &(leaf->get_payload(idx));
    }

    // Edge case: Need to check previous data node(s)
    while (true) {
      if (leaf->prev_leaf_ == nullptr) {
        return &(leaf->get_payload(leaf->first_pos()));
      }
      leaf = leaf->prev_leaf_;
      if (leaf->num_keys_ > 0) {
        return &(leaf->get_payload(leaf->last_pos()));
      }
    }
  }

  typename self_type::Iterator begin() {
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      cur = static_cast<model_node_type*>(cur)->children_[0];
    }
    return Iterator(static_cast<data_node_type*>(cur), 0);
  }

  typename self_type::Iterator end() {
    Iterator it = Iterator();
    it.cur_leaf_ = nullptr;
    it.cur_idx_ = 0;
    return it;
  }

  typename self_type::ConstIterator cbegin() const {
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      cur = static_cast<model_node_type*>(cur)->children_[0];
    }
    return ConstIterator(static_cast<data_node_type*>(cur), 0);
  }

  typename self_type::ConstIterator cend() const {
    ConstIterator it = ConstIterator();
    it.cur_leaf_ = nullptr;
    it.cur_idx_ = 0;
    return it;
  }

  typename self_type::ReverseIterator rbegin() {
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      auto model_node = static_cast<model_node_type*>(cur);
      cur = model_node->children_[model_node->num_children_ - 1];
    }
    auto data_node = static_cast<data_node_type*>(cur);
    return ReverseIterator(data_node, data_node->data_capacity_ - 1);
  }

  typename self_type::ReverseIterator rend() {
    ReverseIterator it = ReverseIterator();
    it.cur_leaf_ = nullptr;
    it.cur_idx_ = 0;
    return it;
  }

  typename self_type::ConstReverseIterator crbegin() const {
    AlexNode<T, P>* cur = root_node_;

    while (!cur->is_leaf_) {
      auto model_node = static_cast<model_node_type*>(cur);
      cur = model_node->children_[model_node->num_children_ - 1];
    }
    auto data_node = static_cast<data_node_type*>(cur);
    return ConstReverseIterator(data_node, data_node->data_capacity_ - 1);
  }

  typename self_type::ConstReverseIterator crend() const {
    ConstReverseIterator it = ConstReverseIterator();
    it.cur_leaf_ = nullptr;
    it.cur_idx_ = 0;
    return it;
  }

  /*** Insert ***/

 public:
  std::pair<Iterator, bool> insert(const V& value) {
    return insert(value.first, value.second);
  }

  template <class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    for (auto it = first; it != last; ++it) {
      insert(*it);
    }
  }


  bool should_expand_right_disk() {
    return (metanode.is_leaf == 0 &&
            ((istats_.num_keys_above_key_domain >= kMinOutOfDomainKeys &&
              istats_.num_keys_above_key_domain >=
                  kOutOfDomainToleranceFactor *
                      (stats_.num_keys /
                           istats_.num_keys_at_last_right_domain_resize -
                       1)) ||
             istats_.num_keys_above_key_domain >= kMaxOutOfDomainKeys));
  }

  bool should_expand_left_disk() {
    return (metanode.is_leaf == 0 &&
            ((istats_.num_keys_below_key_domain >= kMinOutOfDomainKeys &&
              istats_.num_keys_below_key_domain >=
                  kOutOfDomainToleranceFactor *
                      (stats_.num_keys /
                           istats_.num_keys_at_last_left_domain_resize -
                       1)) ||
             istats_.num_keys_below_key_domain >= kMaxOutOfDomainKeys));
  }

  void expand_root_disk(T key, bool expand_left) {
    char data[BlockSize];
    sm->get_block(metanode.root_block_id, data);
    ModelHeaderOnDisk *root_disk = (ModelHeaderOnDisk *) (data + metanode.root_offset);
    ModelHeaderOnDisk new_root_disk_;

    T domain_size = istats_.key_domain_max_ - istats_.key_domain_min_;
    int expansion_factor;
    T new_domain_min = istats_.key_domain_min_;
    T new_domain_max = istats_.key_domain_max_;
    int block_outer = -1;
    int offset_outer = -1;
    if (expand_left) {
      auto key_difference = static_cast<double>(istats_.key_domain_min_ -
                                                std::min(key, istats_.min_key_in_left));
      expansion_factor = pow_2_round_up(static_cast<int>(
          std::ceil((key_difference + domain_size) / domain_size)));
      T half_expandable_domain =
          istats_.key_domain_max_ / 2 - std::numeric_limits<T>::lowest() / 2;
      T half_expanded_domain_size = expansion_factor / 2 * domain_size;
      if (half_expandable_domain < half_expanded_domain_size) {
        new_domain_min = std::numeric_limits<T>::lowest();
      } else {
        new_domain_min = istats_.key_domain_max_;
        new_domain_min -= half_expanded_domain_size;
        new_domain_min -= half_expanded_domain_size;
      }
      istats_.num_keys_at_last_left_domain_resize = stats_.num_keys;
      istats_.num_keys_below_key_domain = 0;
      block_outer = metanode.left_most_block;
      offset_outer = metanode.left_most_offset;
    } else {
      auto key_difference = static_cast<double>(std::max(key, istats_.max_key_in_right) -
                                                istats_.key_domain_max_);
      expansion_factor = pow_2_round_up(static_cast<int>(
          std::ceil((key_difference + domain_size) / domain_size)));
      // Check for overflow. To avoid overflow on signed types while doing
      // this check, we do comparisons using half of the relevant quantities.
      T half_expandable_domain =
          std::numeric_limits<T>::max() / 2 - istats_.key_domain_min_ / 2;
      T half_expanded_domain_size = expansion_factor / 2 * domain_size;
      if (half_expandable_domain < half_expanded_domain_size) {
        new_domain_max = std::numeric_limits<T>::max();
      } else {
        new_domain_max = istats_.key_domain_min_;
        new_domain_max += half_expanded_domain_size;
        new_domain_max += half_expanded_domain_size;
      }
      istats_.num_keys_at_last_right_domain_resize = stats_.num_keys;
      istats_.num_keys_above_key_domain = 0;
      block_outer = metanode.right_most_block;
      offset_outer = metanode.right_most_offset;
    }

    int new_nodes_start;
    int new_nodes_end;
    PhyscialAddr *new_root_children = nullptr;
    if (static_cast<size_t>(root_disk->children_number) * expansion_factor <=
        static_cast<size_t>(derived_params_.max_fanout)) {
      // Expand root node
      stats_.num_model_node_expansions++;
      stats_.num_model_node_expansion_pointers += root_disk->children_number;

      int new_num_children = root_disk->children_number * expansion_factor;

      new_root_children = new PhyscialAddr[new_num_children];
      int copy_start;
      if (expand_left) {
        copy_start = new_num_children - root_disk->children_number;
        new_nodes_start = 0;
        new_nodes_end = copy_start;
        root_disk->intercept += new_num_children - root_disk->children_number;
      } else {
        copy_start = 0;
        new_nodes_start = root_disk->children_number;
        new_nodes_end = new_num_children;
      }
      // load old children pointers
      PhyscialAddr *old_children = new PhyscialAddr[root_disk->children_number];
      {
        long offset = metanode.root_block_id * BlockSize + metanode.root_offset
                + ModelHeaderOnDiskSize + root_disk->empty_byte_count;
        sm->read_arbitrary(old_children, offset, root_disk->children_number * PhyscialAddrSize);
      }
      for (int i = 0; i < root_disk->children_number; i++) {
        new_root_children[copy_start + i] =old_children[i];
      }
      delete []old_children;
      root_disk->children_number = new_num_children;
      // flush on disk & update root node addr
      PhyscialAddr new_root_addr;
      ModelHeaderOnDisk r = write_model_node(*root_disk, new_root_children, &new_root_addr);
      metanode.root_block_id = new_root_addr.block;
      metanode.root_offset = new_root_addr.offset;
      metanode.is_leaf = new_root_addr.flag;
      metanode.dup_root = new_root_addr.duplication_factor_;
      memcpy(&new_root_disk_, &r, ModelHeaderOnDiskSize);
    } else {
      ModelHeaderOnDisk new_root_disk;
      new_root_disk.slope = root_disk->slope / root_disk->children_number;
      new_root_disk.intercept = root_disk->intercept / root_disk->children_number;
      new_root_disk.level = root_disk->level - 1;
      if (expand_left) new_root_disk.intercept += expansion_factor - 1;

      new_root_disk.children_number = expansion_factor;
      new_root_children = new PhyscialAddr[expansion_factor];
      if (expand_left) {
        new_root_children[expansion_factor - 1].block = metanode.root_block_id;
        new_root_children[expansion_factor - 1].offset = metanode.root_offset;
        new_nodes_start = 0;
      } else {
        new_root_children[0].block = metanode.root_block_id;
        new_root_children[0].offset = metanode.root_offset;
        new_nodes_start = 1;
      }
      new_nodes_end = new_nodes_start + expansion_factor - 1;
      PhyscialAddr new_root_addr;
      ModelHeaderOnDisk r = write_model_node(new_root_disk, new_root_children, &new_root_addr);
      metanode.root_block_id = new_root_addr.block;
      metanode.root_offset = new_root_addr.offset;
      metanode.is_leaf = new_root_addr.flag;
      metanode.dup_root = new_root_addr.duplication_factor_;
      memcpy(&new_root_disk_, &r, ModelHeaderOnDiskSize);
    }

    int in_bounds_new_nodes_start = new_nodes_start;
    int in_bounds_new_nodes_end = new_nodes_end;
    //static_cast<int>(new_root_disk_.slope * static_cast<double>(new_domain_min) + new_root_disk_.intercept);
    if (expand_left) {
      in_bounds_new_nodes_start =
          std::max(new_nodes_start, static_cast<int>(new_root_disk_.slope * static_cast<double>(new_domain_min) + new_root_disk_.intercept));
    } else {
      in_bounds_new_nodes_end =
      std::min(new_nodes_end, static_cast<int>(new_root_disk_.slope * static_cast<double>(new_domain_max) + new_root_disk_.intercept) + 1);
    }

    int n = new_root_disk_.children_number - (new_nodes_end - new_nodes_start);
    auto new_node_duplication_factor =
        static_cast<uint8_t>(log_2_round_down(n));
    // build old_data_node in main memory

    //////////// build old node ///////////
      DataHeaderOnDisk old_dh;
      long _start_offset_ = BlockSize * block_outer + offset_outer;
      data_node_sm->read_arbitrary(&old_dh, _start_offset_, DataHeaderOnDiskSize);
      _start_offset_ += (DataHeaderOnDiskSize + old_dh.empty_byte_count_1);
      uint64_t bitmaps_[old_dh.bitmap_size];
      data_node_sm->read_arbitrary(bitmaps_, _start_offset_, old_dh.bitmap_size * sizeof(uint64_t));
      _start_offset_ += old_dh.bitmap_size * sizeof(uint64_t) + old_dh.empty_byte_count_2;
      ItemOnDisk *items = new ItemOnDisk [old_dh.data_capacity];
      data_node_sm->read_arbitrary(items, _start_offset_, old_dh.data_capacity*ItemOnDiskSize);
      auto old_data_node = new (data_node_allocator().allocate(1)) data_node_type(old_dh.level, derived_params_.max_data_node_slots, key_less_, allocator_);
      if (expand_left) {
        old_dh.duplication_factor_ = metanode.dup_left;
      } else {
        old_dh.duplication_factor_ = metanode.dup_right;
      }
      init_from_disk(old_data_node, old_dh, bitmaps_, items);
   ///////////////////////////////////

    if (expand_left) {
      T left_boundary_value = istats_.key_domain_min_;
      int left_boundary = old_data_node->lower_bound(left_boundary_value);//lower_bound_disk(block_outer, offset_outer, left_boundary_value); //outermost_node->lower_bound(left_boundary_value);
      int next_block = block_outer;
      int next_offset = offset_outer;

      for (int i = new_nodes_end; i > new_nodes_start; i -= n) {
        if (i <= in_bounds_new_nodes_start) {
          // Do not initialize nodes that fall outside the key type's domain
          break;
        }
        int right_boundary = left_boundary;
        if (i - n <= in_bounds_new_nodes_start) {
          left_boundary = 0;
        } else {
          left_boundary_value -= domain_size;
          left_boundary = old_data_node->lower_bound(left_boundary_value);///lower_bound_disk(block_outer, offset_outer, left_boundary_value);/// outermost_node->lower_bound(left_boundary_value);
        }
        data_node_type* new_node = bulk_load_leaf_node_from_existing(
            old_data_node, left_boundary, right_boundary, true);
        new_node->level_ = static_cast<short>(new_root_disk_.level + 1);
        new_node->duplication_factor_ = new_node_duplication_factor;
        PhyscialAddr _addr;
        metanode_data_node.last_data_node_block = -1;
        metanode_data_node.last_data_node_offset = -1;
        DataHeaderOnDisk new_dnh = write_data_node(new_node, &_addr, new_node->first_key(), new_node->last_key(), next_block, next_offset);
        // update in write_data_node
        // if (next_block != -1) {
        //   // next->prev_leaf_ = new_node;
        //   cur_dh.pre_leaf_block = _addr.block;
        //   cur_dh.pre_leaf_offset = _addr.offset;
        //   // flush disk
        //   sm->write_arbitrary(next_block * BlockSize + next_offset, &cur_dh, DataHeaderOnDiskSize);
        // }
        // memcpy(&cur_dh, &new_dnh, DataHeaderOnDiskSize);
        next_block = _addr.block;
        next_offset = _addr.offset;
        delete_node(new_node);
        for (int j = i - 1; j >= i - n; j--) {
          new_root_children[j].block = _addr.block;
          new_root_children[j].offset = _addr.offset;
          new_root_children[j].flag = _addr.flag;
          new_root_children[j].duplication_factor_ = _addr.duplication_factor_;
        }
      }
    }
    else {
      T right_boundary_value = istats_.key_domain_max_;
      int right_boundary = old_data_node->lower_bound(right_boundary_value);
      int pre_block = block_outer;
      int pre_offset = offset_outer;
      for (int i = new_nodes_start; i < new_nodes_end; i += n) {
        if (i >= in_bounds_new_nodes_end) {
          // Do not initialize nodes that fall outside the key type's domain
          break;
        }
        int left_boundary = right_boundary;
        if (i + n >= in_bounds_new_nodes_end) {
          right_boundary = old_data_node->data_capacity_;
        } else {
          right_boundary_value += domain_size;
          right_boundary = old_data_node->lower_bound(right_boundary_value);
        }
        data_node_type* new_node = bulk_load_leaf_node_from_existing(
            old_data_node, left_boundary, right_boundary, true);
        new_node->level_ = static_cast<short>(new_root_disk_.level + 1);
        new_node->duplication_factor_ = new_node_duplication_factor;
        PhyscialAddr _addr;
        metanode_data_node.last_data_node_block = pre_block;
        metanode_data_node.last_data_node_offset = pre_offset;
        DataHeaderOnDisk new_dnh = write_data_node(new_node, &_addr, false, new_node->first_key(), new_node->last_key());
        // memcpy(&cur_dh, &new_dnh, DataHeaderOnDiskSize);
        pre_block = _addr.block;
        pre_offset = _addr.offset;
        delete_node(new_node);
        for (int j = i; j < i + n; j++) {
          new_root_children[j].block = _addr.block;
          new_root_children[j].offset = _addr.offset;
          new_root_children[j].flag = _addr.flag;
          new_root_children[j].duplication_factor_ = _addr.duplication_factor_;
        }
      }
    }

    if (expand_left) {
      old_data_node->erase_range(new_domain_min, istats_.key_domain_min_);
      PhyscialAddr last_new_leaf =  new_root_children[new_nodes_end - 1];
      metanode_data_node.last_data_node_block = last_new_leaf.block;
      metanode_data_node.last_data_node_offset = last_new_leaf.offset;
      PhyscialAddr old_addr;
      old_addr.block = block_outer;
      old_addr.offset = offset_outer;
      // we update the its pre_node's next node in write_data_node function
      write_data_node(old_data_node, &old_addr, true, old_data_node->first_key(), old_data_node->last_key(), old_dh.next_leaf_block, old_dh.next_leaf_offset);
    }
    else {
      old_data_node->erase_range(istats_.key_domain_max_, new_domain_max,
                                  true);
      PhyscialAddr old_addr;
      old_addr.block = block_outer;
      old_addr.offset = offset_outer;

      PhyscialAddr first_new_leaf = new_root_children[new_nodes_start];
      metanode_data_node.last_data_node_block = old_dh.pre_leaf_block;
      metanode_data_node.last_data_node_offset = old_dh.pre_leaf_offset;
      write_data_node(old_data_node, &old_addr, true, old_data_node->first_key(), old_data_node->last_key(), first_new_leaf.block, first_new_leaf.offset);
    }
    delete []items;
    PhyscialAddr new_root_addr;
    write_model_node(new_root_disk_, new_root_children, &new_root_addr);
    metanode.root_block_id = new_root_addr.block;
    metanode.root_offset = new_root_addr.offset;
    metanode.is_leaf = new_root_addr.flag;
    metanode.dup_root = new_root_addr.duplication_factor_;
    // flush metanode
    sync_metanode();
    istats_.key_domain_min_ = new_domain_min;
    istats_.key_domain_max_ = new_domain_max;
    delete_node(old_data_node);
  }
  typedef struct {
    PhyscialAddr addr; //node start position
    PhyscialAddr addr_in_p;
    ModelHeaderOnDisk header; // node header
    int child_offset;
  } TraversalPathDisk;
  PhyscialAddr EMPTY_ADDR;
  ModelHeaderOnDisk EMPTY_M_HEADER;


  void correct_traversal_path_on_disk() {

  }
  PhyscialAddr get_leaf_addr(T key, char *block_data2, DataHeaderOnDisk *LEAF_NODE, int *r_block_count, int *p_block, int *p_offset, std::vector<TraversalPathDisk> *traversal_path = nullptr) {
    if (traversal_path) {
      // traversal_path->push_back({-1,0,0}); // specical flag;
      traversal_path->push_back({EMPTY_ADDR,EMPTY_ADDR, EMPTY_M_HEADER,-1});
    }
    PhyscialAddr addr;


    char block_data[BlockSize];
    int block = metanode.root_block_id;
    int offset = metanode.root_offset;
    sm->get_block(block, block_data2);
    *r_block_count += 1;
    uint8_t duplication_factor_ = metanode.dup_root;
    ModelHeaderOnDisk *cur;
    int block_in_p = 0;
    int offset_in_p = 0;
    cur = (ModelHeaderOnDisk *)(block_data2 + offset);
    while (true) {
      double pred_pos_d = cur->slope * static_cast<double>(key) + cur->intercept;
      int pred_pos_i = static_cast<int>(pred_pos_d);
      pred_pos_i = std::min<int>(std::max<int>(pred_pos_i, 0), cur->children_number - 1);

      if (traversal_path) {
        //current node's address, current node's pointer's address in parent, current node, next_child
        traversal_path->push_back({{ALEXModelType, block, offset, duplication_factor_}, {ALEXModelType, block_in_p, offset_in_p, duplication_factor_}, *cur, pred_pos_i});
      }

      long _offset = (block * BlockSize + offset + ModelHeaderOnDiskSize + cur->empty_byte_count) % BlockSize;
      int c1 = (BlockSize - _offset)/PhyscialAddrSize;
      int block_in_inode  = 0;
      int offset_in_inode = 0;
      if (_offset == 0) c1 = 0;
      if (pred_pos_i < c1) {
        block_in_inode = block;
        offset_in_inode = _offset + pred_pos_i * PhyscialAddrSize;
      } else {
        block_in_inode = block;
        block_in_inode += 1;
        block_in_inode += (pred_pos_i - c1) / MaxModeCount;
        offset_in_inode = ((pred_pos_i - c1) % MaxModeCount) * PhyscialAddrSize;
      }
      if (block_in_inode != block) {
          sm->get_block(block_in_inode, block_data2);
          *r_block_count += 1;
          block = block_in_inode;
      }
      block_in_p = block_in_inode;
      offset_in_p = offset_in_inode;
      PhyscialAddr child_addr = *(PhyscialAddr *)(block_data2 + offset_in_inode);
      duplication_factor_ = child_addr.duplication_factor_;
      if (child_addr.flag == ALEXModelType) {
          if (child_addr.block != block) {
            sm->get_block(child_addr.block, block_data2);
            *r_block_count += 1;
            block = child_addr.block;
          }
          offset = child_addr.offset;
          cur = (ModelHeaderOnDisk *)(block_data2 + offset);
          continue;
      }

    // find leaf!!!
      *p_block = block_in_p;
      *p_offset = offset_in_p;
      block = child_addr.block;
      offset = child_addr.offset;
      uint8_t duplication_factor_ = child_addr.duplication_factor_;
      data_node_sm->get_block(child_addr.block, block_data2);
      *r_block_count += 1;
      // data_node = (_meta_data_node*)(block_data + _block_if.offset);
      DataHeaderOnDisk *data_node = (DataHeaderOnDisk *)( block_data2 + child_addr.offset);
      // memcpy(&data_node, block_data + _block_if.offset, data_meta_size);
      int pred_pos_i_rounded = static_cast<int>(pred_pos_d + 0.5);
      double tolerance = 10 * std::numeric_limits<double>::epsilon() * pred_pos_d;
      if (std::abs(pred_pos_d - pred_pos_i_rounded) <= tolerance) {
          if (pred_pos_i_rounded <= pred_pos_d) {
              if (data_node->pre_leaf_block != -1) {
                  // char pre_block_data[BLOCK_SIZE];
                  data_node_sm->get_block(data_node->pre_leaf_block, block_data);
                  *r_block_count += 1;
                  DataHeaderOnDisk *pre_data_node = (DataHeaderOnDisk *) (block_data + data_node->pre_leaf_offset);
                  if (pre_data_node->last_key >= key) {
                      memcpy(LEAF_NODE, pre_data_node, DataHeaderOnDiskSize);
                      printf("correct me!!!\n");
                      correct_traversal_path_on_disk();
                      addr.block = data_node->pre_leaf_block;
                      addr.offset = data_node->pre_leaf_offset;
                      memcpy(block_data2, block_data, BlockSize);
                      return addr;
                  }
              }
          } else {
              if (data_node->next_leaf_block != -1) {
                  data_node_sm->get_block(data_node->next_leaf_block, block_data);
                  *r_block_count += 1;
                  DataHeaderOnDisk *next_data_node = (DataHeaderOnDisk *)(block_data + data_node->next_leaf_offset);
                  // memcpy(&next_data_node, block_data + data_node.next_leaf_offset, data_meta_size);
                  if (next_data_node->first_key <= key) {
                      memcpy(LEAF_NODE, next_data_node, DataHeaderOnDiskSize);
                      printf("correct me!!!\n");
                      correct_traversal_path_on_disk();
                      // get_payload_on_disk(key, data_node->next_leaf_block, data_node->next_leaf_offset, block_data, r_block_count, POS);
                      // *l_block = data_node->next_leaf_block;
                      // *l_offset = data_node->next_leaf_offset;
                      addr.block = data_node->next_leaf_block;
                      addr.offset = data_node->next_leaf_offset;
                      memcpy(block_data2, block_data, BlockSize);
                      return addr;
                  }
              }
          }
      }
      addr.block = block;
      addr.offset = offset;
      addr.duplication_factor_ = duplication_factor_;
      memcpy(LEAF_NODE, data_node, DataHeaderOnDiskSize);
      return addr;
    }
  }
  // Const confs.
  static constexpr double kAppendMostlyThreshold_Disk = 0.9;
  static constexpr double kMaxDensity_Disk = 0.8;
  static constexpr double kInitDensity_Disk = 0.7;
  static constexpr double kMinDensity_Disk = 0.6;

  // Empirical average number of exponential search iterations per operation
  // (either lookup or insert)
  double exp_search_iterations_per_operation(DataHeaderOnDisk leaf_node) {
    if (leaf_node.num_inserts_ + leaf_node.num_lookups_ == 0) {
      return 0;
    }
    return leaf_node.num_exp_search_iterations_ /
           static_cast<double>(leaf_node.num_inserts_ + leaf_node.num_lookups_);
  }

  double shifts_per_insert(DataHeaderOnDisk leaf_node) {
    if (leaf_node.num_inserts_ == 0) {
      return 0;
    }
    return leaf_node.num_shifts_ / static_cast<double>(leaf_node.num_inserts_);
  }

  double empirical_cost(DataHeaderOnDisk leaf_node) {
    if (leaf_node.num_inserts_ + leaf_node.num_lookups_ == 0) {
      return 0;
    }
    double frac_inserts =
        static_cast<double>(leaf_node.num_inserts_) / (leaf_node.num_inserts_ + leaf_node.num_lookups_);
    return kExpSearchIterationsWeight * exp_search_iterations_per_operation(leaf_node) +
           kShiftsWeight * shifts_per_insert(leaf_node) * frac_inserts;
  }

  // Empirical fraction of operations (either lookup or insert) that are inserts
  double frac_inserts(DataHeaderOnDisk leaf_node) const {
    int num_ops = leaf_node.num_inserts_ + leaf_node.num_lookups_;
    if (num_ops == 0) {
      return 0;  // if no operations, assume no inserts
    }
    return static_cast<double>(leaf_node.num_inserts_) / (leaf_node.num_inserts_ + leaf_node.num_lookups_);
  }



  bool significant_cost_deviation(DataHeaderOnDisk leaf_node) {
    double emp_cost = empirical_cost(leaf_node);
    return emp_cost > kNodeLookupsWeight && emp_cost > 1.5 * leaf_node.cost;
  }

  bool catastrophic_cost(DataHeaderOnDisk leaf_node) {
    return shifts_per_insert(leaf_node) > 100 || leaf_node.expected_avg_shifts_ > 100;
  }

  inline bool is_append_mostly_right(DataHeaderOnDisk leaf_node) const {
    return static_cast<double>(leaf_node.num_right_out_of_bounds_inserts_) /
               leaf_node.num_inserts_ >
           kAppendMostlyThreshold_Disk;
  }

  inline bool is_append_mostly_left(DataHeaderOnDisk leaf_node) const {
    return static_cast<double>(leaf_node.num_left_out_of_bounds_inserts_) / leaf_node.num_inserts_ >
           kAppendMostlyThreshold_Disk;
  }


  bool check_exist(int pos, int block, int offset, DataHeaderOnDisk leaf_node, char *DATA, int *last_block) {
    int bitmap_pos = pos >> 6;
    int bit_pos = pos - (bitmap_pos << 6);
    long _start_offset_ = BlockSize * block + offset;
    _start_offset_ += (DataHeaderOnDiskSize + leaf_node.empty_byte_count_1 + bitmap_pos*sizeof(uint64_t));
    int _block = _start_offset_ / BlockSize;
    int _offset = _start_offset_ % BlockSize;
    if (*last_block != _block) {
        data_node_sm->get_block(_block, DATA);
        *last_block = _block;
    }
    uint64_t _bt = *((uint64_t *)(DATA + _offset));
//    data_node_sm->read_arbitrary(&_bt, _start_offset_, sizeof(uint64_t));
    return static_cast<bool>(_bt & (1ULL << bit_pos));
  }


#if Profiling
  long long INSERT_TIME_IN_LEAF_NODE = 0;
  long long SEARCH_TIME_IN_LEAF_NODE = 0;
  long long SMO_TIME_IN_LEAF_NODE = 0;
  long long MAINTAIN_TIME_IN_LEAF_NODE = 0;
  int SMO_COUNT  = 0;
  int SMO_TUPLES = 0;
#endif

  std::pair<int, int> insert_into_leaf_node(T key, P value, char *data, DataHeaderOnDisk *leaf_node, PhyscialAddr *addr,
                                            int p_block, int p_offset,
                                            std::vector<TraversalPathDisk> traversal_path,
                                            std::vector<TraversalNode> traversal_path2) {
#if Profiling
      INSERT_TIME_IN_LEAF_NODE = 0;
      SEARCH_TIME_IN_LEAF_NODE = 0;
      SMO_TIME_IN_LEAF_NODE = 0;
      MAINTAIN_TIME_IN_LEAF_NODE = 0;
#endif

    if (leaf_node->num_inserts_ % 64 == 0 && catastrophic_cost(*leaf_node)) {
      return {2, -1};
    }

#if Profiling
    auto smo_start = std::chrono::high_resolution_clock::now();
#endif
    // Check if node is full (based on expansion_threshold)
    if (leaf_node->num_keys >=leaf_node->expansion_threshold_) {
      if (significant_cost_deviation(*leaf_node)) {
        return {1, -1};
      }
      if (catastrophic_cost(*leaf_node)) {
        return {2, -1};
      }
      if (leaf_node->num_keys > leaf_node->max_slots_ * kMinDensity_Disk) {
        return {3, -1};
      }
      // Expand
#if Profiling
        SMO_COUNT += 1;
#endif
      bool keep_left = is_append_mostly_right(*leaf_node);
      bool keep_right = is_append_mostly_left(*leaf_node);
      // resize(kMinDensity_, false, keep_left, keep_right);
      {
        long _start_offset_ = BlockSize * addr->block + addr->offset;
        // data_node_sm->read_arbitrary(&old_dh, _start_offset_, DataHeaderOnDiskSize);
        _start_offset_ += (DataHeaderOnDiskSize + leaf_node->empty_byte_count_1);
        uint64_t *bitmaps_ = new uint64_t[leaf_node->bitmap_size];
        data_node_sm->read_arbitrary(bitmaps_, _start_offset_, leaf_node->bitmap_size * sizeof(uint64_t));
        _start_offset_ += leaf_node->bitmap_size * sizeof(uint64_t) + leaf_node->empty_byte_count_2;
        ItemOnDisk *items = new ItemOnDisk[leaf_node->data_capacity];
        data_node_sm->read_arbitrary(items, _start_offset_, leaf_node->data_capacity*ItemOnDiskSize);
        auto old_data_node = new (data_node_allocator().allocate(1)) data_node_type(leaf_node->level, derived_params_.max_data_node_slots, key_less_, allocator_);

        leaf_node->duplication_factor_ = addr->duplication_factor_;

        init_from_disk(old_data_node, *leaf_node, bitmaps_, items);
#if _Debug
          printf("###insert_into_leaf_node-datanode resize\n");
#endif
        old_data_node->resize(kMinDensity_Disk, false, keep_left, keep_right);
        old_data_node->num_resizes_++;
        // flush the new node on disk
        // PhyscialAddr new_addr;
        metanode_data_node.last_data_node_block = leaf_node->pre_leaf_block;
        metanode_data_node.last_data_node_offset = leaf_node->pre_leaf_offset;
        DataHeaderOnDisk r = write_data_node(old_data_node, addr, false,
                                             old_data_node->first_key(), old_data_node->last_key(), leaf_node->next_leaf_block, leaf_node->next_leaf_offset);
        delete_node(old_data_node);
        delete []bitmaps_;
        delete []items;
        memcpy(leaf_node, &r, DataHeaderOnDiskSize);
            // update its address on its parent // bug
        // re-construct the main memory parent
        if (leaf_node->duplication_factor_ > 0) {
            if (hybrid_mode == ALL_DISK) {
                TraversalPathDisk parent = traversal_path.back();
                int repeats = 1 << leaf_node->duplication_factor_;
                int start_bucketID = parent.child_offset - (parent.child_offset % repeats);
                int end_bucketID = start_bucketID + repeats;
                int start_offset = parent.addr.block * BlockSize + parent.addr.offset + ModelHeaderOnDiskSize + parent.header.empty_byte_count + start_bucketID*PhyscialAddrSize;
                // PhyscialAddr items[parent.header.children_number];
                int item_size = end_bucketID - start_bucketID;
                PhyscialAddr *items = new PhyscialAddr[item_size];
                for (int i = 0; i < item_size; i++) {
                    items[i] = *addr;
                }
                sm->write_arbitrary(start_offset, items, item_size * PhyscialAddrSize);
                delete []items;
            } else {//update the addr in main memory
                TraversalNode parent = traversal_path2.back();
                int repeats = 1 << leaf_node->duplication_factor_;
                int start_bucketID = parent.bucketID - (parent.bucketID % repeats);
                int end_bucketID = start_bucketID + repeats;
                for (int i = start_bucketID; i < end_bucketID; i++) {
                    parent.node->childrenAddrs[i] = *addr;
                }
                parent.node->childrenAddrs[parent.bucketID] = *addr;
            }

        }
        else {
            if (hybrid_mode == ALL_DISK)
                sm->write_arbitrary(p_block * BlockSize + p_offset, addr, PhyscialAddrSize);
            else { // update main memory only
                TraversalNode parent = traversal_path2.back();
                parent.node->childrenAddrs[parent.bucketID] = *addr;
            }
        }
      }
    }

#if Profiling
      auto smo_end = std::chrono::high_resolution_clock::now();
      auto smo_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
              smo_end - smo_start)
              .count();
      SMO_TIME_IN_LEAF_NODE += smo_b;
      // find the position to insert
      auto search_start = std::chrono::high_resolution_clock::now();
#endif
      double pred_pos_d = leaf_node->slope * static_cast<double>(key) + leaf_node->intercept;
      int pred_pos_i = static_cast<int>(pred_pos_d);
      pred_pos_i = std::min<int>(std::max<int>(pred_pos_i, 0), leaf_node->data_capacity - 1);
      char data2[BlockSize];
      data_node_sm->_read_block(data2, addr->block);
      int c = 0;
      int pos;
      int upper_bound_pos;
      int insertion_position;
      int last_block = addr->block;
      exp_search(leaf_node, pred_pos_i, key, addr->block, addr->offset, data2, &pos, &c, &last_block, true);
      if (pred_pos_i <= pos || check_exist(pos, addr->block, addr->offset, *leaf_node, data2, &last_block)) {
        upper_bound_pos = pos;
        insertion_position = pos;
      } else {
        //
        insertion_position = std::min(pred_pos_i, get_next_filled_position(addr->block, addr->offset, pos, true, *leaf_node) - 1);
        upper_bound_pos = pos;
      }

#if Profiling
      auto search_end = std::chrono::high_resolution_clock::now();
      auto search_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
              search_end - search_start)
              .count();
      SEARCH_TIME_IN_LEAF_NODE += search_b;

//      char *BitmapCache = new char[BlockSize];
//      int last_bitmap_block = -1;
      auto insert_start = std::chrono::high_resolution_clock::now();
#endif
      if (insertion_position < leaf_node->data_capacity && !check_exist(insertion_position, addr->block, addr->offset, *leaf_node, data2, &last_block)) {
        // insert at
        insert_element_at(key, value, insertion_position, addr->block, addr->offset, *leaf_node);
      } else {
        insert_using_shifts(key, value, insertion_position, addr->block, addr->offset, leaf_node);
      }

#if Profiling
      auto insert_end = std::chrono::high_resolution_clock::now();
      auto insert_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
              insert_end - insert_start)
              .count();
      INSERT_TIME_IN_LEAF_NODE += insert_b;
      // update stats & refresh header on disk
      auto maintain_start = std::chrono::high_resolution_clock::now();
#endif
      leaf_node->num_keys++;
      leaf_node->num_inserts_++;
      if (key > leaf_node->max_key_) {
        leaf_node->max_key_ = key;
        leaf_node->num_right_out_of_bounds_inserts_++;
      }
      if (key < leaf_node->min_key_) {
        leaf_node->min_key_ = key;
        leaf_node->num_left_out_of_bounds_inserts_++;
      }
      data_node_sm->write_arbitrary(addr->block * BlockSize + addr->offset, leaf_node, DataHeaderOnDiskSize);

#if Profiling
      auto maintain_end = std::chrono::high_resolution_clock::now();
      auto maintain_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
              maintain_end - maintain_start)
              .count();
      MAINTAIN_TIME_IN_LEAF_NODE += maintain_b;
#endif

      return {0, -1};
  }

  void set_bit(int pos, int block, int offset, DataHeaderOnDisk leaf_node) {
    // read the byte and then write
    int bitmap_pos = pos >> 6;
    int bit_pos = pos - (bitmap_pos << 6);
    long _start_offset_ = BlockSize * block + offset;
    _start_offset_ += (DataHeaderOnDiskSize + leaf_node.empty_byte_count_1 + bitmap_pos*sizeof(uint64_t));
    uint64_t _bt;
    data_node_sm->read_arbitrary(&_bt, _start_offset_, sizeof(uint64_t));
    _bt |= (1ULL << bit_pos);
    data_node_sm->write_arbitrary(_start_offset_, &_bt, sizeof(uint64_t));
  }

  // update as a block-based manner
  void set_key(T key, int pos, int block, int offset, DataHeaderOnDisk leaf_node) {
    ItemOnDisk item;
    item.key = key;
    long _start_offset_ = BlockSize * block + offset;
    _start_offset_ += (DataHeaderOnDiskSize + leaf_node.empty_byte_count_1 + leaf_node.bitmap_size*sizeof(uint64_t)) + leaf_node.empty_byte_count_2 + pos * ItemOnDiskSize;
    data_node_sm->write_arbitrary(_start_offset_, &item, ItemOnDiskSize);
  }

  void insert_element_at(T key, P value, int pos, int block, int offset, DataHeaderOnDisk leaf_node) {
    ItemOnDisk item;
    item.key = key;
    item.value = value;
    long _start_offset_ = BlockSize * block + offset;
    _start_offset_ += (DataHeaderOnDiskSize + leaf_node.empty_byte_count_1 + leaf_node.bitmap_size*sizeof(uint64_t)) + leaf_node.empty_byte_count_2 + pos * ItemOnDiskSize;
    data_node_sm->write_arbitrary(_start_offset_, &item, ItemOnDiskSize);
    set_bit(pos, block, offset, leaf_node);
    pos--;
    // overhead
    int last_block = -1;
    char DATA[BlockSize];
    while (pos >= 0 && !check_exist(pos, block, offset, leaf_node, DATA, &last_block)) {
      set_key(key, pos, block, offset, leaf_node);
      pos--;
    }
  }

  int closest_gap(int pos, int block, int offset, DataHeaderOnDisk leaf_node) {
    int max_left_offset = pos;
    int max_right_offset = leaf_node.data_capacity - pos - 1;
    int max_bidirectional_offset =
        std::min<int>(max_left_offset, max_right_offset);
    int distance = 1;
    int last_block = -1;
    char DATA[BlockSize];
    while (distance <= max_bidirectional_offset) {
      if (!check_exist(pos - distance, block, offset, leaf_node, DATA, &last_block)) {
        return pos - distance;
      }
      if (!check_exist(pos + distance, block, offset, leaf_node, DATA,&last_block)) {
        return pos + distance;
      }
      distance++;
    }
    if (max_left_offset > max_right_offset) {
      for (int i = pos - distance; i >= 0; i--) {
        if (!check_exist(i, block, offset, leaf_node, DATA,&last_block)) return i;
      }
    } else {
      for (int i = pos + distance; i < leaf_node.data_capacity; i++) {
        if (!check_exist(i, block, offset, leaf_node, DATA,&last_block)) return i;
      }
    }
    return -1;
  }
  void insert_using_shifts(T key, P value, int pos, int block ,int offset, DataHeaderOnDisk *leaf_node) {
    int gap_pos = closest_gap(pos, block, offset, *leaf_node);
    set_bit(gap_pos, block, offset, *leaf_node);
    if (gap_pos >= pos) {
      int read_size = (gap_pos-pos) * ItemOnDiskSize;
      ItemOnDisk items[gap_pos-pos];
      long start_offset = block * BlockSize + offset + DataHeaderOnDiskSize + leaf_node->empty_byte_count_1 + leaf_node->bitmap_size*sizeof(uint64_t)
        + leaf_node->empty_byte_count_2 + pos * ItemOnDiskSize;
      data_node_sm->read_arbitrary(items, start_offset, read_size);
      // items[0].key = key;
      // items[0].value = value;
      data_node_sm->write_arbitrary(start_offset + ItemOnDiskSize, items, read_size);
      insert_element_at(key, value, pos, block, offset, *leaf_node);
      //read from gap to pos
      leaf_node->num_shifts_ += gap_pos - pos;
    } else {
      long start_offset = block * BlockSize + offset + DataHeaderOnDiskSize + leaf_node->empty_byte_count_1 + leaf_node->bitmap_size*sizeof(uint64_t)
        + leaf_node->empty_byte_count_2 + (gap_pos+1)* ItemOnDiskSize;
      long read_size = (pos - gap_pos-1) * ItemOnDiskSize;
      ItemOnDisk items[pos-gap_pos-1];
      data_node_sm->read_arbitrary(items, start_offset, read_size);
      // items[pos-gap_pos-1].key = key;
      // items[pos-gap_pos-1].value = value;
      data_node_sm->write_arbitrary(start_offset - ItemOnDiskSize, items, read_size);
      insert_element_at(key, value, pos-1, block, offset, *leaf_node);
      leaf_node->num_shifts_ += pos - gap_pos - 1;
    }
  }

  int get_next_filled_position(int block, int offset, int pos, bool exclusive, DataHeaderOnDisk leaf_node) {
    if (exclusive) {
      pos++;
      if (pos == leaf_node.data_capacity) {
        return leaf_node.data_capacity;
      }
    }

    int curBitmapIdx = pos >> 6;
    uint64_t curBitmapData;
    long _start_offset_ = BlockSize * block + offset;

    char DATA[BlockSize];
    _start_offset_ += (DataHeaderOnDiskSize + leaf_node.empty_byte_count_1 + curBitmapIdx*sizeof(uint64_t));
    int _block = _start_offset_ / BlockSize;
    int _offset = _start_offset_ % BlockSize;
    int last_block  = -1;
    data_node_sm->get_block(_block, DATA);
    last_block = _block;
      curBitmapData = *((uint64_t *)(DATA + _offset));
//    data_node_sm->read_arbitrary(&curBitmapData, _start_offset_, sizeof(uint64_t));
      int bit_pos = pos - (curBitmapIdx << 6);
      curBitmapData &= ~((1ULL << (bit_pos)) - 1);
    while (curBitmapData == 0) {
      curBitmapIdx++;
      if (curBitmapIdx >= leaf_node.bitmap_size) {
        return leaf_node.data_capacity;
      }
      _start_offset_ = BlockSize * block + offset + (DataHeaderOnDiskSize + leaf_node.empty_byte_count_1 + curBitmapIdx*sizeof(uint64_t));
      _block = _start_offset_ / BlockSize;
      _offset = _start_offset_ % BlockSize;
      if (_block != last_block) {
          data_node_sm->get_block(_block, DATA);
          last_block = _block;
      }
      curBitmapData = *((uint64_t *)(DATA + _offset));
//      data_node_sm->read_arbitrary(&curBitmapData, _start_offset_, sizeof(uint64_t));
    }
    uint64_t bit = extract_rightmost_one(curBitmapData);
    return get_offset(curBitmapIdx, bit);
  }

  void expand_root_hybrid(T key, bool expand_left) {
      auto root = static_cast<model_node_type*>(root_node_);

      // Find the new bounds of the key domain.
      // Need to be careful to avoid overflows in the key type.
      T domain_size = istats_.key_domain_max_ - istats_.key_domain_min_;
      int expansion_factor;
      T new_domain_min = istats_.key_domain_min_;
      T new_domain_max = istats_.key_domain_max_;
//      data_node_type* outermost_node;
      int block_outer = -1;
      int offset_outer = -1;
      if (expand_left) {
          auto key_difference = static_cast<double>(istats_.key_domain_min_ -
                                                    std::min(key, istats_.min_key_in_left));
          expansion_factor = pow_2_round_up(static_cast<int>(
                                                    std::ceil((key_difference + domain_size) / domain_size)));
          // Check for overflow. To avoid overflow on signed types while doing
          // this check, we do comparisons using half of the relevant quantities.
          T half_expandable_domain =
                      istats_.key_domain_max_ / 2 - std::numeric_limits<T>::lowest() / 2;
          T half_expanded_domain_size = expansion_factor / 2 * domain_size;
          if (half_expandable_domain < half_expanded_domain_size) {
              new_domain_min = std::numeric_limits<T>::lowest();
          } else {
              new_domain_min = istats_.key_domain_max_;
              new_domain_min -= half_expanded_domain_size;
              new_domain_min -= half_expanded_domain_size;
          }
          istats_.num_keys_at_last_left_domain_resize = stats_.num_keys;
          istats_.num_keys_below_key_domain = 0;
          block_outer = metanode.left_most_block;
          offset_outer = metanode.left_most_offset;
      }
      else {
          auto key_difference = static_cast<double>(std::max(key, istats_.max_key_in_right) -
                                                    istats_.key_domain_max_);
          expansion_factor = pow_2_round_up(static_cast<int>(
                                                    std::ceil((key_difference + domain_size) / domain_size)));
          // Check for overflow. To avoid overflow on signed types while doing
          // this check, we do comparisons using half of the relevant quantities.
          T half_expandable_domain =
                  std::numeric_limits<T>::max() / 2 - istats_.key_domain_min_ / 2;
          T half_expanded_domain_size = expansion_factor / 2 * domain_size;
          if (half_expandable_domain < half_expanded_domain_size) {
              new_domain_max = std::numeric_limits<T>::max();
          } else {
              new_domain_max = istats_.key_domain_min_;
              new_domain_max += half_expanded_domain_size;
              new_domain_max += half_expanded_domain_size;
          }
          istats_.num_keys_at_last_right_domain_resize = stats_.num_keys;
          istats_.num_keys_above_key_domain = 0;
          block_outer = metanode.right_most_block;
          offset_outer = metanode.right_most_offset;
      }
      // Modify the root node appropriately
      int new_nodes_start;  // index of first pointer to a new node
      int new_nodes_end;    // exclusive
      if (static_cast<size_t>(root->num_children_) * expansion_factor <=
          static_cast<size_t>(derived_params_.max_fanout)) {
          // Expand root node
          stats_.num_model_node_expansions++;
          stats_.num_model_node_expansion_pointers += root->num_children_;

          int new_num_children = root->num_children_ * expansion_factor;
          auto new_children = new (pointer_allocator().allocate(new_num_children))
                  AlexNode<T, P>*[new_num_children];
          auto new_addr = new PhyscialAddr [new_num_children];;
          auto new_is_leaf = new int [new_num_children];;
          int copy_start;
          if (expand_left) {
              copy_start = new_num_children - root->num_children_;
              new_nodes_start = 0;
              new_nodes_end = copy_start;
              root->model_.b_ += new_num_children - root->num_children_;
          } else {
              copy_start = 0;
              new_nodes_start = root->num_children_;
              new_nodes_end = new_num_children;
          }
          for (int i = 0; i < root->num_children_; i++) {
              new_children[copy_start + i] = root->children_[i];
              new_addr[copy_start + i] = root->childrenAddrs[i];
              new_is_leaf[copy_start + i] = root->is_leaf[i];
          }
          pointer_allocator().deallocate(root->children_, root->num_children_);
//          pointer_allocator().deallocate(root->childrenAddrs, root->num_children_);
//          pointer_allocator().deallocate(root->is_leaf, root->num_children_);
          delete []root->childrenAddrs;
          delete []root->is_leaf;
          root->childrenAddrs = new_addr;
          root->children_ = new_children;
          root->is_leaf = new_is_leaf;
          root->num_children_ = new_num_children;
      }
      else {
          // Create new root node
          auto new_root = new (model_node_allocator().allocate(1))
                  model_node_type(static_cast<short>(root->level_ - 1), allocator_);
          new_root->model_.a_ = root->model_.a_ / root->num_children_;
          new_root->model_.b_ = root->model_.b_ / root->num_children_;
          if (expand_left) {
              new_root->model_.b_ += expansion_factor - 1;
          }
          new_root->num_children_ = expansion_factor;
          new_root->children_ = new (pointer_allocator().allocate(expansion_factor))
                  AlexNode<T, P>*[expansion_factor];
          new_root->is_leaf = new int[expansion_factor];
          new_root->childrenAddrs = new PhyscialAddr[expansion_factor];
          if (expand_left) {
              new_root->children_[expansion_factor - 1] = root;
              new_root->is_leaf[expansion_factor - 1] = 0;
              new_nodes_start = 0;
          } else {
              new_root->children_[0] = root;
              new_root->is_leaf[0] = 0;
              new_nodes_start = 1;
          }
          new_nodes_end = new_nodes_start + expansion_factor - 1;
          root_node_ = new_root;
          update_superroot_pointer();
          root = new_root;
      }
      // Determine if new nodes represent a range outside the key type's domain.
      // This happens when we're preventing overflows.
      int in_bounds_new_nodes_start = new_nodes_start;
      int in_bounds_new_nodes_end = new_nodes_end;
      if (expand_left) {
          in_bounds_new_nodes_start =
                  std::max(new_nodes_start, root->model_.predict(new_domain_min));
      } else {
          in_bounds_new_nodes_end =
                  std::min(new_nodes_end, root->model_.predict(new_domain_max) + 1);
      }
      // Fill newly created child pointers of the root node with new data nodes.
      // To minimize empty new data nodes, we create a new data node per n child
      // pointers, where n is the number of pointers to existing nodes.
      // Requires reassigning some keys from the outermost pre-existing data node
      // to the new data nodes.
      int n = root->num_children_ - (new_nodes_end - new_nodes_start);
      assert(root->num_children_ % n == 0);
      auto new_node_duplication_factor =
              static_cast<uint8_t>(log_2_round_down(n));

      // build the old node
      DataHeaderOnDisk old_dh;
      long _start_offset_ = BlockSize * block_outer + offset_outer;
      data_node_sm->read_arbitrary(&old_dh, _start_offset_, DataHeaderOnDiskSize);
      _start_offset_ += (DataHeaderOnDiskSize + old_dh.empty_byte_count_1);
      uint64_t bitmaps_[old_dh.bitmap_size];
      data_node_sm->read_arbitrary(bitmaps_, _start_offset_, old_dh.bitmap_size * sizeof(uint64_t));
      _start_offset_ += old_dh.bitmap_size * sizeof(uint64_t) + old_dh.empty_byte_count_2;
      ItemOnDisk *items = new ItemOnDisk [old_dh.data_capacity];
      data_node_sm->read_arbitrary(items, _start_offset_, old_dh.data_capacity*ItemOnDiskSize);
      auto old_data_node = new (data_node_allocator().allocate(1)) data_node_type(old_dh.level, derived_params_.max_data_node_slots, key_less_, allocator_);
      if (expand_left) {
          old_dh.duplication_factor_ = metanode.dup_left;
      }
      else {
          old_dh.duplication_factor_ = metanode.dup_right;
      }
      init_from_disk(old_data_node, old_dh, bitmaps_, items);

      if (expand_left) {
          T left_boundary_value = istats_.key_domain_min_;
          int left_boundary = old_data_node->lower_bound(left_boundary_value);//lower_bound_disk(block_outer, offset_outer, left_boundary_value); //outermost_node->lower_bound(left_boundary_value);
          int next_block = block_outer;
          int next_offset = offset_outer;

          for (int i = new_nodes_end; i > new_nodes_start; i -= n) {
              if (i <= in_bounds_new_nodes_start) {
                  // Do not initialize nodes that fall outside the key type's domain
                  break;
              }
              int right_boundary = left_boundary;
              if (i - n <= in_bounds_new_nodes_start) {
                  left_boundary = 0;
              } else {
                  left_boundary_value -= domain_size;
                  left_boundary = old_data_node->lower_bound(left_boundary_value);///lower_bound_disk(block_outer, offset_outer, left_boundary_value);/// outermost_node->lower_bound(left_boundary_value);
              }
              data_node_type* new_node = bulk_load_leaf_node_from_existing(
                      old_data_node, left_boundary, right_boundary, true);
              new_node->level_ = static_cast<short>(root->level_ + 1);
              new_node->duplication_factor_ = new_node_duplication_factor;
              PhyscialAddr _addr;
              metanode_data_node.last_data_node_block = -1;
              metanode_data_node.last_data_node_offset = -1;
              DataHeaderOnDisk new_dnh = write_data_node(new_node, &_addr, new_node->first_key(), new_node->last_key(), next_block, next_offset);
              next_block = _addr.block;
              next_offset = _addr.offset;
              delete_node(new_node);
              for (int j = i - 1; j >= i - n; j--) {
                  root->is_leaf[j] = 1;
                  root->childrenAddrs[j] = _addr;
              }
          }
      }
      else {
          T right_boundary_value = istats_.key_domain_max_;
          int right_boundary = old_data_node->lower_bound(right_boundary_value);
          int pre_block = block_outer;
          int pre_offset = offset_outer;
          for (int i = new_nodes_start; i < new_nodes_end; i += n) {
              if (i >= in_bounds_new_nodes_end) {
                  // Do not initialize nodes that fall outside the key type's domain
                  break;
              }
              int left_boundary = right_boundary;
              if (i + n >= in_bounds_new_nodes_end) {
                  right_boundary = old_data_node->data_capacity_;
              } else {
                  right_boundary_value += domain_size;
                  right_boundary = old_data_node->lower_bound(right_boundary_value);
              }
              data_node_type* new_node = bulk_load_leaf_node_from_existing(
                      old_data_node, left_boundary, right_boundary, true);
              new_node->level_ = static_cast<short>(root->level_ + 1);
              new_node->duplication_factor_ = new_node_duplication_factor;
              PhyscialAddr _addr;
              metanode_data_node.last_data_node_block = pre_block;
              metanode_data_node.last_data_node_offset = pre_offset;
              DataHeaderOnDisk new_dnh = write_data_node(new_node, &_addr, false, new_node->first_key(), new_node->last_key());
              // memcpy(&cur_dh, &new_dnh, DataHeaderOnDiskSize);
              delete_node(new_node);
              pre_block = _addr.block;
              pre_offset = _addr.offset;
              for (int j = i; j < i + n; j++) {
                  root->is_leaf[j] = 1;
                  root->childrenAddrs[j] = _addr;
              }
          }
      }

      if (expand_left) {
          old_data_node->erase_range(new_domain_min, istats_.key_domain_min_);
          PhyscialAddr last_new_leaf =  root->childrenAddrs[new_nodes_end - 1];
          metanode_data_node.last_data_node_block = last_new_leaf.block;
          metanode_data_node.last_data_node_offset = last_new_leaf.offset;
          PhyscialAddr old_addr;
          old_addr.block = block_outer;
          old_addr.offset = offset_outer;
          // we update the its pre_node's next node in write_data_node function
          write_data_node(old_data_node, &old_addr, true, old_data_node->first_key(), old_data_node->last_key(), old_dh.next_leaf_block, old_dh.next_leaf_offset);
      }
      else {
          old_data_node->erase_range(istats_.key_domain_max_, new_domain_max,
                                     true);
          PhyscialAddr old_addr;
          old_addr.block = block_outer;
          old_addr.offset = offset_outer;

          PhyscialAddr first_new_leaf = root->childrenAddrs[new_nodes_start];
          metanode_data_node.last_data_node_block = old_dh.pre_leaf_block;
          metanode_data_node.last_data_node_offset = old_dh.pre_leaf_offset;
          write_data_node(old_data_node, &old_addr, true, old_data_node->first_key(), old_data_node->last_key(), first_new_leaf.block, first_new_leaf.offset);
      }
      delete []items;
      metanode.is_leaf = ALEXModelType;
      // flush metanode
      sync_metanode();
      istats_.key_domain_min_ = new_domain_min;
      istats_.key_domain_max_ = new_domain_max;
      delete_node(old_data_node);

  }

  PhyscialAddr get_leaf_addr_hybrid(int *r_block_count, char *block_data2, DataHeaderOnDisk *LEAF_NODE,
          T key, std::vector<TraversalNode>* traversal_path = nullptr) {
      PhyscialAddr addr;
      if (metanode.is_leaf == ALEXDataType) {
          addr.block = metanode.root_block_id;
          addr.offset = metanode.root_offset;
          addr.duplication_factor_ = metanode.dup_root;
          addr.flag = metanode.is_leaf;
          return addr;
      }
      AlexNode<T, P>* cur = root_node_;
      char *block_data = new char[BlockSize];

      while (true) {
          auto node = static_cast<model_node_type*>(cur);
          double bucketID_prediction = node->model_.predict_double(key);
          int bucketID = static_cast<int>(bucketID_prediction);
          bucketID =
                  std::min<int>(std::max<int>(bucketID, 0), node->num_children_ - 1);
          if (traversal_path) {
              traversal_path->push_back({node, bucketID});
          }
          if (node->is_leaf[bucketID] == 0) {cur = node->children_[bucketID]; continue;}
          PhyscialAddr child_addr = node->childrenAddrs[bucketID];
          // meet leaf node
          data_node_sm->get_block(child_addr.block, block_data2);
          *r_block_count += 1;
          DataHeaderOnDisk *data_node = (DataHeaderOnDisk *)( block_data2 + child_addr.offset);

          int pred_pos_i_rounded = static_cast<int>(bucketID + 0.5);
          double tolerance = 10 * std::numeric_limits<double>::epsilon() * bucketID_prediction;
          if (std::abs(bucketID_prediction - pred_pos_i_rounded) <= tolerance) {
              if (pred_pos_i_rounded <= bucketID_prediction) {
                  if (data_node->pre_leaf_block != -1) {
                      // char pre_block_data[BLOCK_SIZE];
                      data_node_sm->get_block(data_node->pre_leaf_block, block_data);
                      *r_block_count += 1;
                      DataHeaderOnDisk *pre_data_node = (DataHeaderOnDisk *) (block_data + data_node->pre_leaf_offset);
                      if (pre_data_node->last_key >= key) {
                          memcpy(LEAF_NODE, pre_data_node, DataHeaderOnDiskSize);
                          if (traversal_path) correct_traversal_path(nullptr, *traversal_path, true, LEAF_NODE->duplication_factor_);
                          addr.block = data_node->pre_leaf_block;
                          addr.offset = data_node->pre_leaf_offset;
                          addr.duplication_factor_ = pre_data_node->duplication_factor_;
                          addr.flag = ALEXDataType;
                          memcpy(block_data2, block_data, BlockSize);
                          delete []block_data;
                          return addr;
                      }
                  }
              } else {
                  if (data_node->next_leaf_block != -1) {
                      data_node_sm->get_block(data_node->next_leaf_block, block_data);
                      *r_block_count += 1;
                      DataHeaderOnDisk *next_data_node = (DataHeaderOnDisk *)(block_data + data_node->next_leaf_offset);
                      // memcpy(&next_data_node, block_data + data_node.next_leaf_offset, data_meta_size);
                      if (next_data_node->first_key <= key) {
                          memcpy(LEAF_NODE, next_data_node, DataHeaderOnDiskSize);
                          if (traversal_path) correct_traversal_path(nullptr, *traversal_path, false, LEAF_NODE->duplication_factor_);
                          addr.block = data_node->next_leaf_block;
                          addr.offset = data_node->next_leaf_offset;
                          addr.duplication_factor_ = next_data_node->duplication_factor_;
                          addr.flag = ALEXDataType;
                          memcpy(block_data2, block_data, BlockSize);
                          delete []block_data;
                          return addr;
                      }
                  }
              }
          }
          addr.block = child_addr.block;
          addr.offset = child_addr.offset;
          addr.duplication_factor_ = child_addr.duplication_factor_;
          addr.flag = ALEXDataType;
          memcpy(LEAF_NODE, data_node, DataHeaderOnDiskSize);
          delete []block_data;
          return addr;

      }
  }

  void insert_disk_hybrid(T key, P value) {
      if (key > istats_.key_domain_max_) {
          istats_.num_keys_above_key_domain++;
          istats_.max_key_in_right = std::max(istats_.max_key_in_right, key);
          if (should_expand_right()) {
              printf("expand_root_disk-right\n");
              expand_root_hybrid(key, false);  // expand to the right
          }
      } else if (key < istats_.key_domain_min_) {
          istats_.num_keys_below_key_domain++;
          istats_.min_key_in_left = std::min(istats_.min_key_in_left, key);
          if (should_expand_left()) {
              printf("expand_root_disk-left\n");
              expand_root_hybrid(key, true);  // expand to the left
          }
      }

      int count = 0;
      char *block_data = new char[BlockSize];
      DataHeaderOnDisk LEAF_NODE;
      std::vector<TraversalNode> traversal_path;
      PhyscialAddr leaf_addr = get_leaf_addr_hybrid(&count, block_data, &LEAF_NODE, key, &traversal_path);
      std::vector<TraversalPathDisk> _tp;
      std::pair<int, int> ret = insert_into_leaf_node(key, value, block_data, &LEAF_NODE, &leaf_addr, 0, 0, _tp, traversal_path);

      int fail = ret.first;
      if (fail == -1) { printf("repeated\n"); delete []block_data; return;}

      bool is_first = true;
      if (fail) {
          model_node_type* parent = traversal_path.back().node;
          TraversalNode parent_path = traversal_path.back();
          while (fail) {
              if (parent == superroot_) {
                  update_superroot_key_domain();
              }
              int bucketID = parent->model_.predict(key);
              bucketID = std::min<int>(std::max<int>(bucketID, 0),
                                       parent->num_children_ - 1);
              std::vector<fanout_tree::FTNode> used_fanout_tree_nodes;
              int _pos;
              leaf_addr = parent->get_child_node(key, &_pos);
              DataHeaderOnDisk dh;
              long _start_offset_ = BlockSize * leaf_addr.block + leaf_addr.offset;
              data_node_sm->read_arbitrary(&dh, _start_offset_, DataHeaderOnDiskSize);
              _start_offset_ += (DataHeaderOnDiskSize + dh.empty_byte_count_1);
              uint64_t *bitmaps_ = new uint64_t [dh.bitmap_size];
              data_node_sm->read_arbitrary(bitmaps_, _start_offset_, dh.bitmap_size * sizeof(uint64_t));
              _start_offset_ += dh.bitmap_size * sizeof(uint64_t) + dh.empty_byte_count_2;
              ItemOnDisk *items2 = new ItemOnDisk[dh.data_capacity];
              data_node_sm->read_arbitrary(items2, _start_offset_, dh.data_capacity*ItemOnDiskSize);
              auto node = new (data_node_allocator().allocate(1)) data_node_type(dh.level, derived_params_.max_data_node_slots, key_less_, allocator_);
              dh.duplication_factor_ = leaf_addr.duplication_factor_;
              init_from_disk(node, dh, bitmaps_, items2);
              delete []items2;
              delete []bitmaps_;

              int fanout_tree_depth = 1;
              if (experimental_params_.splitting_policy_method == 0 || fail >= 2) {
                  // always split in 2. No extra work required here
              }
              else if (experimental_params_.splitting_policy_method == 1) {
                  fanout_tree_depth = find_best_fanout_existing_node(parent, {}, {}, {},
                                                                     bucketID, stats_.num_keys, used_fanout_tree_nodes, 2, node);
              }
              else if (experimental_params_.splitting_policy_method == 2) {
                  fanout_tree_depth = find_best_fanout_existing_node(parent, {}, {}, {},
                                                                     bucketID, stats_.num_keys, used_fanout_tree_nodes, derived_params_.max_fanout, node);
              }
              int best_fanout = 1 << fanout_tree_depth;
              if (fanout_tree_depth == 0) {
#if _Debug
                  printf("###insert_into_leaf_node-datanode resize fanout 0\n");
#endif
                  node->resize(data_node_type::kMinDensity_, true,
                               node->is_append_mostly_right(),
                               node->is_append_mostly_left());
                  fanout_tree::FTNode& tree_node = used_fanout_tree_nodes[0];
                  node->cost_ = tree_node.cost;
                  node->expected_avg_exp_search_iterations_ =
                          tree_node.expected_avg_search_iterations;
                  node->expected_avg_shifts_ = tree_node.expected_avg_shifts;
                  node->reset_stats();
                  stats_.num_expand_and_retrains++;
                  // TODO: update the root node as the suc case
                  PhyscialAddr l_addr;
                  metanode_data_node.last_data_node_block = dh.pre_leaf_block;
                  metanode_data_node.last_data_node_offset = dh.pre_leaf_offset;
                  DataHeaderOnDisk r =  write_data_node(node, &l_addr, false, node->first_key(), node->last_key(), dh.next_leaf_block, dh.next_leaf_offset);
                  if (node->duplication_factor_ > 0) {
//                      TraversalPathDisk parent = traversal_path.back();
                      int repeats = 1 << node->duplication_factor_;
                      int start_bucketID = parent_path.bucketID - (parent_path.bucketID % repeats);
                      int end_bucketID = start_bucketID + repeats;
                      for (int i = start_bucketID; i < end_bucketID; i++) {
                          parent_path.node->childrenAddrs[i] = l_addr;
                      }
                  }
                  parent_path.node->childrenAddrs[parent_path.bucketID] = l_addr;
                  memcpy(&LEAF_NODE, &r, DataHeaderOnDiskSize);
                  std::vector<TraversalPathDisk> _tp;
                  std::pair<int, int> ret = insert_into_leaf_node(key, value, block_data, &LEAF_NODE, &l_addr, 0, 0, _tp, traversal_path);
                  fail = ret.first;
                  if (fail == -1) {
                      delete []block_data;
                      printf("repeated!!\n");
                      return;
                  }
              }
              else {

                  bool reuse_model = (fail == 3);
                  // either split sideways or downwards
                  bool should_split_downwards =
                          (parent->num_children_ * best_fanout /
                           (1 << node->duplication_factor_) >
                           derived_params_.max_fanout ||
                           parent->level_ == superroot_->level_);
                  int empty_size = 0;
                  if (should_split_downwards) {
#if _Debug
                      printf("###split_down\n");
#endif
                      parent = split_downwards_disk(parent, bucketID, fanout_tree_depth,
                                                        used_fanout_tree_nodes, node, dh, nullptr, reuse_model, &empty_size);
                  } else {
                      // FIXME: pp_block
#if _Debug
                      printf("###split_sideway\n");
#endif
                      split_sideways_hybrid_disk(parent, bucketID, fanout_tree_depth, dh, used_fanout_tree_nodes, reuse_model, node);
                  }

                  int pos = 0;
                  leaf_addr = parent->get_child_node(key, &pos);
                  data_node_sm->read_arbitrary(&LEAF_NODE, leaf_addr.block*BlockSize + leaf_addr.offset, DataHeaderOnDiskSize);
                  // the traversal_path here is tale while we do not need it.
                  std::vector<TraversalPathDisk> _tp;
                  std::pair<int, int> ret = insert_into_leaf_node(key, value, block_data, &LEAF_NODE,
                                                                  &leaf_addr, 0, 0,
                                                                  _tp, traversal_path);
                  fail = ret.first;
                  if (fail == -1) {
                      delete []block_data;
                      printf("repeated!!\n");
                      return;
                  }
              }
              delete_node(node);
          }
      }
      std::vector<TraversalNode>().swap(traversal_path);
      delete []block_data;
      stats_.num_inserts++;
      stats_.num_keys++;
  }

  void insert_disk(T key, P value, long long *search_time, long long *insert_time,
                   long long *smo_time, long long *maintain_time, int *smo_count) {
#if Profiling
      *search_time = 0, *insert_time = 0, *smo_time = 0; *maintain_time = 0;
      *smo_count = 0; SMO_COUNT = 0;
#endif

      if (hybrid_mode == ALL_DISK) insert_disk_all(key, value, search_time, insert_time, smo_time, maintain_time, smo_count);
      else if (hybrid_mode == LEAF_DISK) insert_disk_hybrid(key, value);
      else {throw std::invalid_argument("not support in insert disk...");}

#if Profiling
      *smo_count = SMO_COUNT;
#endif

      return;
  }
  void insert_disk_all(T key, P value, long long* search_time, long long *insert_time,
                       long long *smo_time, long long *maintain_time, int *smo_count) {
#if Profiling
    auto smo_start = std::chrono::high_resolution_clock::now();
#endif
    if (key > istats_.key_domain_max_) {
        istats_.max_key_in_right = std::max(istats_.max_key_in_right, key);
      istats_.num_keys_above_key_domain++;
      if (should_expand_right_disk()) {
#if _Debug
          printf("###expand_root_disk-right\n");
#endif
        expand_root_disk(key, false);

#if Profiling
          SMO_COUNT += 1;
#endif
      }
    }else if (key < istats_.key_domain_min_) {
      istats_.num_keys_below_key_domain++;
        istats_.min_key_in_left = std::min(istats_.min_key_in_left, key);
      if (should_expand_left()) {
#if _Debug
          printf("###expand_root_disk-left\n");
#endif
        expand_root_disk(key, true);  // expand to the left
#if Profiling
          SMO_COUNT += 1;
#endif
      }
    }

#if Profiling
    auto smo_end = std::chrono::high_resolution_clock::now();
    auto smo_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
              smo_end - smo_start)
              .count();
    *smo_time += smo_b;

    auto search_start = std::chrono::high_resolution_clock::now();
#endif

    char block_data[BlockSize];
    DataHeaderOnDisk LEAF_NODE;
    int count = 0;
    int p_block;
    int p_offset;
    std::vector<TraversalPathDisk> traversal_path;
    std::vector<TraversalNode> traversal_path2;
    PhyscialAddr leaf_addr = get_leaf_addr(key, block_data, &LEAF_NODE, &count, &p_block, &p_offset, &traversal_path);

#if Profiling
    auto search_end = std::chrono::high_resolution_clock::now();
    auto search_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
              search_end - search_start)
              .count();
    *search_time += search_b;
#endif

    std::pair<int, int> ret = insert_into_leaf_node(key, value, block_data, &LEAF_NODE, &leaf_addr, p_block, p_offset, traversal_path,
                                                    traversal_path2);

#if Profiling
    *insert_time += INSERT_TIME_IN_LEAF_NODE;
    *search_time += SEARCH_TIME_IN_LEAF_NODE;
    *smo_time += SMO_TIME_IN_LEAF_NODE;
    *maintain_time += MAINTAIN_TIME_IN_LEAF_NODE;
#endif

    int fail = ret.first;
    if (fail == -1) {
      printf("repeated!!\n");
      return;
    }

    bool is_first = true;
    while (fail) {
#if Profiling
        search_start = std::chrono::high_resolution_clock::now();
#endif
        if (!is_first) {
            leaf_addr = get_leaf_addr(key, block_data, &LEAF_NODE, &count, &p_block, &p_offset, &traversal_path);
        } else {
            is_first = false;
        }
      TraversalPathDisk parent = traversal_path.back();
      if (parent.child_offset == -1) { // supperroot
        printf("child_offset is -1\n");
      }
      // int pred_pos = parent.header.slope
      double pred_pos_d = parent.header.slope * static_cast<double>(key) +  parent.header.intercept;
      int pred_pos_i = static_cast<int>(pred_pos_d);
      pred_pos_i = std::min<int>(std::max<int>(pred_pos_i, 0),  parent.header.children_number - 1);

#if Profiling
      search_end = std::chrono::high_resolution_clock::now();
      search_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                search_end - search_start)
                .count();
      *search_time += search_b;
#endif

      std::vector<fanout_tree::FTNode> used_fanout_tree_nodes;


      // Hai LAN this can be optimized to reduce one read...
#if Profiling
      SMO_COUNT += 1;
      smo_start = std::chrono::high_resolution_clock::now();
#endif

      DataHeaderOnDisk dh;
      long _start_offset_ = BlockSize * leaf_addr.block + leaf_addr.offset;
      data_node_sm->read_arbitrary(&dh, _start_offset_, DataHeaderOnDiskSize);
      _start_offset_ += (DataHeaderOnDiskSize + dh.empty_byte_count_1);
      uint64_t *bitmaps_ = new uint64_t [dh.bitmap_size];
      data_node_sm->read_arbitrary(bitmaps_, _start_offset_, dh.bitmap_size * sizeof(uint64_t));
      _start_offset_ += dh.bitmap_size * sizeof(uint64_t) + dh.empty_byte_count_2;
      ItemOnDisk *items2 = new ItemOnDisk[dh.data_capacity];
      data_node_sm->read_arbitrary(items2, _start_offset_, dh.data_capacity*ItemOnDiskSize);
      auto node = new (data_node_allocator().allocate(1)) data_node_type(dh.level, derived_params_.max_data_node_slots, key_less_, allocator_);
      dh.duplication_factor_ = leaf_addr.duplication_factor_;
      init_from_disk(node, dh, bitmaps_, items2);
      delete []items2;
      delete []bitmaps_;

      int fanout_tree_depth = 1;
      if (experimental_params_.splitting_policy_method == 0 || fail >= 2) {
          // always split in 2. No extra work required here
      }
      else if (experimental_params_.splitting_policy_method == 1) {
          fanout_tree_depth = find_best_fanout_existing_node(nullptr, parent.header, parent.addr, leaf_addr,
                                                            pred_pos_i, stats_.num_keys, used_fanout_tree_nodes, 2, node);
      }
      else if (experimental_params_.splitting_policy_method == 2) {
          fanout_tree_depth = find_best_fanout_existing_node(nullptr, parent.header, parent.addr, leaf_addr,
                                                            pred_pos_i, stats_.num_keys, used_fanout_tree_nodes, derived_params_.max_fanout, node);
      }
      int best_fanout = 1 << fanout_tree_depth;
      if (fanout_tree_depth == 0) {
          // expand existing data node and retrain model
#if _Debug
          printf("###insert_into_leaf_node-datanode resize fanout 0\n");
#endif
          node->resize(data_node_type::kMinDensity_, true,
                       node->is_append_mostly_right(),
                       node->is_append_mostly_left());
          fanout_tree::FTNode& tree_node = used_fanout_tree_nodes[0];
          node->cost_ = tree_node.cost;
          node->expected_avg_exp_search_iterations_ =
              tree_node.expected_avg_search_iterations;
          node->expected_avg_shifts_ = tree_node.expected_avg_shifts;
          node->reset_stats();
          stats_.num_expand_and_retrains++;
          PhyscialAddr l_addr;
          metanode_data_node.last_data_node_block = dh.pre_leaf_block;
          metanode_data_node.last_data_node_offset = dh.pre_leaf_offset;
          DataHeaderOnDisk r =  write_data_node(node, &l_addr, false, node->first_key(), node->last_key(), dh.next_leaf_block, dh.next_leaf_offset);
          if (node->duplication_factor_ > 0) {
              TraversalPathDisk parent = traversal_path.back();
              int repeats = 1 << node->duplication_factor_;
              int start_bucketID = parent.child_offset - (parent.child_offset % repeats);
              int end_bucketID = start_bucketID + repeats;
              int start_offset = parent.addr.block * BlockSize + parent.addr.offset + ModelHeaderOnDiskSize + parent.header.empty_byte_count + start_bucketID*PhyscialAddrSize;
              // PhyscialAddr items[parent.header.children_number];
              int item_size = end_bucketID - start_bucketID;
              PhyscialAddr *items = new PhyscialAddr[item_size];
              for (int i = 0; i < item_size; i++) {
                  items[i] = l_addr;
              }
              sm->write_arbitrary(start_offset, items, item_size * PhyscialAddrSize);
              delete []items;
          }
          else sm->write_arbitrary(p_block * BlockSize + p_offset, &l_addr, PhyscialAddrSize);

#if Profiling
          smo_end = std::chrono::high_resolution_clock::now();
          smo_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                  smo_end - smo_start)
                  .count();
          *smo_time += smo_b;
#endif
          memcpy(&LEAF_NODE, &r, DataHeaderOnDiskSize);
          std::vector<TraversalNode> traversal_path2;

          std::pair<int, int> ret = insert_into_leaf_node(key, value, block_data, &LEAF_NODE, &l_addr, p_block, p_offset, traversal_path,
                                                          traversal_path2);

#if Profiling
        *insert_time += INSERT_TIME_IN_LEAF_NODE;
        *search_time += SEARCH_TIME_IN_LEAF_NODE;
        *smo_time += SEARCH_TIME_IN_LEAF_NODE;
        *maintain_time += MAINTAIN_TIME_IN_LEAF_NODE;
#endif
          fail = ret.first;

          if (fail == -1) {
            printf("repeated!!\n");
            return;
          }
      }
      else {
          /// build the parent in main memory data structure/////
          // ModelHeaderOnDisk ph;
          int start_offset = parent.addr.block * BlockSize + parent.addr.offset + ModelHeaderOnDiskSize + parent.header.empty_byte_count;
          // PhyscialAddr items[parent.header.children_number];
          PhyscialAddr *items = new PhyscialAddr[parent.header.children_number];
          sm->read_arbitrary(items, start_offset, parent.header.children_number * PhyscialAddrSize);
          auto model_node = new (model_node_allocator().allocate(1))
          model_node_type(parent.header.level, allocator_);
          parent.header.duplication_factor_ = parent.addr.duplication_factor_;
          init_model_node_from_disk(model_node, parent.header, items);

          bool reuse_model = (fail == 3);
          // either split sideways or downwards
          bool should_split_downwards =
              (model_node->num_children_ * best_fanout /
                       (1 << node->duplication_factor_) >
                   derived_params_.max_fanout ||
              //FIXME: check superroot....
               model_node->level_ == -1);
          int empty_size = 0;
          if (should_split_downwards) {
#if _Debug
              printf("###split_down\n");
#endif
            model_node = split_downwards_disk(model_node, pred_pos_i, fanout_tree_depth,
                                     used_fanout_tree_nodes, node, dh, &parent.addr, reuse_model, &empty_size);
          } else {
#if _Debug
              printf("###split_sideway\n");
#endif
            split_sideways_disk(model_node, pred_pos_i, fanout_tree_depth, dh,
                           used_fanout_tree_nodes, reuse_model, node, parent.addr_in_p.block, parent.addr_in_p.offset, &parent.addr, &empty_size, traversal_path);
          }

#if Profiling
          smo_end = std::chrono::high_resolution_clock::now(); // align with if branch
          smo_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                  smo_end - smo_start)
                  .count();
          *smo_time += smo_b;

          search_start = std::chrono::high_resolution_clock::now();
#endif
          long start = parent.addr.block * BlockSize + parent.addr.offset + ModelHeaderOnDiskSize + empty_size;
          int pos = 0;
          PhyscialAddr leaf_addr = model_node->get_child_node(key, &pos);
          start += pos * PhyscialAddrSize;
          p_block = start / BlockSize;
          p_offset = start % BlockSize;
          data_node_sm->read_arbitrary(&LEAF_NODE, leaf_addr.block*BlockSize + leaf_addr.offset, DataHeaderOnDiskSize);
          // the traversal_path here is tale while we do not need it.
          std::vector<TraversalNode> traversal_path2;

#if Profiling
          search_end = std::chrono::high_resolution_clock::now();
          search_b = std::chrono::duration_cast<std::chrono::nanoseconds>(
                  search_end - search_start)
                  .count();
          *search_time += search_b;
#endif
          std::pair<int, int> ret = insert_into_leaf_node(key, value, block_data, &LEAF_NODE, &leaf_addr, p_block, p_offset, traversal_path,
                                                          traversal_path2);

#if Profiling
          *insert_time += INSERT_TIME_IN_LEAF_NODE;
          *search_time += SEARCH_TIME_IN_LEAF_NODE;
          *smo_time += SEARCH_TIME_IN_LEAF_NODE;
          *maintain_time += MAINTAIN_TIME_IN_LEAF_NODE;
#endif
          fail = ret.first;
          delete []items;
          if (fail == -1) {
            printf("repeated!!\n");
            return;
          }
          delete []model_node->childrenAddrs;
          delete_node(model_node);
      }
      delete_node(node);
    } // end if (fail)
    std::vector<TraversalPathDisk>().swap(traversal_path);
    stats_.num_inserts++;
    stats_.num_keys++;
  }

  // With the bad implementation before, we implement fanout function here, we'd better set the Node as template
  int find_best_fanout_existing_node(const AlexModelNode<T, P> *parent, ModelHeaderOnDisk header, PhyscialAddr parent_addr,
                                     PhyscialAddr child_addr, int pos, int total_keys,
                                    std::vector<fanout_tree::FTNode> &used_fanout_tree_nodes, int max_fanout, data_node_type *node) {
    // int start_offset = addr.block * BlockSize + addr.offset + ModelHeaderOnDiskSize
    //                   + header.empty_byte_count + pos * PhyscialAddrSize;

    int num_keys = node->num_keys_;
    int repeats = 1 << node->duplication_factor_;
    int best_level = 0;
    double best_cost = std::numeric_limits<double>::max();
    std::vector<double> fanout_costs;
    std::vector<std::vector<fanout_tree::FTNode>> fanout_tree;

    // int repeats = 1 << dh.duplication_factor_;

    int start_pos = pos - (pos % repeats);
    int end_pos = start_pos + repeats;

    double left_boundary_value = (start_pos - (hybrid_mode == ALL_DISK ? header.intercept: parent->model_.b_))
            / (hybrid_mode == ALL_DISK ? header.slope : parent->model_.a_);
    double right_boundary_value = (end_pos - (hybrid_mode == ALL_DISK ? header.intercept: parent->model_.b_))
            / (hybrid_mode == ALL_DISK ? header.slope : parent->model_.a_);
    LinearModel<T> base_model;
    base_model.a_ = 1.0 / (right_boundary_value - left_boundary_value);
    base_model.b_ = -1.0 * base_model.a_ * left_boundary_value;

    for (int fanout = 1, fanout_tree_level = 0; fanout <= max_fanout;
          fanout *= 2, fanout_tree_level++) {
      std::vector<fanout_tree::FTNode> new_level;
      double cost = 0.0;
      double a = base_model.a_ * fanout;
      double b = base_model.b_ * fanout;
      int left_boundary = 0;
      int right_boundary = 0;
      for (int i = 0; i < fanout; i++) {
        left_boundary = right_boundary;
        right_boundary = i == fanout - 1 ? node->data_capacity_
                                        : node->lower_bound(((i + 1) - b) / a);
        if (left_boundary == right_boundary) {
          new_level.push_back({fanout_tree_level, i, 0, left_boundary,
                              right_boundary, false, 0, 0, 0, 0, 0});
          continue;
        }
        int num_actual_keys = 0;
        LinearModel<T> model;
        typename AlexDataNode<T, P>::const_iterator_type it(node, left_boundary);
        LinearModelBuilder<T> builder(&model);
        for (int j = 0; it.cur_idx_ < right_boundary && !it.is_end(); it++, j++) {
          builder.add(it.key(), j);
          num_actual_keys++;
        }
        builder.build();

        double empirical_insert_frac = node->frac_inserts();
        DataNodeStats stats;
        double node_cost =
            AlexDataNode<T, P>::compute_expected_cost_from_existing(
                node, left_boundary, right_boundary,
                AlexDataNode<T, P>::kInitDensity_, empirical_insert_frac, &model,
                &stats);

        cost += node_cost * num_actual_keys / num_keys;

        new_level.push_back({fanout_tree_level, i, node_cost, left_boundary,
                            right_boundary, false, stats.num_search_iterations,
                            stats.num_shifts, model.a_, model.b_,
                            num_actual_keys});
      }
      // model weight reflects that it has global effect, not local effect
      double traversal_cost =
          kNodeLookupsWeight +
          (kModelSizeWeight * fanout *
          (sizeof(AlexDataNode<T, P>) + sizeof(void*)) * total_keys / num_keys);
      cost += traversal_cost;
      fanout_costs.push_back(cost);
      // stop after expanding fanout increases cost twice in a row
      if (fanout_costs.size() >= 3 &&
          fanout_costs[fanout_costs.size() - 1] >
              fanout_costs[fanout_costs.size() - 2] &&
          fanout_costs[fanout_costs.size() - 2] >
              fanout_costs[fanout_costs.size() - 3]) {
      break;
      }
      if (cost < best_cost) {
        best_cost = cost;
        best_level = fanout_tree_level;
      }
      fanout_tree.push_back(new_level);
    }
    for (fanout_tree::FTNode& tree_node : fanout_tree[best_level]) {
      tree_node.use = true;
    }
    fanout_tree::merge_nodes_upwards<T, P>(best_level, best_cost, num_keys, total_keys,
                            fanout_tree);
    fanout_tree::collect_used_nodes(fanout_tree, best_level, used_fanout_tree_nodes);
    return best_level;
    // return 1;
  }


  // This will NOT do an update of an existing key.
  // To perform an update or read-modify-write, do a lookup and modify the
  // payload's value.
  // Returns iterator to inserted element, and whether the insert happened or
  // not.
  // Insert does not happen if duplicates are not allowed and duplicate is
  // found.
  std::pair<Iterator, bool> insert(const T& key, const P& payload) {
    // If enough keys fall outside the key domain, expand the root to expand the
    // key domain
    if (key > istats_.key_domain_max_) {
      istats_.num_keys_above_key_domain++;
      if (should_expand_right()) {
          printf("expand_root_disk-right\n");
        expand_root(key, false);  // expand to the right
      }
    } else if (key < istats_.key_domain_min_) {
      istats_.num_keys_below_key_domain++;
      if (should_expand_left()) {
          printf("expand_root_disk-left\n");
        expand_root(key, true);  // expand to the left
      }
    }

    data_node_type* leaf = get_leaf(key);

    // Nonzero fail flag means that the insert did not happen
    std::pair<int, int> ret = leaf->insert(key, payload);
    int fail = ret.first;
    int insert_pos = ret.second;
    if (fail == -1) {
      // Duplicate found and duplicates not allowed
        printf("repeated!!\n");
      return {Iterator(leaf, insert_pos), false};
    }

    // If no insert, figure out what to do with the data node to decrease the
    // cost
    if (fail) {
      std::vector<TraversalNode> traversal_path;
      get_leaf(key, &traversal_path);
      model_node_type* parent = traversal_path.back().node;

      while (fail) {
        auto start_time = std::chrono::high_resolution_clock::now();
        stats_.num_expand_and_scales += leaf->num_resizes_;

        if (parent == superroot_) {
          update_superroot_key_domain();
        }
        int bucketID = parent->model_.predict(key);
        bucketID = std::min<int>(std::max<int>(bucketID, 0),
                                 parent->num_children_ - 1);
        std::vector<fanout_tree::FTNode> used_fanout_tree_nodes;

        int fanout_tree_depth = 1;
        if (experimental_params_.splitting_policy_method == 0 || fail >= 2) {
          // always split in 2. No extra work required here
        } else if (experimental_params_.splitting_policy_method == 1) {
          // decide between no split (i.e., expand and retrain) or splitting in
          // 2
          fanout_tree_depth = fanout_tree::find_best_fanout_existing_node<T, P>(
              parent, bucketID, stats_.num_keys, used_fanout_tree_nodes, 2);
        } else if (experimental_params_.splitting_policy_method == 2) {
          // use full fanout tree to decide fanout
          fanout_tree_depth = fanout_tree::find_best_fanout_existing_node<T, P>(
              parent, bucketID, stats_.num_keys, used_fanout_tree_nodes,
              derived_params_.max_fanout);
        }
        int best_fanout = 1 << fanout_tree_depth;
        stats_.cost_computation_time +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - start_time)
                .count();

        if (fanout_tree_depth == 0) {
          // expand existing data node and retrain model
#if _Debug
            printf("###insert_into_leaf_node-datanode resize fanout 0\n");
#endif
          leaf->resize(data_node_type::kMinDensity_, true,
                       leaf->is_append_mostly_right(),
                       leaf->is_append_mostly_left());
          fanout_tree::FTNode& tree_node = used_fanout_tree_nodes[0];
          leaf->cost_ = tree_node.cost;
          leaf->expected_avg_exp_search_iterations_ =
              tree_node.expected_avg_search_iterations;
          leaf->expected_avg_shifts_ = tree_node.expected_avg_shifts;
          leaf->reset_stats();
          stats_.num_expand_and_retrains++;
        } else {
          // split data node: always try to split sideways/upwards, only split
          // downwards if necessary
          bool reuse_model = (fail == 3);
          if (experimental_params_.allow_splitting_upwards) {
            // allow splitting upwards
            assert(experimental_params_.splitting_policy_method != 2);
            int stop_propagation_level = best_split_propagation(traversal_path);
            if (stop_propagation_level <= superroot_->level_) {
              parent = split_downwards(parent, bucketID, fanout_tree_depth,
                                       used_fanout_tree_nodes, reuse_model);
            } else {
              split_upwards(key, stop_propagation_level, traversal_path,
                            reuse_model, &parent);
            }
          } else {
            // either split sideways or downwards
            bool should_split_downwards =
                (parent->num_children_ * best_fanout /
                         (1 << leaf->duplication_factor_) >
                     derived_params_.max_fanout ||
                 parent->level_ == superroot_->level_);
            if (should_split_downwards) {
#if _Debug
                printf("###split_down\n");
#endif
              parent = split_downwards(parent, bucketID, fanout_tree_depth,
                                       used_fanout_tree_nodes, reuse_model);
            } else {
#if _Debug
                printf("###split_sideway\n");
#endif
              split_sideways(parent, bucketID, fanout_tree_depth,
                             used_fanout_tree_nodes, reuse_model);
            }
          }
          leaf = static_cast<data_node_type*>(parent->get_child_node(key));
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = end_time - start_time;
        stats_.splitting_time +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(duration)
                .count();

        // Try again to insert the key
        ret = leaf->insert(key, payload);
        fail = ret.first;
        insert_pos = ret.second;
        if (fail == -1) {
          // Duplicate found and duplicates not allowed
          return {Iterator(leaf, insert_pos), false};
        }
      }
    }
    stats_.num_inserts++;
    stats_.num_keys++;
    return {Iterator(leaf, insert_pos), true};
  }

 private:
  // Our criteria for when to expand the root, thereby expanding the key domain.
  // We want to strike a balance between expanding too aggressively and too
  // slowly.
  // Specifically, the number of inserted keys falling to the right of the key
  // domain must have one of two properties: (1) above some maximum threshold,
  // or
  // (2) above some minimum threshold and the number is much more than we would
  // expect from randomness alone.
  bool should_expand_right() const {
    return (!root_node_->is_leaf_ &&
            ((istats_.num_keys_above_key_domain >= kMinOutOfDomainKeys &&
              istats_.num_keys_above_key_domain >=
                  kOutOfDomainToleranceFactor *
                      (stats_.num_keys /
                           istats_.num_keys_at_last_right_domain_resize -
                       1)) ||
             istats_.num_keys_above_key_domain >= kMaxOutOfDomainKeys));
  }

  // Similar to should_expand_right, but for insertions to the left of the key
  // domain.
  bool should_expand_left() const {
    return (!root_node_->is_leaf_ &&
            ((istats_.num_keys_below_key_domain >= kMinOutOfDomainKeys &&
              istats_.num_keys_below_key_domain >=
                  kOutOfDomainToleranceFactor *
                      (stats_.num_keys /
                           istats_.num_keys_at_last_left_domain_resize -
                       1)) ||
             istats_.num_keys_below_key_domain >= kMaxOutOfDomainKeys));
  }

  int best_split_propagation(std::vector<TraversalPathDisk>& traversal_path, DataHeaderOnDisk leaf_node) {
    // if (root_node_->is_leaf_) {
    //   return superroot_->level_;
    // }
    std::vector<SplitDecisionCosts> traversal_costs;
    for (int i = 0; i < traversal_path.size(); i++) {
      double stop_cost;
      uint8_t next_duplication_factor_ = i == traversal_path.size() - 1 ? leaf_node.duplication_factor_ : traversal_path.at(i+1).duplication_factor_;
      int num_children = i == 0 ? 0 : traversal_path.at(i).children_number;
      if (next_duplication_factor_ > 0) {
        stop_cost = 0;
      } else {
        stop_cost = num_children > derived_params_.max_fanout ? std::numeric_limits<double>::max() : num_children * SplitDecisionCosts::base_cost;
      }
      traversal_costs.push_back({stop_cost, num_children <= 2 ? 0 : num_children / 2 + SplitDecisionCosts::base_cost});
    }
    double cumulative_cost = 0;
    double best_cost = std::numeric_limits<double>::max();
    int best_path_level = traversal_path.at(1).header.level - 1;
    for (int i = traversal_costs.size() - 1; i >= 0; i--) {
      SplitDecisionCosts& c = traversal_costs[i];
      if (c.stop_cost != std::numeric_limits<double>::max() &&
          cumulative_cost + c.stop_cost < best_cost) {
        best_cost = cumulative_cost + c.stop_cost;
        best_path_level = traversal_path[i].node->level_;
      }
      cumulative_cost += c.split_cost;
    }
    return best_path_level;
  }

  // When splitting upwards, find best internal node to propagate upwards to.
  // Returns the level of that node. Returns superroot's level if splitting
  // sideways not possible.
  int best_split_propagation(const std::vector<TraversalNode>& traversal_path,
                             bool verbose = false) const {
    if (root_node_->is_leaf_) {
      return superroot_->level_;
    }

    // Find costs on the path down to data node
    std::vector<SplitDecisionCosts> traversal_costs;
    for (const TraversalNode& tn : traversal_path) {
      double stop_cost;
      AlexNode<T, P>* next = tn.node->children_[tn.bucketID];
      if (next->duplication_factor_ > 0) {
        stop_cost = 0;
      } else {
        stop_cost =
            tn.node->num_children_ >= derived_params_.max_fanout
                ? std::numeric_limits<double>::max()
                : tn.node->num_children_ + SplitDecisionCosts::base_cost;
      }
      traversal_costs.push_back(
          {stop_cost,
           tn.node->num_children_ <= 2
               ? 0
               : tn.node->num_children_ / 2 + SplitDecisionCosts::base_cost});
    }

    // Compute back upwards to find the optimal node to stop propagation.
    // Ignore the superroot (the first node in the traversal path).
    double cumulative_cost = 0;
    double best_cost = std::numeric_limits<double>::max();
    int best_path_level = superroot_->level_;
    for (int i = traversal_costs.size() - 1; i >= 0; i--) {
      SplitDecisionCosts& c = traversal_costs[i];
      if (c.stop_cost != std::numeric_limits<double>::max() &&
          cumulative_cost + c.stop_cost < best_cost) {
        best_cost = cumulative_cost + c.stop_cost;
        best_path_level = traversal_path[i].node->level_;
      }
      cumulative_cost += c.split_cost;
    }

    if (verbose) {
      std::cout << "[Best split propagation] best level: " << best_path_level
                << ", parent level: " << traversal_path.back().node->level_
                << ", best cost: " << best_cost
                << ", traversal path (level/addr/n_children): ";
      for (const TraversalNode& tn : traversal_path) {
        std::cout << tn.node->level_ << "/" << tn.node << "/"
                  << tn.node->num_children_ << " ";
      }
      std::cout << std::endl;
    }
    return best_path_level;
  }

  // Expand the key value space that is covered by the index.
  // Expands the root node (which is a model node).
  // If the root node is at the max node size, then we split the root and create
  // a new root node.
  void expand_root(T key, bool expand_left) {
    auto root = static_cast<model_node_type*>(root_node_);

    // Find the new bounds of the key domain.
    // Need to be careful to avoid overflows in the key type.
    T domain_size = istats_.key_domain_max_ - istats_.key_domain_min_;
    int expansion_factor;
    T new_domain_min = istats_.key_domain_min_;
    T new_domain_max = istats_.key_domain_max_;
    data_node_type* outermost_node;
    if (expand_left) {
      auto key_difference = static_cast<double>(istats_.key_domain_min_ -
                                                std::min(key, get_min_key()));
      expansion_factor = pow_2_round_up(static_cast<int>(
          std::ceil((key_difference + domain_size) / domain_size)));
      // Check for overflow. To avoid overflow on signed types while doing
      // this check, we do comparisons using half of the relevant quantities.
      T half_expandable_domain =
          istats_.key_domain_max_ / 2 - std::numeric_limits<T>::lowest() / 2;
      T half_expanded_domain_size = expansion_factor / 2 * domain_size;
      if (half_expandable_domain < half_expanded_domain_size) {
        new_domain_min = std::numeric_limits<T>::lowest();
      } else {
        new_domain_min = istats_.key_domain_max_;
        new_domain_min -= half_expanded_domain_size;
        new_domain_min -= half_expanded_domain_size;
      }
      istats_.num_keys_at_last_left_domain_resize = stats_.num_keys;
      istats_.num_keys_below_key_domain = 0;
      outermost_node = first_data_node();
    } else {
      auto key_difference = static_cast<double>(std::max(key, get_max_key()) -
                                                istats_.key_domain_max_);
      expansion_factor = pow_2_round_up(static_cast<int>(
          std::ceil((key_difference + domain_size) / domain_size)));
      // Check for overflow. To avoid overflow on signed types while doing
      // this check, we do comparisons using half of the relevant quantities.
      T half_expandable_domain =
          std::numeric_limits<T>::max() / 2 - istats_.key_domain_min_ / 2;
      T half_expanded_domain_size = expansion_factor / 2 * domain_size;
      if (half_expandable_domain < half_expanded_domain_size) {
        new_domain_max = std::numeric_limits<T>::max();
      } else {
        new_domain_max = istats_.key_domain_min_;
        new_domain_max += half_expanded_domain_size;
        new_domain_max += half_expanded_domain_size;
      }
      istats_.num_keys_at_last_right_domain_resize = stats_.num_keys;
      istats_.num_keys_above_key_domain = 0;
      outermost_node = last_data_node();
    }
    assert(expansion_factor > 1);

    // Modify the root node appropriately
    int new_nodes_start;  // index of first pointer to a new node
    int new_nodes_end;    // exclusive
    if (static_cast<size_t>(root->num_children_) * expansion_factor <=
        static_cast<size_t>(derived_params_.max_fanout)) {
      // Expand root node
      stats_.num_model_node_expansions++;
      stats_.num_model_node_expansion_pointers += root->num_children_;

      int new_num_children = root->num_children_ * expansion_factor;
      auto new_children = new (pointer_allocator().allocate(new_num_children))
          AlexNode<T, P>*[new_num_children];
      int copy_start;
      if (expand_left) {
        copy_start = new_num_children - root->num_children_;
        new_nodes_start = 0;
        new_nodes_end = copy_start;
        root->model_.b_ += new_num_children - root->num_children_;
      } else {
        copy_start = 0;
        new_nodes_start = root->num_children_;
        new_nodes_end = new_num_children;
      }
      for (int i = 0; i < root->num_children_; i++) {
        new_children[copy_start + i] = root->children_[i];
      }
      pointer_allocator().deallocate(root->children_, root->num_children_);
      root->children_ = new_children;
      root->num_children_ = new_num_children;
    } else {
      // Create new root node
      auto new_root = new (model_node_allocator().allocate(1))
          model_node_type(static_cast<short>(root->level_ - 1), allocator_);
      new_root->model_.a_ = root->model_.a_ / root->num_children_;
      new_root->model_.b_ = root->model_.b_ / root->num_children_;
      if (expand_left) {
        new_root->model_.b_ += expansion_factor - 1;
      }
      new_root->num_children_ = expansion_factor;
      new_root->children_ = new (pointer_allocator().allocate(expansion_factor))
          AlexNode<T, P>*[expansion_factor];
      if (expand_left) {
        new_root->children_[expansion_factor - 1] = root;
        new_nodes_start = 0;
      } else {
        new_root->children_[0] = root;
        new_nodes_start = 1;
      }
      new_nodes_end = new_nodes_start + expansion_factor - 1;
      root_node_ = new_root;
      update_superroot_pointer();
      root = new_root;
    }
    // Determine if new nodes represent a range outside the key type's domain.
    // This happens when we're preventing overflows.
    int in_bounds_new_nodes_start = new_nodes_start;
    int in_bounds_new_nodes_end = new_nodes_end;
    if (expand_left) {
      in_bounds_new_nodes_start =
          std::max(new_nodes_start, root->model_.predict(new_domain_min));
    } else {
      in_bounds_new_nodes_end =
          std::min(new_nodes_end, root->model_.predict(new_domain_max) + 1);
    }

    // Fill newly created child pointers of the root node with new data nodes.
    // To minimize empty new data nodes, we create a new data node per n child
    // pointers, where n is the number of pointers to existing nodes.
    // Requires reassigning some keys from the outermost pre-existing data node
    // to the new data nodes.
    int n = root->num_children_ - (new_nodes_end - new_nodes_start);
    assert(root->num_children_ % n == 0);
    auto new_node_duplication_factor =
        static_cast<uint8_t>(log_2_round_down(n));

    if (expand_left) {
      T left_boundary_value = istats_.key_domain_min_;
      int left_boundary = outermost_node->lower_bound(left_boundary_value);
      data_node_type* next = outermost_node;
      for (int i = new_nodes_end; i > new_nodes_start; i -= n) {
        if (i <= in_bounds_new_nodes_start) {
          // Do not initialize nodes that fall outside the key type's domain
          break;
        }
        int right_boundary = left_boundary;
        if (i - n <= in_bounds_new_nodes_start) {
          left_boundary = 0;
        } else {
          left_boundary_value -= domain_size;
          left_boundary = outermost_node->lower_bound(left_boundary_value);
        }
        data_node_type* new_node = bulk_load_leaf_node_from_existing(
            outermost_node, left_boundary, right_boundary, true);
        new_node->level_ = static_cast<short>(root->level_ + 1);
        new_node->duplication_factor_ = new_node_duplication_factor;
        if (next) {
          next->prev_leaf_ = new_node;
        }
        new_node->next_leaf_ = next;
        next = new_node;
        for (int j = i - 1; j >= i - n; j--) {
          root->children_[j] = new_node;
        }
      }
    }
    else {
      T right_boundary_value = istats_.key_domain_max_;
      int right_boundary = outermost_node->lower_bound(right_boundary_value);
      data_node_type* prev = nullptr;
      for (int i = new_nodes_start; i < new_nodes_end; i += n) {
        if (i >= in_bounds_new_nodes_end) {
          // Do not initialize nodes that fall outside the key type's domain
          break;
        }
        int left_boundary = right_boundary;
        if (i + n >= in_bounds_new_nodes_end) {
          right_boundary = outermost_node->data_capacity_;
        } else {
          right_boundary_value += domain_size;
          right_boundary = outermost_node->lower_bound(right_boundary_value);
        }
        data_node_type* new_node = bulk_load_leaf_node_from_existing(
            outermost_node, left_boundary, right_boundary, true);
        new_node->level_ = static_cast<short>(root->level_ + 1);
        new_node->duplication_factor_ = new_node_duplication_factor;
        if (prev) {
          prev->next_leaf_ = new_node;
        }
        new_node->prev_leaf_ = prev;
        prev = new_node;
        for (int j = i; j < i + n; j++) {
          root->children_[j] = new_node;
        }
      }
    }

    // Connect leaf nodes and remove reassigned keys from outermost pre-existing
    // node.
    if (expand_left) {
      outermost_node->erase_range(new_domain_min, istats_.key_domain_min_);
      auto last_new_leaf =
          static_cast<data_node_type*>(root->children_[new_nodes_end - 1]);
      outermost_node->prev_leaf_ = last_new_leaf;
      last_new_leaf->next_leaf_ = outermost_node;
    }
    else {
      outermost_node->erase_range(istats_.key_domain_max_, new_domain_max,
                                  true);
      auto first_new_leaf =
          static_cast<data_node_type*>(root->children_[new_nodes_start]);
      outermost_node->next_leaf_ = first_new_leaf;
      first_new_leaf->prev_leaf_ = outermost_node;
    }

    istats_.key_domain_min_ = new_domain_min;
    istats_.key_domain_max_ = new_domain_max;
  }

  model_node_type* split_downwards_disk(model_node_type* parent, int bucketID, int fanout_tree_depth,
                      std::vector<fanout_tree::FTNode>& used_fanout_tree_nodes,
                      data_node_type *leaf, DataHeaderOnDisk dh, PhyscialAddr *parent_addr,
                      bool reuse_model, int *empty_size) {
    stats_.num_downward_splits++;
    stats_.num_downward_split_keys += leaf->num_keys_;
    int fanout = 1 << fanout_tree_depth;
    auto new_node = new (model_node_allocator().allocate(1))
        model_node_type(leaf->level_, allocator_);
    new_node->duplication_factor_ = leaf->duplication_factor_;
    new_node->num_children_ = fanout;
    if (hybrid_mode == LEAF_DISK)
        new_node->children_ =
            new (pointer_allocator().allocate(fanout)) AlexNode<T, P>*[fanout];
    new_node->childrenAddrs = new PhyscialAddr[fanout];
    new_node->is_leaf = new int[fanout];
    int repeats = 1 << leaf->duplication_factor_;
    int start_bucketID =
        bucketID - (bucketID % repeats);  // first bucket with same child
    int end_bucketID =
        start_bucketID + repeats;  // first bucket with different child
    double left_boundary_value =
        (start_bucketID - parent->model_.b_) / parent->model_.a_;
    double right_boundary_value =
        (end_bucketID - parent->model_.b_) / parent->model_.a_;
    new_node->model_.a_ =
        1.0 / (right_boundary_value - left_boundary_value) * fanout;
    new_node->model_.b_ = -new_node->model_.a_ * left_boundary_value;

    if (used_fanout_tree_nodes.empty()) {
      assert(fanout_tree_depth == 1);
      create_two_new_data_nodes(leaf, new_node, dh, fanout_tree_depth, reuse_model);
    } else {
      create_new_data_nodes(leaf, new_node, dh, fanout_tree_depth,
                            used_fanout_tree_nodes);
    }
    stats_.num_data_nodes--;
    stats_.num_model_nodes++;



    if (hybrid_mode == ALL_DISK) {
        // write new node on disk
        PhyscialAddr n_addr;
        ModelHeaderOnDisk _mhd = write_model_node(new_node, &n_addr);
        *empty_size = _mhd.empty_byte_count;
        for (int i = start_bucketID; i < end_bucketID; i++) {
            parent->childrenAddrs[i] = n_addr;
        }
        //update its parent in place
        write_model_node(parent, parent_addr, true);
        parent_addr->block = n_addr.block;
        parent_addr->offset = n_addr.offset;
    } else {
        for (int i = start_bucketID; i < end_bucketID; i++) {
            parent->children_[i] = new_node;
            parent->is_leaf[i] = 0;
        }
    }
    return new_node;

  }

  // Splits downwards in the manner determined by the fanout tree and updates
  // the pointers of the parent.
  // If no fanout tree is provided, then splits downward in two. Returns the
  // newly created model node.
  model_node_type* split_downwards(
      model_node_type* parent, int bucketID, int fanout_tree_depth,
      std::vector<fanout_tree::FTNode>& used_fanout_tree_nodes,
      bool reuse_model) {
    auto leaf = static_cast<data_node_type*>(parent->children_[bucketID]);
    stats_.num_downward_splits++;
    stats_.num_downward_split_keys += leaf->num_keys_;

    // Create the new model node that will replace the current data node
    int fanout = 1 << fanout_tree_depth;
    auto new_node = new (model_node_allocator().allocate(1))
        model_node_type(leaf->level_, allocator_);
    new_node->duplication_factor_ = leaf->duplication_factor_;
    new_node->num_children_ = fanout;
    new_node->children_ =
        new (pointer_allocator().allocate(fanout)) AlexNode<T, P>*[fanout];

    int repeats = 1 << leaf->duplication_factor_;
    int start_bucketID =
        bucketID - (bucketID % repeats);  // first bucket with same child
    int end_bucketID =
        start_bucketID + repeats;  // first bucket with different child
    double left_boundary_value =
        (start_bucketID - parent->model_.b_) / parent->model_.a_;
    double right_boundary_value =
        (end_bucketID - parent->model_.b_) / parent->model_.a_;
    new_node->model_.a_ =
        1.0 / (right_boundary_value - left_boundary_value) * fanout;
    new_node->model_.b_ = -new_node->model_.a_ * left_boundary_value;

    // Create new data nodes
    if (used_fanout_tree_nodes.empty()) {
      assert(fanout_tree_depth == 1);
      // just for patameter
      DataHeaderOnDisk old_header;
      create_two_new_data_nodes(leaf, new_node, old_header, fanout_tree_depth, reuse_model);
    } else {
      DataHeaderOnDisk old_header;
      create_new_data_nodes(leaf, new_node, old_header, fanout_tree_depth,
                            used_fanout_tree_nodes);
    }

    delete_node(leaf);
    stats_.num_data_nodes--;
    stats_.num_model_nodes++;
    for (int i = start_bucketID; i < end_bucketID; i++) {
      parent->children_[i] = new_node;
    }
    if (parent == superroot_) {
      root_node_ = new_node;
      update_superroot_pointer();
    }

    return new_node;
  }

    void split_sideways_hybrid_disk(model_node_type* parent, int bucketID,
                             int fanout_tree_depth, DataHeaderOnDisk old_header,
                             std::vector<fanout_tree::FTNode>& used_fanout_tree_nodes,
                             bool reuse_model, data_node_type *leaf) {
        stats_.num_sideways_splits++;
        stats_.num_sideways_split_keys += leaf->num_keys_;
        int fanout = 1 << fanout_tree_depth;
        int repeats = 1 << leaf->duplication_factor_;
        bool has_expand = false;
        if (fanout > repeats) {
            // Expand the pointer array in the parent model node if there are not
            // enough redundant pointers
            stats_.num_model_node_expansions++;
            stats_.num_model_node_expansion_pointers += parent->num_children_;
            int expansion_factor =
                    parent->expand(fanout_tree_depth - leaf->duplication_factor_, 1);
            repeats *= expansion_factor;
            bucketID *= expansion_factor;
            has_expand = true;
        }
        int start_bucketID =
                bucketID - (bucketID % repeats);  // first bucket with same child

        if (used_fanout_tree_nodes.empty()) {
            assert(fanout_tree_depth == 1);
#if _Debug
            printf("create_two\n");
#endif
            create_two_new_data_nodes(
                    leaf, parent, old_header,
                    std::max(fanout_tree_depth,
                             static_cast<int>(leaf->duplication_factor_)),
                    reuse_model, start_bucketID);
        } else {
            // Extra duplication factor is required when there are more redundant
            // pointers than necessary
#if _Debug
            printf("create_new\n");
#endif
            int extra_duplication_factor =
                    std::max(0, leaf->duplication_factor_ - fanout_tree_depth);
            create_new_data_nodes(leaf, parent, old_header, fanout_tree_depth,
                                  used_fanout_tree_nodes, start_bucketID,
                                  extra_duplication_factor);
        }
        // re-write parent
        // delete_node(leaf);
        stats_.num_data_nodes--;
    }

  void split_sideways_disk(model_node_type* parent, int bucketID,
                      int fanout_tree_depth, DataHeaderOnDisk old_header,
                      std::vector<fanout_tree::FTNode>& used_fanout_tree_nodes,
                      bool reuse_model, data_node_type *leaf, int pp_block,
                      int pp_offset, PhyscialAddr *parent_addr, int *empty_size, std::vector<TraversalPathDisk> traversal_path) {
    stats_.num_sideways_splits++;
    stats_.num_sideways_split_keys += leaf->num_keys_;
    int fanout = 1 << fanout_tree_depth;
    int repeats = 1 << leaf->duplication_factor_;
    bool has_expand = false;
    if (fanout > repeats) {
      // Expand the pointer array in the parent model node if there are not
      // enough redundant pointers
      stats_.num_model_node_expansions++;
      stats_.num_model_node_expansion_pointers += parent->num_children_;
      int expansion_factor =
          parent->expand_disk(fanout_tree_depth - leaf->duplication_factor_);
      repeats *= expansion_factor;
      bucketID *= expansion_factor;
      has_expand = true;
    }
    int start_bucketID =
        bucketID - (bucketID % repeats);  // first bucket with same child

    if (used_fanout_tree_nodes.empty()) {
      assert(fanout_tree_depth == 1);
      #if _Debug
      printf("create_two\n");
      #endif
      create_two_new_data_nodes(
          leaf, parent, old_header,
          std::max(fanout_tree_depth,
                   static_cast<int>(leaf->duplication_factor_)),
          reuse_model, start_bucketID);
    } else {
      // Extra duplication factor is required when there are more redundant
      // pointers than necessary
      #if _Debug
      printf("create_new\n");
      #endif
      int extra_duplication_factor =
          std::max(0, leaf->duplication_factor_ - fanout_tree_depth);
      create_new_data_nodes(leaf, parent, old_header, fanout_tree_depth,
                            used_fanout_tree_nodes, start_bucketID,
                            extra_duplication_factor);
    }
    // re-write parent
    if (has_expand) {
      // PhyscialAddr n_addr;
      ModelHeaderOnDisk _mhd = write_model_node(parent, parent_addr);
      *empty_size = _mhd.empty_byte_count;
      if (pp_block == 0) {
        metanode.root_block_id = parent_addr->block;
        metanode.root_offset = parent_addr->offset;
        metanode.is_leaf = parent_addr->flag;
        metanode.dup_root = parent_addr->duplication_factor_;
      }
      else
      {
          if (parent_addr->duplication_factor_ == 0)
            sm->write_arbitrary(pp_block * BlockSize + pp_offset, parent_addr, PhyscialAddrSize);
          else {
              int _size = traversal_path.size();
              assert(_size > 2);
              TraversalPathDisk pparent = traversal_path.at(_size - 2);
              int repeats = 1 << parent_addr->duplication_factor_;
              int start_bucketID = pparent.child_offset - (pparent.child_offset % repeats);
              int end_bucketID = start_bucketID + repeats;
              int start_offset = pparent.addr.block * BlockSize + pparent.addr.offset + ModelHeaderOnDiskSize + pparent.header.empty_byte_count;
              // PhyscialAddr items[parent.header.children_number];
              PhyscialAddr *items = new PhyscialAddr[pparent.header.children_number];
              sm->read_arbitrary(items, start_offset, pparent.header.children_number * PhyscialAddrSize);
              auto model_node = new (model_node_allocator().allocate(1))
                      model_node_type(pparent.header.level, allocator_);
              pparent.header.duplication_factor_ = pparent.addr.duplication_factor_;
              init_model_node_from_disk(model_node, pparent.header, items);
              for (int i = start_bucketID; i < end_bucketID; i++) {
                  model_node->childrenAddrs[i] = *parent_addr;
              }
              write_model_node(model_node, &(pparent.addr), true);
              delete []items;
              delete []model_node->childrenAddrs;
              delete_node(model_node);
          }
      }
    } else {
      ModelHeaderOnDisk _mhd = write_model_node(parent, parent_addr, true);
      *empty_size = _mhd.empty_byte_count;
    }
    // delete_node(leaf);
    stats_.num_data_nodes--;
  }
  // Splits data node sideways in the manner determined by the fanout tree.
  // If no fanout tree is provided, then splits sideways in two.
  void split_sideways(model_node_type* parent, int bucketID,
                      int fanout_tree_depth,
                      std::vector<fanout_tree::FTNode>& used_fanout_tree_nodes,
                      bool reuse_model) {
    auto leaf = static_cast<data_node_type*>(parent->children_[bucketID]);
    stats_.num_sideways_splits++;
    stats_.num_sideways_split_keys += leaf->num_keys_;

    int fanout = 1 << fanout_tree_depth;
    int repeats = 1 << leaf->duplication_factor_;
    if (fanout > repeats) {
      // Expand the pointer array in the parent model node if there are not
      // enough redundant pointers
      stats_.num_model_node_expansions++;
      stats_.num_model_node_expansion_pointers += parent->num_children_;
      int expansion_factor =
          parent->expand(fanout_tree_depth - leaf->duplication_factor_, 0);
      repeats *= expansion_factor;
      bucketID *= expansion_factor;
    }
    int start_bucketID =
        bucketID - (bucketID % repeats);  // first bucket with same child

    if (used_fanout_tree_nodes.empty()) {
      assert(fanout_tree_depth == 1);
      DataHeaderOnDisk old_header;
      #if _Debug
        printf("create_two\n");
        #endif
      create_two_new_data_nodes(
          leaf, parent, old_header,
          std::max(fanout_tree_depth,
                   static_cast<int>(leaf->duplication_factor_)),
          reuse_model, start_bucketID);
    } else {
      // Extra duplication factor is required when there are more redundant
      // pointers than necessary
      DataHeaderOnDisk old_header;
      #if _Debug
        printf("create_new\n");
        #endif
      int extra_duplication_factor =
          std::max(0, leaf->duplication_factor_ - fanout_tree_depth);
      create_new_data_nodes(leaf, parent, old_header, fanout_tree_depth,
                            used_fanout_tree_nodes, start_bucketID,
                            extra_duplication_factor);
    }

    delete_node(leaf);
    stats_.num_data_nodes--;
  }

  // Create two new data nodes by equally dividing the key space of the old data
  // node, insert the new
  // nodes as children of the parent model node starting from a given position,
  // and link the new data nodes together.
  // duplication_factor denotes how many child pointer slots were assigned to
  // the old data node.
  void create_two_new_data_nodes(data_node_type* old_node,
                                 model_node_type* parent, DataHeaderOnDisk old_header,
                                 int duplication_factor, bool reuse_model,
                                 int start_bucketID = 0) {
    assert(duplication_factor >= 1);
    int num_buckets = 1 << duplication_factor;
    int end_bucketID = start_bucketID + num_buckets;
    int mid_bucketID = start_bucketID + num_buckets / 2;

    bool append_mostly_right = old_node->is_append_mostly_right();
    int appending_right_bucketID = std::min<int>(
        std::max<int>(parent->model_.predict(old_node->max_key_), 0),
        parent->num_children_ - 1);
    bool append_mostly_left = old_node->is_append_mostly_left();
    int appending_left_bucketID = std::min<int>(
        std::max<int>(parent->model_.predict(old_node->min_key_), 0),
        parent->num_children_ - 1);

    int right_boundary = old_node->lower_bound(
        (mid_bucketID - parent->model_.b_) / parent->model_.a_);
    // Account for off-by-one errors due to floating-point precision issues.
    while (right_boundary < old_node->data_capacity_ &&
           old_node->get_key(right_boundary) != data_node_type::kEndSentinel_ &&
           parent->model_.predict(old_node->get_key(right_boundary)) <
               mid_bucketID) {
      right_boundary = std::min(
          old_node->get_next_filled_position(right_boundary, false) + 1,
          old_node->data_capacity_);
    }
    data_node_type* left_leaf = bulk_load_leaf_node_from_existing(
        old_node, 0, right_boundary, true, nullptr, reuse_model,
        append_mostly_right && start_bucketID <= appending_right_bucketID &&
            appending_right_bucketID < mid_bucketID,
        append_mostly_left && start_bucketID <= appending_left_bucketID &&
            appending_left_bucketID < mid_bucketID);
    data_node_type* right_leaf = bulk_load_leaf_node_from_existing(
        old_node, right_boundary, old_node->data_capacity_, true, nullptr,
        reuse_model,
        append_mostly_right && mid_bucketID <= appending_right_bucketID &&
            appending_right_bucketID < end_bucketID,
        append_mostly_left && mid_bucketID <= appending_left_bucketID &&
            appending_left_bucketID < end_bucketID);
    left_leaf->level_ = static_cast<short>(parent->level_ + 1);
    right_leaf->level_ = static_cast<short>(parent->level_ + 1);
    left_leaf->duplication_factor_ =
        static_cast<uint8_t>(duplication_factor - 1);
    right_leaf->duplication_factor_ =
        static_cast<uint8_t>(duplication_factor - 1);

    #if DiskSetting
    PhyscialAddr l_addr;
    metanode_data_node.last_data_node_block = old_header.pre_leaf_block;
    metanode_data_node.last_data_node_offset = old_header.pre_leaf_offset;
    write_data_node(left_leaf, &l_addr, false, left_leaf->first_key(), left_leaf->last_key(), 0, 0);
    PhyscialAddr r_addr;
    metanode_data_node.last_data_node_block = l_addr.block;
    metanode_data_node.last_data_node_offset = l_addr.offset;
    write_data_node(right_leaf, &r_addr, false, right_leaf->first_key(), right_leaf->last_key(), old_header.next_leaf_block, old_header.next_leaf_offset);
    #endif

    for (int i = start_bucketID; i < mid_bucketID; i++) {
      #if DiskSetting
        parent->childrenAddrs[i] = l_addr;
        if (hybrid_mode) parent->is_leaf[i] = 1;
      #else
      parent->children_[i] = left_leaf;
      #endif
    }
    for (int i = mid_bucketID; i < end_bucketID; i++) {
      #if DiskSetting
        parent->childrenAddrs[i] = r_addr;
        if (hybrid_mode) parent->is_leaf[i] = 1;
      #else
      parent->children_[i] = right_leaf;
      #endif
    }
    #if !DiskSetting
    link_data_nodes(old_node, left_leaf, right_leaf);
    #else
      delete_node(left_leaf);
      delete_node(right_leaf);
    #endif
  }

  // Create new data nodes from the keys in the old data node according to the
  // fanout tree, insert the new
  // nodes as children of the parent model node starting from a given position,
  // and link the new data nodes together.
  // Helper for splitting when using a fanout tree.
  void create_new_data_nodes(
      const data_node_type* old_node, model_node_type* parent, DataHeaderOnDisk old_header,
      int fanout_tree_depth,
      std::vector<fanout_tree::FTNode>& used_fanout_tree_nodes,
      int start_bucketID = 0, int extra_duplication_factor = 0) {
    bool append_mostly_right = old_node->is_append_mostly_right();
    int appending_right_bucketID = std::min<int>(
        std::max<int>(parent->model_.predict(old_node->max_key_), 0),
        parent->num_children_ - 1);
    bool append_mostly_left = old_node->is_append_mostly_left();
    int appending_left_bucketID = std::min<int>(
        std::max<int>(parent->model_.predict(old_node->min_key_), 0),
        parent->num_children_ - 1);

    // Create the new data nodes
    int cur = start_bucketID;  // first bucket with same child
    int left_boundary = 0;
    int right_boundary = 0;
    // Keys may be re-assigned to an adjacent fanout tree node due to off-by-one
    // errors
    #if DiskSetting
    metanode_data_node.last_data_node_block = old_header.pre_leaf_block;
    metanode_data_node.last_data_node_offset = old_header.pre_leaf_offset;
    #else
     data_node_type* prev_leaf =
        old_node->prev_leaf_;  // used for linking the new data nodes
    #endif
    int num_reassigned_keys = 0;
    int I = 0;
    for (fanout_tree::FTNode& tree_node : used_fanout_tree_nodes) {
      left_boundary = right_boundary;
      auto duplication_factor = static_cast<uint8_t>(
          fanout_tree_depth - tree_node.level + extra_duplication_factor);
      int child_node_repeats = 1 << duplication_factor;
      bool keep_left = append_mostly_right && cur <= appending_right_bucketID &&
                       appending_right_bucketID < cur + child_node_repeats;
      bool keep_right = append_mostly_left && cur <= appending_left_bucketID &&
                        appending_left_bucketID < cur + child_node_repeats;
      right_boundary = tree_node.right_boundary;
      // Account for off-by-one errors due to floating-point precision issues.
      tree_node.num_keys -= num_reassigned_keys;
      num_reassigned_keys = 0;
      while (right_boundary < old_node->data_capacity_ &&
             old_node->get_key(right_boundary) !=
                 data_node_type::kEndSentinel_ &&
             parent->model_.predict(old_node->get_key(right_boundary)) <
                 cur + child_node_repeats) {
        num_reassigned_keys++;
        right_boundary = std::min(
            old_node->get_next_filled_position(right_boundary, false) + 1,
            old_node->data_capacity_);
      }
      tree_node.num_keys += num_reassigned_keys;
      data_node_type* child_node = bulk_load_leaf_node_from_existing(
          old_node, left_boundary, right_boundary, false, &tree_node, false,
          keep_left, keep_right);
      child_node->level_ = static_cast<short>(parent->level_ + 1);
      child_node->cost_ = tree_node.cost;
      child_node->duplication_factor_ = duplication_factor;
      child_node->expected_avg_exp_search_iterations_ =
          tree_node.expected_avg_search_iterations;
      child_node->expected_avg_shifts_ = tree_node.expected_avg_shifts;
      #if DiskSetting
      PhyscialAddr c_addr;
      int next_block = 0;
      int next_offset = 0;
      if (I == used_fanout_tree_nodes.size() - 1) {
        next_block = old_header.next_leaf_block;
        next_offset = old_header.next_leaf_offset;
      }

//        DataHeaderOnDisk r =  write_data_node(node, &l_addr, false, node->first_key(), node->last_key(), dh.next_leaf_block, dh.next_leaf_offset);
      write_data_node(child_node, &c_addr, false, child_node->first_key(), child_node->last_key(), next_block, next_offset);
//      metanode_data_node.last_data_node_block = c_addr.block;
//      metanode_data_node.last_data_node_offset = c_addr.offset;
      #else
      child_node->prev_leaf_ = prev_leaf;
      if (prev_leaf != nullptr) {
        prev_leaf->next_leaf_ = child_node;
      }
      #endif
      for (int i = cur; i < cur + child_node_repeats; i++) {
        #if DiskSetting
        parent->childrenAddrs[i] = c_addr;
          if (hybrid_mode) parent->is_leaf[i] = 1;
        #else
        parent->children_[i] = child_node;
        #endif
      }
      cur += child_node_repeats;
      #if DiskSetting
      metanode_data_node.last_data_node_block = c_addr.block;
      metanode_data_node.last_data_node_offset = c_addr.offset;
      delete_node(child_node);
      #else
      prev_leaf = child_node;
      #endif
      I += 1;
    }
    #if !DiskSetting
    prev_leaf->next_leaf_ = old_node->next_leaf_;
    if (old_node->next_leaf_ != nullptr) {
      old_node->next_leaf_->prev_leaf_ = prev_leaf;
    }
    #endif
  }

  data_node_type* split_upwards(
      T key, int stop_propagation_level,
      const std::vector<TraversalNode>& traversal_path, bool reuse_model,
      model_node_type** new_parent, bool verbose = false) {
    assert(stop_propagation_level >= root_node_->level_);
    std::vector<AlexNode<T, P>*> to_delete;  // nodes that need to be deleted

    // Split the data node into two new data nodes
    const TraversalNode& parent_path_node = traversal_path.back();
    model_node_type* parent = parent_path_node.node;
    auto leaf = static_cast<data_node_type*>(
        parent->children_[parent_path_node.bucketID]);
    int leaf_repeats = 1 << (leaf->duplication_factor_);
    int leaf_start_bucketID =
        parent_path_node.bucketID - (parent_path_node.bucketID % leaf_repeats);
    double leaf_mid_bucketID = leaf_start_bucketID + leaf_repeats / 2.0;
    int leaf_end_bucketID =
        leaf_start_bucketID + leaf_repeats;  // first bucket with next child
    stats_.num_sideways_splits++;
    stats_.num_sideways_split_keys += leaf->num_keys_;

    // Determine if either of the two new data nodes will need to adapt to
    // append-mostly behavior
    bool append_mostly_right = leaf->is_append_mostly_right();
    bool left_half_appending_right = false, right_half_appending_right = false;
    if (append_mostly_right) {
      double appending_right_bucketID =
          parent->model_.predict_double(leaf->max_key_);
      if (appending_right_bucketID >= leaf_start_bucketID &&
          appending_right_bucketID < leaf_mid_bucketID) {
        left_half_appending_right = true;
      } else if (appending_right_bucketID >= leaf_mid_bucketID &&
                 appending_right_bucketID < leaf_end_bucketID) {
        right_half_appending_right = true;
      }
    }
    bool append_mostly_left = leaf->is_append_mostly_left();
    bool left_half_appending_left = false, right_half_appending_left = false;
    if (append_mostly_left) {
      double appending_left_bucketID =
          parent->model_.predict_double(leaf->min_key_);
      if (appending_left_bucketID >= leaf_start_bucketID &&
          appending_left_bucketID < leaf_mid_bucketID) {
        left_half_appending_left = true;
      } else if (appending_left_bucketID >= leaf_mid_bucketID &&
                 appending_left_bucketID < leaf_end_bucketID) {
        right_half_appending_left = true;
      }
    }

    int mid_boundary = leaf->lower_bound(
        (leaf_mid_bucketID - parent->model_.b_) / parent->model_.a_);
    data_node_type* left_leaf = bulk_load_leaf_node_from_existing(
        leaf, 0, mid_boundary, true, nullptr, reuse_model,
        append_mostly_right && left_half_appending_right,
        append_mostly_left && left_half_appending_left);
    data_node_type* right_leaf = bulk_load_leaf_node_from_existing(
        leaf, mid_boundary, leaf->data_capacity_, true, nullptr, reuse_model,
        append_mostly_right && right_half_appending_right,
        append_mostly_left && right_half_appending_left);
    // This is the expected duplication factor; it will be correct once we
    // split/expand the parent
    left_leaf->duplication_factor_ = leaf->duplication_factor_;
    right_leaf->duplication_factor_ = leaf->duplication_factor_;
    left_leaf->level_ = leaf->level_;
    right_leaf->level_ = leaf->level_;
    link_data_nodes(leaf, left_leaf, right_leaf);
    to_delete.push_back(leaf);
    stats_.num_data_nodes--;

    if (verbose) {
      std::cout << "[Splitting upwards data node] level " << leaf->level_
                << ", node addr: " << leaf
                << ", node repeats in parent: " << leaf_repeats
                << ", node indexes in parent: [" << leaf_start_bucketID << ", "
                << leaf_end_bucketID << ")"
                << ", left leaf indexes: [0, " << mid_boundary << ")"
                << ", right leaf indexes: [" << mid_boundary << ", "
                << leaf->data_capacity_ << ")"
                << ", new nodes addr: " << left_leaf << "," << right_leaf
                << std::endl;
    }

    // The new data node that the key falls into is the one we return
    data_node_type* new_data_node;
    if (parent->model_.predict_double(key) < leaf_mid_bucketID) {
      new_data_node = left_leaf;
    } else {
      new_data_node = right_leaf;
    }

    // Split all internal nodes from the parent up to the highest node along the
    // traversal path.
    // As this happens, the entries of the traversal path will go stale, which
    // is fine because we no longer use them.
    // Splitting an internal node involves dividing the child pointers into two
    // halves, and doubling the relevant half.
    AlexNode<T, P>* prev_left_split = left_leaf;
    AlexNode<T, P>* prev_right_split = right_leaf;
    int path_idx = static_cast<int>(traversal_path.size()) - 1;
    while (traversal_path[path_idx].node->level_ > stop_propagation_level) {
      // Decide which half to double
      const TraversalNode& path_node = traversal_path[path_idx];
      model_node_type* cur_node = path_node.node;
      stats_.num_model_node_splits++;
      stats_.num_model_node_split_pointers += cur_node->num_children_;
      bool double_left_half = path_node.bucketID < cur_node->num_children_ / 2;
      model_node_type* left_split = nullptr;
      model_node_type* right_split = nullptr;

      // If one of the resulting halves will only have one child pointer, we
      // should "pull up" that child
      bool pull_up_left_child = false, pull_up_right_child = false;
      AlexNode<T, P>* left_half_first_child = cur_node->children_[0];
      AlexNode<T, P>* right_half_first_child =
          cur_node->children_[cur_node->num_children_ / 2];
      if (double_left_half &&
          (1 << right_half_first_child->duplication_factor_) ==
              cur_node->num_children_ / 2) {
        // pull up right child if all children in the right half are the same
        pull_up_right_child = true;
        left_split = new (model_node_allocator().allocate(1))
            model_node_type(cur_node->level_, allocator_);
      } else if (!double_left_half &&
                 (1 << left_half_first_child->duplication_factor_) ==
                     cur_node->num_children_ / 2) {
        // pull up left child if all children in the left half are the same
        pull_up_left_child = true;
        right_split = new (model_node_allocator().allocate(1))
            model_node_type(cur_node->level_, allocator_);
      } else {
        left_split = new (model_node_allocator().allocate(1))
            model_node_type(cur_node->level_, allocator_);
        right_split = new (model_node_allocator().allocate(1))
            model_node_type(cur_node->level_, allocator_);
      }

      // Do the split
      AlexNode<T, P>* next_left_split = nullptr;
      AlexNode<T, P>* next_right_split = nullptr;
      if (double_left_half) {
        // double left half
        assert(left_split != nullptr);
        if (path_idx == static_cast<int>(traversal_path.size()) - 1) {
          *new_parent = left_split;
        }
        left_split->num_children_ = cur_node->num_children_;
        left_split->children_ =
            new (pointer_allocator().allocate(left_split->num_children_))
                AlexNode<T, P>*[left_split->num_children_];
        left_split->model_.a_ = cur_node->model_.a_ * 2;
        left_split->model_.b_ = cur_node->model_.b_ * 2;
        int cur = 0;
        while (cur < cur_node->num_children_ / 2) {
          AlexNode<T, P>* cur_child = cur_node->children_[cur];
          int cur_child_repeats = 1 << cur_child->duplication_factor_;
          for (int i = 2 * cur; i < 2 * (cur + cur_child_repeats); i++) {
            left_split->children_[i] = cur_child;
          }
          cur_child->duplication_factor_++;
          cur += cur_child_repeats;
        }
        assert(cur == cur_node->num_children_ / 2);

        if (pull_up_right_child) {
          next_right_split = cur_node->children_[cur_node->num_children_ / 2];
          next_right_split->level_ = cur_node->level_;
        } else {
          right_split->num_children_ = cur_node->num_children_ / 2;
          right_split->children_ =
              new (pointer_allocator().allocate(right_split->num_children_))
                  AlexNode<T, P>*[right_split->num_children_];
          right_split->model_.a_ = cur_node->model_.a_;
          right_split->model_.b_ =
              cur_node->model_.b_ - cur_node->num_children_ / 2;
          int j = 0;
          for (int i = cur_node->num_children_ / 2; i < cur_node->num_children_;
               i++) {
            right_split->children_[j] = cur_node->children_[i];
            j++;
          }
          next_right_split = right_split;
        }

        int new_bucketID = path_node.bucketID * 2;
        int repeats = 1 << (prev_left_split->duplication_factor_ + 1);
        int start_bucketID =
            new_bucketID -
            (new_bucketID % repeats);  // first bucket with same child
        int mid_bucketID = start_bucketID + repeats / 2;
        int end_bucketID =
            start_bucketID + repeats;  // first bucket with next child
        for (int i = start_bucketID; i < mid_bucketID; i++) {
          left_split->children_[i] = prev_left_split;
        }
        for (int i = mid_bucketID; i < end_bucketID; i++) {
          left_split->children_[i] = prev_right_split;
        }
        next_left_split = left_split;
      } else {
        // double right half
        assert(right_split != nullptr);
        if (path_idx == static_cast<int>(traversal_path.size()) - 1) {
          *new_parent = right_split;
        }
        if (pull_up_left_child) {
          next_left_split = cur_node->children_[0];
          next_left_split->level_ = cur_node->level_;
        } else {
          left_split->num_children_ = cur_node->num_children_ / 2;
          left_split->children_ =
              new (pointer_allocator().allocate(left_split->num_children_))
                  AlexNode<T, P>*[left_split->num_children_];
          left_split->model_.a_ = cur_node->model_.a_;
          left_split->model_.b_ = cur_node->model_.b_;
          int j = 0;
          for (int i = 0; i < cur_node->num_children_ / 2; i++) {
            left_split->children_[j] = cur_node->children_[i];
            j++;
          }
          next_left_split = left_split;
        }

        right_split->num_children_ = cur_node->num_children_;
        right_split->children_ =
            new (pointer_allocator().allocate(right_split->num_children_))
                AlexNode<T, P>*[right_split->num_children_];
        right_split->model_.a_ = cur_node->model_.a_ * 2;
        right_split->model_.b_ =
            (cur_node->model_.b_ - cur_node->num_children_ / 2) * 2;
        int cur = cur_node->num_children_ / 2;
        while (cur < cur_node->num_children_) {
          AlexNode<T, P>* cur_child = cur_node->children_[cur];
          int cur_child_repeats = 1 << cur_child->duplication_factor_;
          int right_child_idx = cur - cur_node->num_children_ / 2;
          for (int i = 2 * right_child_idx;
               i < 2 * (right_child_idx + cur_child_repeats); i++) {
            right_split->children_[i] = cur_child;
          }
          cur_child->duplication_factor_++;
          cur += cur_child_repeats;
        }
        assert(cur == cur_node->num_children_);

        int new_bucketID =
            (path_node.bucketID - cur_node->num_children_ / 2) * 2;
        int repeats = 1 << (prev_left_split->duplication_factor_ + 1);
        int start_bucketID =
            new_bucketID -
            (new_bucketID % repeats);  // first bucket with same child
        int mid_bucketID = start_bucketID + repeats / 2;
        int end_bucketID =
            start_bucketID + repeats;  // first bucket with next child
        for (int i = start_bucketID; i < mid_bucketID; i++) {
          right_split->children_[i] = prev_left_split;
        }
        for (int i = mid_bucketID; i < end_bucketID; i++) {
          right_split->children_[i] = prev_right_split;
        }
        next_right_split = right_split;
      }
      assert(next_left_split != nullptr && next_right_split != nullptr);
      if (verbose) {
        std::cout << "[Splitting upwards through-node] level "
                  << cur_node->level_ << ", node addr: " << path_node.node
                  << ", node children: " << path_node.node->num_children_
                  << ", child index: " << path_node.bucketID
                  << ", child repeats in node: "
                  << (1 << prev_left_split->duplication_factor_)
                  << ", node repeats in parent: "
                  << (1 << path_node.node->duplication_factor_)
                  << ", new nodes addr: " << left_split << "," << right_split
                  << std::endl;
      }
      to_delete.push_back(cur_node);
      if (!pull_up_left_child && !pull_up_right_child) {
        stats_.num_model_nodes++;
      }
      // This is the expected duplication factor; it will be correct once we
      // split/expand the parent
      next_left_split->duplication_factor_ = cur_node->duplication_factor_;
      next_right_split->duplication_factor_ = cur_node->duplication_factor_;
      prev_left_split = next_left_split;
      prev_right_split = next_right_split;
      path_idx--;
    }

    // Insert into the top node
    const TraversalNode& top_path_node = traversal_path[path_idx];
    model_node_type* top_node = top_path_node.node;
    assert(top_node->level_ == stop_propagation_level);
    if (path_idx == static_cast<int>(traversal_path.size()) - 1) {
      *new_parent = top_node;
    }
    int top_bucketID = top_path_node.bucketID;
    int repeats =
        1 << prev_left_split->duplication_factor_;  // this was the duplication
                                                    // factor of the child that
                                                    // was deleted
    if (verbose) {
      std::cout << "[Splitting upwards top node] level "
                << stop_propagation_level << ", node addr: " << top_node
                << ", node children: " << top_node->num_children_
                << ", child index: " << top_bucketID
                << ", child repeats in node: " << repeats
                << ", node repeats in parent: "
                << (1 << top_node->duplication_factor_) << std::endl;
    }

    // Expand the top node if necessary
    if (repeats == 1) {
      stats_.num_model_node_expansions++;
      stats_.num_model_node_expansion_pointers += top_node->num_children_;
      top_node->expand(1, 0);  // double size of top node
      top_bucketID *= 2;
      repeats *= 2;
    } else {
      prev_left_split->duplication_factor_--;
      prev_right_split->duplication_factor_--;
    }

    int start_bucketID =
        top_bucketID -
        (top_bucketID % repeats);  // first bucket with same child
    int mid_bucketID = start_bucketID + repeats / 2;
    int end_bucketID =
        start_bucketID + repeats;  // first bucket with next child
    for (int i = start_bucketID; i < mid_bucketID; i++) {
      top_node->children_[i] = prev_left_split;
    }
    for (int i = mid_bucketID; i < end_bucketID; i++) {
      top_node->children_[i] = prev_right_split;
    }

    for (auto node : to_delete) {
      delete_node(node);
    }

    return new_data_node;
  }

  /*** Delete ***/

 public:
  // Erases the left-most key with the given key value
  int erase_one(const T& key) {
    data_node_type* leaf = get_leaf(key);
    int num_erased = leaf->erase_one(key);
    stats_.num_keys -= num_erased;
    if (leaf->num_keys_ == 0) {
      merge(leaf, key);
    }
    if (key > istats_.key_domain_max_) {
      istats_.num_keys_above_key_domain -= num_erased;
    } else if (key < istats_.key_domain_min_) {
      istats_.num_keys_below_key_domain -= num_erased;
    }
    return num_erased;
  }

  // Erases all keys with a certain key value
  int erase(const T& key) {
    data_node_type* leaf = get_leaf(key);
    int num_erased = leaf->erase(key);
    stats_.num_keys -= num_erased;
    if (leaf->num_keys_ == 0) {
      merge(leaf, key);
    }
    if (key > istats_.key_domain_max_) {
      istats_.num_keys_above_key_domain -= num_erased;
    } else if (key < istats_.key_domain_min_) {
      istats_.num_keys_below_key_domain -= num_erased;
    }
    return num_erased;
  }

  // Erases element pointed to by iterator
  void erase(Iterator it) {
    if (it.is_end()) {
      return;
    }
    T key = it.key();
    it.cur_leaf_->erase_one_at(it.cur_idx_);
    stats_.num_keys--;
    if (it.cur_leaf_->num_keys_ == 0) {
      merge(it.cur_leaf_, key);
    }
    if (key > istats_.key_domain_max_) {
      istats_.num_keys_above_key_domain--;
    } else if (key < istats_.key_domain_min_) {
      istats_.num_keys_below_key_domain--;
    }
  }

  // Removes all elements
  void clear() {
    for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
         node_it.next()) {
      delete_node(node_it.current());
    }
    auto empty_data_node = new (data_node_allocator().allocate(1))
        data_node_type(key_less_, allocator_);
    empty_data_node->bulk_load(nullptr, 0);
    root_node_ = empty_data_node;
    create_superroot();
    stats_.num_keys = 0;
  }

 private:
  // Try to merge empty leaf, which can be traversed to by looking up key
  // This may cause the parent node to merge up into its own parent
  void merge(data_node_type* leaf, T key) {
    // first save the complete path down to data node
    std::vector<TraversalNode> traversal_path;
    auto leaf_dup = get_leaf(key, &traversal_path);
    // We might need to correct the traversal path in edge cases
    if (leaf_dup != leaf) {
      if (leaf_dup->prev_leaf_ == leaf) {
        correct_traversal_path(leaf, traversal_path, true);
      } else if (leaf_dup->next_leaf_ == leaf) {
        correct_traversal_path(leaf, traversal_path, false);
      } else {
        assert(false);
        return;
      }
    }
    if (traversal_path.size() == 1) {
      return;
    }
    int path_pos = static_cast<int>(traversal_path.size()) - 1;
    TraversalNode tn = traversal_path[path_pos];
    model_node_type* parent = tn.node;
    int bucketID = tn.bucketID;
    int repeats = 1 << leaf->duplication_factor_;

    while (path_pos >= 0) {
      // repeatedly merge leaf with "sibling" leaf by redirecting pointers in
      // the parent
      while (leaf->num_keys_ == 0 && repeats < parent->num_children_) {
        int start_bucketID = bucketID - (bucketID % repeats);
        int end_bucketID = start_bucketID + repeats;
        // determine if the potential sibling leaf is adjacent to the right or
        // left
        bool adjacent_to_right =
            (bucketID % (repeats << 1) == bucketID % repeats);
        data_node_type* adjacent_leaf = nullptr;

        // check if adjacent node is a leaf
        if (adjacent_to_right && parent->children_[end_bucketID]->is_leaf_) {
          adjacent_leaf =
              static_cast<data_node_type*>(parent->children_[end_bucketID]);
        } else if (!adjacent_to_right &&
                   parent->children_[start_bucketID - 1]->is_leaf_) {
          adjacent_leaf = static_cast<data_node_type*>(
              parent->children_[start_bucketID - 1]);
        } else {
          break;  // unable to merge with sibling leaf
        }

        // check if adjacent node is a sibling
        if (leaf->duplication_factor_ != adjacent_leaf->duplication_factor_) {
          break;  // unable to merge with sibling leaf
        }

        // merge with adjacent leaf
        for (int i = start_bucketID; i < end_bucketID; i++) {
          parent->children_[i] = adjacent_leaf;
        }
        if (adjacent_to_right) {
          adjacent_leaf->prev_leaf_ = leaf->prev_leaf_;
          if (leaf->prev_leaf_) {
            leaf->prev_leaf_->next_leaf_ = adjacent_leaf;
          }
        } else {
          adjacent_leaf->next_leaf_ = leaf->next_leaf_;
          if (leaf->next_leaf_) {
            leaf->next_leaf_->prev_leaf_ = adjacent_leaf;
          }
        }
        adjacent_leaf->duplication_factor_++;
        delete_node(leaf);
        stats_.num_data_nodes--;
        leaf = adjacent_leaf;
        repeats = 1 << leaf->duplication_factor_;
      }

      // try to merge up by removing parent and replacing pointers to parent
      // with pointers to leaf in grandparent
      if (repeats == parent->num_children_) {
        leaf->duplication_factor_ = parent->duplication_factor_;
        repeats = 1 << leaf->duplication_factor_;
        bool is_root_node = (parent == root_node_);
        delete_node(parent);
        stats_.num_model_nodes--;

        if (is_root_node) {
          root_node_ = leaf;
          update_superroot_pointer();
          break;
        }

        path_pos--;
        tn = traversal_path[path_pos];
        parent = tn.node;
        bucketID = tn.bucketID;
        int start_bucketID = bucketID - (bucketID % repeats);
        int end_bucketID = start_bucketID + repeats;
        for (int i = start_bucketID; i < end_bucketID; i++) {
          parent->children_[i] = leaf;
        }
      } else {
        break;  // unable to merge up
      }
    }
  }

  /*** Stats ***/

 public:
  // Number of elements
  size_t size() const { return static_cast<size_t>(stats_.num_keys); }

  // True if there are no elements
  bool empty() const { return (size() == 0); }

  // This is just a function required by the STL standard. ALEX can hold more
  // items.
  size_t max_size() const { return size_t(-1); }

  // Size in bytes of all the keys, payloads, and bitmaps stored in this index
  long long data_size() const {
    long long size = 0;
    for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
         node_it.next()) {
      AlexNode<T, P>* cur = node_it.current();
      if (cur->is_leaf_) {
        size += static_cast<data_node_type*>(cur)->data_size();
      }
    }
    return size;
  }

  // Size in bytes of all the model nodes (including pointers) and metadata in
  // data nodes
  long long model_size() const {
    long long size = 0;
    for (NodeIterator node_it = NodeIterator(this); !node_it.is_end();
         node_it.next()) {
      size += node_it.current()->node_size();
    }
    return size;
  }

  // Total number of nodes in the RMI
  int num_nodes() const {
    return stats_.num_data_nodes + stats_.num_model_nodes;
  };

  // Number of data nodes in the RMI
  int num_leaves() const { return stats_.num_data_nodes; };

  // Return a const reference to the current statistics
  const struct Stats& get_stats() const { return stats_; }

  /*** Debugging ***/

 public:
  // If short_circuit is true, then we stop validating after detecting the first
  // invalid property
  bool validate_structure(bool validate_data_nodes = false,
                          bool short_circuit = false) const {
    bool is_valid = true;
    std::stack<AlexNode<T, P>*> node_stack;
    AlexNode<T, P>* cur;
    node_stack.push(root_node_);

    while (!node_stack.empty()) {
      cur = node_stack.top();
      node_stack.pop();

      if (!cur->is_leaf_) {
        auto node = static_cast<model_node_type*>(cur);
        if (!node->validate_structure(true)) {
          std::cout << "[Model node invalid structure]"
                    << " node addr: " << node
                    << ", node level: " << node->level_ << std::endl;
          if (short_circuit) {
            return false;
          } else {
            is_valid = false;
          }
        }

        node_stack.push(node->children_[node->num_children_ - 1]);
        for (int i = node->num_children_ - 2; i >= 0; i--) {
          if (node->children_[i] != node->children_[i + 1]) {
            node_stack.push(node->children_[i]);
          }
        }
      } else {
        if (validate_data_nodes) {
          auto node = static_cast<data_node_type*>(cur);
          if (!node->validate_structure(true)) {
            std::cout << "[Data node invalid structure]"
                      << " node addr: " << node
                      << ", node level: " << node->level_ << std::endl;
            if (short_circuit) {
              return false;
            } else {
              is_valid = false;
            }
          }
          if (node->num_keys_ > 0) {
            data_node_type* prev_nonempty_leaf = node->prev_leaf_;
            while (prev_nonempty_leaf != nullptr &&
                   prev_nonempty_leaf->num_keys_ == 0) {
              prev_nonempty_leaf = prev_nonempty_leaf->prev_leaf_;
            }
            if (prev_nonempty_leaf) {
              T last_in_prev_leaf = prev_nonempty_leaf->last_key();
              T first_in_cur_leaf = node->first_key();
              if (last_in_prev_leaf >= first_in_cur_leaf) {
                std::cout
                    << "[Data node keys not in sorted order with prev node]"
                    << " node addr: " << node
                    << ", node level: " << node->level_
                    << ", last in prev leaf: " << last_in_prev_leaf
                    << ", first in cur leaf: " << first_in_cur_leaf
                    << std::endl;
                if (short_circuit) {
                  return false;
                } else {
                  is_valid = false;
                }
              }
            }
            data_node_type* next_nonempty_leaf = node->next_leaf_;
            while (next_nonempty_leaf != nullptr &&
                   next_nonempty_leaf->num_keys_ == 0) {
              next_nonempty_leaf = next_nonempty_leaf->next_leaf_;
            }
            if (next_nonempty_leaf) {
              T first_in_next_leaf = next_nonempty_leaf->first_key();
              T last_in_cur_leaf = node->last_key();
              if (last_in_cur_leaf >= first_in_next_leaf) {
                std::cout
                    << "[Data node keys not in sorted order with next node]"
                    << " node addr: " << node
                    << ", node level: " << node->level_
                    << ", last in cur leaf: " << last_in_cur_leaf
                    << ", first in next leaf: " << first_in_next_leaf
                    << std::endl;
                if (short_circuit) {
                  return false;
                } else {
                  is_valid = false;
                }
              }
            }
          }
        }
      }
    }
    return is_valid;
  }

  /*** Iterators ***/

 public:
  class Iterator {
   public:
    data_node_type* cur_leaf_ = nullptr;  // current data node
    int cur_idx_ = 0;         // current position in key/data_slots of data node
    int cur_bitmap_idx_ = 0;  // current position in bitmap
    uint64_t cur_bitmap_data_ = 0;  // caches the relevant data in the current
                                    // bitmap position

    Iterator() {}

    Iterator(data_node_type* leaf, int idx) : cur_leaf_(leaf), cur_idx_(idx) {
      initialize();
    }

    Iterator(const Iterator& other)
        : cur_leaf_(other.cur_leaf_),
          cur_idx_(other.cur_idx_),
          cur_bitmap_idx_(other.cur_bitmap_idx_),
          cur_bitmap_data_(other.cur_bitmap_data_) {}

    Iterator(const ReverseIterator& other)
        : cur_leaf_(other.cur_leaf_), cur_idx_(other.cur_idx_) {
      initialize();
    }

    Iterator& operator=(const Iterator& other) {
      if (this != &other) {
        cur_idx_ = other.cur_idx_;
        cur_leaf_ = other.cur_leaf_;
        cur_bitmap_idx_ = other.cur_bitmap_idx_;
        cur_bitmap_data_ = other.cur_bitmap_data_;
      }
      return *this;
    }

    Iterator& operator++() {
      advance();
      return *this;
    }

    Iterator operator++(int) {
      Iterator tmp = *this;
      advance();
      return tmp;
    }

#if ALEX_DATA_NODE_SEP_ARRAYS
    // Does not return a reference because keys and payloads are stored
    // separately.
    // If possible, use key() and payload() instead.
    V operator*() const {
      return std::make_pair(cur_leaf_->key_slots_[cur_idx_],
                            cur_leaf_->payload_slots_[cur_idx_]);
    }
#else
    // If data node stores key-payload pairs contiguously, return reference to V
    V& operator*() const { return cur_leaf_->data_slots_[cur_idx_]; }
#endif

    const T& key() const { return cur_leaf_->get_key(cur_idx_); }

    P& payload() const { return cur_leaf_->get_payload(cur_idx_); }

    bool is_end() const { return cur_leaf_ == nullptr; }

    bool operator==(const Iterator& rhs) const {
      return cur_idx_ == rhs.cur_idx_ && cur_leaf_ == rhs.cur_leaf_;
    }

    bool operator!=(const Iterator& rhs) const { return !(*this == rhs); };

   private:
    void initialize() {
      if (!cur_leaf_) return;
      assert(cur_idx_ >= 0);
      if (cur_idx_ >= cur_leaf_->data_capacity_) {
        cur_leaf_ = cur_leaf_->next_leaf_;
        cur_idx_ = 0;
        if (!cur_leaf_) return;
      }

      cur_bitmap_idx_ = cur_idx_ >> 6;
      cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];

      // Zero out extra bits
      int bit_pos = cur_idx_ - (cur_bitmap_idx_ << 6);
      cur_bitmap_data_ &= ~((1ULL << bit_pos) - 1);

      (*this)++;
    }

    forceinline void advance() {
      while (cur_bitmap_data_ == 0) {
        cur_bitmap_idx_++;
        if (cur_bitmap_idx_ >= cur_leaf_->bitmap_size_) {
          cur_leaf_ = cur_leaf_->next_leaf_;
          cur_idx_ = 0;
          if (cur_leaf_ == nullptr) {
            return;
          }
          cur_bitmap_idx_ = 0;
        }
        cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];
      }
      uint64_t bit = extract_rightmost_one(cur_bitmap_data_);
      cur_idx_ = get_offset(cur_bitmap_idx_, bit);
      cur_bitmap_data_ = remove_rightmost_one(cur_bitmap_data_);
    }
  };

  class ConstIterator {
   public:
    const data_node_type* cur_leaf_ = nullptr;  // current data node
    int cur_idx_ = 0;         // current position in key/data_slots of data node
    int cur_bitmap_idx_ = 0;  // current position in bitmap
    uint64_t cur_bitmap_data_ = 0;  // caches the relevant data in the current
                                    // bitmap position

    ConstIterator() {}

    ConstIterator(const data_node_type* leaf, int idx)
        : cur_leaf_(leaf), cur_idx_(idx) {
      initialize();
    }

    ConstIterator(const Iterator& other)
        : cur_leaf_(other.cur_leaf_),
          cur_idx_(other.cur_idx_),
          cur_bitmap_idx_(other.cur_bitmap_idx_),
          cur_bitmap_data_(other.cur_bitmap_data_) {}

    ConstIterator(const ConstIterator& other)
        : cur_leaf_(other.cur_leaf_),
          cur_idx_(other.cur_idx_),
          cur_bitmap_idx_(other.cur_bitmap_idx_),
          cur_bitmap_data_(other.cur_bitmap_data_) {}

    ConstIterator(const ReverseIterator& other)
        : cur_leaf_(other.cur_leaf_), cur_idx_(other.cur_idx_) {
      initialize();
    }

    ConstIterator(const ConstReverseIterator& other)
        : cur_leaf_(other.cur_leaf_), cur_idx_(other.cur_idx_) {
      initialize();
    }

    ConstIterator& operator=(const ConstIterator& other) {
      if (this != &other) {
        cur_idx_ = other.cur_idx_;
        cur_leaf_ = other.cur_leaf_;
        cur_bitmap_idx_ = other.cur_bitmap_idx_;
        cur_bitmap_data_ = other.cur_bitmap_data_;
      }
      return *this;
    }

    ConstIterator& operator++() {
      advance();
      return *this;
    }

    ConstIterator operator++(int) {
      ConstIterator tmp = *this;
      advance();
      return tmp;
    }

#if ALEX_DATA_NODE_SEP_ARRAYS
    // Does not return a reference because keys and payloads are stored
    // separately.
    // If possible, use key() and payload() instead.
    V operator*() const {
      return std::make_pair(cur_leaf_->key_slots_[cur_idx_],
                            cur_leaf_->payload_slots_[cur_idx_]);
    }
#else
    // If data node stores key-payload pairs contiguously, return reference to V
    const V& operator*() const { return cur_leaf_->data_slots_[cur_idx_]; }
#endif

    const T& key() const { return cur_leaf_->get_key(cur_idx_); }

    const P& payload() const { return cur_leaf_->get_payload(cur_idx_); }

    bool is_end() const { return cur_leaf_ == nullptr; }

    bool operator==(const ConstIterator& rhs) const {
      return cur_idx_ == rhs.cur_idx_ && cur_leaf_ == rhs.cur_leaf_;
    }

    bool operator!=(const ConstIterator& rhs) const { return !(*this == rhs); };

   private:
    void initialize() {
      if (!cur_leaf_) return;
      assert(cur_idx_ >= 0);
      if (cur_idx_ >= cur_leaf_->data_capacity_) {
        cur_leaf_ = cur_leaf_->next_leaf_;
        cur_idx_ = 0;
        if (!cur_leaf_) return;
      }

      cur_bitmap_idx_ = cur_idx_ >> 6;
      cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];

      // Zero out extra bits
      int bit_pos = cur_idx_ - (cur_bitmap_idx_ << 6);
      cur_bitmap_data_ &= ~((1ULL << bit_pos) - 1);

      (*this)++;
    }

    forceinline void advance() {
      while (cur_bitmap_data_ == 0) {
        cur_bitmap_idx_++;
        if (cur_bitmap_idx_ >= cur_leaf_->bitmap_size_) {
          cur_leaf_ = cur_leaf_->next_leaf_;
          cur_idx_ = 0;
          if (cur_leaf_ == nullptr) {
            return;
          }
          cur_bitmap_idx_ = 0;
        }
        cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];
      }
      uint64_t bit = extract_rightmost_one(cur_bitmap_data_);
      cur_idx_ = get_offset(cur_bitmap_idx_, bit);
      cur_bitmap_data_ = remove_rightmost_one(cur_bitmap_data_);
    }
  };

  class ReverseIterator {
   public:
    data_node_type* cur_leaf_ = nullptr;  // current data node
    int cur_idx_ = 0;         // current position in key/data_slots of data node
    int cur_bitmap_idx_ = 0;  // current position in bitmap
    uint64_t cur_bitmap_data_ = 0;  // caches the relevant data in the current
                                    // bitmap position

    ReverseIterator() {}

    ReverseIterator(data_node_type* leaf, int idx)
        : cur_leaf_(leaf), cur_idx_(idx) {
      initialize();
    }

    ReverseIterator(const ReverseIterator& other)
        : cur_leaf_(other.cur_leaf_),
          cur_idx_(other.cur_idx_),
          cur_bitmap_idx_(other.cur_bitmap_idx_),
          cur_bitmap_data_(other.cur_bitmap_data_) {}

    ReverseIterator(const Iterator& other)
        : cur_leaf_(other.cur_leaf_), cur_idx_(other.cur_idx_) {
      initialize();
    }

    ReverseIterator& operator=(const ReverseIterator& other) {
      if (this != &other) {
        cur_idx_ = other.cur_idx_;
        cur_leaf_ = other.cur_leaf_;
        cur_bitmap_idx_ = other.cur_bitmap_idx_;
        cur_bitmap_data_ = other.cur_bitmap_data_;
      }
      return *this;
    }

    ReverseIterator& operator++() {
      advance();
      return *this;
    }

    ReverseIterator operator++(int) {
      ReverseIterator tmp = *this;
      advance();
      return tmp;
    }

#if ALEX_DATA_NODE_SEP_ARRAYS
    // Does not return a reference because keys and payloads are stored
    // separately.
    // If possible, use key() and payload() instead.
    V operator*() const {
      return std::make_pair(cur_leaf_->key_slots_[cur_idx_],
                            cur_leaf_->payload_slots_[cur_idx_]);
    }
#else
    // If data node stores key-payload pairs contiguously, return reference to V
    V& operator*() const { return cur_leaf_->data_slots_[cur_idx_]; }
#endif

    const T& key() const { return cur_leaf_->get_key(cur_idx_); }

    P& payload() const { return cur_leaf_->get_payload(cur_idx_); }

    bool is_end() const { return cur_leaf_ == nullptr; }

    bool operator==(const ReverseIterator& rhs) const {
      return cur_idx_ == rhs.cur_idx_ && cur_leaf_ == rhs.cur_leaf_;
    }

    bool operator!=(const ReverseIterator& rhs) const {
      return !(*this == rhs);
    };

   private:
    void initialize() {
      if (!cur_leaf_) return;
      assert(cur_idx_ >= 0);
      if (cur_idx_ >= cur_leaf_->data_capacity_) {
        cur_leaf_ = cur_leaf_->next_leaf_;
        cur_idx_ = 0;
        if (!cur_leaf_) return;
      }

      cur_bitmap_idx_ = cur_idx_ >> 6;
      cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];

      // Zero out extra bits
      int bit_pos = cur_idx_ - (cur_bitmap_idx_ << 6);
      cur_bitmap_data_ &= (1ULL << bit_pos) | ((1ULL << bit_pos) - 1);

      advance();
    }

    forceinline void advance() {
      while (cur_bitmap_data_ == 0) {
        cur_bitmap_idx_--;
        if (cur_bitmap_idx_ < 0) {
          cur_leaf_ = cur_leaf_->prev_leaf_;
          if (cur_leaf_ == nullptr) {
            cur_idx_ = 0;
            return;
          }
          cur_idx_ = cur_leaf_->data_capacity_ - 1;
          cur_bitmap_idx_ = cur_leaf_->bitmap_size_ - 1;
        }
        cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];
      }
      assert(cpu_supports_bmi());
      #if M1_Chip
      int bit_pos = static_cast<int>(63 - _nlz(cur_bitmap_data_));
      #else
      int bit_pos = static_cast<int>(63 - _lzcnt_u64(cur_bitmap_data_));
      #endif
      cur_idx_ = (cur_bitmap_idx_ << 6) + bit_pos;
      cur_bitmap_data_ &= ~(1ULL << bit_pos);
    }
  };

  class ConstReverseIterator {
   public:
    const data_node_type* cur_leaf_ = nullptr;  // current data node
    int cur_idx_ = 0;         // current position in key/data_slots of data node
    int cur_bitmap_idx_ = 0;  // current position in bitmap
    uint64_t cur_bitmap_data_ = 0;  // caches the relevant data in the current
                                    // bitmap position

    ConstReverseIterator() {}

    ConstReverseIterator(const data_node_type* leaf, int idx)
        : cur_leaf_(leaf), cur_idx_(idx) {
      initialize();
    }

    ConstReverseIterator(const ConstReverseIterator& other)
        : cur_leaf_(other.cur_leaf_),
          cur_idx_(other.cur_idx_),
          cur_bitmap_idx_(other.cur_bitmap_idx_),
          cur_bitmap_data_(other.cur_bitmap_data_) {}

    ConstReverseIterator(const ReverseIterator& other)
        : cur_leaf_(other.cur_leaf_),
          cur_idx_(other.cur_idx_),
          cur_bitmap_idx_(other.cur_bitmap_idx_),
          cur_bitmap_data_(other.cur_bitmap_data_) {}

    ConstReverseIterator(const Iterator& other)
        : cur_leaf_(other.cur_leaf_), cur_idx_(other.cur_idx_) {
      initialize();
    }

    ConstReverseIterator(const ConstIterator& other)
        : cur_leaf_(other.cur_leaf_), cur_idx_(other.cur_idx_) {
      initialize();
    }

    ConstReverseIterator& operator=(const ConstReverseIterator& other) {
      if (this != &other) {
        cur_idx_ = other.cur_idx_;
        cur_leaf_ = other.cur_leaf_;
        cur_bitmap_idx_ = other.cur_bitmap_idx_;
        cur_bitmap_data_ = other.cur_bitmap_data_;
      }
      return *this;
    }

    ConstReverseIterator& operator++() {
      advance();
      return *this;
    }

    ConstReverseIterator operator++(int) {
      ConstReverseIterator tmp = *this;
      advance();
      return tmp;
    }

#if ALEX_DATA_NODE_SEP_ARRAYS
    // Does not return a reference because keys and payloads are stored
    // separately.
    // If possible, use key() and payload() instead.
    V operator*() const {
      return std::make_pair(cur_leaf_->key_slots_[cur_idx_],
                            cur_leaf_->payload_slots_[cur_idx_]);
    }
#else
    // If data node stores key-payload pairs contiguously, return reference to V
    const V& operator*() const { return cur_leaf_->data_slots_[cur_idx_]; }
#endif

    const T& key() const { return cur_leaf_->get_key(cur_idx_); }

    const P& payload() const { return cur_leaf_->get_payload(cur_idx_); }

    bool is_end() const { return cur_leaf_ == nullptr; }

    bool operator==(const ConstReverseIterator& rhs) const {
      return cur_idx_ == rhs.cur_idx_ && cur_leaf_ == rhs.cur_leaf_;
    }

    bool operator!=(const ConstReverseIterator& rhs) const {
      return !(*this == rhs);
    };

   private:
    void initialize() {
      if (!cur_leaf_) return;
      assert(cur_idx_ >= 0);
      if (cur_idx_ >= cur_leaf_->data_capacity_) {
        cur_leaf_ = cur_leaf_->next_leaf_;
        cur_idx_ = 0;
        if (!cur_leaf_) return;
      }

      cur_bitmap_idx_ = cur_idx_ >> 6;
      cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];

      // Zero out extra bits
      int bit_pos = cur_idx_ - (cur_bitmap_idx_ << 6);
      cur_bitmap_data_ &= (1ULL << bit_pos) | ((1ULL << bit_pos) - 1);

      advance();
    }

    forceinline void advance() {
      while (cur_bitmap_data_ == 0) {
        cur_bitmap_idx_--;
        if (cur_bitmap_idx_ < 0) {
          cur_leaf_ = cur_leaf_->prev_leaf_;
          if (cur_leaf_ == nullptr) {
            cur_idx_ = 0;
            return;
          }
          cur_idx_ = cur_leaf_->data_capacity_ - 1;
          cur_bitmap_idx_ = cur_leaf_->bitmap_size_ - 1;
        }
        cur_bitmap_data_ = cur_leaf_->bitmap_[cur_bitmap_idx_];
      }
      assert(cpu_supports_bmi());
      #if M1_Chip
      int bit_pos = static_cast<int>(63 - _nlz(cur_bitmap_data_));
      #else
      int bit_pos = static_cast<int>(63 - _lzcnt_u64(cur_bitmap_data_));
      #endif
      cur_idx_ = (cur_bitmap_idx_ << 6) + bit_pos;
      cur_bitmap_data_ &= ~(1ULL << bit_pos);
    }
  };

  // Iterates through all nodes with pre-order traversal
  class NodeIterator {
   public:
    const self_type* index_;
    AlexNode<T, P>* cur_node_;
    std::stack<AlexNode<T, P>*> node_stack_;  // helps with traversal
    bool is_leaf_disk = false;

    // Start with root as cur and all children of root in stack
    explicit NodeIterator(const self_type* index)
        : index_(index), cur_node_(index->root_node_) {
      if (cur_node_ && !cur_node_->is_leaf_) {
        auto node = static_cast<model_node_type*>(cur_node_);
          if (!is_leaf_disk || is_leaf_disk && node->is_leaf[node->num_children_ - 1] == 0)
            node_stack_.push(node->children_[node->num_children_ - 1]);
        for (int i = node->num_children_ - 2; i >= 0; i--) {
            if (!is_leaf_disk || is_leaf_disk && node->is_leaf[i] == 0)
              if (node->children_[i] != node->children_[i + 1]) {
                node_stack_.push(node->children_[i]);
          }
        }
      }
    }

    AlexNode<T, P>* current() const { return cur_node_; }

    AlexNode<T, P>* next() {
      if (node_stack_.empty()) {
        cur_node_ = nullptr;
        return nullptr;
      }

      cur_node_ = node_stack_.top();
      node_stack_.pop();

      if (!cur_node_->is_leaf_) {
        auto node = static_cast<model_node_type*>(cur_node_);
          if (!is_leaf_disk || is_leaf_disk && node->is_leaf[node->num_children_ - 1] == 0)
            node_stack_.push(node->children_[node->num_children_ - 1]);
        for (int i = node->num_children_ - 2; i >= 0; i--) {
            if (!is_leaf_disk || is_leaf_disk && node->is_leaf[i] == 0)
              if (node->children_[i] != node->children_[i + 1]) {
                node_stack_.push(node->children_[i]);
              }
        }
      }

      return cur_node_;
    }

    bool is_end() const { return cur_node_ == nullptr; }
  };
};
}  // namespace alex
