/**
 * @file utility.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-02-12
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <stdint.h>
const long BlockSize =  long(8192/2);
#define KeyType uint64_t
#define ValueType uint64_t

#define AddNewEntry 1
#define NoOperation 2

#define SplitInnerNode 1
#define SplitLeafNode 2
#define AppendInnerNode 3
#define AppendLeafNode 4

#define InnerNodeType 0
#define LeafNodeType 1

// typedef struct NodePointer {
//     char flag;
//     short offset;
// } NodePointer;
// #define NodePointerSize sizeof(NodePointer)

typedef struct InnerNodeIterm {
    KeyType key;
    int block_id;
} InnerNodeIterm;
#define InnerNodeItermSize sizeof(InnerNodeIterm)

typedef struct LeafNodeIterm{
    KeyType key;
    ValueType value;
} LeafNodeIterm;
#define LeafNodeItemSize sizeof(LeafNodeIterm)

typedef struct InnerNodeHeader {
    char node_type;
    short item_count;
    int next_block_id;
    int level;
} InnerNodeHeader;
#define InnerNodeHeaderSize sizeof(InnerNodeHeader)

typedef struct LeaftNodeHeader {
    char node_type;
    short item_count;
    int next_block_id;
    int level;
} LeaftNodeHeader;
#define LeaftNodeHeaderSize sizeof(LeaftNodeHeader)

typedef struct MetaNode {
    int block_count;
    int root_block_id;
    int level;
    int offset;
    // use for hybrid case
    int last_block;
} MetaNode;
#define MetaNodeSize sizeof(MetaNode)

typedef struct {
    char data[BlockSize];
} Block;

typedef struct {
    int status;
    InnerNodeIterm ini;
    int level;
    int added_block;
    KeyType update_key;
} BuildStatus;

/*
ItemCount
p1 to item1
p2 to item2
...

...
item2
item1
*/
#define MaxItemInInnerNode ((BlockSize - InnerNodeHeaderSize) / InnerNodeItermSize )
#define MaxItemInLeafNode ((BlockSize - LeaftNodeHeaderSize) / LeafNodeItemSize)

#define SearchThreshold 40


