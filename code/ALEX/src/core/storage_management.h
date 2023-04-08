/**
 * @file storage_management.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-02-12
 *
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include<stdio.h>
#include<stdlib.h>
#include <map>
#include<iostream>
#include <cstring>
#define Caching 0
const long BlockSize = 8192L;
typedef struct {
    int next_block;
    int next_offset;
    int root_block_id;
    int root_offset;
    // used in bulk node
    int last_data_node_block;
    int last_data_node_offset;
    //
    int left_most_block;
    int left_most_offset;
    int right_most_block;
    int right_most_offset;
    int is_leaf;
    uint8_t dup_left;
    uint8_t dup_right;
    uint8_t dup_root;
} MetaNode;

typedef struct {
    char data[BlockSize];
} Block;
#define MetaNodeSize sizeof(MetaNode)
namespace alex {
    class StorageManager {
        private:
        char *file_name = nullptr;
        FILE *fp = nullptr;

        void _write_block(void *data, int block_id) {
            fseek(fp, block_id * BlockSize, SEEK_SET);
            fwrite(data, BlockSize, 1, fp);
            return;
        }

        char* _allocate_new_block() {
            char* ptr = nullptr;
            ptr = (char *)malloc(BlockSize * sizeof(char));
            return ptr;
        }

        void _get_file_handle() {
            fp = fopen(file_name,"r+b");
        }

        void _create_file() {
            fp = fopen(file_name,"wb");
            char empty_block[BlockSize];
            for (int i = 0; i < BlockSize; i++) {
                empty_block[i] = 0;
            }
            MetaNode mn;
            mn.next_block = 1;
            mn.next_offset = 0;
            mn.root_block_id = -1;
            mn.root_offset = -1;
            mn.last_data_node_block = -1;
            mn.last_data_node_offset = -1;
            // mn.level = 1;
            memcpy(empty_block, &mn, MetaNodeSize);
            _write_block(empty_block, 0);
            _close_file_handle();
            _get_file_handle();
        }

        void _close_file_handle() {
            fclose(fp);
        }

        public:
            StorageManager(char *fn, bool first = false) {
                file_name = fn;
                if (first) {
                    _create_file();
                } else {
                    _get_file_handle();
                }
                
            }

            StorageManager() = default;

            ~StorageManager() {
                if (fp != nullptr) _close_file_handle();
            }  

            size_t get_file_size() {
                fseek(fp,0,SEEK_END);
                return ftell(fp);
            }

            void _read_block(void *data, int block_id) {
                fseek(fp, block_id * BlockSize, SEEK_SET);
                fread(data, BlockSize, 1, fp);
            return;
        }

            Block get_block(int block_id) {
                Block block;
                //char data[BlockSize];
                _read_block(block.data, block_id);
                //memcpy(block.data, data, BlockSize);
                return block;
            }

            void get_block(int block_id, char *data) {
                _read_block(data, block_id);
            }

            void write_block(int block_id, Block block) {
                _write_block(&block, block_id);
            }

            void write_with_size(int block_id, void *data, long size) {
                fseek(fp, block_id * BlockSize, SEEK_SET);
                fwrite(data, size, 1, fp);
                return;
            }

            void write_arbitrary(long offset, void *data, long size) {
                fseek(fp, offset, SEEK_SET);
                fwrite(data, size, 1, fp);
                return;
            }

            void read_block_arbitrary(void *data, long offset) {
                fseek(fp, offset, SEEK_SET);
                fread(data, BlockSize, 1, fp);
                return;
            }

            void read_arbitrary(void *data, long offset, long len) {
                fseek(fp, offset, SEEK_SET);
                fread(data, len, 1, fp);
                return;
            }

    };
}