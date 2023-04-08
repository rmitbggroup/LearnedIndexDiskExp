namespace pgm {
    const long BLOCK_SIZE = 8192/2;

    void write_data(FILE *fp, void *data, long offset, int len) {
        fseek(fp, offset, SEEK_SET);
        fwrite(data, len, 1, fp);
        return;
    }

    void read_block(FILE *fp, void *data, int block_id) {
        fseek(fp, block_id * BLOCK_SIZE, SEEK_SET);
        fread(data, BLOCK_SIZE, 1, fp);
        return;
    }

    void write_block(FILE *fp, void *data, int block_id) {
        fseek(fp, block_id * BLOCK_SIZE, SEEK_SET);
        fwrite(data, BLOCK_SIZE, 1, fp);
        return;
    }

    void read_data(FILE *fp, void *data, long offset, int len) {
        fseek(fp, offset, SEEK_SET);
        fread(data, len, 1, fp);
        return;
    }
}
