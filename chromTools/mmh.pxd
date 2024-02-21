# mmh.pxd
    
cdef extern from "murmurhash.c":
    int murmurhash_32(char *key, int len, int seed)