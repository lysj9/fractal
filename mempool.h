#ifndef _MEMPOOL_H_
#define _MEMPOOL_H_

void mempool_init(size_t data_size);
void* mempool_malloc();
void mempool_free(void *p);
void mempool_destroy();

#endif
