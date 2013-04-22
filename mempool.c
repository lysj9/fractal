#include <stdio.h>
#include <stdlib.h>

typedef union FreeNode{
	union FreeNode *next;
}FreeNode;

typedef struct MemBlock{
	struct MemBlock *next;
	int length;
	void *buf;
}MemBlock;

static FreeNode *FreeNodeHeader=NULL;
static MemBlock *MemBlockHeader=NULL;
static size_t nodesize=sizeof(FreeNode);
static int ITEM_NUM=8192;

void mempool_init(size_t data_size)
{
	if ( data_size>nodesize )
		nodesize = data_size;
	int i;
	MemBlockHeader = (MemBlock*) malloc( sizeof(MemBlock) );
	MemBlockHeader->next = NULL;
	MemBlockHeader->length = ITEM_NUM;
	MemBlockHeader->buf = malloc(ITEM_NUM*nodesize);

	FreeNodeHeader = (FreeNode*) (MemBlockHeader->buf);
	for (i=1;i<ITEM_NUM;++i){
//		FreeNodeHeader->next = (FreeNode*) (MemBlockHeader->buf+i*nodesize);
//	void* -> char*, otherwise icc would have compiling error
//	"expression must be pointer to a complete object type" (not to void)
		FreeNodeHeader->next = (FreeNode*) ((char*)MemBlockHeader->buf+i*nodesize);
		FreeNodeHeader = FreeNodeHeader->next;
	}
	FreeNodeHeader = (FreeNode*) (MemBlockHeader->buf);
	
	return;
}

void* mempool_malloc()
{
	int i;
//	printf("mempool_malloc in\n");
	if ( FreeNodeHeader==NULL ){
//		printf("generate a memblock!\n");
		ITEM_NUM *= 2;
		MemBlock *newblock = (MemBlock*) malloc( sizeof(MemBlock) );
		newblock->buf = malloc(ITEM_NUM*nodesize);

		newblock->next = MemBlockHeader;
		MemBlockHeader = newblock;

		FreeNodeHeader = (FreeNode*) (MemBlockHeader->buf);
		for (i=1;i<ITEM_NUM;++i){
//			FreeNodeHeader->next = (FreeNode*) (MemBlockHeader->buf+i*nodesize);
//	void* -> char*, otherwise icc would have compiling error
//	"expression must be pointer to a complete object type" (not to void)
			FreeNodeHeader->next = (FreeNode*) ((char*)MemBlockHeader->buf+i*nodesize);
			FreeNodeHeader = FreeNodeHeader->next;
		}
		FreeNodeHeader = (FreeNode*) (MemBlockHeader->buf);
	}
	void *p = FreeNodeHeader;
	FreeNodeHeader = FreeNodeHeader->next;
	return p;
}

void mempool_free(void *p)
{
	FreeNode *pNode = (FreeNode*) p;
	pNode->next = FreeNodeHeader;
	FreeNodeHeader = pNode;
	return;
}

void mempool_destroy()
{
	while(!MemBlockHeader){
		MemBlock *pBlock = MemBlockHeader->next;
		free(MemBlockHeader);
		MemBlockHeader = pBlock;
	}
	return;
}
