CC 		= 	gcc

EXEC 	= 	frac.exe

C_SOURCES 	= 	main.c frac.c randomz.c randqueue.c mempool.c getopt.c generate_binary.c generate_mass.c make_mass.c scale.c quick_sort_widx.c output.c
#C_SOURCES 	= 	main.c frac.c randomz.c randqueue.c getopt.c generate_binary.c generate_mass.c make_mass.c scale.c quick_sort_widx.c output.c
C_OBJS 		= 	$(patsubst %.c, %.o, $(C_SOURCES))
OBJS 		= 	$(C_OBJS)
INCL 		= 	type.h func.h constant.h Makefile

CFLAGS 		= 	-W -Wall -O3 -ffast-math
LDFLAGS 	= 	-lm

$(EXEC) 	: 	$(OBJS) $(INCL)
	$(CC) -o $(EXEC) $(OBJS) $(LDFLAGS)

$(OBJS) 	: 	$(C_SOURCES) $(INCL)
	$(CC) $(CFLAGS) -c $(C_SOURCES)

.PHONY 			: 	clean
clean:
	-rm -f $(OBJS) $(EXEC)
