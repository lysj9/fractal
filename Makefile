CC 		= 	gcc

EXEC 	= 	frac.exe

C_SOURCES 	= 	main.c frac.c randomz.c randqueue.c mempool.c getopt.c generate_binary.c generate_mass.c make_mass.c scale.c quick_sort_widx.c output.c quick_sort.c quick_select.c
#C_OBJS 		= 	$(patsubst %.c, %.o, $(C_SOURCES))
C_OBJS 		= 	$(C_SOURCES:%.c=%.o)
C_DEPS 		= 	$(C_SOURCES:%.c=%.d)

OBJS 		= 	$(C_OBJS)
DEPS 		= 	$(C_DEPS)
INCL 		= 	type.h func.h constant.h Makefile

CFLAGS 		= 	-W -Wall -O3 -ffast-math -fopenmp
LDFLAGS 	= 	-lm -fopenmp

$(EXEC) 	: 	$(OBJS) Makefile
	$(CC) -o $(EXEC) $(OBJS) $(LDFLAGS)

$(C_OBJS) : %.o : %.c Makefile
	$(CC) $(CFLAGS) -c -o $@ $<

-include $(DEPS)

$(C_DEPS) : %.d : %.c
	@set -e; rm -f $@; \
		$(CC) -MM $< > $@.$$$$; \
		sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' < $@.$$$$ > $@; \
		rm -f $@.$$$$

.PHONY 			: 	clean
clean:
	-rm -f $(OBJS) $(DEPS) $(EXEC)
