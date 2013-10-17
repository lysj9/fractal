CC 		= 	gcc

EXEC 	= 	frac.exe

C_SOURCES 	= 	main.c binary_pairing.c frac.c generate_binary.c generate_mass.c get_wtime.c getopt.c make_mass.c mempool.c output.c quick_sort.c quick_select.c randomz.c randqueue.c scale.c
#C_OBJS 		= 	$(patsubst %.c, %.o, $(C_SOURCES))
C_OBJS 		= 	$(C_SOURCES:%.c=%.o)
C_DEPS 		= 	$(C_SOURCES:%.c=%.d)

OBJS 		= 	$(C_OBJS)
DEPS 		= 	$(C_DEPS)
INCL 		=

CFLAGS 		= 	-W -Wall -O3 -ffast-math -fopenmp
LDFLAGS 	= 	-lm -fopenmp

$(EXEC) 	: 	$(OBJS) Makefile
	$(CC) -o $(EXEC) $(OBJS) $(LDFLAGS)

$(C_OBJS) : %.o : %.c Makefile
	$(CC) $(CFLAGS) -c -o $@ $<

-include $(DEPS)

#$(C_DEPS) : %.d : %.c
#	@set -e; rm -f $@; \
#		$(CC) -MM $< > $@.$$$$; \
#		sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' < $@.$$$$ > $@; \
#		rm -f $@.$$$$
$(C_DEPS) : %.d : %.c
	@set -e; rm -f $@; \
		$(CC) -MM $< | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@

.PHONY 			: 	clean
clean:
	-rm -f $(OBJS) $(DEPS) $(EXEC)
