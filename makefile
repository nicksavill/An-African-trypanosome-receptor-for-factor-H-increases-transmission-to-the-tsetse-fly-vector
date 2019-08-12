DEBUG =
SRCSTEM = model
EXESTEM = ba
CC = gcc
WARNINGS = -Wall -Wno-format -Wno-unused-result -Wno-maybe-uninitialized \
-Wno-misleading-indentation -fno-stack-protector
ifdef DEBUG
	OPT = -g
else
	OPT = -O3
endif
TARGETS = $(patsubst $(SRCSTEM)%.c, $(EXESTEM)%, $(wildcard $(SRCSTEM)*.c))
BAYES = $(HOME)/Bayes/Dynamical
FLAGS = $(OPT) $(WARNINGS) -I$(BAYES) -fopenmp
LIBS = -lm -lgsl -lgslcblas
OBJECTS = $(BAYES)/bayes.o $(BAYES)/func.o $(BAYES)/output.o

all: $(TARGETS)
$(EXESTEM)%: $(SRCSTEM)%.o $(OBJECTS) $(FOBJS)
	$(CC) $(FLAGS) $+ -o $@ $(LIBS) 
$(SRCSTEM)%.o: $(SRCSTEM)%.c $(BAYES)/model.h $(HEADERS)
	$(CC) -c $(FLAGS) $<
$(BAYES)/bayes.o: $(BAYES)/bayes.c $(BAYES)/bayes.h $(BAYES)/model.h
	$(CC) -c $(FLAGS) -o $(BAYES)/bayes.o $(BAYES)/bayes.c
$(BAYES)/func.o: $(BAYES)/func.c $(BAYES)/func.h
	$(CC) -c $(FLAGS) -o $(BAYES)/func.o $(BAYES)/func.c
$(BAYES)/output.o: $(BAYES)/output.c $(BAYES)/bayes.h $(BAYES)/model.h
	$(CC) -c $(FLAGS) -o $(BAYES)/output.o $(BAYES)/output.c


clean:
	rm -f *.o
	rm -f $(BAYES)/*.o
	rm -f $(TARGETS)
