
CXXFLAGS = -W -Wall -O3

.PHONY : clean default

default : dumb-code-counter

dumb-code-counter : dumb-code-counter.cc

clean :
	$(RM) *~ dumb-code-counter
