
CXXFLAGS = -W -Wall -O3

LDLIBS=-lgmpxx

.PHONY : clean default

default : dumb-code-counter smart-code-counter

dumb-code-counter : dumb-code-counter.cc

smart-code-counter : smart-code-counter.cc

clean :
	$(RM) *~ dumb-code-counter smart-code-counter
