
CXXFLAGS = -W -Wall -O3

LDLIBS=-lgmpxx -lgmp

.PHONY : clean default

default : dumb-code-counter smart-code-counter

dumb-code-counter : dumb-code-counter.cc

smart-code-counter : smart-code-counter.cc

paper.pdf : paper.tex
	pdflatex paper.tex
	pdflatex paper.tex

clean :
	$(RM) *~ dumb-code-counter smart-code-counter paper.aux paper.log paper.pdf
