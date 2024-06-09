
CXXFLAGS = -W -Wall -O3

LDLIBS=-lgmpxx -lgmp

.PHONY : clean default

default : dumb_code_counter smart_code_counter

dumb_code_counter : dumb_code_counter.cc

smart_code_counter : smart_code_counter.cc

paper.pdf : paper.tex
	pdflatex paper.tex
	pdflatex paper.tex

clean :
	$(RM) *~ dumb_code_counter smart_code_counter paper.aux paper.log paper.pdf
