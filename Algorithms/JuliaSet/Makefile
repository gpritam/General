CCXX := g++
IDIR := include
LDIR := lib
ODIR := build
SRCDIR := src
DOCDIR := documents

CXXFLAGS := -I${IDIR}

TARGET := MandlebrotSet.cpp

all: 
	${CCXX} -std=c++17 -O3 -o run ${TARGET}

run:
	./run

.PHONY: clean run document

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o *~ run ${IDIR}/*~ Output/*.tec Output/*.vtk Output/*.jpeg Output/*.png
