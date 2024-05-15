CCXX := g++
IDIR := include
LDIR := lib
ODIR := build
SRCDIR := src
DOCDIR := documents

CXXFLAGS := -I${IDIR}

TARGET := TestHilbertCurve.cpp

_DEPENDENCIES := SpaceFillingCurves.h
DEPENDENCIES  := ${patsubst %,${IDIR}/%,${_DEPENDENCIES}}

_OBJECTFILES := ${patsubst %.h,%.o,${_DEPENDENCIES}} ${patsubst %.cpp,%.o,${TARGET}}
OBJECTFILES  := ${patsubst %,${ODIR}/%,${_OBJECTFILES}}

all: ${OBJECTFILES}
	${CCXX} -std=c++17 -O3 -o run $^

${ODIR}/%.o: %.cpp ${DEPENDENCIES}
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

${ODIR}/%.o: ${SRCDIR}/%.cpp ${DEPENDENCIES}
	${CCXX} -std=c++17 -O3 -c -o $@ $< ${CXXFLAGS}

run:
	./run

.PHONY: clean run document

document:
	@xelatex -synctex=1 -interaction=nonstopmode ${DOCDIR}/Document.tex >> /dev/null
	@rm -rf ${DOCDIR}/*.aux ${DOCDIR}/*.log ${DOCDIR}/*.gz ${DOCDIR}/*.out
	@rm -rf *.aux *.log *.gz *.out
	@mv *.pdf ${DOCDIR}/

clean:
	rm -rf ${ODIR}/*.o *~ run ${IDIR}/*~
	rm -rf Output/*.js Output/*.tec
