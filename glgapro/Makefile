# /*****************************************************************
# //
# //  Author:             Ehsan Kourkchi
# //
# //  DATE:             August 20, 2016
# //
# //  FILE:             Makefile        
# //
# //  DESCRIPTION:      Makefile for GLGA binary files
# //
# //****************************************************************/
# 
# Prerequisites: To create pixval, you need to have cfitsio installed
#

all: bin/pixval
	

bin/pixval:
	c++ -o bin/pixval cpp/pixval.cpp -Wno-deprecated -lcfitsio



.PHONY: clean

clean:
	rm -f bin/*
