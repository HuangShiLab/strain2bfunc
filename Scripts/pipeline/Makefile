CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
        exist = $(shell if [ -e '/usr/local/bin/g++-10' ]; then echo "exist"; else echo "notexist"; fi;)
        ifeq ($(exist),exist)
                CC:=g++-10
        else
                exist = $(shell if [ -e '/usr/local/bin/g++-9' ]; then echo "exist"; else echo "notexist"; fi;)
                ifeq ($(exist),exist)
                        CC:=g++-9
                else
                        CC:=g++-8
                endif
        endif
endif
OMPFLG=-fopenmp
EXE_TAX=Strain2bfunc-pipeline
tax:$(OBJ_TAX) pipeline.cpp
	$(CC) -o $(EXE_TAX) pipeline.cpp $(OMPFLG)
