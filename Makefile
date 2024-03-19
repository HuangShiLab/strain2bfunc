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
BUILDFLG=-Wno-attributes
EXE_TAX=bin/Strain2bfunc-pipeline
EXE_FUN=Scripts/func/calculate_ko_abd
tax:$(OBJ_TAX) Scripts/pipeline/pipeline.cpp
	$(CC) -o $(EXE_TAX) Scripts/pipeline/pipeline.cpp $(OMPFLG)
	$(CC) -o $(EXE_FUN) Scripts/func/calculate_ko_abd.cpp $(BUILDFLG)
