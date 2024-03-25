CC:=g++

# Check if the operating system is macOS
ifeq ($(shell uname), Darwin)
    # Check if g++-13 is installed on macOS
    GPP_13 := $(shell command -v g++-13 2> /dev/null)
    ifeq ($(GPP_13),)
        # If g++-13 is not installed, check if g++-12 is installed
        GPP_12 := $(shell command -v g++-12 2> /dev/null)
        ifeq ($(GPP_12),)
            # If g++-12 is not installed, check if g++-11 is installed
            GPP_11 := $(shell command -v g++-11 2> /dev/null)
            ifeq ($(GPP_11),)
                # If g++-11 is not installed, check if g++-10 is installed
                GPP_10 := $(shell command -v g++-10 2> /dev/null)
                ifeq ($(GPP_10),)
                    # If none of g++-13, g++-12, g++-11, and g++-10 is installed, set CC to g++
                    CC := g++
                else
                    # If g++-10 is installed, set CC to g++-10
                    CC := g++-10
                endif
            else
                # If g++-11 is installed, set CC to g++-11
                CC := g++-11
            endif
        else
            # If g++-12 is installed, set CC to g++-12
            CC := g++-12
        endif
    else
        # If g++-13 is installed, set CC to g++-13
        CC := g++-13
    endif
endif

OMPFLG=-fopenmp
BUILDFLG=-Wno-attributes
EXE_TAX=bin/Strain2bFunc-pipeline
EXE_FUN=Scripts/func/calculate_ko_abd
tax:$(OBJ_TAX) Scripts/pipeline/pipeline.cpp
	$(CC) -o $(EXE_TAX) Scripts/pipeline/pipeline.cpp $(OMPFLG)
	$(CC) -o $(EXE_FUN) Scripts/func/calculate_ko_abd.cpp $(BUILDFLG)
