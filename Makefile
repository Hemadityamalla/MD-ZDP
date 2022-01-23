MAKE_DIR := $(PWD)
SRC_DIR := $(MAKE_DIR)/src
PROG := zdPlasma




.PHONY: all 

all: $(PROG)
	

$(PROG): $(SRC_DIR)/main.f90
	@$(MAKE) -C $(SRC_DIR) -f Makefile.makefile
	cp src/main .

clean:
	@$(MAKE) -C $(SRC_DIR) -f Makefile.makefile clean
