NAME := cheetah

CC := g++
SRCDIR := src
BUILDDIR := build
MODULES := app lib solvers
MKBDIR := $(addprefix build/, $(MODULES))
TARGET := bin/$(NAME)

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -Wall -Wextra -pedantic -std=c++11 -O3 -ffast-math
LIB := -fopenmp
INC := -I src

TESTDIR := tests
TESTSRCS := $(shell find $(TESTDIR) -type f -name *.$(SRCEXT))
TESTOBJS := $(patsubst $(TESTDIR)/%,$(BUILDDIR)/$(TESTDIR)/%,$(TESTSRCS:.$(SRCEXT)=.o))
TESTTARGET := bin/tests
TESTLIB := -lgtest -lgtest_main -lpthread
TESTFLAG :=

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	$(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(MKBDIR)
	$(CC) $(CFLAGS) $(TESTFLAG) $(INC) -c -o $@ $<

$(BUILDDIR)/$(TESTDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)/$(TESTDIR)
	$(eval TESTFLAG := -DTEST_MAIN)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	rm -rf $(BUILDDIR) $(TARGET)

run: $(TARGET)
	$(TARGET)

test: $(TESTOBJS) $(OBJECTS)
	@echo " Running tests..."
	$(CC) $^ -o $(TESTTARGET) $(LIB) $(TESTLIB)
	$(TESTTARGET)

.PHONY: clean
