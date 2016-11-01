NAME := cheetah
LIBNAME := libcheetah.a

CC := g++
SRCDIR := src
BUILDDIR := build
MODULES := app approximators cheetah lib solvers
MKBDIR := $(addprefix build/, $(MODULES))
TARGET := bin/$(NAME)
LIBTARGET := bin/$(LIBNAME)
LIBINCL := bin/include
LIBMODULES := $(addprefix $(LIBINCL)/,$(MODULES))

SRCEXT := cpp
HDREXT := h
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
HEADERSRC := $(shell find $(SRCDIR) -type f -name *.$(HDREXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
HEADERS := $(patsubst $(SRCDIR)/%,$(LIBINCL)/%,$(HEADERSRC))
CFLAGS := -g -Wall -Wextra -pedantic -std=c++11 -O3 -ffast-math -march=native -fopenmp
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

$(LIBMODULES):
	mkdir -p $@

$(HEADERS): $(LIBMODULES)
	cp $(patsubst $(LIBINCL)/%,$(SRCDIR)/%,$@) $@

linklib: $(OBJECTS)
	@echo " Linking to static library..."
	ar rcs $(LIBTARGET) $^

lib: linklib $(HEADERS)

clean:
	@echo " Cleaning...";
	rm -rf $(BUILDDIR) $(TARGET) $(LIBTARGET) $(LIBINCL)

run: $(TARGET)
	$(TARGET)

test: $(TESTOBJS) $(OBJECTS)
	@echo " Running tests..."
	$(CC) $^ -o $(TESTTARGET) $(LIB) $(TESTLIB)
	$(TESTTARGET)

.PHONY: clean
