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
CFLAGS := -g -Wall -fopenmp -pedantic -std=c++11
LIB :=
INC := -I src

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	$(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(MKBDIR)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	rm -rf $(BUILDDIR) $(TARGET)

run:
	bin/$(NAME)

.PHONY: clean
