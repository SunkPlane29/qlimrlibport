CXX = g++
TARGET = libqlimr.so
OBJDIR = build
SRCDIR = src

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
$(info Compiling sources: $(SOURCES))

HEADERS = $(wildcard $(SRCDIR)/*.h)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
DEPFILES = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.d,$(SOURCES))

CXXFLAGS = -I./ -Wall -Wextra -g3 -fPIC
LDFLAGS = -L/usr/local/lib -g3 -shared
STD = -std=gnu++11
LIBS = -lgsl -lgslcblas -lm
OPTIMIZATION = -O3

.SECONDARY: $(DEPFILES)

all: $(TARGET)

$(TARGET): $(OBJECTS) | build 					 
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJDIR)/%.d: $(SRCDIR)/%.cpp | build
	$(CXX) $(CXXFLAGS) $(OPTIMIZATION) $(STD) -MM -MT '$(OBJDIR)/$*.o' -MF $@ $<

$(OBJDIR)/%.o:  $(SRCDIR)/%.cpp $(OBJDIR)/%.d | build
	$(CXX) $(CXXFLAGS) $(OPTIMIZATION) $(STD) -c $< -o $@

.PHONY: build
build:
	@mkdir -p $(OBJDIR)

.PHONY: clean
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(OBJDIR) *.o *.so $(TARGET) $(DEPFILES)

.PHONY: build-main
build-main:
	$(CXX) -g main.cpp -L. -l:libqlimr.so -I./src -o main.out