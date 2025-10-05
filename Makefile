CXX       = g++
CXXFLAGS  = -std=c++17 -O2 -Wall -Wextra -mbmi2
LDFLAGS   =
TARGETS   = merge_cli

# adjust paths if your sources are elsewhere
SRCS_MERGE_CLI = test_merge_cli.cpp louds-trie.cpp
OBJS_MERGE_CLI = $(SRCS_MERGE_CLI:.cpp=.o)

all: $(TARGETS)

merge_cli: $(OBJS_MERGE_CLI)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS_MERGE_CLI) $(TARGETS)
