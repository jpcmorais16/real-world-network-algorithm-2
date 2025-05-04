CXX = g++
CXXFLAGS = -Wall -std=c++11 -fopenmp
TARGET = network_metrics
SRCS = main.cpp networkMetrics.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
	
run:
	./$(TARGET) 