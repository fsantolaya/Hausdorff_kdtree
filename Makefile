.SILENT:
CPP=g++

BINS=main hausdorff

CPPFLAGS=-O3 -DNDEBUG -w -fcompare-debug-second
DEST=.

%.o: %.c
	$(CPP) $(CPPFLAGS) -c $< -o $@

all: clean bin

bin: $(BINS)

main:
	g++ $(CPPFLAGS) -o $(DEST)/program main.cpp TimeMesure.c Build_cBiK.cpp NodeSKQ.cpp Node.cpp MinHeap.cpp MaxHeap.cpp -lm

hausdorff:
	g++ $(CPPFLAGS) -o $(DEST)/hausdorff getTime.cpp utils.cpp ryu-kamata.cpp TimeMesure.c Build_cBiK.cpp NodeSKQ.cpp Node.cpp MinHeap.cpp MaxHeap.cpp -lm
clean:
	rm -f  $(BINS)
	cd $(DEST); rm -f *.a $(BINS)
