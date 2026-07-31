BASILISK ?= $(HOME)/basilisk/src
QCC       = qcc
MPICC     = mpicc
CFLAGS    = -O3 -Wall -D_GNU_SOURCE -disable-dimensions
GLLIBS    = -L$(BASILISK)/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lGL
LDLIBS    = -lm

BIN  = bubble
RBIN = render

.PHONY: all movie render serial clean

all: $(BIN)

$(BIN): src/main.c
	CC99='$(MPICC) -std=c99' $(QCC) -D_MPI=1 $(CFLAGS) src/main.c -o $(BIN) $(LDLIBS)

movie: src/main.c
	CC99='$(MPICC) -std=c99' $(QCC) -D_MPI=1 -DMOVIE=1 $(CFLAGS) src/main.c -o $(BIN) $(GLLIBS) $(LDLIBS)

render: src/render.c
	$(QCC) $(CFLAGS) src/render.c -o $(RBIN) $(GLLIBS) $(LDLIBS)

serial: src/main.c
	$(QCC) -g -O2 -Wall -D_GNU_SOURCE -disable-dimensions src/main.c -o $(BIN)_serial $(LDLIBS)

clean:
	rm -f $(BIN) $(BIN)_serial $(RBIN) *.o _*.c
