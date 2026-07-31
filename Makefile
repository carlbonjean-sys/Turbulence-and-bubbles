# ============================================================================
#  Makefile -- DNS bulle en HIT (Basilisk)
# ----------------------------------------------------------------------------
#  Cibles :
#    make           -> build production MPI, sans graphique  (bin/bubble)
#    make movie      -> build MPI + film inline (necessite OSMesa/GL)
#    make render     -> outil de rendu offline               (bin/render)
#    make serial     -> build sequentiel pour debug
#    make clean
#
#  Prerequis : Basilisk installe, $BASILISK pointant vers .../basilisk/src,
#  qcc et mpicc dans le PATH.
# ============================================================================

BASILISK ?= $(HOME)/basilisk/src
QCC       = qcc
MPICC     = mpicc
# -D_GNU_SOURCE       : expose madvise()/MADV_DONTNEED (grilles Basilisk).
# -disable-dimensions : desactive l'analyse dimensionnelle de qcc (incompatible
#                       avec le melange reduced.h + forcage maison).
CFLAGS    = -O3 -Wall -D_GNU_SOURCE -disable-dimensions
GLLIBS    = -L$(BASILISK)/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lGL
LDLIBS    = -lm

BIN  = bubble
RBIN = render

.PHONY: all movie render serial clean

# --- Production MPI (sans graphique) : c'est ce que lance run_pipeline.sh ----
all: $(BIN)

$(BIN): src/main.c
	CC99='$(MPICC) -std=c99' $(QCC) -D_MPI=1 $(CFLAGS) src/main.c -o $(BIN) $(LDLIBS)

# --- Production MPI + film inline (MOVIE) -----------------------------------
movie: src/main.c
	CC99='$(MPICC) -std=c99' $(QCC) -D_MPI=1 -DMOVIE=1 $(CFLAGS) \
	    src/main.c -o $(BIN) $(GLLIBS) $(LDLIBS)

# --- Outil de rendu offline (sequentiel, avec GL) ---------------------------
render: src/render.c
	$(QCC) $(CFLAGS) src/render.c -o $(RBIN) $(GLLIBS) $(LDLIBS)

# --- Build sequentiel pour debug --------------------------------------------
serial: src/main.c
	$(QCC) -g -O2 -Wall -D_GNU_SOURCE -disable-dimensions src/main.c -o $(BIN)_serial $(LDLIBS)

clean:
	rm -f $(BIN) $(BIN)_serial $(RBIN) *.o _*.c
