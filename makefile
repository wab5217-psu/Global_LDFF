CFLAGS= -D _GNU_SOURCE -D __USE_XOPEN  -I$(RSTPATH)/include -O3 -ipo -mcmodel medium -fPIC -Wall -D_GNU_SOURCE -D_LINUX -Wno-unused-but-set-variable


# -D_Float32=float -D_Float64=float -D_Float32x=float -D_Float64x=float

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn -I$(IPATH)/analysis -I/home/wab5217/src/lib

LIBS=-L /home/wab5217/src/rst/lib -lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lrpos.1 -laacgm.1 -lrfile.1 -lcfit.1 -ligrf.1 -lmlt_v2.1


# LIBS=-lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lrpos.1  -lcnvmap.1 -loldcnvmap.1 -lshf.1 -lgrd.1 -loldgrd.1 -laacgm.1 -ldmap.1 -lrfile.1 -lopt.1 -lcfit.1 -ligrf.1 -lmlt_v2.1

# LIBS= -lcnvmap.1 -lshf.1 -loldcnvmap.1 -lgrd.1 -loldgrd.1 -lradar.1 \
# 	-ldmap.1 -lrfile.1 -lrtime.1 -lopt.1 -lrcnv.1 -laacgm.1 \
# 	-laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lmlt_v2.1 -lrtime.1 \
# 	-lfit.1 -lrscan.1


MYLIB=/home/wab5217/src/lib/


.c.o:
	icx $(CFLAGS) $(INCLUDE) -c -o $@ $<

OBJS=ml_df_fit.o sub_sphazm.o sub_sphcal.o sub_gcp.o cgls.o

ml_df_fit_ng:	$(OBJS)
	icx -o ~/bin/ml_df_fit_ng $(CFLAGS) $(INCLUDE) $(OBJS)  -lgsl -lgslcblas -lm  -lz -qmkl=parallel -qopenmp $(LIBS)
	cp ~/bin/ml_df_fit_ng ~/bin/ml_df_fit_ng_s
