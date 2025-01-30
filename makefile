CFLAGS= -D _GNU_SOURCE -D __USE_XOPEN  -I$(RSTPATH)/include -O3 -fPIC -Wall -D_GNU_SOURCE -D_LINUX

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn -I$(IPATH)/analysis

LIBS=-L$(RSTPATH)/lib -lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lrpos.1  -lcnvmap.1 -loldcnvmap.1 -lshf.1 -lgrd.1 -loldgrd.1 -laacgm.1 -ldmap.1 -lrfile.1 -lopt.1 -lcfit.1 -ligrf.1

.c.o:
	icx $(CFLAGS) $(INCLUDE) -c -o $@ $<


v_image_merge:	v_image_merge.o 
	icx -o ~/bin/v_image_merge $(CFLAGS) $(INCLUDE) v_image_merge.o -lm $(LIBS)
