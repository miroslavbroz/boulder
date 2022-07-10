# Makefile
# Makefile for simplex and chi2.
# Miroslav Broz (miroslav.broz@email.cz), Mar 13th 2011

f77 = gfortran
cc = gcc

#opt = -O -pg -fbacktrace -Wconversion -Wsurprising -Wunderflow
opt = -O3

libs = -L.

obj = \
  boulder/coll_rate.o \
  boulder/coll_rate_interp.o \
  boulder/cumul.o \
  boulder/cumul4.o \
  boulder/cumul5.o \
  boulder/focusing_factor.o \
  boulder/iso_bodies.o \
  boulder/piscaling.o \
  boulder/pop_change.o \
  boulder/pop_check.o \
  boulder/pop_decay.o \
  boulder/pop_decay_interp.o \
  boulder/read_collprob.o \
  boulder/read_collprob_time.o \
  boulder/read_ib.o \
  boulder/read_param.o \
  boulder/read_pop_decay.o \
  boulder/read_phys_par.o \
  boulder/read_trans.o \
  boulder/read_yarko.o \
  boulder/read_yarko_tau.o \
  boulder/skip.o \
  boulder/testchar.o \
  boulder/ucrm_ejecta_mass.o \
  boulder/ucrm_partition.o \
  boulder/ucrm_round.o \
  boulder/ucrm_step_ij.o \
  boulder/ucrm_update_arrs.o \
  boulder/update_pops.o \
  boulder/write_ob.o \
  boulder/write_families.o \
  boulder/write_craters.o \
  boulder/yarko_decay.o \
  misc/ran3.o \
  misc/linterp.o \
  misc/loginterp.o \
  simplex/add_ceres_vesta.o \
  simplex/amoeba.o \
  simplex/amotry.o \
  simplex/boulder_call.o \
  simplex/chi2_func.o \
  simplex/gammln.o \
  simplex/gammp.o \
  simplex/gammq.o \
  simplex/gcf.o \
  simplex/gen_sfd_ic.o \
  simplex/gser.o \
  simplex/length.o \
  simplex/read_ascii.o \
  simplex/save_families.o \
  simplex/save_ob.o \
  simplex/sort2.o \
  simplex/srtidx.o \

objc = \

inc = \
  boulder/ucrm3.4.inc \
  boulder/yarko.inc \
  boulder/trans.inc \
  simplex/dependent.inc \
  simplex/simplex.inc \
  simplex/cb_fam.inc \
  simplex/cb_sfd.inc \
  simplex/cb_itmax.inc \

all: main/boulder main/gen_sfd3 main/gen_sfd4 main/gen_sfd5 main/gen_ic main/simplex main/chi2

main/boulder: main/boulder.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

main/gen_sfd3: main/gen_sfd3.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

main/gen_sfd4: main/gen_sfd4.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

main/gen_sfd5: main/gen_sfd5.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

main/gen_ic: main/gen_ic.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

main/simplex: main/simplex.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

main/chi2: main/chi2.f $(obj) $(objc) $(inc)
	$(f77) $(opt) $(obj) $(objc) -o $@ $< $(libs)

$(obj) : %.o:%.f $(inc)
	$(f77) $(opt) -c -o $@ $<

$(objc) : %.o:%.c
	$(cc) $(opt) -c -o $@ $<

clean : FORCE
	rm -f $(obj) $(objc)

FORCE :


