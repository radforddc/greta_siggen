# Makefile for signal generation from PPC detectors
#   - uses .c library codes by Karin Lagergren, heavily modified by David Radford
#
# [-lreadline option required for readline, addhistory...]

CC = gcc
CPP = g++
CFLAGS = -O3 -Wall
RM = rm -f

# common files and headers
mk_signal_files = calc_signal.c cyl_point.c detector_geometry.c fields.c point.c \
                  read_config.c signal_calc_util.c
mk_signal_headers = calc_signal.h cyl_point.h detector_geometry.h fields.h greta_siggen.h \
                    point.h signal_calc_util.h

All: stester greta_fieldgen
# mass drift_time_map

# interactive interface for signal calculation code
stester: $(mk_signal_files) $(mk_signal_headers) signal_tester.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) signal_tester.c -lm -lreadline

greta_fieldgen: fieldgen3d.c field_init.c point.c read_config.c signal_calc_util.c greta_siggen.h point.h
	$(CC) $(CFLAGS) -o $@ fieldgen3d.c field_init.c point.c read_config.c signal_calc_util.c -lm

mass: mass.c read_config.c greta_siggen.h
	$(CC) $(CFLAGS) -o $@ mass.c read_config.c -lm

sgrid: $(mk_signal_files) $(mk_signal_headers) sgrid.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) sgrid.c -lm

FORCE:

clean: 
	$(RM) *.o core* *[~%] *.trace
	$(RM) stester greta_fieldgen mass
