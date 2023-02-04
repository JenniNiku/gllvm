all:	gllvmVA.so	gllvmLA.so
		# Comment here preserves the prior tab
gllvmVA.so:	gllvmVA.cpp
		"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript"	-e	"TMB::compile('gllvmVA.cpp', safebounds = FALSE, safeunload = FALSE)"
gllvmLA.so:	gllvmLA.cpp
		"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript"	-e	"TMB::compile('gllvmLA.cpp', safebounds = FALSE, safeunload = FALSE)"

clean:
		rm	-rf	*o