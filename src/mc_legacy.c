#include "mc_globals.h"
#include <errno.h>
#include <stdarg.h>
#include <sys/wait.h>

static int checked_snprintf(char *buffer, size_t size, const char *fmt, ...)
{
	va_list args;
	int written;

	va_start(args, fmt);
	written = vsnprintf(buffer, size, fmt, args);
	va_end(args);

	if (written < 0 || (size_t)written >= size) {
		fprintf(stderr, "Formatted path exceeded its destination buffer.\n");
		return 0;
	}
	return 1;
}

static int run_child(char *const argv[])
{
	pid_t child = fork();
	if (child < 0) {
		fprintf(stderr, "Failed to fork process: %s\n", strerror(errno));
		return 1;
	}
	if (child == 0) {
		execvp(argv[0], argv);
		fprintf(stderr, "Failed to run %s: %s\n", argv[0], strerror(errno));
		_exit(127);
	}

	int status = 0;
	if (waitpid(child, &status, 0) < 0) {
		fprintf(stderr, "Failed waiting for %s: %s\n", argv[0], strerror(errno));
		return 1;
	}
	if (WIFEXITED(status)) return WEXITSTATUS(status);
	if (WIFSIGNALED(status)) {
		fprintf(stderr, "%s exited from signal %d.\n", argv[0], WTERMSIG(status));
	}
	return 1;
}

static void rewrite_legacy_include(FILE *out, char *line)
{
	char *p = line;
	while (*p == ' ' || *p == '\t') p++;
	if (strncmp(p, "#include", 8) != 0) {
		fputs(line, out);
		return;
	}
	const char *roots[] = {
		"../topology/",
		"../sections/",
		"../utils/",
		"../analysis/statistics/"
	};
	const char *prefixes[] = {
		"../../deps/archive/topology/",
		"../../deps/archive/sections/",
		"../../deps/archive/utils/",
		"../../deps/archive/analysis/statistics/"
	};
	for (int i = 0; i < 4; i++) {
		char *root = strstr(line, roots[i]);
		if (root != NULL) {
			char *name = root + strlen(roots[i]);
			char *end = strchr(name, '"');
			if (end != NULL) *end = '\0';
			fprintf(out, "#include \"%s%s\"\n", prefixes[i], name);
			return;
		}
	}
	if (strstr(line, "#include \"pw_meth_") != NULL) {
		char *name = strchr(line, '"') + 1;
		char *end = strchr(name, '"');
		if (end != NULL) *end = '\0';
		fprintf(out, "#include \"../../deps/archive/transfer_matrix/%s\"\n", name);
		return;
	}
	fputs(line, out);
}

static int write_legacy_2sap_source(const char *src_path, const char *dst_path, const char *data_prefix)
{
	FILE *in = fopen(src_path, "r");
	if (in == NULL) {
		fprintf(stderr, "Could not open legacy 2SAP sampler source at %s\n", src_path);
		return 0;
	}
	FILE *out = fopen(dst_path, "w");
	if (out == NULL) {
		fprintf(stderr, "Could not write generated 2SAP sampler source at %s\n", dst_path);
		fclose(in);
		return 0;
	}

	char line[4096];
	int skip_evector_file_read = 0;
	while (fgets(line, sizeof(line), in) != NULL) {
		if (skip_evector_file_read) {
			if (strncmp(line, "/*", 2) == 0 || strncmp(line, "\t/*", 3) == 0) {
				fputs(line, out);
				skip_evector_file_read = 0;
			}
			continue;
		}

		char define_name[128];
		if (sscanf(line, " # define %127s", define_name) == 1 || sscanf(line, "#define %127s", define_name) == 1) {
			if (strcmp(define_name, "L") == 0) {
				fprintf(out, "#define\tL %d\n", L);
				continue;
			}
			if (strcmp(define_name, "M") == 0) {
				fprintf(out, "#define\tM %d\n", M);
				continue;
			}
			if (strcmp(define_name, "totalspan") == 0) {
				fprintf(out, "#define\ttotalspan %d\n", totalspan);
				continue;
			}
			if (strcmp(define_name, "samplesize") == 0) {
				fprintf(out, "#define\tsamplesize %d\n", samplesize);
				continue;
			}
			if (strcmp(define_name, "runnum") == 0) {
				fprintf(out, "#define\trunnum %d\n", runnum);
				continue;
			}
			if (strcmp(define_name, "seednum") == 0) {
				fprintf(out, "#define\tseednum %u\n", seednum);
				continue;
			}
			if (strcmp(define_name, "maxpolys") == 0) {
				fprintf(out, "#define\tmaxpolys %d\n", maxpolys);
				continue;
			}
		}

		if (strstr(line, "#include \"../../include/marsaglia.h\"") != NULL) {
			fputs(line, out);
			fprintf(out, "\n");
			fprintf(out, "#include <sys/stat.h>\n");
			fprintf(out, "double *MC_L_Evector[2];\n");
			fprintf(out, "double *MC_R_Evector[2];\n");
			fprintf(out, "unsigned long int **MC_tspans_edges;\n");
			fprintf(out, "double fval = 0.0;\n");
			fprintf(out, "double max_eval_LRvec(double fugacity);\n");
			continue;
		}

		if (strstr(line, "sprintf(filename2, \"2SAP_R_Evector_TS_L%dM%d.txt\"") != NULL) {
			fprintf(out, "\tsprintf(filename2, \"%sdata/TMresults/2SAP_R_Evector_TS_L%%dM%%d.txt\", L, M);\n", data_prefix);
			continue;
		}

		if (strstr(line, "double* R_Evector;") != NULL) {
			fprintf(out, "\tMC_L_Evector[0] = (double*)malloc(sizeof(double)*(max_tspans+1));\n");
			fprintf(out, "\tMC_L_Evector[1] = (double*)malloc(sizeof(double)*(max_tspans+1));\n");
			fprintf(out, "\tMC_R_Evector[0] = (double*)malloc(sizeof(double)*(max_tspans+1));\n");
			fprintf(out, "\tMC_R_Evector[1] = (double*)malloc(sizeof(double)*(max_tspans+1));\n");
			fprintf(out, "\tif(MC_L_Evector[0]==NULL || MC_L_Evector[1]==NULL || MC_R_Evector[0]==NULL || MC_R_Evector[1]==NULL){\n");
			fprintf(out, "\t\tfprintf(stderr, \"Out of memory\");\n");
			fprintf(out, "\t\texit(0);\n");
			fprintf(out, "\t}\n");
			fprintf(out, "\tMC_tspans_edges = (unsigned long int**)malloc(sizeof(unsigned long int*)*(max_keynum+1));\n");
			fprintf(out, "\tif(MC_tspans_edges==NULL){ fprintf(stderr, \"Out of memory\"); exit(0); }\n");
			fprintf(out, "\tfor(i=1; i<=max_keynum; i++){\n");
			fprintf(out, "\t\tMC_tspans_edges[i] = unsgnlong_vecalloc(1, num_outsections[i]);\n");
			fprintf(out, "\t\tfor(j=1; j<=num_outsections[i]; j++) MC_tspans_edges[i][j] = 0;\n");
			fprintf(out, "\t}\n");
			fprintf(out, "\tdouble calculated_dom_evalue = max_eval_LRvec(1.0) + 1.0;\n");
			fprintf(out, "\tprintf(\"Calculated in-process TS eigenvalue=%%.15f (Expected: %%.15f)\\n\", calculated_dom_evalue, dom_evalue);\n");
			fprintf(out, "\tmkdir(\"data/MC_Evectors\", 0775);\n");
			fprintf(out, "\tchar export_fn[128];\n");
			fprintf(out, "\tsprintf(export_fn, \"data/MC_Evectors/2SAP_R_Evector%s_TS_L%%dM%%d.txt\", L, M);\n", (strstr(src_path, "_Ham.c") != NULL ? "Ham" : ""));
			fprintf(out, "\tFILE *export_fp = fopen(export_fn, \"w\");\n");
			fprintf(out, "\tif (export_fp != NULL) {\n");
			fprintf(out, "\t\tfor (int i = 1; i <= max_tspans; i++) {\n");
			fprintf(out, "\t\t\tfprintf(export_fp, \"%%.15f\\n\", MC_R_Evector[0][i]);\n");
			fprintf(out, "\t\t}\n");
			fprintf(out, "\t\tfclose(export_fp);\n");
			fprintf(out, "\t\tprintf(\"Calculated eigenvectors exported to %%s\\n\", export_fn);\n");
			fprintf(out, "\t}\n");
			fprintf(out, "\tdouble* R_Evector = MC_R_Evector[0];\n");
			skip_evector_file_read = 1;
			continue;
		}
		if (strstr(line, "sprintf(filename2, \"2SAP_R_EvectorHam_TS_L%dM%d.txt\"") != NULL) {
			fprintf(out, "\tsprintf(filename2, \"%sdata/TMresults/2SAP_R_EvectorHam_TS_L%%dM%%d.txt\", L, M);\n", data_prefix);
			continue;
		}
		if (strstr(line, "sprintf(filename, \"MC2SAPsL%dM%dspan%drun%dnum%lu.txt\"") != NULL) {
			fprintf(out, "\t\tsprintf(filename, \"%sdata/2SAPs/MC2SAPsL%%dM%%dspan%%drun%%dnum%%lu.txt\", L, M, totalspan, runnum, filenum);\n", data_prefix);
			continue;
		}
		if (strstr(line, "sprintf(filename, \"MC2SAPsHamL%dM%dspan%drun%dnum%lu.txt\"") != NULL) {
			fprintf(out, "\t\tsprintf(filename, \"%sdata/Ham2SAPs/MC2SAPsHamL%%dM%%dspan%%drun%%dnum%%lu.txt\", L, M, totalspan, runnum, filenum);\n", data_prefix);
			continue;
		}
		char *include_line = line;
		while (*include_line == ' ' || *include_line == '\t') include_line++;
		if (strncmp(include_line, "#include", 8) == 0 &&
			(strstr(line, "\"../") != NULL || strstr(line, "\"pw_meth_") != NULL)) {
			rewrite_legacy_include(out, line);
			continue;
		}
		if (strstr(line, "#include \"Num_section_12V_endhinge_nonordered2.c\"") != NULL) {
		fprintf(out, "#include \"../../deps/archive/sections/Num_section_12V_endhinge_nonordered2.c\"\n");
		}

		fputs(line, out);
	}

	fprintf(out, "\n#define L_Evector MC_L_Evector\n");
	fprintf(out, "#define R_Evector MC_R_Evector\n");
	fprintf(out, "#define tspans_outsection t_outsection\n");
	fprintf(out, "#define tspans_nrr t_nrr\n");
	fprintf(out, "#define tspans_edges MC_tspans_edges\n");
	if (strstr(src_path, "_Ham.c") != NULL) {
		fprintf(out, "#include \"../../deps/archive/transfer_matrix/pw_meth_ts_LRvec_fcheck_2SAP_HAM.c\"\n");
	} else {
		fprintf(out, "#include \"../../deps/archive/transfer_matrix/pw_meth_ts_LRvec_fcheck_2SAP.c\"\n");
	}
	fprintf(out, "#undef tspans_edges\n");
	fprintf(out, "#undef tspans_nrr\n");
	fprintf(out, "#undef tspans_outsection\n");
	fprintf(out, "#undef L_Evector\n");
	fprintf(out, "#undef R_Evector\n");

	fclose(out);
	fclose(in);
	return 1;
}

int run_legacy_2sap_sampler(void)
{
	int ham_mode = (mode == 3);
	const char *legacy_src = ham_mode ? "deps/archive/monte_carlo/2SAP_MCsample_Ham.c" : "deps/archive/monte_carlo/2SAP_MCsample.c";
	FILE *in = fopen(legacy_src, "r");
	if (!in) {
		fprintf(stderr, "Fatal: Could not open legacy source %s\n", legacy_src);
		exit(1);
	}
	fclose(in);

	const char *prefix = "";
	char data_dir[1024];
	char mc_dir[1024];
	char generated_src[1024];
	char exe_path[1024];

	if (!checked_snprintf(data_dir, sizeof(data_dir), "%sdata", prefix)) return 1;
	if (!checked_snprintf(mc_dir, sizeof(mc_dir), "%sbuild/output", prefix)) return 1;
	mkdir(data_dir, 0775);
	if (!checked_snprintf(data_dir, sizeof(data_dir), "%sdata/2SAPs", prefix)) return 1;
	mkdir(data_dir, 0775);
	if (!checked_snprintf(data_dir, sizeof(data_dir), "%sdata/Ham2SAPs", prefix)) return 1;
	mkdir(data_dir, 0775);
	mkdir(mc_dir, 0775);
	if (!checked_snprintf(generated_src, sizeof(generated_src), "%s/%s.generated.c", mc_dir, ham_mode ? "2SAP_MCsample_Ham" : "2SAP_MCsample")) return 1;

	if (!write_legacy_2sap_source(legacy_src, generated_src, prefix)) {
		return 1;
	}

	if (!checked_snprintf(exe_path, sizeof(exe_path), "%s/%s_L%dM%d", mc_dir, ham_mode ? "2sap_sampler_ham" : "2sap_sampler", L, M)) return 1;
	char *compile_args[] = {"gcc", "-O3", "-Wno-unused-result", "-o", exe_path, generated_src, "-lm", NULL};
	printf("Compiling legacy %s2SAP sampler: gcc -O3 -Wno-unused-result -o %s %s -lm\n", ham_mode ? "Hamiltonian " : "", exe_path, generated_src);
	if (run_child(compile_args) != 0) {
		fprintf(stderr, "Failed to compile generated 2SAP sampler.\n");
		return 1;
	}

	printf("Running legacy %s2SAP sampler with L=%d M=%d span=%d samples=%d run=%d seed=%u\n",
		ham_mode ? "Hamiltonian " : "", L, M, totalspan, samplesize, runnum, seednum);
	char *run_args[] = {exe_path, NULL};
	return run_child(run_args);
}
