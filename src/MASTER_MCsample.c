#include "mc_globals.h"
#include "../include/marsaglia.h"
#include <errno.h>
#include <stdarg.h>

static void checked_snprintf(char *buffer, size_t size, const char *fmt, ...) {
	va_list args;
	int written;

	va_start(args, fmt);
	written = vsnprintf(buffer, size, fmt, args);
	va_end(args);

	if (written < 0 || (size_t)written >= size) {
		fprintf(stderr, "Formatted path exceeded its destination buffer.\n");
		exit(EXIT_FAILURE);
	}
}

static void ensure_directory(const char *path) {
	if (mkdir(path, 0775) != 0 && errno != EEXIST) {
		fprintf(stderr, "Could not create directory '%s': %s\n", path, strerror(errno));
		exit(EXIT_FAILURE);
	}
}

void generate_evectors() {
	printf("Starting integrated spectral solver...\n");
	double fugacity = 1.0; // Uniform sampling uses eigenvectors at z=1.0
	
	// The integrated power method uses L_Evector and R_Evector_solve
	// We need to provide them as expected by pw_meth_ts_LRvec_fcheck.c
	extern double max_eval_LRvec(double fugacity);
	double calculated_evalue = max_eval_LRvec(fugacity) + 1.0;
	printf("Calculated dominant eigenvalue: %.15f (Expected: %.15f)\n", calculated_evalue, dom_evalue);
	
	// VERIFICATION: Check against archival file if it exists
	char r_filename[100];
	if (ham_check) {
		checked_snprintf(r_filename, sizeof(r_filename), "data/MC_Evectors/R_EvectorHam_TS_L%dM%d.txt", L, M);
	} else {
		checked_snprintf(r_filename, sizeof(r_filename), "data/MC_Evectors/R_Evector_TS_L%dM%d.txt", L, M);
	}

	FILE *r_fp = fopen(r_filename, "r");
	if (r_fp != NULL) {
		printf("Verifying calculated eigenvectors against archival data in %s...\n", r_filename);
		double max_diff = 0.0;
		int diff_count = 0;
		for (int i = 1; i <= (int)num_tspans; i++) {
			double file_val;
			if (fscanf(r_fp, "%lf", &file_val) == 1) {
				double diff = fabs(R_Evector_solve[0][i] - file_val);
				if (diff > max_diff) max_diff = diff;
				if (diff > 1e-9) diff_count++;
			}
		}
		printf("Verification complete. Max difference: %e, Significant differences (>1e-9): %d\n", max_diff, diff_count);
		fclose(r_fp);
	} else {
		printf("No archival file found for verification (%s). Proceeding with calculated values.\n", r_filename);
	}

	// EXPORT: Save calculated eigenvectors for later use/audit
	ensure_directory("data/MC_Evectors");
	char export_fn[100];
	if (ham_check) {
		checked_snprintf(export_fn, sizeof(export_fn), "data/MC_Evectors/R_EvectorHam_TS_L%dM%d.txt", L, M);
	} else {
		checked_snprintf(export_fn, sizeof(export_fn), "data/MC_Evectors/R_Evector_TS_L%dM%d.txt", L, M);
	}
	
	FILE *export_fp = fopen(export_fn, "w");
	if (export_fp != NULL) {
		for (int i = 1; i <= (int)num_tspans; i++) {
			fprintf(export_fp, "%.15f\n", R_Evector_solve[0][i]);
		}
		fclose(export_fp);
		printf("Calculated eigenvectors exported to %s\n", export_fn);
	}
}

void parse_args(int argc, char *argv[]) {
	int opt;
	while ((opt = getopt(argc, argv, "L:M:s:n:r:S:m:")) != -1) {
		switch (opt) {
			case 'L': L = atoi(optarg); break;
			case 'M': M = atoi(optarg); break;
			case 's': totalspan = atoi(optarg); break;
			case 'n': samplesize = atoi(optarg); break;
			case 'r': runnum = atoi(optarg); break;
			case 'S': seednum = (unsigned int)atoi(optarg); break;
			case 'm': 
				mode = atoi(optarg);
				if (mode == 1 || mode == 3) ham_check = 1;
				break;
			default:
				fprintf(stderr, "Usage: %s [-L L] [-M M] [-s totalspan] [-n samplesize] [-r runnum] [-S seednum] [-m mode]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}
}

void build_transfer_matrix_graph() {
	int ordNum[3], side = 0;
	int i, j;

	clock();	/* Start clock to time program */

	ordNum[0] = 1;	/* first edge in section 0 gets numbered 1 (left side of 2-span)	*/
	ordNum[1] = 1;	/* first edge in section 1 gets numbered 1 (right side of 2-span)	*/
	ordNum[2] = 0;	/* to start there are no edges in the 2-span (edges in hinge)		*/

	for (i=1; i<=max_sections; i++){
		current_hinge_span[i] = newhinge(1);		/* free up space of hinge structure */
		first_hinge_span[i] = current_hinge_span[i];	/* make first_hinge_span[i] point to this first one */
		currentendhinge[i] = newendhinge(1);		/* free up space of endhinge structure */
		firstendhinge[i] = currentendhinge[i];		/* make firstendhinge[i] point to this first one */
		sectionkey[i] = 0;
	}
	printf("space created for hinge structures, endhinge structures, and sectionkey initialized\n");

	//initializations
	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			colhingeedges[i][j]=0;
		}
	}
	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			rowhingeedges[i][j]=0;
		}
	}

	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
		}
	}

	//create 2-spans
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){
				printf("2-span: entering at i=%d,j=%d\n", i, j);
				enterhinge(i, j, side, &ordNum, 0);
			}
			alreadyentered[i][j]=1;
		}
	}

	printf("%lu 2-spans recorded\n", num_tspans);
	printf("%lu duplicate 2spans\n", num_duplicate_tspans);

	//generate the endhinges
	num_walks=0;
	for(i=0; i<=vM*vL-1; i++){
		for(j=0; j<=2; j++){
			curstart[j][i]=-1;
			curend[j][i]=-1;
		}
		for(j=0; j<=vM*vL; j++){
			curwalks[j][i]=0;
		}
	}

	for(i=0; i<=M; i++){
		for(j=0; j<=L; j++){
			hingestatus[i][j]=0;
			alreadyentered[i][j]=0;
		}
	}

	for(i=0; i<=M-1; i++){
		for(j=0; j<=L; j++){
			colhingeedges[i][j]=0;
		}
	}
	for(i=0; i<=M; i++){
		for(j=0; j<=L-1; j++){
			rowhingeedges[i][j]=0;
		}
	}

	EndOrdNum[0]=1;
	EndOrdNum[1]=0;

	//create endhinges
	for (i = 0; i <= M; i++) {
		for (j = 0; j <= L; j++) {
			if( !(i==M && j==L) ){
				printf("Endhinge: entering at i=%d,j=%d\n", i, j);
				enterendhinge(i, j, 0);
			}
			alreadyentered[i][j]=1;
		}
	}

	printf("%lu right endhinges recorded\n", tot_right_endhinges);
	printf("%lu left endhinges recorded\n", tot_left_endhinges);
	printf("%lu duplicate right endhinges\n", num_duplicate_right_endhinges);
	printf("%lu duplicate left endhinges\n", num_duplicate_left_endhinges);

	printf("\nDone in %f seconds.\n\n",(double)clock()/(double)CLOCKS_PER_SEC);

	printf("converting hinge_span linked lists to arrays...\n");
	conv_to_array();
	printf("done converting hinge_spans to arrays.\n\n");

	printf("converting endhinge linked lists to arrays...\n");
	conv_endhinges_to_array();
	printf("done converting endhinges to arrays.\n\n");

	unsigned long int totalsecs = 0;
	for (i = 1; i <= (unsigned long int)max_sections; i++) {
		if (num_outsections[i] > 0) {
			totalsecs++;
		}
	}
	printf("Final: Total sections=%lu, total 2spans=%lu\n", totalsecs, num_tspans);
	printf("max_sections was set = %d\n", max_sections);
	printf("max_tspans was set = %d\n", max_tspans);
	printf("There are a total of %lu potential two-span patterns.\n", num_2_spans - 1);
	printf("There are a total of %lu valid two-span patterns.\n", num_tspans);
	printf("Number of duplicate 2-spans generated that were not recorded=%lu\n", dupcounter);
	printf("Dominant Eigenvalue: %.15f\n", dom_evalue);

	if (ham_check) {
		check_reachability();
	}

	printf("Finished generating all 2-spans and endhinges.\n");
}

void run_rejection_sampler() {
	if(totalspan < 2){
		printf("totalspan=%d. This program only works for span>=2. Exitting\n", totalspan);
		exit(1);
	}

	printf("\nNOW SAMPLING: %d samples from L=%d, M=%d, span=%d%s\n", samplesize, L, M, totalspan, (ham_check ? " (Hamiltonian)" : ""));

	const char *file_prefix = (ham_check) ? "MCpolysHam" : "MCpolys";
	ensure_directory("data");
	ensure_directory(output_dir);
	checked_snprintf(filename, sizeof(filename), "%s/%sL%dM%dspan%drun%dnum%lu.txt", output_dir, file_prefix, L, M, totalspan, runnum, filenum);
	fp = fopen(filename, "w");

	if(fp != NULL){
		fprintf(fp, "UofS\n");
		printf("printed UofS in file '%s'\n", filename);
	} else {
		fprintf(stderr, "Could not open '%s' for writing: %s\n", filename, strerror(errno));
		exit(1);
	}

	unsigned long int secnum, sec1, sec2;
	int nth_endhinge, nth_walk;
	int curspan;
	int nth_tspan;

	int chosenLEH;
	int curLEH;
	int chosenREH;

	unsigned long int curSample=0;

	double prob;
	double sumofprobs;
	unsigned long int cur_tspan_num;

	// Always calculate eigenvectors at runtime
	generate_evectors();

	double* R_Evector;
	R_Evector = R_Evector_solve[0]; // Point directly to the calculated global vector

	double maxt_one=0.0;
	double t_one[tot_left_endhinges];
	curLEH=0;
	int i, j, k;
	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_left_endhinges[i]; j++){
			t_one[curLEH]=0;
			for(k=1; k<=num_outsections[i]; k++){
				t_one[curLEH] = t_one[curLEH] + R_Evector[t_nrr[i][k]];
			}
			if(t_one[curLEH]>maxt_one){
				maxt_one = t_one[curLEH];
			}
			curLEH++;
		}
	}
	unsigned long int reject_one=0;

	double maxt_two=0.0;
	double* t_two;
	t_two = (double*)malloc(sizeof(double)*(max_tspans+1));
	if(t_two==NULL){
		fprintf(stderr, "Out of memory");
		exit(0);
	}

	for(i=1; i<=max_sections; i++){
		for(j=1; j<=num_outsections[i]; j++){
			t_two[t_nrr[i][j]] = num_right_endhinges[t_outsection[i][j]] / R_Evector[t_nrr[i][j]];
			if(t_two[t_nrr[i][j]] > maxt_two){
				maxt_two = t_two[t_nrr[i][j]];
			}
		}
	}
	unsigned long int reject_two=0;

	while(curSample <= samplesize-1){
		for(i=0; i<=vM*vL/2-1; i++){
			for(j=0; j<=2; j++){
				built_walks_start[i][j]=-1;
				built_walks_end[i][j]=-1;
			}
			for(j=0; j<=vM*vL*(totalspan+1)-1; j++){
				built_walks_direcs[i][j]=0;
			}
		}

		chosenLEH = floor(ran1real_()*tot_left_endhinges);

		curLEH=0;
		for(secnum=1; secnum<=max_sections; secnum++){
			for(nth_endhinge=1; nth_endhinge<=num_left_endhinges[secnum]; nth_endhinge++){
				if(curLEH==chosenLEH){
					if( ran1real_() < t_one[curLEH] / maxt_one ){
						num_built_walks = Lend_num_walks[secnum][nth_endhinge];
						for(nth_walk=0; nth_walk<=num_built_walks-1; nth_walk++){
							for(i=0; i<=2; i++){
								built_walks_start[nth_walk][i] = Lend_start[secnum][nth_endhinge][i][nth_walk];
								built_walks_end[nth_walk][i] = Lend_end[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=0; i<=vM*vL; i++){
								built_walks_direcs[nth_walk][i] = Lend_walks[secnum][nth_endhinge][i][nth_walk];
							}
							for(i=vM*vL+1; i<=vM*vL*(totalspan+1)-1; i++){
								built_walks_direcs[nth_walk][i]=0;
							}
						}
						curspan=1;
						sec1=secnum;

						prob = ran1real_();
						nth_tspan=1;
						sumofprobs = R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH];
						while(sumofprobs<prob){
							nth_tspan++;
							sumofprobs = sumofprobs + R_Evector[t_nrr[sec1][nth_tspan]] / t_one[curLEH];
						}
						add_to_built_walks(sec1, nth_tspan);
						sec2 = t_outsection[sec1][nth_tspan];
						cur_tspan_num = t_nrr[sec1][nth_tspan];
						curspan++;

						while(curspan < totalspan){
							sumofprobs=0.0;
							for(i=1; i<=num_outsections[sec2]; i++){
								sumofprobs = sumofprobs + R_Evector[t_nrr[sec2][i]] / R_Evector[cur_tspan_num] / dom_evalue;
							}

							if(sumofprobs<0.9 || sumofprobs>1.1){
								printf("second sumofprobs=%f\n", sumofprobs);
								printf("problem with second sumofprobs when leaving tspan number %lu, (made up of %lu and %lu)\n", cur_tspan_num, sec1, sec2);
								printf("This was on sample %lu. Exitting\n", curSample);
								exit(1);
							}

							prob = ran1real_();
							prob = prob * sumofprobs;
							nth_tspan=1;
							sumofprobs = R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
							while(sumofprobs<prob){
								nth_tspan++;
								sumofprobs = sumofprobs + R_Evector[t_nrr[sec2][nth_tspan]] / R_Evector[cur_tspan_num] / dom_evalue;
							}
							sec1=sec2;
							add_to_built_walks(sec1, nth_tspan);
							sec2 = t_outsection[sec1][nth_tspan];
							cur_tspan_num = t_nrr[sec1][nth_tspan];
							curspan++;
						}

						if( ran1real_() < t_two[cur_tspan_num] / maxt_two ){
							chosenREH = 1 + floor(ran1real_()*num_right_endhinges[sec2]);
							add_right_endhinge(sec2, chosenREH);
							curSample++;
							printtofile();
						}
						else{
							reject_two++;
						}
					}
					else{
						reject_one++;
					}
				}
				else if(curLEH>chosenLEH){
					break;
				}
				curLEH++;
			}
		}
	}

	fprintf(fp, "-999\n");
	fclose(fp);
	printf("File %lu contains %lu polys.\n", filenum, filetotal);

	printf("\nSampling Complete.\n");

	printf("\nRESULTS: %lu samples created from L=%d, M=%d, span=%d\n", curSample, L, M, totalspan);
	printf("Rejected at first step %lu times\n", reject_one);
	printf("Rejected at last step %lu times\n", reject_two);
	
	free(t_two);
}

void free_sampler_memory() {
	// Free global dynamic arrays allocated in mc_globals.c
	int i;
	if (built_walks_start) {
		for (i = 0; i < vM * vL / 2 + 1; i++) free(built_walks_start[i]);
		free(built_walks_start);
	}
	if (built_walks_end) {
		for (i = 0; i < vM * vL / 2 + 1; i++) free(built_walks_end[i]);
		free(built_walks_end);
	}
	if (built_walks_direcs) {
		for (i = 0; i < vM * vL / 2 + 1; i++) free(built_walks_direcs[i]);
		free(built_walks_direcs);
	}
	
	// Also free the R_Evector matrix
	for(i = 0; i <= 1; i++) {
		if (R_Evector_solve[i]) free(R_Evector_solve[i]);
		if (L_Evector[i]) free(L_Evector[i]);
	}
}

int main(int argc, char *argv[]) {
	parse_args(argc, argv);

	if (ham_check && ((L + 1) * (M + 1) * (totalspan + 1)) % 2 != 0) {
		fprintf(stderr, "\nError: Impossible to construct a Hamiltonian polygon with an odd number of vertices on a bipartite lattice (L=%d, M=%d, span=%d).\n\n", L, M, totalspan);
		exit(EXIT_FAILURE);
	}

	if (mode == 2 || mode == 3) {
		return run_legacy_2sap_sampler();
	}

	set_system_params();
	allocate_globals();

	if (mode == 1) output_dir = "data/HamSAPs";
	else if (mode == 2) output_dir = "data/2SAPs";
	else if (mode == 3) output_dir = "data/Ham2SAPs";
	else output_dir = "data/SAPs";

	unsigned int seed=seednum;
	initran_(&seed);

	if(system("clear")){}
	if(system("date")){}

	printf("section-section version\n");
	printf("\nThe variables are as follows M: %i L: %i\n", M, L);

	build_transfer_matrix_graph();
	run_rejection_sampler();

	free_sampler_memory();

	if(system("date")){}
	printf("\nProgram Complete.\n\n");

	return 0;
}
