// Select one particular run from the genetic algorithm (ranked)
//Compilation and execution:
//gcc map-get-evolution.c -lm -o map-get-evolution; ./map-get-evolution N
//where N is the chosen rank


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <time.h>

#define PERTURB_STEP 5

#define FILEIN1 "evols.dat"
#define FILEIN2 "fitness.dat"
#define FILEOUT1 "evols_sort.dat"
#define FILEOUT2 "params.dat"
#define FILEGNP1 "map-plot_results.gnp"
#define FILEGNP2 "map-fitness.gnp"

#define PI 3.141592654

/*------------------------------------------------------------*/
typedef struct {
	double fitness;
	unsigned int loop;
} par_fit_loop_t;


/*------------------------------------------------------------*/
int greater_fun (const void *aa, const void *bb) {
	const par_fit_loop_t *a_struct = (const par_fit_loop_t *)aa;
	const par_fit_loop_t *b_struct = (const par_fit_loop_t *)bb;
	double a,b;
	a = a_struct->fitness;
	b = b_struct->fitness;

	if (isnan(a) || isnan(b) || isinf(a) || isinf(b)) {
		if ((isnan(a) || isinf(a)) && (isnan(b) || isinf(b))) {
			return 1;
		}
		else if (isnan(a) || isinf(a)) {
			return 1;
		}
		else {
			return -1;
		}
	}
	else {
		return a>=b?-1:1;
	}

}


/*------------------------------------------------------------*/
void main(int argc, char *argv[]){
	FILE *in1,*in2,*out1,*out2,*gnp1,*gnp2;
	unsigned int i,j;
	double *fitness;
	double *a,*b,*c,*d;
	double *alfa1,*beta1,*gamma1,*delta1,*epsilon1,*dseta1,*eta1;
	double *alfa2,*beta2,*gamma2,*delta2,*epsilon2,*dseta2,*eta2;
	double discr;
	double eigenvector1_x,eigenvector1_y;
	double eigenvector2_x,eigenvector2_y;
	double lambda1,lambda2,ang1,ang2;
	double slope1,slope2;
	unsigned int *loop,n_loops,select_rank,select_rank_indx;
	par_fit_loop_t *fitness_sort;
	char garbage[500];

	select_rank = atoi(argv[1]);
	select_rank_indx = select_rank - 1;

	// read parameter values from file
	in1 = fopen(FILEIN1,"r");
	fgets(garbage,500,in1);
	sscanf(garbage,"%d\n",&n_loops);
	loop = (unsigned int*)calloc(n_loops,sizeof(unsigned int));
	fitness = (double*)calloc(n_loops,sizeof(double));
	fitness_sort = (par_fit_loop_t*)calloc(n_loops,sizeof(par_fit_loop_t));
	a = (double*)calloc(n_loops,sizeof(double));
	b = (double*)calloc(n_loops,sizeof(double));
	c = (double*)calloc(n_loops,sizeof(double));
	d = (double*)calloc(n_loops,sizeof(double));
	alfa1 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	beta1 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	gamma1 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	delta1 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: alpha
	epsilon1 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	dseta1= (double*)calloc(n_loops,sizeof(double)); //in the manuscript: beta
	eta1= (double*)calloc(n_loops,sizeof(double)); //in the manuscript: gamma
	alfa2 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: delta
	beta2 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	gamma2 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	delta2 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	epsilon2 = (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	dseta2= (double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0
	eta2=(double*)calloc(n_loops,sizeof(double)); //in the manuscript: 0


	// load evolutions
	for (i = 0; i < n_loops; i++) {
		fgets(garbage,500,in1);
		sscanf(garbage,"%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg \
			\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
			&loop[i],&fitness[i],&a[i],&b[i],&c[i],&d[i],
			&alfa1[i],&beta1[i],&gamma1[i],&delta1[i],&epsilon1[i],
			&dseta1[i],&eta1[i],
			&alfa2[i],&beta2[i],&gamma2[i],&delta2[i],&epsilon2[i],
			&dseta2[i],&eta2[i]);
		fitness_sort[i].fitness = fitness[i];
		fitness_sort[i].loop = loop[i];
	}
	fclose(in1);

	// save all evolutions, sorted
	qsort(fitness_sort,n_loops,sizeof(par_fit_loop_t),greater_fun);
	out1 = fopen(FILEOUT1,"w");
	for (i = 0; i < n_loops; i++) {
		j = fitness_sort[i].loop - 1;
		fprintf(out1,"%d\t%lg\t\%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg \
			\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
			fitness_sort[i].loop,fitness_sort[i].fitness,a[j],b[j],c[j],d[j],
			alfa1[j],beta1[j],gamma1[j],delta1[j],epsilon1[j],
			dseta1[j],eta1[j],
			alfa2[j],beta2[j],delta2[j],gamma2[j],epsilon2[j],
			dseta2[j],eta2[j]);
	}
	fclose(out1);

	// save selected evolution
	i = fitness_sort[select_rank_indx].loop - 1;
	printf("select_rank = %d, fitness = %g, loop = %d\n",
		select_rank,fitness_sort[select_rank_indx].fitness,i+1);
	printf("a=%1.4g, b=%1.4g, alfa1=%1.4g, beta1=%1.4g, gamma1=%1.4g\ndelta1=%1.4g, epsilon1=%1.4g, dseta1=%1.4g, eta1=%1.4g\nc=%1.4g, d=%1.4g, alfa2=%1.4g, beta2=%1.4g, gamma2=%1.4g\ndelta2=%1.4g, epsilon2=%1.4g dseta2=%1.4g, eta2=%1.4g\n",
		a[i],b[i],alfa1[i],beta1[i],gamma1[i],delta1[i],epsilon1[i],
		dseta1[i],eta1[i],
		c[i],d[i],alfa2[i],beta2[i],gamma2[i],delta2[i],epsilon2[i],
		dseta2[i],eta2[i]);

	discr = a[i]*a[i] + 4*b[i]*c[i] - 2*a[i]*d[i] + d[i]*d[i];
	lambda1 = 0.5*(a[i] + d[i] - sqrt(discr));
	lambda2 = 0.5*(a[i] + d[i] + sqrt(discr));
	eigenvector1_x = lambda1 - d[i];
	eigenvector1_y = c[i];
	eigenvector2_x = lambda2 - d[i];
	eigenvector2_y = c[i];
	slope1 = eigenvector1_y/eigenvector1_x;
	slope2 = eigenvector2_y/eigenvector2_x;
	ang1 = (180/PI)*atan(slope1);
	ang2 = (180/PI)*atan(slope2);

	printf("lambda1=%.2g, lambda2=%.2g\n",lambda1,lambda2);
	printf("eigenvectors (degrees): ang1=%.2g, ang2=%.2g\n\n",ang1,ang2);

	out2 = fopen(FILEOUT2,"w");
	fprintf(out2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf \
		\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		a[i],b[i],c[i],d[i],alfa1[i],beta1[i],gamma1[i],delta1[i],
		epsilon1[i],dseta1[i],eta1[i],
		alfa2[i],beta2[i],gamma2[i],delta2[i],epsilon2[i],dseta2[i],eta2[i]);
	fclose(out2);


	// gnuplot script for plotting time series and trajectories
	gnp1 = fopen(FILEGNP1,"w");
	fprintf(gnp1,"#gnuplot script\n\n\nlwidth = 2\npsize = 1\n");
	fprintf(gnp1,"T = %d\na = %g\nb = %g\nc = %g\nd = %g\nalfa1 = %g\nbeta1 = %g\ngamma1 = %g\ndelta1 = %g\nepsilon1 = %g\ndseta1= %g\neta1 = %g\nalfa2 = %g\nbeta2= %g\ngamma2 = %g\ndelta2 = %g\nepsilon2 = %g\ndseta2 = %g\neta2 = %g\npstep = %d\n\n\n",500,a[i],b[i],c[i],d[i],alfa1[i],beta1[i],gamma1[i],delta1[i],epsilon1[i],dseta1[i],eta1[i],alfa2[i],beta2[i],gamma2[i],delta2[i],epsilon2[i],dseta2[i],eta2[i],PERTURB_STEP);
	fprintf(gnp1,"discr = a*a + 4*b*c - 2*a*d + d*d\n");
	fprintf(gnp1,"trace = a + d\n");
	fprintf(gnp1,"lambda1 = 0.5*trace + sqrt(discr)\n");
	fprintf(gnp1,"lambda2 = 0.5*trace - sqrt(discr)\n");
	fprintf(gnp1,"eigenvector1_x = lambda1 - d\n");
	fprintf(gnp1,"eigenvector1_y = c\n");
	fprintf(gnp1,"eigenvector2_x = lambda2 - d\n");
	fprintf(gnp1,"eigenvector2_y = c\n");
	fprintf(gnp1,"slope1 = eigenvector1_y/eigenvector1_x\n");
	fprintf(gnp1,"slope2 = eigenvector2_y/eigenvector2_x\n");
	fprintf(gnp1,"eigenvector1(x) = slope1*x\n");
	fprintf(gnp1,"eigenvector2(x) = slope2*x\n\n\n");

	fprintf(gnp1,"set auto\n");
	fprintf(gnp1,"set size 1,1\n");
	fprintf(gnp1,"set origin 0,0\n");
	fprintf(gnp1,"set multiplot layout 2,3\n");

	fprintf(gnp1,"set xrange [-5:20]\n");
	fprintf(gnp1,"set yrange [-60:60]\n");
	fprintf(gnp1,"set zeroaxis lt 3\n");
	fprintf(gnp1,"set key top right noreverse\n");
	fprintf(gnp1,"set xlabel \"step n\"\n");
	fprintf(gnp1,"set ylabel \"synchrony error e_n (ms)\"\n");
	fprintf(gnp1,"set xlabel \"step n\"\nset ylabel \"synchrony error e_n (ms)\"\n");
	fprintf(gnp1,"plot 'perturb_exp_plus50.dat' u 1:2 t \"exp +/-50\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_exp_minus50.dat' u 1:2 t \"\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:2 t \"fit +/-50\" w lp pt 5 lt 1 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:20 t \"\" w lp pt 5 lt 1 lw 2 ps psize\n\n\n");

	fprintf(gnp1,"set xrange [-5:20]\n");
	fprintf(gnp1,"set yrange [-60:60]\n");
	fprintf(gnp1,"set zeroaxis lt 3\n");
	fprintf(gnp1,"set key top right noreverse\n");
	fprintf(gnp1,"set xlabel \"step n\"\n");
	fprintf(gnp1,"set ylabel \"synchrony error e_n (ms)\"\n");
	fprintf(gnp1,"plot 'perturb_exp_plus40.dat' u 1:2 t \"exp +/-40\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_exp_minus40.dat' u 1:2 t \"\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:4 t \"fit +/-40\" w lp pt 5 lt 1 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:18 t \"\" w lp pt 5 lt 1 lw 2 ps psize\n\n\n");

	fprintf(gnp1,"set xrange [-5:20]\n");
	fprintf(gnp1,"set yrange [-60:60]\n");
	fprintf(gnp1,"set zeroaxis lt 3\n");
	fprintf(gnp1,"set key top right noreverse\n");
	fprintf(gnp1,"set xlabel \"step n\"\n");
	fprintf(gnp1,"set ylabel \"synchrony error e_n (ms)\"\n");
	fprintf(gnp1,"plot 'perturb_exp_plus30.dat' u 1:2 t \"exp +/-30\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_exp_minus30.dat' u 1:2 t \"\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:6 t \"fit +/-30\" w lp pt 5 lt 1 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:16 t \"\" w lp pt 5 lt 1 lw 2 ps psize\n\n\n");

	fprintf(gnp1,"set xrange [-5:20]\n");
	fprintf(gnp1,"set yrange [-60:60]\n");
	fprintf(gnp1,"set zeroaxis lt 3\n");
	fprintf(gnp1,"set key top right noreverse\n");
	fprintf(gnp1,"set xlabel \"step n\"\n");
	fprintf(gnp1,"set ylabel \"synchrony error e_n (ms)\"\n");
	fprintf(gnp1,"plot 'perturb_exp_plus20.dat' u 1:2 t \"exp +/-20\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_exp_minus20.dat' u 1:2 t \"\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:8 t \"fit +/-20\" w lp pt 5 lt 1 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:14 t \"\" w lp pt 5 lt 1 lw 2 ps psize\n\n\n");

	fprintf(gnp1,"set xrange [-5:20]\n");
	fprintf(gnp1,"set yrange [-60:60]\n");
	fprintf(gnp1,"set zeroaxis lt 3\n");
	fprintf(gnp1,"set key top right noreverse\n");
	fprintf(gnp1,"set xlabel \"step n\"\n");
	fprintf(gnp1,"set ylabel \"synchrony error e_n (ms)\"\n");
	fprintf(gnp1,"plot 'perturb_exp_plus10.dat' u 1:2 t \"exp +/-10\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_exp_minus10.dat' u 1:2 t \"\" w lp pt 7 lt 3 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:10 t \"fit +/-10\" w lp pt 5 lt 1 lw 2 ps psize, \\\n");
	fprintf(gnp1,"'perturb_model_withbl.dat' u 1:12 t \"\" w lp pt 5 lt 1 lw 2 ps psize\n\n\n");

	fprintf(gnp1,"set xrange [-60:60]\n");
	fprintf(gnp1,"set yrange [-60:60]\n");
	fprintf(gnp1,"set zeroaxis lt 3\n");
	fprintf(gnp1,"set key top left Left reverse\n");
	fprintf(gnp1,"set xlabel \"e_n (ms)\"\n");
	fprintf(gnp1,"set ylabel \"x_n - T_{post} (ms)\"\n");
	fprintf(gnp1,"plot 'perturb_model_nobl.dat' every ::pstep u 2:($3-(T-50)) t \"fit ({/Symbol D}<0, nobaseline)\" w lp pt 5 lt 1 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 4:($5-(T-40)) t \"\" w lp pt 5 lt 1 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 6:($7-(T-30)) t \"\" w lp pt 5 lt 1 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 8:($9-(T-20)) t \"\" w lp pt 5 lt 1 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 10:($11-(T-10)) t \"\" w lp pt 5 lt 1 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 12:($13-(T+10)) t \"\" w lp pt 5 lt 2 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 14:($15-(T+20)) t \"\" w lp pt 5 lt 2 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 16:($17-(T+30)) t \"\" w lp pt 5 lt 2 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 18:($19-(T+40)) t \"\" w lp pt 5 lt 2 lw 2 ps psize,\\\n");
	fprintf(gnp1,"'perturb_model_nobl.dat' every ::pstep u 20:($21-(T+50)) t \"fit ({/Symbol D}>0, nobaseline)\" w lp pt 5 lt 2 lw 2 ps psize,\\\n");
	fprintf(gnp1,"eigenvector1(x) t \"eigenvectors\" w l lt -1 lw 1,\\\n");
	fprintf(gnp1,"eigenvector2(x) t \"\" w l lt -1 lw 1\n");

	fprintf(gnp1,"unset multiplot\n");
	fprintf(gnp1,"set auto\n");
	fprintf(gnp1,"set xlabel \"\"\n");
	fprintf(gnp1,"set ylabel \"\"\n");
	fclose(gnp1);



	// gnuplot script for plotting fitness
	gnp2 = fopen(FILEGNP2,"w");
	fprintf(gnp2,"#gnuplot script\n\n\n");
	fprintf(gnp2,"set key right bottom\n");
	fprintf(gnp2,"plot 'fitness.dat' u 1 t \"\" w l, \\\n");
	for (j=2; j<=n_loops; j++) {
		fprintf(gnp2,"'fitness.dat' u %u t \"\" w l, \\\n",j);
	}
	fprintf(gnp2,"'fitness.dat' u %u t \"best\" w l lw 4\n",i+1);
	fclose(gnp2);



	free(loop);
	free(fitness);
	free(fitness_sort);
	free(a);
	free(b);
	free(c);
	free(d);
	free(alfa1);
	free(beta1);
	free(gamma1);
	free(delta1);
	free(epsilon1);
	free(dseta1);
	free(eta1);
	free(alfa2);
	free(beta2);
	free(gamma2);
	free(delta2);
	free(epsilon2);
	free(dseta2);
	free(eta2);

	return;
}

