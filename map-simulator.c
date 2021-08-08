// Generate model time series from selected evolution
//Compilation and execution:
//gcc map-simulator.c -lm -o map-simulator; ./map-simulator


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <time.h>

#define N_STIMS 26
#define T_BASE 500
#define N_PERTS 10
#define PERTURB_BIP 5

#define FILEIN1 "params.dat"
#define FILEIN2 "list_expdatafiles.dat"
#define FILEOUT1 "perturb_model_withbl.dat"
#define FILEOUT2 "perturb_model_nobl.dat"

#define PI 3.141592654


//--------------------------------------------------------------------
//parameters
typedef struct {
	double a;
	double b;
	double c;
	double d;
	double alfa1;
	double alfa2;
	double beta1;
	double beta2;
	double gamma1;
	double gamma2;
	double delta1;
	double delta2;
	double epsilon1;
	double epsilon2;
	double dseta1;
	double dseta2;
	double eta1;
	double eta2;
} params_t;



//--------------------------------------------------------------------
void load_params(FILE *in1, params_t *p) {
	char garbage[500];

	fgets(garbage,500,in1);
	sscanf(garbage,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg \
		\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
		&(p->a),&(p->b),&(p->c),&(p->d),
		&(p->alfa1),&(p->beta1),&(p->gamma1),&(p->delta1),
		&(p->epsilon1),&(p->dseta1),&(p->eta1),
		&(p->alfa2),&(p->beta2),&(p->gamma2),&(p->delta2),
		&(p->epsilon2),&(p->dseta2),&(p->eta2));
	return;
}

//--------------------------------------------------------------------
// Numerical simulation
void trial_loop(params_t *params, double *model_asyn, double *model_y, int perturb, double postbaseline) {
	unsigned int i,j;
	double p[N_STIMS],y[N_STIMS],t[N_STIMS],s[N_STIMS];
	double baseline;
	double asyn,f,g,h,k;
	double a,b,c,d,alfa1,alfa2,beta1,beta2,gamma1,gamma2;
	double delta1,delta2,epsilon1,epsilon2;
	double dseta1,dseta2,eta1,eta2;

	a=params->a; b=params->b; c=params->c; d=params->d;
	alfa1=params->alfa1; beta1=params->beta1; gamma1=params->gamma1;
	delta1=params->delta1; epsilon1=params->epsilon1;
	dseta1=params->dseta1; eta1=params->eta1;
	alfa2=params->alfa2; beta2=params->beta2; gamma2=params->gamma2;
	delta2=params->delta2; epsilon2=params->epsilon2;
	dseta2=params->dseta2; eta2=params->eta2;
	
	p[0] = 0;
	y[0] = T_BASE;
	t[0] = T_BASE;
	s[0] = T_BASE;
	baseline = 0;

	// loop for stimuli
	for (i=0; i<N_STIMS-1; i++) {
		if (i==PERTURB_BIP) {
			//the only thing we change "manually" during a perturbation:
			//period or interstimulus interval
			t[i] += perturb;
			baseline = postbaseline;
		}

		asyn = p[i] - (t[i] - s[i]);
		model_asyn[i] = asyn;
		model_y[i] = y[i];

		f = a*(asyn-baseline) + b*(y[i]-t[i])
			+ alfa1*(asyn-baseline)*(asyn-baseline)
			+ beta1*(asyn-baseline)*(y[i]-t[i])
			+ gamma1*(y[i]-t[i])*(y[i]-t[i])
			+ delta1*(asyn-baseline)*(asyn-baseline)*(asyn-baseline)
			+ epsilon1*(asyn-baseline)*(asyn-baseline)*(y[i]-t[i])
			+ dseta1*(asyn-baseline)*(y[i]-t[i])*(y[i]-t[i])
			+ eta1*(y[i]-t[i])*(y[i]-t[i])*(y[i]-t[i])
			+ baseline;

		g = c*(asyn-baseline) + d*(y[i]-t[i])
			+ alfa2*(asyn-baseline)*(asyn-baseline)
			+ beta2*(asyn-baseline)*(y[i]-t[i])
			+ gamma2*(y[i]-t[i])*(y[i]-t[i])
			+ delta2*(asyn-baseline)*(asyn-baseline)*(asyn-baseline)
			+ epsilon2*(asyn-baseline)*(asyn-baseline)*(y[i]-t[i])
			+ dseta2*(asyn-baseline)*(y[i]-t[i])*(y[i]-t[i])
			+ eta2*(y[i]-t[i])*(y[i]-t[i])*(y[i]-t[i])
			+ t[i];

		h = t[i];
		k = t[i];


		p[i+1] = f;
		y[i+1] = g;
		t[i+1] = h;
		s[i+1] = k;
	}
	model_asyn[i] = p[i] - (t[i] - s[i]);
	model_y[i] = y[i];

	return;
}


//--------------------------------------------------------------------
void main(int argc, char *argv[]){
	FILE *in1,*in2,*out1,*out2;
	unsigned int i,j;
	int perturb;
	int pert_sizes_withbl[N_PERTS]={-50,-40,-30,-20,-10,10,20,30,40,50};
	int pert_sizes_nobl[N_PERTS]={-50,-40,-30,-20,-10,10,20,30,40,50};
	double postbaseline,postbl_list[N_PERTS];
	double model_asyn_withbl[N_PERTS][N_STIMS],model_y_withbl[N_PERTS][N_STIMS];
	double model_asyn_nobl[N_PERTS][N_STIMS],model_y_nobl[N_PERTS][N_STIMS];
	double a,b,c,d;
	double alfa1,beta1,gamma1,delta1,epsilon1,dseta1,eta1;
	double alfa2,beta2,gamma2,delta2,epsilon2,dseta2,eta2;
	params_t params;
	char garbage[500],aux_char[256];
	int aux_int;

	// read parameter values from file
	in1 = fopen(FILEIN1,"r");
	load_params(in1,&params);
	fclose(in1);

	// read baseline values
	in2 = fopen(FILEIN2,"r");
	for (j=0; j<N_PERTS; j++) {
		fgets(garbage,500,in2);
		sscanf(garbage,"%s\t%d\t%lf",aux_char,&aux_int,&(postbl_list[j]));
	}
	fclose(in2);


	// simulation (with baseline)
	for (j=0; j<N_PERTS; j++) {
		perturb = pert_sizes_withbl[j];
		postbaseline = postbl_list[j];
		trial_loop(&params,model_asyn_withbl[j],model_y_withbl[j],perturb,postbaseline);
	}
	// save simulation
	out1 = fopen(FILEOUT1,"w");
	for (i=0; i<N_STIMS; i++) {
		fprintf(out1,"%d\t",i-PERTURB_BIP);
		for (j=0; j<N_PERTS-1; j++) {
			fprintf(out1,"%lg\t%lg\t",model_asyn_withbl[j][i],model_y_withbl[j][i]);
		}
		fprintf(out1,"%lg\t%lg\n",model_asyn_withbl[j][i],model_y_withbl[j][i]);
	}
	fclose(out1);


	// simulation (no baseline)
	postbaseline = 0;
	for (j=0; j<N_PERTS; j++) {
		perturb = pert_sizes_nobl[j];
		trial_loop(&params,model_asyn_nobl[j],model_y_nobl[j],perturb,postbaseline);
	}
	// save simulation
	out2 = fopen(FILEOUT2,"w");
	for (i=0; i<N_STIMS; i++) {
		fprintf(out2,"%d\t",i-PERTURB_BIP);
		for (j=0; j<N_PERTS-1; j++) {
			fprintf(out2,"%lg\t%lg\t",model_asyn_nobl[j][i],model_y_nobl[j][i]);
		}
		fprintf(out2,"%lg\t%lg\n",model_asyn_nobl[j][i],model_y_nobl[j][i]);
	}
	fclose(out2);


	return;
}

