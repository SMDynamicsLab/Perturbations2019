// Genetic algorithm for fitting model to experimental time series.
// Compilation and execution:
// gcc map-fitting.c -lm -I/usr/local/include -L/usr/local/lib -lgaul -lgaul_util -o map-fitting; ./map-fitting


#include <gaul.h>
#include <stdio.h>

#define N_LOOPS 200
#define N_GENERATIONS 200
#define POP_SIZE 2000
#define CROSSOVER_RATE 0.9
#define MUTATION_RATE 0.1
#define LINEAR_INIT_RANGE 1.0		//parameter initialization
#define QUAD_INIT_RANGE 0.01
#define CUBIC_INIT_RANGE 0.0001
#define PERTURB_BIP 5

#define COEFF_X_1 1
#define COEFF_Y_1 1
#define COEFF_X_2 1
#define COEFF_Y_2 1

#define COEFF_X2_1 0
#define COEFF_XY_1 0
#define COEFF_Y2_1 0
#define COEFF_X2_2 1	//in the manuscript: delta
#define COEFF_XY_2 0
#define COEFF_Y2_2 0

#define COEFF_X3_1 1	//in the manuscript: alpha
#define COEFF_X2Y_1 0
#define COEFF_XY2_1 1	//in the manuscript: beta
#define COEFF_Y3_1 1	//in the manuscript: gamma
#define COEFF_X3_2 0
#define COEFF_X2Y_2 0
#define COEFF_XY2_2 0
#define COEFF_Y3_2 0


#define PENALTY 100000
#define FILEIN_LIST "list_expdatafiles.dat"
#define FILEOUT_HYPERPARAMS "hyperparams.dat"
#define FILEOUT_EVOLUTIONS "evols.dat"
#define FILEOUT_FITNESS "fitness.dat"



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
//training data
typedef struct {
	int *num_data;
	int max_data;
	int n_datafiles;
	int *perturb;
	double *postbaseline;
	double *x;
	double **y;
	unsigned int **weight;
	unsigned int n_params;
	float **fitness;
	unsigned int loop;
} exp_data_t;


//--------------------------------------------------------------------
//compute eigenvalues and eigenvectors
void eigen(float a, float b, float c, float d) {
	float discr;
	float real,imag;
	float ang1,ang2;
	float lambda1,lambda2;

	//check whether eigenvalues are real or complex
	discr = a*a + 4*b*c - 2*a*d + d*d;
	if (discr < 0) {
		printf("eigenvalues are complex.\n");
		real = 0.5*(a + d);
		imag = 0.5*sqrt(-discr);
		printf("lambda1,2 = %.2g +/- %.2g i\n",real,imag);
		//rotation
		if (c > 0)
			printf("rotation counterclockwise.\n");
		else
			printf("rotation clockwise!\n");
	}
	else {
		printf("eigenvalues are real:");
		lambda1 = 0.5*(a + d - sqrt(discr));
		lambda2 = 0.5*(a + d + sqrt(discr));
		printf("lambda1=%.2g, lambda2=%.2g\n",lambda1,lambda2);
		ang1 = (180/PI)*atan(1/((a - d - sqrt(discr))/(2*c)));
		ang2 = (180/PI)*atan(1/((a - d + sqrt(discr))/(2*c)));
		printf("eigenvectors (degrees):");
		printf("ang1=%.2g, ang2=%.2g\n",ang1,ang2);
	}
	return;
}


//--------------------------------------------------------------------
//translation between gene number and parameter name
//it only translates parameters defined in the chromosome; all other
//parameters are set to zero
void param_translation(entity *entity, unsigned int n_params, params_t *p) {
	unsigned int i=0;

	if (COEFF_X_1==1 && i<n_params) {
		p->a = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y_1==1 && i<n_params) {
		p->b = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X_2==1 && i<n_params) {
		p->c = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y_2==1 && i<n_params) {
		p->d = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2_1==1 && i<n_params) {
		p->alfa1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY_1==1 && i<n_params) {
		p->beta1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y2_1==1 && i<n_params) {
		p->gamma1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2_2==1 && i<n_params) {
		p->alfa2 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY_2==1 && i<n_params) {
		p->beta2 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y2_2==1 && i<n_params) {
		p->gamma2 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X3_1==1 && i<n_params) {
		p->delta1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2Y_1==1 && i<n_params) {
		p->epsilon1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY2_1==1 && i<n_params) {
		p->dseta1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y3_1==1 && i<n_params) {
		p->eta1 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X3_2==1 && i<n_params) {
		p->delta2 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2Y_2==1 && i<n_params) {
		p->epsilon2 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY2_2==1 && i<n_params) {
		p->dseta2 = ((double *)entity->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y3_2==1 && i<n_params) {
		p->eta2 = ((double *)entity->chromosome[0])[i];
		i++;
	}

	return;
}


//--------------------------------------------------------------------
//model
void map2D(population *pop, entity *entity, params_t *params, double *model, int j) {
	FILE		*in;
	int			i;
	double		*p,*y,*t,*s;
	double		f,g,h,k;
	int			tau_base,perturb;
	exp_data_t	*data;
	double		baseline,postbaseline,aux1,aux2;
	char		garbage[200];
	double		a,b,c,d,alfa1,alfa2,beta1,beta2,gamma1,gamma2;
	double		delta1,delta2,epsilon1,epsilon2;
	double		dseta1,dseta2,eta1,eta2;

	data = (exp_data_t *)pop->data;
	p = (double*)calloc(data->num_data[j],sizeof(double));
	y = (double*)calloc(data->num_data[j],sizeof(double));
	t = (double*)calloc(data->num_data[j],sizeof(double));
	s = (double*)calloc(data->num_data[j],sizeof(double));
	a=params->a; b=params->b; c=params->c; d=params->d;
	alfa1=params->alfa1; beta1=params->beta1; gamma1=params->gamma1;
	delta1=params->delta1; epsilon1=params->epsilon1;
	dseta1=params->dseta1; eta1=params->eta1;
	alfa2=params->alfa2; beta2=params->beta2; gamma2=params->gamma2;
	delta2=params->delta2; epsilon2=params->epsilon2;
	dseta2=params->dseta2; eta2=params->eta2;


	tau_base = 500;
	perturb = data->perturb[j];
	baseline = 0;
	postbaseline = data->postbaseline[j];
	p[0] = 0;
	y[0] = tau_base;
	t[0] = tau_base;
	s[0] = tau_base;


	for (i=0; i<data->num_data[j]-1; i++) {
		if (i==PERTURB_BIP) {
			//the only thing we change "manually" during a perturbation:
			//period or interstimulus interval
			t[i] += perturb;
			baseline = postbaseline;
		}

		model[i] = p[i] - (t[i] - s[i]);

		//parameters not subject to fitting are set to zero
		f = a*(model[i]-baseline) + b*(y[i]-t[i])
			+ alfa1*(model[i]-baseline)*(model[i]-baseline)
			+ beta1*(model[i]-baseline)*(y[i]-t[i])
			+ gamma1*(y[i]-t[i])*(y[i]-t[i])
			+ delta1*(model[i]-baseline)*(model[i]-baseline)*(model[i]-baseline)
			+ epsilon1*(model[i]-baseline)*(model[i]-baseline)*(y[i]-t[i])
			+ dseta1*(model[i]-baseline)*(y[i]-t[i])*(y[i]-t[i])
			+ eta1*(y[i]-t[i])*(y[i]-t[i])*(y[i]-t[i])
			+ baseline;

		g = c*(model[i]-baseline) + d*(y[i]-t[i])
			+ alfa2*(model[i]-baseline)*(model[i]-baseline)
			+ beta2*(model[i]-baseline)*(y[i]-t[i])
			+ gamma2*(y[i]-t[i])*(y[i]-t[i])
			+ delta2*(model[i]-baseline)*(model[i]-baseline)*(model[i]-baseline)
			+ epsilon2*(model[i]-baseline)*(model[i]-baseline)*(y[i]-t[i])
			+ dseta2*(model[i]-baseline)*(y[i]-t[i])*(y[i]-t[i])
			+ eta2*(y[i]-t[i])*(y[i]-t[i])*(y[i]-t[i])
			+ t[i];

		h = t[i];
		k = t[i];


		p[i+1] = f;
		y[i+1] = g;
		t[i+1] = h;
		s[i+1] = k;
	}
	model[i] = p[i] - (t[i] - s[i]);

	free(p);
	free(y);
	free(t);
	free(s);

	return;
}


//--------------------------------------------------------------------
//fitness function
boolean fitting_score(population *pop, entity *entity) {
	int			i,j;
	double		fitness,score=0,penalty=0;
	double		**model;
	exp_data_t	*data;
	params_t	params = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double		a,b,c,d;
	double		lambda1,lambda2,realpart,ang1,ang2,discriminant;
	FILE		*out1;

	entity->fitness = 0;
	data = (exp_data_t *)pop->data;
	model = (double**)calloc(data->n_datafiles,sizeof(double*));

	//copy parameter values from chromosome to params struct for easier manip
	param_translation(entity, data->n_params, &params);


	//fitness weights
	for (j=0; j<data->n_datafiles; j++) {
		for (i=0; i<data->num_data[0]; i++) {
			(data->weight[j])[i] = 1;
		}
	}


	//run model for each input datafile
	for (j=0; j<data->n_datafiles; j++) {
		model[j] = (double*)calloc(data->num_data[j],sizeof(double));
		map2D(pop,entity,&params,model[j],j);
	}

	//compute fitness
	for (j=0; j<data->n_datafiles; j++) {
		for (i=0; i<data->num_data[j]; i++) {
			score += (data->weight[j])[i]*SQU((data->y[j])[i] - model[j][i]);
		}
	}

	//constraint penalization
	a = params.a;
	b = params.b;
	c = params.c;
	d = params.d;

	discriminant = a*a + 4*b*c - 2*a*d + d*d;
	if (discriminant >= 0) {
		// eigenvalues are real
		lambda1 = 0.5*(a + d + sqrt(discriminant));
		lambda2 = 0.5*(a + d - sqrt(discriminant));
		ang1 = (180/PI)*atan(1/((a - d - sqrt(discriminant))/(2*c)));
		ang2 = (180/PI)*atan(1/((a - d + sqrt(discriminant))/(2*c)));

		penalty += (lambda1<0||lambda1>=1)||(lambda2<0||lambda2>=1)?PENALTY:0;
	}
	else {
		// eigenvalues are complex
		penalty += PENALTY;
	}


	fitness = -(sqrt(score/(float)(data->n_datafiles*data->num_data[0])) + penalty);  
	if (isnan(fitness) || abs(isinf(fitness))) {
		// keep it finite
		fitness = -1E10;
	}
	entity->fitness = fitness;

	for (j=0; j<data->n_datafiles; j++)
		free(model[j]);
	free(model);
	return TRUE;
}


//--------------------------------------------------------------------
//generation callback
//user-defined function
boolean fitting_generation_callback(int generation, population *pop) {
	exp_data_t	*data;

	data = (exp_data_t *)pop->data;
	(data->fitness[data->loop-1])[generation] = pop->entity_iarray[0]->fitness;

	return TRUE;
}


//--------------------------------------------------------------------
//initialize genetic data
boolean fitting_seed(population *pop, entity *adam) {

	unsigned int i=0,n_params;
	exp_data_t	*data;
	float coeff1=2*LINEAR_INIT_RANGE;
	float coeff2=2*QUAD_INIT_RANGE;
	float coeff3=2*CUBIC_INIT_RANGE; 

	// Checks.
	if (!pop) die("Null pointer to population structure passed.");
	if (!adam) die("Null pointer to entity structure passed.");

	data = (exp_data_t *)pop->data;
	n_params = data->n_params;
	
	if (COEFF_X_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff1*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_Y_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff1*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff1*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_Y_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff1*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X2_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff2*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_XY_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff2*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_Y2_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff2*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X2_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff2*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_XY_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff2*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_Y2_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff2*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X3_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X2Y_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_XY2_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_Y3_1==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X3_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_X2Y_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_XY2_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}
	if (COEFF_Y3_2==1 && i<n_params) {
		((double *)adam->chromosome[0])[i] = coeff3*(random_double(1.0)-0.5);
		i++;
	}

	return TRUE;
}


//--------------------------------------------------------------------
//read training data from file
void get_data(exp_data_t *data, char *infile) {
	unsigned int j;
	int line_count=0;					// Number of lines read from datafile
	int file_count;
	char buffer[MAX_LINE_LEN], *line;	// Buffer for input.
	char datafilename[256];
	FILE *listfile,*datafile;

	if (!data) die("Null pointer to data structure passed.");

	data->n_params = COEFF_X_1 + COEFF_Y_1 + COEFF_X_2 + COEFF_Y_2 +
		COEFF_X2_1 + COEFF_XY_1 + COEFF_Y2_1 + COEFF_X2_2 + COEFF_XY_2 +
		COEFF_Y2_2 + COEFF_X3_1 + COEFF_X2Y_1 + COEFF_XY2_1 + COEFF_Y3_1 +
		COEFF_X3_2 + COEFF_X2Y_2 + COEFF_XY2_2 + COEFF_Y3_2;

	//get number of input datafiles
	listfile = fopen(FILEIN_LIST,"r");
	if (!listfile) die("No input listfile.");
	data->n_datafiles = 0;
	while (fgets(buffer,MAX_LINE_LEN,listfile) != NULL)
		data->n_datafiles++;
	fclose(listfile);


	data->y = (double**)calloc(data->n_datafiles,sizeof(double*));
	data->perturb = (int*)calloc(data->n_datafiles,sizeof(int));
	data->postbaseline = (double*)calloc(data->n_datafiles,sizeof(double));
	data->num_data = (int*)calloc(data->n_datafiles,sizeof(int));
	data->weight = (unsigned int**)calloc(data->n_datafiles,sizeof(unsigned int*));
	//get name of input datafiles
	listfile = fopen(FILEIN_LIST,"r");
	file_count = 0;
	while (!feof(listfile) && fgets(buffer,MAX_LINE_LEN,listfile)!=NULL) {
		line = buffer;
		// Skip leading whitespace.
		while (*line == ' ' || *line == '\t') line++;
		// Ignore commented or empty lines
		if (*line == '#' || *line == '!' || *line == '\n') {}

		sscanf(line,"%s\t%d\t%lf",datafilename,&(data->perturb[file_count]),
			&(data->postbaseline[file_count]));
		datafile = fopen(datafilename,"r");
		if (!datafile) die("No input datafile.");


// Read lines.  Each specifies one x,y pair except those starting with '#'
// or '!' which are comment lines and are ignored.  Don't bother parsing
// blank lines either.

		data->num_data[file_count] = 0;
		data->max_data = 0;
		while (!feof(datafile) && fgets(buffer,MAX_LINE_LEN,datafile)!=NULL) {
			line = buffer;

			// Skip leading whitespace.
			while (*line == ' ' || *line == '\t') line++;
			// Ignore this line
			if (*line == '#' || *line == '!' || *line == '\n') {}
			else {
				// Ensure sufficient memory is available.
				if (data->num_data[file_count] == data->max_data) {
					data->max_data += 256;
					data->x = s_realloc(data->x, sizeof(double)*data->max_data);
					data->y[file_count] = s_realloc(data->y[file_count],sizeof(double)*data->max_data);
				}
				sscanf(line, "%lf %lf", &(data->x[data->num_data[file_count]]),
					&((data->y[file_count])[data->num_data[file_count]]));
				data->num_data[file_count]++;
			}
			line_count++;
		}
		for (j=0; j<data->n_datafiles; j++)
			data->weight[j] = (unsigned int*)calloc(data->num_data[j],sizeof(unsigned int));

		file_count++;
		fclose(datafile);
	}
	fclose(listfile);
	return;
}


//--------------------------------------------------------------------
//save parameters to file
void param_save(population *pop, unsigned int loop) {
	unsigned int i,n_params;
	exp_data_t *data;
	double fitness;
	double a=0,b=0,c=0,d=0;
	double alfa1=0,alfa2=0,beta1=0,beta2=0,gamma1=0,gamma2=0;
	double delta1=0,delta2=0,epsilon1=0,epsilon2=0;
	double dseta1=0,dseta2=0,eta1=0,eta2=0;
	FILE *out;


	data = (exp_data_t *)pop->data;
	n_params = data->n_params;

	i = 0;
	if (COEFF_X_1==1 && i<n_params) {
		a=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y_1==1 && i<n_params) {
		b=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X_2==1 && i<n_params) {
		c=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y_2==1 && i<n_params) {
		d=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2_1==1 && i<n_params) {
		alfa1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY_1==1 && i<n_params) {
		beta1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y2_1==1 && i<n_params) {
		gamma1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2_2==1 && i<n_params) {
		alfa2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY_2==1 && i<n_params) {
		beta2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y2_2==1 && i<n_params) {
		gamma2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X3_1==1 && i<n_params) {
		delta1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2Y_1==1 && i<n_params) {
		epsilon1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY2_1==1 && i<n_params) {
		dseta1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y3_1==1 && i<n_params) {
		eta1=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X3_2==1 && i<n_params) {
		delta2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_X2Y_2==1 && i<n_params) {
		epsilon2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_XY2_2==1 && i<n_params) {
		dseta2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}
	if (COEFF_Y3_2==1 && i<n_params) {
		eta2=((double *)ga_get_entity_from_rank(pop,0)->chromosome[0])[i];
		i++;
	}

	out = fopen(FILEOUT_EVOLUTIONS,"a");
	fitness = ga_get_entity_from_rank(pop,0)->fitness;
	fprintf(out,"%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g \
			\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			loop,fitness,a,b,c,d,alfa1,beta1,gamma1,delta1,epsilon1,dseta1,eta1,
			alfa2,beta2,gamma2,delta2,epsilon2,dseta2,eta2);
	fclose(out);

	return;
}


//--------------------------------------------------------------------
//save hyperparameters to file
void hyperparam_save(population *pop) {
	FILE *out;
	unsigned int i,j;
	exp_data_t *data;

	data = (exp_data_t *)pop->data;

	out = fopen(FILEOUT_HYPERPARAMS,"w");
	fprintf(out,"N_LOOPS %u\n",N_LOOPS);
	fprintf(out,"N_GENERATIONS %u\n",N_GENERATIONS);
	fprintf(out,"POP_SIZE %u\n",POP_SIZE);
	fprintf(out,"CROSSOVER_RATE %f\n",CROSSOVER_RATE);
	fprintf(out,"MUTATION_RATE %f\n",MUTATION_RATE);
	fprintf(out,"LINEAR_INIT_RANGE %f\n",LINEAR_INIT_RANGE);
	fprintf(out,"QUAD_INIT_RANGE %f\n",QUAD_INIT_RANGE);
	fprintf(out,"CUBIC_INIT_RANGE %f\n",CUBIC_INIT_RANGE);
	fprintf(out,"PERTURB_BIP %u\n",PERTURB_BIP);
	fprintf(out,"COEFF_X_1 %u\n",COEFF_X_1);
	fprintf(out,"COEFF_Y_1 %u\n",COEFF_Y_1);
	fprintf(out,"COEFF_X_2 %u\n",COEFF_X_2);
	fprintf(out,"COEFF_Y_2 %u\n",COEFF_Y_2);
	fprintf(out,"COEFF_X2_1 %u\n",COEFF_X2_1);
	fprintf(out,"COEFF_XY_1 %u\n",COEFF_XY_1);
	fprintf(out,"COEFF_Y2_1 %u\n",COEFF_Y2_1);
	fprintf(out,"COEFF_X2_2 %u\n",COEFF_X2_2);
	fprintf(out,"COEFF_XY_2 %u\n",COEFF_XY_2);
	fprintf(out,"COEFF_Y2_2 %u\n",COEFF_Y2_2);
	fprintf(out,"COEFF_X3_1 %u\n",COEFF_X3_1);
	fprintf(out,"COEFF_X2Y_1 %u\n",COEFF_X2Y_1);
	fprintf(out,"COEFF_XY2_1 %u\n",COEFF_XY2_1);
	fprintf(out,"COEFF_Y3_1 %u\n",COEFF_Y3_1);
	fprintf(out,"COEFF_X3_2 %u\n",COEFF_X3_2);
	fprintf(out,"COEFF_X2Y_2 %u\n",COEFF_X2Y_2);
	fprintf(out,"COEFF_XY2_2 %u\n",COEFF_XY2_2);
	fprintf(out,"COEFF_Y3_2 %u\n",COEFF_Y3_2);
	for (i=0; i<data->num_data[0]; i++)
		fprintf(out,"%u\t%u\n",(data->weight[0])[i],(data->weight[1])[i]);

	fclose(out);
	return;
}


//--------------------------------------------------------------------
//save fitness, column wise
void fitness_save(population *pop) {
	FILE *out;
	unsigned int loop,i;
	exp_data_t *data;

	data = (exp_data_t *)pop->data;

	out = fopen(FILEOUT_FITNESS,"w");
	fprintf(out,"#loop");
	for(loop=1; loop<N_LOOPS; loop++)
		fprintf(out,"%u\tloop",loop);
	fprintf(out,"%u\n",N_LOOPS);
	for (i=0; i<N_GENERATIONS; i++) {
		for(loop=1; loop<N_LOOPS; loop++)
			fprintf(out,"%f\t",(data->fitness[loop-1])[i]);
		fprintf(out,"%f\n",(data->fitness[N_LOOPS-1])[i]);
	}
	fclose(out);
	return;
}


//--------------------------------------------------------------------
//main function
int main(int argc, char **argv) {
	unsigned int loop;
	population *pop;			/* Population of solutions. */
	exp_data_t data={NULL,0,0,NULL,NULL,NULL,0,NULL,0};	/* Training data. */
	FILE *out;
	

	out = fopen(FILEOUT_EVOLUTIONS,"w");
	fprintf(out,"%d\n",N_LOOPS);
	fclose(out);

	// load training data
	get_data(&data,FILEIN_LIST);
	data.fitness = (float**)calloc(N_LOOPS,sizeof(float*));

	random_seed(12); 
	for(loop = 1; loop <= N_LOOPS; loop++) {
		data.fitness[loop-1] = (float*)calloc(N_GENERATIONS,sizeof(float));
		data.loop = loop;
		pop = ga_genesis_double(
			POP_SIZE,			/* const int              population_size */
			1,					/* const int              num_chromo */
			data.n_params,		/* const int              len_chromo */
			fitting_generation_callback,/* GAgeneration_hook generation_hook */
			NULL,	/* GAiteration_hook       iteration_hook */
			NULL,	/* GAdata_destructor      data_destructor */
			NULL,	/* GAdata_ref_incrementor data_ref_incrementor */
			fitting_score,		/* GAevaluate             evaluate */
			fitting_seed,		/* GAseed                 seed */
			NULL,				/* GAadapt                adapt */
			ga_select_one_aggressive,
			ga_select_two_aggressive,
			ga_mutate_double_singlepoint_drift,	/* GAmutate   mutate */
			ga_crossover_double_singlepoints,	/* GAcrossover    crossover */
			NULL,				/* GAreplace              replace */
			&data				/* vpointer	User data */
		);

		printf("\rloop = %.3d/%d",loop,N_LOOPS);
		fflush(stdout);

		ga_population_set_parameters(
			pop,				/* population      *pop */
			GA_SCHEME_DARWIN,		/* const ga_scheme_type     scheme */
			GA_ELITISM_ONE_PARENT_SURVIVES,	/* const ga_elitism_type  elitism */
			CROSSOVER_RATE,				/* double  crossover */
			MUTATION_RATE,				/* double  mutation */
			0.0      		        /* double  migration */
		);


		//rutina principal
		ga_evolution_threaded(pop,N_GENERATIONS);

		//save parameters to file
		param_save(pop,loop);

		//save hyperparameters to file
		if (loop==1)
			hyperparam_save(pop);

		//end
		ga_extinction(pop);
	}
	//save fitness, column wise
	fitness_save(pop);

	printf("\n");
	exit(EXIT_SUCCESS);
}

