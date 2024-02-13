#include "gkyl_read_radiation.h"


struct all_radiation_states* gkyl_read_rad_fit_params(){

  char *filepath="./data/adas/radiation_fit_parameters.txt";
  FILE *fptr = fopen(filepath,"r");
  if (fptr == NULL){
    printf("Error opening radiation fit file\n");
    exit(EXIT_FAILURE);
  }
  
  // Read header: Max z of elements, number of elements
  int number_elements, max_atomic_number;
  read_two_numbers(fptr, &max_atomic_number, &number_elements);
  int max_charge_state = max_atomic_number;

  struct all_radiation_states *rad_data = malloc(sizeof(struct all_radiation_states));
  rad_data->max_atomic_number = max_atomic_number;
  
  /* allocate pointers for each species. Unphysical states of charge_state>atomic_number
   * are allocated but not used for simplicity 
   */
  rad_data->all_states = (struct radiating_state*)malloc(max_atomic_number*max_charge_state*sizeof(struct radiating_state));
  for (int i=0; i<max_atomic_number; i++){  
    for (int j=0; j<max_charge_state; j++){
      rad_data->all_states[i*max_charge_state+j].atomic_number = i+1;
      rad_data->all_states[i*max_charge_state+j].charge_state = j;
      rad_data->all_states[i*max_charge_state+j].state_exists = false;
      rad_data->all_states[i*max_charge_state+j].number_of_densities = 0;
      rad_data->all_states[i*max_charge_state+j].rad_fits = NULL;
      rad_data->all_states[i*max_charge_state+j].electron_densities = NULL;
    }
  }
  int atomic_number;
  int num_of_charge_states;
  int density_intervals;
  int charge_state;
  for (int i=0; i<number_elements; i++) {
    // For each element, read atomic number and # of charge states
    num_of_charge_states = 0;
    while (num_of_charge_states == 0)
      read_two_numbers(fptr, &atomic_number, &num_of_charge_states);
    
    atomic_number = atomic_number-1;  // convert to based 0 array
    for (int j=0; j<num_of_charge_states; j++) {
      // For each charge state, read # of density intervals
      read_two_numbers(fptr, &charge_state, &density_intervals);
      charge_state = charge_state-1;  // convert to based 0 array
      if (density_intervals>0) {
	int index = atomic_number*max_charge_state+charge_state;
	rad_data->all_states[index].number_of_densities = density_intervals;
	rad_data->all_states[index].state_exists = true;
	rad_data->all_states[index].electron_densities = (double*)malloc(density_intervals*sizeof(double));
	rad_data->all_states[index].rad_fits = (struct rad_fit_parameters*)malloc(density_intervals*sizeof(struct rad_fit_parameters));
	char buf[1000*15];
	for (int k=0; k<density_intervals; k++) {
	  if (fgets(buf,1000*15, fptr) == NULL) {
	    rad_data->all_states[index].number_of_densities = rad_data->all_states[index].number_of_densities-1;
	    rad_data->all_states[index].electron_densities[k] = 0;
	  } else {
	    sscanf(buf, "%lf%lf%lf%lf%lf%lf%d",&rad_data->all_states[index].electron_densities[k],
		   &rad_data->all_states[index].rad_fits[k].A,
		   &rad_data->all_states[index].rad_fits[k].alpha,
		   &rad_data->all_states[index].rad_fits[k].beta,
		   &rad_data->all_states[index].rad_fits[k].V0,
		   &rad_data->all_states[index].rad_fits[k].gamma,
		   &rad_data->all_states[index].rad_fits[k].te_intervals);
	    int te_intervals = rad_data->all_states[index].rad_fits[k].te_intervals;
	    rad_data->all_states[index].rad_fits[k].te = (double*)malloc(te_intervals*sizeof(double));
	    rad_data->all_states[index].rad_fits[k].Lz = (double*)malloc(te_intervals*sizeof(double));
	    int offset=0;
	    if (fgets(buf, 1000*15, fptr)!=NULL) {
	      char *data = buf;	
	      for (int l=0; l<te_intervals; l++) {
		sscanf(data, " %lf%n", &rad_data->all_states[index].rad_fits[k].te[l], &offset);
		data += offset;
	      }
	    }
	    if (fgets(buf, 1000*15, fptr) !=NULL) {
	      char *data2 = buf;
	      for (int l=0; l<te_intervals; l++) {
		sscanf(data2, " %lf%n", &rad_data->all_states[index].rad_fits[k].Lz[l], &offset);
		data2 += offset;
	      }
	    }
	  }
	}
      }
    }
  }
  fclose(fptr);
  return rad_data;
}

int gkyl_get_fit_params(const struct all_radiation_states rad_data, int atomic_z, int charge_state, double *a, double *alpha, double *beta, double *gamma, double *V0, int num_densities){
  int location = 0;
  double ref_dens = 19;
  atomic_z = atomic_z-1;
  int index = atomic_z*rad_data.max_atomic_number+charge_state;
  if (!rad_data.all_states[index].state_exists)
    return 1;
  if (num_densities>1)
    return 2;
  
  if (num_densities==1) {
    for (int i = 0; i<rad_data.all_states[index].number_of_densities; i++){
      if ( fabs(rad_data.all_states[index].electron_densities[i]-ref_dens)<
	   fabs(rad_data.all_states[index].electron_densities[location]-ref_dens)) {
	location = i;
      }
    }
    a[0] = rad_data.all_states[index].rad_fits[location].A;
    alpha[0] = rad_data.all_states[index].rad_fits[location].alpha;
    beta[0] = rad_data.all_states[index].rad_fits[location].beta;
    gamma[0] = rad_data.all_states[index].rad_fits[location].gamma;
    V0[0] = rad_data.all_states[index].rad_fits[location].V0;
  } else { // More than 1 density untested
    int count = 0;
    for (int i=0; i<rad_data.all_states[index].number_of_densities;
	 i=i+rad_data.all_states[index].number_of_densities/num_densities) {
      a[count] = rad_data.all_states[index].rad_fits[i].A;
      alpha[count] = rad_data.all_states[index].rad_fits[i].alpha;
      beta[count] = rad_data.all_states[index].rad_fits[i].beta;
      gamma[count] = rad_data.all_states[index].rad_fits[i].gamma;
      V0[count] = rad_data.all_states[index].rad_fits[i].V0;
    }
  }
  return 0;
}

int gkyl_get_fit_lz(const struct all_radiation_states rad_data, int atomic_z, int charge_state, double ne, double* te, double* Lz){
  int location = 0;
  atomic_z = atomic_z-1;
  int index = atomic_z*rad_data.max_atomic_number+charge_state;
  if (!rad_data.all_states[index].state_exists)
    return 1;
  
  for (int i = 0; i<rad_data.all_states[index].number_of_densities; i++){
    if (fabs(rad_data.all_states[index].electron_densities[i]-ne)<
	 fabs(rad_data.all_states[index].electron_densities[location]-ne)) 
      location = i;
  }

  int location2 = 0;
  for (int i=0; i<rad_data.all_states[index].rad_fits[location].te_intervals; i++){
    if (fabs(rad_data.all_states[index].rad_fits[location].te[i]-te[0])<
	fabs(rad_data.all_states[index].rad_fits[location].te[location2]-te[0])) 
      location2 = i;
  }
  te[0]=rad_data.all_states[index].rad_fits[location].te[location2];
  Lz[0]=rad_data.all_states[index].rad_fits[location].Lz[location2];
  return 0;
}

void gkyl_release_fit_params(struct all_radiation_states *rad_data){
  int max_Z = rad_data->max_atomic_number;
  for (int i=0; i<max_Z; i++){
    for (int j=0; j<max_Z; j++){
      if (rad_data->all_states[i*max_Z+j].state_exists) {
	free(rad_data->all_states[i*max_Z+j].electron_densities);
	for (int k=0; k<rad_data->all_states[i*max_Z+j].number_of_densities; k++){
	  free(rad_data->all_states[i*max_Z+j].rad_fits[k].te);
	  free(rad_data->all_states[i*max_Z+j].rad_fits[k].Lz);
	}
      }
      free(rad_data->all_states[i*max_Z+j].rad_fits);
    }
  }
  free(rad_data->all_states);
  free(rad_data);
}

