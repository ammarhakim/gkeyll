typedef struct atomic_data{
  char* reaction_type;
  char* symbol;
  char* name;
  int te_intervals;
  int ne_intervals;
  int Z;
  int* charge_states;
  double* logT;
  double* logNe;
  double* logData;
} atomic_data;
