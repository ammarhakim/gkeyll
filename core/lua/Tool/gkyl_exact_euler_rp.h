// IFC for solving exact Euler RP

struct _ProblemState {
  double dl, ul, pl;
  double dr, ur, pr;
  double lower, upper;
  int ncell;
  double tEnd;
  double gas_gamma;
  double disLoc;
};

void solveRiemannProblem(struct _ProblemState _ps, const char *out_prefix);
