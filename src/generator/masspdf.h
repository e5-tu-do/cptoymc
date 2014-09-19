namespace cptoymc {
namespace generator {
namespace masspdf {

double BK(double ni, double x);
double LnBK(double ni, double x);
double LogEval(double d, double l, double alpha, double beta, double delta);
double diff_eval(double d, double l, double alpha, double beta, double delta);
double EvalIpatia(double x, double l, double zeta, double fb, double sigma, double mu, double a, double n, double a2, double n2);
 
} // namespace masspdf
} // namespace generator
} // namespace cptoymc