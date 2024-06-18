#Some numerical constants and parameters
using GLPK

const DEFAULT_SOLVER = GLPK.Optimizer
const EPSILON = 0.0001

approx_equal(x, y) = abs(x - y) < EPSILON
is_approx_zero(x) = approx_equal(x, 0.0)
