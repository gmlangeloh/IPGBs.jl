#Some numerical constants and useful functions
const EPSILON = 0.0001

approx_equal(x, y) = abs(x - y) < EPSILON
is_approx_zero(x) = approx_equal(x, 0.0)
