using Logging

#Some numerical constants and parameters
const EPSILON = 0.0001

#4ti2 uses 2500 here.
AUTO_REDUCE_FREQ :: Int = 0
DEBUG :: Bool = false
INFO :: Bool = false

function initialize_parameters(;
    auto_reduce_freq = 2500,
    debug = false,
    info = false,
)
    global AUTO_REDUCE_FREQ = auto_reduce_freq
    global DEBUG = debug
    global INFO = info
    if DEBUG || INFO
        loglevel = DEBUG ? Logging.Debug : Logging.Info
        logger = SimpleLogger(stderr, loglevel)
        global_logger(logger)
    end
end

approx_equal(x, y) = abs(x - y) < EPSILON
is_approx_zero(x) = approx_equal(x, 0.0)
should_auto_reduce(i) = i != 0 && i % AUTO_REDUCE_FREQ == 0