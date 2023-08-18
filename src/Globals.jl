using Logging

#Some numerical constants and parameters
const EPSILON = 0.0001

@enum AUTO_REDUCTION_TYPE begin
    NONE
    FIXED_ITERATIONS
    FIXED_ELEMENTS
    FRACTION_ELEMENTS
end

#4ti2 uses 2500 here.
AUTO_RED_FREQ :: Float64 = 0
AUTO_RED_TYPE :: AUTO_REDUCTION_TYPE = NONE
DEBUG :: Bool = false
INFO :: Bool = false
CACHE_TREE_SIZE :: Int = 100

function initialize_parameters(;
    auto_reduce_freq = 2500,
    auto_reduce_type = FIXED_ITERATIONS,
    cache_tree_size = 100,
    debug = false,
    info = false,
)
    global AUTO_RED_FREQ = auto_reduce_freq
    global AUTO_RED_TYPE = auto_reduce_type
    global CACHE_TREE_SIZE = cache_tree_size
    global DEBUG = debug
    global INFO = info
    if DEBUG || INFO
        loglevel = DEBUG ? Logging.Debug : Logging.Info
        logger = SimpleLogger(stderr, loglevel)
        global_logger(logger)
    else
        logger = NullLogger()
        global_logger(logger)
    end
end

approx_equal(x, y) = abs(x - y) < EPSILON
is_approx_zero(x) = approx_equal(x, 0.0)

prev_gb_size = 0
function should_auto_reduce(i, gb_size)
    if AUTO_RED_TYPE == NONE
        return false
    elseif AUTO_RED_TYPE == FIXED_ITERATIONS
        return i != 0 && i % AUTO_RED_FREQ == 0
    elseif AUTO_RED_TYPE == FIXED_ELEMENTS
        delta = gb_size - prev_gb_size
        should_reduce = gb_size != 0 && delta >= AUTO_RED_FREQ
        global prev_gb_size = gb_size
        return should_reduce
    elseif AUTO_RED_TYPE == FRACTION_ELEMENTS
        delta = gb_size - prev_gb_size
        should_reduce = gb_size != 0 && (delta / gb_size >= AUTO_RED_FREQ)
        global prev_gb_size = gb_size
        return should_reduce
    end
    @assert false #Should never reach this
    return false
end