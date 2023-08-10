module GBOptions

auto_reduce_freq :: Int = 0
debug :: Bool = false

function initialize_parameters(
    auto_reduce_freq = 5,
    debug = false
)
    global auto_reduce_freq = auto_reduce_freq
    global debug = debug
end

end