# Higer order function

# Reduce for paired sum

my.sum <- function (x) {Reduce(`+`, x)}

my.sum(c(1:10))

# Filter for get certain values

small.even.numbers <- Filter(function (x) {x %% 2 == 0}, 1:10)
small.odd.numbers <- Filter(function (x) {x %% 2 == 1}, 1:10)

# Find for the first number of Filter

Find(is.prime, 1000:2000)

# Map for apply
is.divisor <- function(a, x) {x %% a == 0}
is.prime <- function (x) {length(Filter(function (a) {is.divisor(a, x)}, 1:x)) == 2}
proper.divisors <- function (x) {Filter(function (a) {is.divisor(a, x)}, 1:(x - 1))}
is.perfect <- function(x) {x == Reduce(`+`, proper.divisors(x))}
small.perfect.numbers <- Filter(is.perfect, 1:1000)
Map(proper.divisors, 1:5)

# Position is index version for Find

Position(is.prime, 1000:2000)

# Negate for flips Boolean

is.composite <- Negate(is.prime)
