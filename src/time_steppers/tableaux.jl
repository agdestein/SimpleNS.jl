"""
Set up Butcher arrays `A`, `b`, and `c`, as well as and SSP coefficient `r`.
For families of methods, optional input `s` is the number of stages.

Original (MATLAB) by David Ketcheson, extended by Benjamin Sanderse.
"""

# Default SSP coefficient
# r = 0

## ================Explicit Methods=========================

"""
    FE11(; kwargs...)

FE11.
"""
function FE11(; kwargs...)
    # Forward Euler
    A = fill(0, 1, 1)
    b = [1]
    c = [0]
    r = 1
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    SSP22(; kwargs...)

SSP22.
"""
function SSP22(; kwargs...)
    A = [0 0; 1 0]
    b = [1 // 2, 1 // 2]
    c = sum(eachcol(A))
    r = 1
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    SSP42(; kwargs...)

SSP42.
"""
function SSP42(; kwargs...)
    s = 4
    A = [0 0 0 0; 1//3 0 0 0; 1//3 1//3 0 0; 1//3 1//3 1//3 0]
    b = fill(1 // 4, s)
    c = sum(eachcol(A))
    r = 3
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    SSP33(; kwargs...)

SSP33.
"""
function SSP33(; kwargs...)
    A = [0 0 0; 1 0 0; 1//4 1//4 0]
    b = [1 // 6, 1 // 6, 2 // 3]
    c = sum(eachcol(A))
    r = 1
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    SSP43(; kwargs...)

SSP43.
"""
function SSP43(; kwargs...)
    A = [0 0 0 0; 1//2 0 0 0; 1//2 1//2 0 0; 1//6 1//6 1//6 0]
    b = [1 // 6, 1 // 6, 1 // 6, 1 // 2]
    c = sum(eachcol(A))
    r = 2
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    SSP104(; kwargs...)

SSP104.
"""
function SSP104(; kwargs...)
    s = 10
    α0 = diagm(-1 => fill(1 // 1, s - 1))
    α0[6, 5] = 2 // 5
    α0[6, 1] = 3 // 5
    β0 = 1 // 6 * diagm(-1 => fill(1 // 1, s - 1))
    β0[6, 5] = 1 // 15
    A = (I(s) - α0) \ β0
    b = fill(1 // 10, s)
    c = sum(eachcol(A))
    r = 6
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    rSSPs2(s = 2; kwargs...)

Rational (optimal, low-storage) `s`-stage 2nd order SSP.
"""
function rSSPs2(s = 2; kwargs...)
    s ≥ 2 || error("Explicit second order SSP family requires s ≥ 2")
    r = s - 1
    α = [fill(0, 1, s); 1 // 1 * I(s)]
    α[s+1, s] = (s - 1) // s
    β = α .// r
    α[s+1, 1] = 1 // s
    A = (I(s) - α[1:s, :]) \ β[1:s, :]
    b = β[s+1, :] + A'α[s+1, :]
    c = sum(eachcol(A))
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    rSSPs3(s = 4; kwargs...)

Rational (optimal, low-storage) `s^2`-stage 3rd order SSP.
"""
function rSSPs3(s = 4; kwargs...)
    if !(round(sqrt(s)) ≈ sqrt(s)) || s < 4
        error("Explicit third order SSP family requires s = n^2, n > 1")
    end
    n = s^2
    r = n - s
    α = [fill(0, 1, n); 1 // 1 * I(n)]
    α[s*(s+1)÷2+1, s*(s+1)÷2] = (s - 1) // (2s - 1)
    β = α .// r
    α[s*(s+1)÷2+1, (s-1)*(s-2)÷2+1] = s // (2s - 1)
    A = (I(n) - α[1:n, :]) \ β[1:n, :]
    b = β[n+1, :] + A'α[n+1, :]
    c = sum(eachcol(A))
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    Wray3(; kwargs...)

Wray's RK3.
"""
function Wray3(; kwargs...)
    A = zeros(Rational, 3, 3)
    A[2, 1] = 8 // 15
    A[3, 1] = (8 // 15) - (17 // 60)
    A[3, 2] = 5 // 12
    b = [(8 // 15) - (17 // 60), 0, 3 // 4]
    c = [0, A[2, 1], A[3, 1] + A[3, 2]]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RK56(; kwargs...)

RK56.
"""
function RK56(; kwargs...)
    r = 0
    A = [
        0 0 0 0 0 0
        1//4 0 0 0 0 0
        1//8 1//8 0 0 0 0
        0 0 1//2 0 0 0
        3//16 -3//8 3//8 9//16 0 0
        -3//7 8//7 6//7 -12//7 8//7 0
    ]
    b = [7 // 90, 0, 16 // 45, 2 // 15, 16 // 45, 7 // 90]
    c = [0, 1 // 4, 1 // 4, 1 // 2, 3 // 4, 1]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    DOPRI6(; kwargs...)

Dormand-Price pair.
"""
function DOPRI6(; kwargs...)
    A = [
        0 0 0 0 0 0
        1//5 0 0 0 0 0
        3//40 9//40 0 0 0 0
        44//45 -56//15 32//9 0 0 0
        19372//6561 -25360//2187 64448//6561 -212//729 0 0
        9017//3168 -355//33 46732//5247 49//176 -5103//18656 0
    ]
    b = [35 // 384, 0, 500 // 1113, 125 // 192, -2187 // 6784, 11 // 84]
    c = sum(eachcol(A))
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end


## ===================Half explicit methods========================

"""
    HEM3(; kwargs...)

Brasey and Hairer.
"""
function HEM3(; kwargs...)
    r = 0
    A = [0 0 0; 1//3 0 0; -1 2 0]
    b = [0, 3 // 4, 1 // 4]
    c = sum(eachcol(A))
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    HEM3BS(; kwargs...)

HEM3BS.
"""
function HEM3BS(; kwargs...)
    r = 0
    A = [0 0 0; 1//2 0 0; -1 2 0]
    b = [1 // 6, 2 // 3, 1 // 6]
    c = sum(eachcol(A))
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    HEM5(; kwargs...)

Brasey and Hairer, 5 stage, 4th order.
"""
function HEM5(; kwargs...)
    r = 0
    A = [
        0 0 0 0 0
        3//10 0 0 0 0
        (1+sqrt(6))/30 (11-4*sqrt(6))/30 0 0 0
        (-79-31*sqrt(6))/150 (-1-4*sqrt(6))/30 (24+11*sqrt(6))/25 0 0
        (14+5*sqrt(6))/6 (-8+7*sqrt(6))/6 (-9-7*sqrt(6))/4 (9-sqrt(6))/4 0
    ]
    b = [0, 0, (16 - sqrt(6)) / 36, (16 + sqrt(6)) / 36, 1 / 9]
    c = sum(eachcol(A))
    runge_kutta_method(A, b, c, r; kwargs...)
end

## ================Classical Methods=========================

# Gauss-Legendre methods -- order 2s

"""
    GL1(; kwargs...)

GL1.
"""
function GL1(; kwargs...)
    r = 2
    A = fill(1 // 2, 1, 1)
    b = [1]
    c = [1 // 2]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    GL2(; kwargs...)

GL2.
"""
function GL2(; kwargs...)
    r = 0
    A = [
        1/4 1/4-sqrt(3)/6
        1/4+sqrt(3)/6 1/4
    ]
    b = [1 / 2, 1 / 2]
    c = [1 / 2 - sqrt(3) / 6, 1 / 2 + sqrt(3) / 6]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    GL3(; kwargs...)

GL3.
"""
function GL3(; kwargs...)
    r = 0
    A = [
        5/36 (80-24*sqrt(15))/360 (50-12*sqrt(15))/360
        (50+15*sqrt(15))/360 2/9 (50-15*sqrt(15))/360
        (50+12*sqrt(15))/360 (80+24*sqrt(15))/360 5/36
    ]
    b = [5 / 18, 4 / 9, 5 / 18]
    c = [(5 - sqrt(15)) / 10, 1 / 2, (5 + sqrt(15)) / 10]
    runge_kutta_method(A, b, c, r; kwargs...)
end

# Radau IA methods -- order 2s-1

"""
    RIA1(; kwargs...)

This is implicit Euler.
"""
function RIA1(; kwargs...)
    r = 1
    A = fill(1, 1, 1)
    b = [1]
    c = [0]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RIA2(; kwargs...)

RIA2.
"""
function RIA2(; kwargs...)
    r = 0
    A = [
        1//4 -1//4
        1//4 5//12
    ]
    b = [1 // 4, 3 // 4]
    c = [0, 2 // 3]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RIA3(; kwargs...)

RIA3.
"""
function RIA3(; kwargs...)
    r = 0
    A = [
        1/9 (-1-sqrt(6))/18 (-1+sqrt(6))/18
        1/9 (88+7*sqrt(6))/360 (88-43*sqrt(6))/360
        1/9 (88+43*sqrt(6))/360 (88-7*sqrt(6))/360
    ]
    b = [1 // 9, (16 + sqrt(6)) / 36, (16 - sqrt(6)) / 36]
    c = [0, (6 - sqrt(6)) / 10, (6 + sqrt(6)) / 10]
    runge_kutta_method(A, b, c, r; kwargs...)
end

# Radau IIA methods -- order 2s-1

"""
    RIIA1(; kwargs...)

RIIA1.
"""
function RIIA1(; kwargs...)
    r = 1
    A = fill(1, 1, 1)
    b = [1]
    c = [1]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RIIA2(; kwargs...)

RIIA2.
"""
function RIIA2(; kwargs...)
    r = 0
    A = [
        5//12 -1//12
        3//4 1//4
    ]
    b = [3 // 4, 1 // 4]
    c = [1 // 3, 1]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RIIA3(; kwargs...)

RIIA3.
"""
function RIIA3(; kwargs...)
    r = 0
    A = [
        (88-7*sqrt(6))/360 (296-169*sqrt(6))/1800 (-2+3*sqrt(6))/225
        (296+169*sqrt(6))/1800 (88+7*sqrt(6))/360 (-2-3*sqrt(6))/225
        (16-sqrt(6))/36 (16+sqrt(6))/36 1/9
    ]
    b = [(16 - sqrt(6)) / 36, (16 + sqrt(6)) / 36, 1 / 9]
    c = [(4 - sqrt(6)) / 10, (4 + sqrt(6)) / 10, 1]

    runge_kutta_method(A, b, c, r; kwargs...)
end

# Lobatto IIIA methods -- order 2s-2

"""
    LIIIA2(; kwargs...)

LIIIA2.
"""
function LIIIA2(; kwargs...)
    r = 0
    A = [0 0; 1//2 1//2]
    b = [1 // 2, 1 // 2]
    c = [0, 1]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    LIIIA3(; kwargs...)

LIIIA3.
"""
function LIIIA3(; kwargs...)
    r = 0
    A = [
        0 0 0
        5//24 1//3 -1//24
        1//6 2//3 1//6
    ]
    b = [1 // 6, 2 // 3, 1 // 6]
    c = [0, 1 // 2, 1]
    runge_kutta_method(A, b, c, r; kwargs...)
end

# Chebyshev methods

"""
    CHDIRK3(; kwargs...)

Chebyshev based DIRK (not algebraically stable).
"""
function CHDIRK3(; kwargs...)
    A = [
        0 0 0
        1//4 1//4 0
        0 1 0
    ]
    b = [1 // 6, 2 // 3, 1 // 6]
    c = [0, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    CHCONS3(; kwargs...)

CHCONS3.
"""
function CHCONS3(; kwargs...)
    A = [
        1//12 -1//6 1//12
        5//24 1//3 -1//24
        1//12 5//6 1//12
    ]
    b = [1 // 6, 2 // 3, 1 // 6]
    c = [0, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    CHC3(; kwargs...)

Chebyshev quadrature and C(3) satisfied. Note this equals Lobatto IIIA.
"""
function CHC3(; kwargs...)
    A = [
        0 0 0
        5//24 1//3 -1//24
        1//6 2//3 1//6
    ]
    b = [1 // 6, 2 // 3, 1 // 6]
    c = [0, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    CHC5(; kwargs...)

CHC5.
"""
function CHC5(; kwargs...)
    A = [
        0 0 0 0 0
        0.059701779686442 0.095031716019062 -0.012132034355964 0.006643368370744 -0.002798220313558
        0.016666666666667 0.310110028629970 0.200000000000000 -0.043443361963304 0.016666666666667
        0.036131553646891 0.260023298295923 0.412132034355964 0.171634950647605 -0.026368446353109
        0.033333333333333 0.266666666666667 0.400000000000000 0.266666666666667 0.033333333333333
    ]
    b = [1 // 30, 4 // 15, 2 // 5, 4 // 15, 1 // 30]
    c = [0, 0.146446609406726, 0.500000000000000, 0.853553390593274, 1.000000000000000]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

## ==================Miscellaneous Methods================

"""
    Mid22(; kwargs...)

Midpoint 22 method.
"""
function Mid22(; kwargs...)
    A = [0 0; 1//2 0]
    b = [0, 1]
    c = [0, 1 // 2]
    r = 1 // 2
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    MTE22(; kwargs...)

Minimal truncation error 22 method (Heun).
"""
function MTE22(; kwargs...)
    A = [0 0; 2//3 0]
    b = [1 // 4, 3 // 4]
    c = [0, 2 // 3]
    r = 1 // 2
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    CN22(; kwargs...)

Crank-Nicholson.
"""
function CN22(; kwargs...)
    A = [0 0; 1//2 1//2]
    b = [1 // 2, 1 // 2]
    c = [0, 1]
    r = 2
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    Heun33(; kwargs...)

Heun33.
"""
function Heun33(; kwargs...)
    A = [0 0 0; 1//3 0 0; 0 2//3 0]
    b = [1 // 4, 0, 3 // 4]
    c = sum(eachcol(A))
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RK33C2(; kwargs...)

RK3 satisfying C(2) for i=3.
"""
function RK33C2(; kwargs...)
    A = [0 0 0; 2//3 0 0; 1//3 1//3 0]
    b = [1 // 4, 0, 3 // 4]
    c = [0, 2 // 3, 2 // 3]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RK33P2(; kwargs...)

RK3 satisfying the second order condition for the pressure.
"""
function RK33P2(; kwargs...)
    A = [0 0 0; 1//3 0 0; -1 2 0]
    b = [0, 3 // 4, 1 // 4]
    c = [0, 1 // 3, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RK44(; kwargs...)

Classical fourth order.
"""
function RK44()
    A = [0 0 0 0; 1//2 0 0 0; 0 1//2 0 0; 0 0 1 0]
    b = [1 // 6, 1 // 3, 1 // 3, 1 // 6]
    c = sum(eachcol(A))
    r = 0
    runge_kutta_method(A, b, c, r)
end

"""
    RK44C2(; kwargs...)

RK4 satisfying C(2) for i=3.
"""
function RK44C2(; kwargs...)
    A = [0 0 0 0; 1//4 0 0 0; 0 1//2 0 0; 1 -2 2 0]
    b = [1 // 6, 0, 2 // 3, 1 // 6]
    c = [0, 1 // 4, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RK44C23(; kwargs...)

RK4 satisfying C(2) for i=3 and c2=c3.
"""
function RK44C23(; kwargs...)
    A = [0 0 0 0; 1//2 0 0 0; 1//4 1//4 0 0; 0 -1 2 0]
    b = [1 // 6, 0, 2 // 3, 1 // 6]
    c = [0, 1 // 2, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    RK44P2(; kwargs...)

RK4 satisfying the second order condition for the pressure (but not third order).
"""
function RK44P2(; kwargs...)
    A = [0 0 0 0; 1 0 0 0; 3//8 1//8 0 0; -1//8 -3//8 3//2 0]
    b = [1 // 6, -1 // 18, 2 // 3, 2 // 9]
    c = [0, 1, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

## ===================DSRK Methods========================

"""
    DSso2(; kwargs...)

CBM's DSRKso2.
"""
function DSso2(; kwargs...)
    A = [3//4 -1//4; 1 0]
    b = [1, 0]
    c = [1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    DSRK2(; kwargs...)

CBM's DSRK2.
"""
function DSRK2(; kwargs...)
    A = [
        1//2 -1//2
        1//2 1//2
    ]
    b = [1 // 2, 1 // 2]
    c = [0, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    DSRK3(; kwargs...)

Zennaro's DSRK3.
"""
function DSRK3(; kwargs...)
    A = [
        5//2 -2 -1//2
        -1 2 -1//2
        1//6 2//3 1//6
    ]
    b = [1 // 6, 2 // 3, 1 // 6]
    c = [0, 1 // 2, 1]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

## ==================="Non-SSP" Methods of Wong & Spiteri========================

"""
    NSSP21(; kwargs...)

NSSP21.
"""
function NSSP21(; kwargs...)
    A = [
        0 0
        3//4 0
    ]
    b = [0, 1]
    c = [0, 3 // 4]
    r = 0
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    NSSP32(; kwargs...)

NSSP32.
"""
function NSSP32(; kwargs...)
    r = 0
    A = [
        0 0 0
        1//3 0 0
        0 1 0
    ]
    b = [1 // 2, 0, 1 // 2]
    c = [0, 1 // 3, 1]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    NSSP33(; kwargs...)

NSSP33.
"""
function NSSP33(; kwargs...)
    r = 0
    A = [
        0 0 0
        -4//9 0 0
        7//6 -1//2 0
    ]
    b = [1 // 4, 0, 3 // 4]
    c = [0, -4 // 9, 2 // 3]
    runge_kutta_method(A, b, c, r; kwargs...)
end

"""
    NSSP53(; kwargs...)

NSSP53.
"""
function NSSP53(; kwargs...)
    r = 0
    A = [
        0 0 0 0 0
        1//7 0 0 0 0
        0 3//16 0 0 0
        0 0 1//3 0 0
        0 0 0 2//3 0
    ]
    b = [1 // 4, 0, 0, 0, 3 // 4]
    c = [0, 1 // 7, 3 // 16, 1 // 3, 2 // 3]
    runge_kutta_method(A, b, c, r; kwargs...)
end
