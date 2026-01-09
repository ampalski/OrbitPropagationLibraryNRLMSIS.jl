function _get_base_mjd(jd::JDate)
    return floor(Int, sum(jdate_to_mjdate(jd).epoch))
end
function _get_base_mjd(jd::MJDate)
    return floor(Int, sum(jd.epoch))
end

# Borrowed from Statistics.jl
function mean(input::AbstractVector)
    y = iterate(input)
    if y === nothing
        error("Empty Vector")
    end
    count = 1
    value, state = y
    total = value
    y = iterate(input, state)
    while y !== nothing
        value, state = y
        total += value
        count += 1
        y = iterate(input, state)
    end
    return total / count
end

"""
    _spline!(y2, x, y, yp1, ypn)
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
N: SIZE OF ARRAYS X,Y
YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
         >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
"""
function _spline!(y2, x, y, yp1::Float64, ypn::Float64)
    u = zeros(100)
    n = length(x)
    if yp1 > 0.99e30
        y2[1] = 0.0
        u[1] = 0.0
    else
        y2[1] = -0.5
        u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1)
    end
    @inbounds for i in 2:(n - 1)
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (
            6.0 * (
                (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                    (y[i] - y[i - 1]) / (x[i] - x[i - 1])
            ) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]
        ) / p
    end

    if ypn > 0.99e30
        qn = 0
        un = 0
    else
        qn = 0.5
        un = (3.0 / (x[n] - x[n - 1])) * (
            ypn - (y[n] - y[n - 1]) /
                (x[n] - x[n - 1])
        )
    end
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0)
    @inbounds for i in (n - 1):-1:1
        y2[i] = y2[i] * y2[i + 1] + u[i]
    end

    return nothing
end

"""
    _splint(xa, ya, y2a, x)
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
CALCULATE CUBIC SPLINE INTERP VALUE
ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
Y2A: ARRAY OF SECOND DERIVATIVES
N: SIZE OF ARRAYS XA,YA,Y2A
X: ABSCISSA FOR INTERPOLATION
Y: OUTPUT VALUE
"""
function _splint(xa, ya, y2a, x)
    n = length(xa)
    k_lo = 1
    k_hi = n
    while k_hi - k_lo > 1
        k = div(k_hi + k_lo, 2)
        if xa[k] > x
            k_hi = k
        else
            k_lo = k
        end
    end
    h = xa[k_hi] - xa[k_lo]
    if h == 0.0
        @info ("Bad xa input to splint")
    end
    a = (xa[k_hi] - x) / h
    b = (x - xa[k_lo]) / h
    y = a * ya[k_lo] + b * ya[k_hi] + (
        (a * a * a - a) * y2a[k_lo] +
            (b * b * b - b) * y2a[k_hi]
    ) * h * h / 6
    return y
end

"""
    _splini(xa, ya, y2a, x)
      SUBROUTINE SPLINI(XA,YA,Y2A,N,X,YI)
INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
 XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 Y2A: ARRAY OF SECOND DERIVATIVES
 N: SIZE OF ARRAYS XA,YA,Y2A
 X: ABSCISSA ENDPOINT FOR INTEGRATION
 Y: OUTPUT VALUE
"""
function _splini(yi, xa, ya, y2a, x)
    yi = 0.0
    n = length(xa)
    k_lo = 1
    k_hi = 2
    while x > xa[k_lo] && k_hi <= n
        xx = x
        if k_hi < n
            xx = min(x, xa[k_hi])
        end
        h = xa[k_hi] - xa[k_lo]
        a = (xa[k_hi] - xx) / h
        b = (xx - xa[k_lo]) / h
        a2 = a * a
        b2 = b * b
        yi += h * (
            (1.0 - a2) * ya[k_lo] / 2.0 + b2 * ya[k_hi] / 2.0 +
                (
                (-(1.0 + a2 * a2) / 4.0 + a2 / 2.0) * y2a[k_lo] +
                    (b2 * b2 / 4.0 - b2 / 2.0) * y2a[k_hi]
            ) * h * h / 6.0
        )
        k_lo += 1
        k_hi += 1
    end

    return yi
end
#=
=#
