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
    _spline!(y2, x, y, n, yp1, ypn)
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
* X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
* N: SIZE OF ARRAYS X,Y
* YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
         >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
* Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
"""
function _spline!(y2, x, y, n::Int, yp1::Float64, ypn::Float64)
    u = zeros(n)
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
    # @info "y2 $y2"
    # @info "u $u"

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
    _splint(xa, ya, y2a, n, x)
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
CALCULATE CUBIC SPLINE INTERP VALUE
ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
* XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
* Y2A: ARRAY OF SECOND DERIVATIVES
* N: SIZE OF ARRAYS XA,YA,Y2A
* X: ABSCISSA FOR INTERPOLATION
* Y: OUTPUT VALUE
"""
function _splint(xa, ya, y2a, n::Int, x)
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
    _splini(xa, ya, y2a, n, x)
      SUBROUTINE SPLINI(XA,YA,Y2A,N,X,YI)
INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
* XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
* Y2A: ARRAY OF SECOND DERIVATIVES
* N: SIZE OF ARRAYS XA,YA,Y2A
* X: ABSCISSA ENDPOINT FOR INTEGRATION
* Y: OUTPUT VALUE
"""
function _splini(xa, ya, y2a, n, x)
    yi = 0.0
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

"""
    _ccor2(alt, r, h1, zh, h2)
      FUNCTION  CCOR2(ALT, R,H1,ZH,H2)

O&O2 CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
"""
function _ccor2(alt, r, h1, zh, h2)
    e1 = (alt - zh) / h1
    e2 = (alt - zh) / h2
    if e1 > 70 || e2 > 70
        return 1.0
    elseif e1 < -70 && e2 < -70
        return exp(r)
    end
    ex1 = exp(e1)
    ex2 = exp(e2)
    return exp(r / (1.0 + 0.5 * (ex1 + ex2)))
end

"""
    _ccor(alt, r, h1, zh)
      FUNCTION  CCOR(ALT, R,H1,ZH)

CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
* ALT - altitude
* R - target ratio
* H1 - transition scale length
* ZH - altitude of 1/2 R
"""
function _ccor(alt, r, h1, zh)
    e = (alt - zh) / h1
    if e > 70
        return 1.0
    elseif e < -70
        return exp(r)
    end

    ex = exp(e)
    return exp(r / (1.0 + ex))
end

"""
    _dnet(dd, dm, zhm, xmm, xm)
      FUNCTION DNET(DD,DM,ZHM,XMM,XM)

TURBOPAUSE CORRECTION FOR MSIS MODELS
  Root mean density
8/20/80
*   DD - diffusive density
*   DM - full mixed density
*   ZHM - transition scale length
*   XMM - full mixed molecular weight
*   XM  - species molecular weight
*   DNET - combined density
"""
function _dnet(dd, dm, zhm, xmm, xm)
    a = zhm / (xmm - xm)
    if dm <= 0 || dd <= 0
        @warn "DNET LOG ERROR $dm $dd $xm"
        if dd == 0 && dm == 0
            dd = 1
            #this might not be 1-1 with the original. if inputs are passed by
            #reference then this should change dd external to the call
        end
        (dm == 0) && (return dd)
        (dd == 0) && (return dm)
    end
    ylog = a * log(dm / dd)
    (ylog < -10) && (return dd)
    (ylog > 10) && (return dm)
    return dd * (1.0 + exp(ylog))^(1.0 / a)
end

"""
    _scalh(alt, xm, temp, gsurf, re)
      FUNCTION SCALH(ALT,XM,TEMP)

Calculate scale height (km)
"""
function _scalh(alt, xm, temp, gsurf, re)
    g = gsurf / (1.0 + alt / re)^2
    return RGAS * temp / (g * xm)
end

"""
    (gv, rEff) = _glatf(lat::Float64)
      SUBROUTINE GLATF(LAT,GV,REFF)

CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
RADIUS (REFF)
"""
function _glatf!(parmb::PARMB, lat::Float64)
    c2 = cos(2.0 * DGTR * lat)
    parmb.GSURF = 980.616 * (1.0 - 0.0026373 * c2)
    parmb.RE = 2.0 * parmb.GSURF / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5
    return nothing
end
