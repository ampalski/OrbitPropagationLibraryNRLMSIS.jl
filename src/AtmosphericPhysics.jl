"""
    (tz, densm) = _densm(alt, d0, xm, tz, mn3, zn3, tn3, tgn3, mn2, zn2, tn2, tgn2, re, gsurf)
      FUNCTION DENSM(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)

Calculate Temperature and Density Profiles for lower atmosphere, where mixing
occurs.
"""
function _densm(
        alt, d0, xm, tz,
        mn3, zn3, tn3, tgn3,
        mn2, zn2, tn2, tgn2,
        re, gsurf
    )
    zeta(zz, zl) = (zz - zl) * (re + zl) / (re + zz)
    xs = MVector{10}(zeros(10))
    ys = MVector{10}(zeros(10))

    densm = d0
    if alt > zn2[1]
        (xm == 0.0) && (return (tz, tz))
        return (tz, densm)
    end
    # STRATOSPHERE/MESOSPHERE TEMPERATURE
    z = max(alt, zn2[mn2])
    mn = mn2
    z1 = zn2[1]
    z2 = zn2[mn]
    t1 = tn2[1]
    t2 = tn2[mn]
    zg = zeta(z, z1)
    zgdif = zeta(z2, z1)
    # Set up spline nodes
    for k in 1:mn
        xs[k] = zeta(zn2[k], z1) / zgdif
        ys[k] = 1.0 / tn2[k]
    end

    yd1 = -tgn2[1] / (t1 * t1) * zgdif
    yd2 = -tgn2[2] / (t2 * t2) * zgdif * ((re + z2) / (re + z1))^2

    # Calculate spline coefficients
    y2out = zeros(length(xs))
    _spline!(y2out, xs, ys, mn, yd1, yd2)
    x = zg / zgdif
    y = _splint(xs, ys, y2out, mn, x)

    # Temperature at altitude
    tz = 1.0 / y

    if xm != 0
        # Calculate stratospheric/mesosphere density
        glb = gsurf / (1.0 + z1 / re)^2
        gamm = xm * glb * zgdif / RGAS

        # Integrate temperature profile
        yi = _splini(xs, ys, y2out, mn, x)
        expl = max(gamm * yi, 50.0)

        # Density at altitude
        densm = densm * (t1 / tz) * exp(-expl)
    end

    if alt > zn3[1]
        (xm == 0.0) && (return (tz, tz))
        return (tz, densm)
    end

    # Troposphere/stratosphere temperature
    z = alt
    mn = mn3
    z1 = zn3[1]
    z2 = zn3[mn]
    t1 = tn3[1]
    t2 = tn3[mn]
    zg = zeta(z, z1)
    zgdif = zeta(z2, z1)

    # Set up spline nodes
    y2out2 = zeros(mn)
    for k in 1:mn
        xs[k] = zeta(zn3[k], z1) / zgdif
        ys[k] = 1.0 / tn3[k]
    end

    yd1 = -tgn3[1] / (t1 * t1) * zgdif
    yd2 = -tgn3[2] / (t2 * t2) * zgdif * ((re + z2) / (re + z1))^2
    _spline!(y2out2, xs, ys, mn, yd1, yd2)
    x = zg / zgdif
    y = _splint(xs, ys, y2out2, mn, x)

    # Temperature at altitude
    tz = 1.0 / y
    if xm == 0.0
        return (tz, tz)
    end

    # Calculate tropospheric/stratosphere density
    glb = gsurf / (1.0 + z1 / re)^2
    gamm = xm * glb * zgdif / RGAS

    # Integrate temperature profile
    yi = _splini(xs, ys, y2out2, mn, x)
    expl = max(gamm * yi, 50.0)

    # Density at altitude
    densm = densm * (t1 / tz) * exp(-expl)
    (xm == 0) && (densm = tz)
    return (tz, densm)
end

"""
    (tz, densu) = _densu(alt, dlb, tinf, tlb, xm, alpha, tz, zlb, s2, mn1, zn1, tn1, tgn1, re, gsurf)
      FUNCTION DENSU(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,
      MN1,ZN1,TN1,TGN1)

Calculate Temperature and Density Profiles for MSIS models

New lower thermo polynomial 10/30/89
"""
function _densu!(
        alt, dlb, tinf, tlb,
        xm, alpha, tz, zlb, s2,
        mn1::Int, zn1, tn1, tgn1,
        re, gsurf
    )
    #NOTE: tn1 and tgn1 are the things that *should* mutate. tz also does, but
    #is a float so I don't think will end up mutating and the current method of
    #using it as a local variable that mutates should be fine
    zeta(zz, zl) = (zz - zl) * (re + zl) / (re + zz)
    xs = MVector{5}(zeros(5))
    ys = MVector{5}(zeros(5))
    y2out = zeros(mn1)
    z1 = 1.0
    zgdif = 1.0
    x = 1.0
    t1 = 1.0
    mn = mn1


    # Joining altitude of Bates and spline
    za = zn1[1]
    z = max(alt, za)

    # Geopotential altitude difference from ZLB
    zg2 = zeta(z, zlb)

    # Bates temperature
    tt = tinf - (tinf - tlb) * exp(-s2 * zg2)
    ta = tt
    tz = tt
    densu = tz

    if alt < za
        # Calculate temperature below ZA
        # Temperature gradient at ZA from Bates profile

        dta = (tinf - ta) * s2 * ((re + zlb) / (re + za))^2
        tgn1[1] = dta
        tn1[1] = ta
        z = max(alt, zn1[mn1])
        mn = mn1
        z1 = zn1[1]
        z2 = zn1[mn]
        t1 = tn1[1]
        t2 = tn1[mn]

        # Geopotential difference from Z1
        zg = zeta(z, z1)
        zgdif = zeta(z2, z1)

        # Set up spline nodes
        for k in 1:mn
            xs[k] = zeta(zn1[k], z1) / zgdif
            ys[k] = 1.0 / tn1[k]
        end

        # End node derivatives
        yd1 = -tgn1[1] / (t1 * t1) * zgdif
        yd2 = -tgn1[2] / (t2 * t2) * zgdif * ((re + z2) / (re + z1))^2

        # Calculate spline coefficients
        _spline!(y2out, xs, ys, mn, yd1, yd2)
        x = zg / zgdif
        y = _splint(xs, ys, y2out, mn, x)

        # Temperature at altitude
        tz = 1.0 / y
        densu = tz
    end

    if xm == 0
        return (tz, densu)
    end

    # Calculate density above za
    glb = gsurf / (1.0 + zlb / re)^2
    gamma = xm * glb / (s2 * RGAS * tinf)
    expl = exp(-s2 * gamma * zg2)
    if expl > 50 || tt <= 0
        expl = 50.0
    end

    # Density at altitude
    densa = dlb * (tlb / tt)^(1.0 + alpha + gamma) * expl
    densu = densa

    if alt >= za
        return (tz, densu)
    end

    # Calculate density below za
    glb = gsurf / (1.0 + z1 / re)^2
    gamm = xm * glb * zgdif / RGAS

    # Integrate spline temperatures
    yi = _splini(xs, ys, y2out, mn, x)
    expl = gamm * yi
    if expl > 50 || tz <= 0
        expl = 50.0
    end

    # Density at altitude
    densu = densu * (t1 / tz)^(1.0 + alpha) * exp(-expl)

    return (tz, densu)
end

"""
    _glob7s(p, input, lpoly)
      FUNCTION GLOB7S(P)

VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
"""
function _glob7s(p, input::OplMsis_Input, lpoly::LPOLY)
    PSET = 2.0

    if p[100] == 0.0
        # this might throw an error, trying to change a static value
        p[100] = PSET
    end
    if p[100] != PSET
        error("Wrong parameter set for GLOB7S")
    end

    t = zeros(14)

    CD14 = cos(DR * (lpoly.DAY - p[14]))
    CD18 = cos(2.0 * DR * (lpoly.DAY - p[18]))
    CD32 = cos(DR * (lpoly.DAY - p[32]))
    CD39 = cos(2.0 * DR * (lpoly.DAY - p[39]))

    # F10.7
    t[1] = p[22] * lpoly.DFA

    # Time independent
    t[2] = p[2] * lpoly.PLG[3, 1] + p[3] * lpoly.PLG[5, 1] +
        p[23] * lpoly.PLG[7, 1] + p[27] * lpoly.PLG[2, 1] +
        p[15] * lpoly.PLG[4, 1] + p[60] * lpoly.PLG[6, 1]

    # Symmetrical annual
    t[3] = (
        p[19] + p[48] * lpoly.PLG[3, 1] +
            p[30] * lpoly.PLG[5, 1]
    ) * CD32

    # Symmetrical semiannual
    t[4] = (
        p[16] + p[17] * lpoly.PLG[3, 1] +
            p[31] * lpoly.PLG[5, 1]
    ) * CD18

    # Asymmetrical annual
    t[5] = (
        p[10] * lpoly.PLG[2, 1] + p[11] * lpoly.PLG[4, 1] +
            p[21] * lpoly.PLG[6, 1]
    ) * CD14

    # Asymmetrical semiannual
    t[6] = (p[38] * lpoly.PLG[2, 1]) * CD39

    # Diurnal
    if input.switches.SW[7] != 0
        T71 = p[12] * lpoly.PLG[3, 2] * CD14 * input.switches.SWC[5]
        T72 = p[13] * lpoly.PLG[3, 2] * CD14 * input.switches.SWC[5]
        t[7] = (
            (p[4] * lpoly.PLG[2, 2] + p[5] * lpoly.PLG[4, 2] + T71) *
                lpoly.CTLOC +
                (p[7] * lpoly.PLG[2, 2] + p[8] * lpoly.PLG[4, 2] + T72) *
                lpoly.STLOC
        )
    end

    # Semidiurnal
    if input.switches.SW[8] != 0
        T81 = (p[24] * lpoly.PLG[4, 3] + p[36] * lpoly.PLG[6, 3]) *
            CD14 * input.switches.SWC[5]
        T82 = (p[34] * lpoly.PLG[4, 3] + p[37] * lpoly.PLG[6, 3]) *
            CD14 * input.switches.SWC[5]
        t[8] = (p[6] * lpoly.PLG[3, 3] + p[42] * lpoly.PLG[5, 3] + T81) *
            lpoly.C2TLOC +
            (p[9] * lpoly.PLG[3, 3] + p[43] * lpoly.PLG[5, 3] + T82) *
            lpoly.S2TLOC

    end

    if input.switches.SW[14] != 0
        t[14] = p[40] * lpoly.PLG[4, 4] * lpoly.S3TLOC +
            p[41] * lpoly.PLG[4, 4] * lpoly.C3TLOC
    end

    # Magnetic activity
    if input.switches.SW == 1
        t[9] = lpoly.APDF * (
            p[33] +
                p[46] * lpoly.PLG[3, 1] * input.switches.SWC[2]
        )
    elseif input.switches.SW == -1
        t[9] = p[51] * lpoly.APT[1] +
            p[97] * lpoly.PLG[3, 1] * lpoly.APT[1] * input.switches.SWC[2]

    end

    if !(input.switches.SW[10] == 0 || input.switches.SW[11] == 0)
        t[11] = (
            1.0 + lpoly.PLG[2, 1] *
                (
                p[81] * input.switches.SWC[5] * cos(DR * (input.iyd - p[82])) +
                    p[86] * input.switches.SWC[6] * cos(2.0 * DR * (input.iyd - p[87]))
            ) +
                p[84] * input.switches.SWC[3] * cos(DR * (input.iyd - p[85])) +
                p[88] * input.switches.SWC[4] * cos(2 * DR * (input.iyd - p[89]))
        ) * (
            (
                p[65] * lpoly.PLG[3, 2] + p[66] * lpoly.PLG[5, 2] + p[67] * lpoly.PLG[7, 2] +
                    p[75] * lpoly.PLG[2, 2] + p[76] * lpoly.PLG[4, 2] + p[77] * lpoly.PLG[6, 2]
            ) *
                cos(DGTR * input.long) + (
                p[91] * lpoly.PLG[3, 2] + p[92] * lpoly.PLG[5, 2] + p[93] * lpoly.PLG[7, 2] +
                    p[78] * lpoly.PLG[2, 2] + p[79] * lpoly.PLG[4, 2] + p[80] * lpoly.PLG[6, 2]
            ) *
                sin(DGTR * input.long)
        )
    end

    tt = 0.0
    for i in 1:14
        tt += abs(input.switches.SW[i]) * t[i]
    end
    return tt
end

function _setLpoly!(lpoly, input)
    HR = 0.2618
    lpoly.DAY = input.iyd
    lpoly.XLONG = input.long

    # Calculate legendre polynomials
    c = sin(input.lat * DGTR)
    s = cos(input.lat * DGTR)
    c2 = c * c
    c4 = c2 * c2
    s2 = s * s
    lpoly.PLG[2, 1] = c
    lpoly.PLG[3, 1] = 0.5 * (3.0 * c2 - 1.0)
    lpoly.PLG[4, 1] = 0.5 * (5.0 * c * c2 - 3.0 * c)
    lpoly.PLG[5, 1] = (35.0 * c4 - 30 * c2 + 3) / 8.0
    lpoly.PLG[6, 1] = (63.0 * c2 * c2 * c - 70 * c2 * c + 15 * c) / 8.0
    lpoly.PLG[7, 1] = (11.0 * c * lpoly.PLG[6, 1] - 5.0 * lpoly.PLG[5, 1]) / 6.0
    #lpoly.PLG[8,1] = (13. * c * lpoly.PLG[7,1] - 6.0 * lpoly.PLG[6,1])/7.0
    lpoly.PLG[2, 2] = s
    lpoly.PLG[3, 2] = 3.0 * c * s
    lpoly.PLG[4, 2] = 1.5 * (5.0 * c2 - 1.0) * s
    lpoly.PLG[5, 2] = 2.5 * (7.0 * c2 * c - 3.0 * c) * s
    lpoly.PLG[6, 2] = 1.875 * (21.0 * c4 - 14.0 * c2 + 1.0) * s
    lpoly.PLG[7, 2] = (11.0 * c * lpoly.PLG[6, 2] - 6.0 * lpoly.PLG[5, 2]) / 5.0
    # lpoly.PLG[8,2) = (13.*c*plg(7,2)-7.*plg(6,2))/6.
    # lpoly.PLG[9,2) = (15.*c*plg(8,2)-8.*plg(7,2))/7.
    lpoly.PLG[3, 3] = 3.0 * s2
    lpoly.PLG[4, 3] = 15.0 * s2 * c
    lpoly.PLG[5, 3] = 7.5 * (7.0 * c2 - 1.0) * s2
    lpoly.PLG[6, 3] = 3.0 * c * lpoly.PLG[5, 3] - 2.0 * lpoly.PLG[4, 3]
    lpoly.PLG[7, 3] = (11.0 * c * lpoly.PLG[6, 3] - 7.0 * lpoly.PLG[5, 3]) / 4.0
    lpoly.PLG[8, 3] = (13.0 * c * lpoly.PLG[7, 3] - 8.0 * lpoly.PLG[6, 3]) / 5.0
    lpoly.PLG[4, 4] = 15.0 * s2 * s
    lpoly.PLG[5, 4] = 105.0 * s2 * s * c
    lpoly.PLG[6, 4] = (9.0 * c * lpoly.PLG[5, 4] - 7.0 * lpoly.PLG[4, 4]) / 2.0
    lpoly.PLG[7, 4] = (11.0 * c * lpoly.PLG[6, 4] - 8.0 * lpoly.PLG[5, 4]) / 3.0

    if !(input.switches.SW[7] == 0 && input.switches.SW[8] == 0 && input.switches.SW[14] == 0)
        lpoly.STLOC = sin(HR * input.stl)
        lpoly.CTLOC = cos(HR * input.stl)
        lpoly.S2TLOC = sin(2 * HR * input.stl)
        lpoly.C2TLOC = cos(2 * HR * input.stl)
        lpoly.S3TLOC = sin(3 * HR * input.stl)
        lpoly.C3TLOC = cos(3 * HR * input.stl)
    end
    return nothing
end

"""
    _globe7!(input, p, lpoly)
      FUNCTION GLOBE7(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)

CALCULATE G(L) FUNCTION 
Upper Thermosphere Parameters
"""
function _globe7!(input::OplMsis_Input, p, lpoly::LPOLY)
    HR = 0.2618
    SR = 7.2722e-5
    NSW = 14
    t = zeros(15)

    # 3 hr Magnetic activity functions
    # EQ A24d
    G0(A) = (
        A - 4.0 + (p[26] - 1.0) *
            (A - 4.0 + (exp(-abs(p[25]) * (A - 4.0)) - 1.0) / abs(p[25]))
    )

    # EQ A24c
    SUMEX(EX) = 1.0 + (1.0 - EX^19) / (1.0 - EX) * EX^(0.5)

    # EQ A24a
    SG0(EX) = (
        G0(input.ap[2]) + (
            G0(input.ap[3]) * EX + G0(input.ap[4]) * EX * EX +
                G0(input.ap[5]) * EX^3 +
                (G0(input.ap[6]) * EX^4 + G0(input.ap[7]) * EX^12)
                * (1 - EX^8) / (1 - EX)
        )
    ) / SUMEX(EX)

    if input.switches.SW[9] > 0
        SW9 = 1
    end
    if input.switches.SW[9] < 0
        SW9 = -1
    end


    CD14 = cos(DR * (lpoly.DAY - p[14]))
    CD18 = cos(2.0 * DR * (lpoly.DAY - p[18]))
    CD32 = cos(DR * (lpoly.DAY - p[32]))
    CD39 = cos(2.0 * DR * (lpoly.DAY - p[39]))

    # F10.7 Effect
    df = input.f107 - input.f107a
    lpoly.DFA = input.f107a - 150.0
    t[1] = p[20] * df * (1.0 + p[60] * lpoly.DFA) + p[21] * df * df + p[22] * lpoly.DFA +
        p[30] * lpoly.DFA^2
    f1 = 1.0 + (p[48] * lpoly.DFA + p[20] * df + p[21] * df * df) * input.switches.SWC[1]
    f2 = 1.0 + (p[50] * lpoly.DFA + p[20] * df + p[21] * df * df) * input.switches.SWC[1]

    # Time Independent
    t[2] = @inbounds (
        p[2] * lpoly.PLG[3, 1] +
            p[3] * lpoly.PLG[5, 1] +
            p[23] * lpoly.PLG[7, 1]
    ) +
        (p[15] * lpoly.PLG[3, 1]) * lpoly.DFA * input.switches.SWC[1] +
        p[27] * lpoly.PLG[2, 1]

    # Symmetrical annual
    t[3] = p[19] * CD32

    # Symmetrical semiannual
    t[4] = (p[16] + p[17] * lpoly.PLG[3, 1]) * CD18

    # Asymmetrical annual
    t[5] = f1 * (p[10] * lpoly.PLG[2, 1] + p[11] * lpoly.PLG[4, 1]) * CD14

    # Asymmetrical semiannual
    t[6] = p[38] * lpoly.PLG[2, 1] * CD39

    # Diurnal
    if input.switches.SW[7] != 0
        T71 = (p[12] * lpoly.PLG[3, 2]) * CD14 * input.switches.SWC[5]
        T72 = (p[13] * lpoly.PLG[3, 2]) * CD14 * input.switches.SWC[5]
        t[7] = f2 * (
            (
                p[4] * lpoly.PLG[2, 2] +
                    p[5] * lpoly.PLG[4, 2] +
                    p[28] * lpoly.PLG[6, 2] + T71
            ) * lpoly.CTLOC + (
                p[7] * lpoly.PLG[2, 2] +
                    p[8] * lpoly.PLG[4, 2] +
                    p[29] * lpoly.PLG[6, 2] + T72
            ) * lpoly.STLOC
        )
    end

    # Semidiurnal
    if input.switches.SW[8] != 0
        T81 = (p[24] * lpoly.PLG[4, 3] + p[36] * lpoly.PLG[6, 3]) *
            CD14 * input.switches.SWC[5]
        T82 = (p[34] * lpoly.PLG[4, 3] + p[37] * lpoly.PLG[6, 3]) *
            CD14 * input.switches.SWC[5]
        t[8] = f2 * (
            (p[6] * lpoly.PLG[3, 3] + p[42] * lpoly.PLG[5, 3] + T81) *
                lpoly.C2TLOC +
                (p[9] * lpoly.PLG[3, 3] + p[43] * lpoly.PLG[5, 3] + T82) *
                lpoly.S2TLOC
        )
    end

    # Terdiurnal
    if input.switches.SW[14] != 0
        t[14] = f2 * (
            (
                p[40] * lpoly.PLG[4, 4] +
                    (p[94] * lpoly.PLG[5, 4] + p[47] * lpoly.PLG[7, 4]) *
                    CD14 * input.switches.SWC[5]
            ) * lpoly.S3TLOC +
                (
                p[41] * lpoly.PLG[4, 4] +
                    (p[95] * lpoly.PLG[5, 4] + p[49] * lpoly.PLG[7, 4]) *
                    CD14 * input.switches.SWC[5]
            ) * lpoly.C3TLOC
        )
    end

    if SW9 != -1
        lpoly.APD = input.ap - 4.0
        p44 = p[44]
        p45 = p[45]
        if p44 < 0
            p44 = 1.0e-5
        end
        lpoly.APDF = lpoly.APD + (p45 - 1.0) *
            (lpoly.APD + (exp(-p44 * lpoly.APD) - 1.0) / p44)
        if input.switches.SW[9] != 0
            t[9] = lpoly.APDF * (
                p[33] + p[46] * lpoly.PLG[3, 1] + p[35] * lpoly.PLG[5, 1] +
                    (
                    p[101] * lpoly.PLG[2, 1] + p[102] * lpoly.PLG[4, 1] +
                        p[103] * lpoly.PLG[6, 1]
                ) * CD14 * input.switches.SWC[5] +
                    (
                    p[122] * lpoly.PLG[2, 2] + p[123] * lpoly.PLG[4, 2] +
                        p[124] * lpoly.PLG[6, 2]
                ) * input.switches.SWC[7] *
                    cos(HR * (input.stl - p[125]))
            )
        end
    else
        if p[52] != 0
            exp1 = min(
                0.99999,
                exp(-10800.0 * abs(p[52]) / (1.0 + p[139] * (45.0 - abs(input.lat))))
            )
            if p[25] < 1.0e-4
                p[25] = 1.0e-4
            end

            lpoly.APT[1] = SG0(exp1)

            if input.switches.SW[9] != 0
                t[9] = lpoly.APT[1] * (
                    p[51] + p[97] * lpoly.PLG[3, 1] + p[55] * lpoly.PLG[5, 1] +
                        (
                        p[126] * lpoly.PLG[2, 1] + p[127] * lpoly.PLG[4, 1] +
                            p[128] * lpoly.PLG[6, 1]
                    ) * CD14 * input.switches.SWC[5] +
                        (
                        p[129] * lpoly.PLG[2, 2] + p[130] * lpoly.PLG[4, 2] +
                            p[131] * lpoly.PLG[6, 2]
                    ) * input.switches.SWC[7] *
                        cos(HR * (input.stl - p[132]))
                )
            end

        end
    end

    if input.switches.SW[10] == 0
        Tinf = p[31]
        for i in 1:NSW
            Tinf += abs(input.switches.SW[i]) * t[i]
        end
        return Tinf
    end

    if input.switches.SW[11] != 0
        t[11] = (1.0 + p[81] * lpoly.DFA * input.switches.SWC[1]) *
            (
            (
                p[65] * lpoly.PLG[3, 2] + p[66] * lpoly.PLG[5, 2] + p[67] * lpoly.PLG[7, 2]
                    + p[104] * lpoly.PLG[2, 2] + p[105] * lpoly.PLG[4, 2] + p[106] * lpoly.PLG[6, 2]
                    + input.switches.SWC[5] * (p[110] * lpoly.PLG[2, 2] + p[111] * lpoly.PLG[4, 2] + p[112] * lpoly.PLG[6, 2]) * CD14
            ) *
                cos(DGTR * input.long)
                + (
                p[91] * lpoly.PLG[3, 2] + p[92] * lpoly.PLG[5, 2] + p[93] * lpoly.PLG[7, 2]
                    + p[107] * lpoly.PLG[2, 2] + p[108] * lpoly.PLG[4, 2] + p[109] * lpoly.PLG[6, 2]
                    + input.switches.SWC[5] * (p[113] * lpoly.PLG[2, 2] + p[114] * lpoly.PLG[4, 2] + p[115] * lpoly.PLG[6, 2]) * CD14
            ) *
                sin(DGTR * input.long)
        )
    end

    # UT and mixed UT, longitude
    if input.switches.SW[12] != 0
        t[12] = (1.0 + p[96] * lpoly.PLG[2, 1]) *
            (1.0 + p[82] * lpoly.DFA * input.switches.SWC[1]) *
            (1.0 + p[120] * lpoly.PLG[2, 1] * CD14 * input.switches.SWC[5]) *
            (
            (p[69] * lpoly.PLG[2, 1] + p[70] * lpoly.PLG[4, 1] + p[71] * lpoly.PLG[6, 1]) *
                cos(SR * (input.sec - p[72]))
        )
        t[12] += input.switches.SWC[11] *
            (p[77] * lpoly.PLG[4, 3] + p[78] * lpoly.PLG[6, 3] + p[79] * lpoly.PLG[8, 3]) *
            cos(
            SR * (input.sec - p[80]) +
                2.0 * DGTR * input.long
        ) * (1.0 + p[138] * lpoly.DFA * input.switches.SWC[1])
    end

    # UT, longitude magnetic activity
    if input.switches.SW[13] != 0
        if SW9 != -1
            t[13] = lpoly.APDF * input.switches.SWC[11] *
                (1.0 + p[121] * lpoly.PLG[2, 1]) *
                (
                (p[61] * lpoly.PLG[3, 2] + p[62] * lpoly.PLG[5, 2] + p[63] * lpoly.PLG[7, 2]) *
                    cos(DGTR * (input.long - p[64]))
            ) +
                lpoly.APDF * input.switches.SWC[11] * input.switches.SWC[5] *
                (p[116] * lpoly.PLG[2, 2] + p[117] * lpoly.PLG[4, 2] + p[118] * lpoly.PLG[6, 2]) *
                CD14 * cos(DGTR * (input.long - p[119])) +
                lpoly.APDF * input.switches.SWC[12] *
                (p[84] * lpoly.PLG[2, 1] + p[85] * lpoly.PLG[4, 1] + p[86] * lpoly.PLG[6, 1]) *
                cos(SR * (input.sec - p[76]))
        else
            if p[52] != 0
                t[13] = lpoly.APT[1] * input.switches.SWC[11] *
                    (1.0 + p[133] * lpoly.PLG[2, 1]) *
                    (
                    (p[53] * lpoly.PLG[3, 2] + p[99] * lpoly.PLG[5, 2] + p[68] * lpoly.PLG[7, 2]) *
                        cos(DGTR * (input.long - p[98]))
                ) +
                    lpoly.APT[1] * input.switches.SWC[11] * input.switches.SWC[5] *
                    (p[134] * lpoly.PLG[2, 2] + p[135] * lpoly.PLG[4, 2] + p[136] * lpoly.PLG[6, 2]) *
                    CD14 * cos(DGTR * (input.long - p[137])) +
                    lpoly.APT[1] * input.switches.SWC[12] *
                    (p[56] * lpoly.PLG[2, 1] + p[57] * lpoly.PLG[4, 1] + p[58] * lpoly.PLG[6, 1]) *
                    cos(SR * (input.sec - p[59]))
            end
        end
    end

    Tinf = p[31]
    for i in 1:NSW
        Tinf += abs(input.switches.SW[i]) * t[i]
    end
    return Tinf
    #=
    =#
end

#=
=#
