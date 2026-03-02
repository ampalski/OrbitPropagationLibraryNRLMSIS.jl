#=
SWITCHES: The following is for test and special purposes:
     
   TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW),
   WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
   FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
   FOR THE FOLLOWING VARIATIONS
          1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
          3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
          5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
          7 - DIURNAL               8 - SEMIDIURNAL
          9 - DAILY AP             10 - ALL UT/LONG EFFECTS
         11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
         13 - MIXED AP/UT/LONG     14 - TERDIURNAL
         15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
         16 - ALL TINF VAR         17 - ALL TLB VAR
         18 - ALL TN1 VAR           19 - ALL S VAR
         20 - ALL TN2 VAR           21 - ALL NLB VAR
         22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
=#
export OplMsis_Switches
mutable struct OplMsis_Switches
    SW::MVector{25, Int}
    SWC::MVector{25, Int}
    IMR::Bool
end

export getMsisSwitches
"""
    switches = getMsisSwitches(; SW=SA[ones(Int, 25)...], SA[ones(Int, 25)...])

Return a switches struct to control execution of the model.
"""
function getMsisSwitches(;
        SW = MVector{25}(ones(Int, 25)),
        SWC = MVector{25}(ones(Int, 25)),
        IMR = true,
    )
    return OplMsis_Switches(SW, SWC, IMR)
end

export OplMsis_Input
mutable struct OplMsis_Input
    iyd::Int # day of the year, integer between 1 and 365(366)
    sec::Float64 # Number of seconds into the current day (UT)
    alt::Float64 # Altitude, in km
    lat::Float64 # Geodetic Latitude, in radians
    long::Float64 # Longitude, in radians
    stl::Float64 # Local apparent satellite time
    f107a::Float64 # 81 day centered average of F10.7 flux data
    f107::Float64 # Daily F10.7 flux for previous day
    ap::Union{Float64, SVector{7, Float64}} # magnetic index (daily)
    # or, when SW[9] = -1.0, the input Ap array
    mass::Int # mass number -- only density for selected gas is calculated.
    # 0 for temperature, 48 for all
    switches::OplMsis_Switches
end

#=
INPUT VARIABLES:
   IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
         (Year ignored in current model)
   SEC - UT(SEC)
   ALT - ALTITUDE(KM)
   GLAT - GEODETIC LATITUDE(DEG)
   GLONG - GEODETIC LONGITUDE(DEG)
   STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
   F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
   F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
   AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
      - ARRAY CONTAINING:
        (1) DAILY AP
        (2) 3 HR AP INDEX FOR CURRENT TIME
        (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
        (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
        (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
        (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
               TO CURRENT TIME
        (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
               TO CURRENT TIME
   MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
            CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
            MASS 17 IS Anomalous O ONLY.)

NOTES ON INPUT VARIABLES: 
   UT, Local Time, and Longitude are used independently in the
   model and are not of equal importance for every situation.  
   For the most physically realistic calculation these three
   variables should be consistent (STL=SEC/3600+GLONG/15).
   The Equation of Time departures from the above formula
   for apparent local time can be included if available but
   are of minor importance.

   F107 and F107A values used to generate the model correspond
   to the 10.7 cm radio flux at the actual distance of the Earth
   from the Sun rather than the radio flux at 1 AU.
=#
export getMsisInput
function getMsisInput(
        jd::JulianDate,
        alt::Float64,
        lat::Float64,
        long::Float64;
        mass::Int = 48,
        switches::OplMsis_Switches = getMsisSwitches(),
    )
    jd = convert_jd(jd, :UT1)
    datevec = jdate2datevec(jd)
    jd0, _ = datevec2jdate([datevec[1], 1, 1, 0, 0, 0], system = :UT1)
    jd_sum = sum(jd.epoch)
    jd0_sum = sum(jd0.epoch)

    iyd = floor(Int, jd_sum - jd0_sum + 1)
    sec = (jd_sum - floor(jd_sum)) * 86400.0

    # Local time as in the MSISE documentation
    stl = sec / 3600 + long / 15

    f107a = getF107a(jd)
    f107 = getF107(jd)

    if switches.SW[9] == -1
        ap = getApInput(jd)
    else
        ap = getAp(jd)
    end

    input = OplMsis_Input(
        iyd, sec,
        alt, lat, long, stl,
        f107a, f107, ap,
        mass, switches
    )

    return input
end

export gtd7d!
"""
    gtd7d!(output, input)

      SUBROUTINE GTD7D(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,
      D,T)

Provide Effective Total Mass Density, to include "anomalous oxygen".

NRLMSISE-00
   This subroutine provides Effective Total Mass Density for
   output D(6) which includes contributions from "anomalous
   oxygen" which can affect satellite drag above 500 km.  This
   subroutine is part of the distribution package for the 
   Neutral Atmosphere Empirical Model from the surface to lower
   exosphere.  See subroutine GTD7 for more extensive comments.

# Inputs
- `input::OplMsis_Input`: Input struct, best constructed with `getMsisInput`

# Outputs (in place)
- `output::OplMsis_Output`: Output struct, best constructed with `getMsisOutput`
"""
function gtd7d!(output::OplMsis_Output, input::OplMsis_Input)

    gtd7!(output, input)
    if input.mass == 48
        output.D[6] = 1.66e-24 *
            (
            4.0 * output.D[1] +
                16.0 * output.D[2] +
                28.0 * output.D[3] +
                32.0 * output.D[4] +
                40.0 * output.D[5] +
                1.0 * output.D[7] +
                14.0 * output.D[8] +
                16.0 * output.D[9]
        )
    end

    if input.switches.IMR
        output.D[6] /= 1000
    end
    return nothing
end

#=
------------------------------------------------------------------
SUBROUTINE GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)

NRLMSISE-00
-----------
   Neutral Atmosphere Empirical Model from the surface to lower
   exosphere

   NEW FEATURES:
     *Extensive satellite drag database used in model generation
     *Revised O2 (and O) in lower thermosphere
     *Additional nonlinear solar activity term
     *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
      At high altitudes (> 500 km), hot atomic oxygen or ionized
      oxygen can become appreciable for some ranges of subroutine
      inputs, thereby affecting drag on satellites and debris. We
      group these species under the term "anomalous oxygen," since
      their individual variations are not presently separable with
      the drag data used to define this model component.

   SUBROUTINES FOR SPECIAL OUTPUTS:
   
   HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY 
   (SUBROUTINE GTD7D, OUTPUT D(6))
      For atmospheric drag calculations at altitudes above 500 km,
      call SUBROUTINE GTD7D to compute the "effective total mass
      density" by including contributions from "anomalous oxygen."
      See "NOTES ON OUTPUT VARIABLES" below on D(6).

   PRESSURE GRID (SUBROUTINE GHP7)
     See subroutine GHP7 to specify outputs at a pressure level
     rather than at an altitude.

   OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)

INPUT VARIABLES:
   IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
         (Year ignored in current model)
   SEC - UT(SEC)
   ALT - ALTITUDE(KM)
   GLAT - GEODETIC LATITUDE(DEG)
   GLONG - GEODETIC LONGITUDE(DEG)
   STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
   F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
   F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
   AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
      - ARRAY CONTAINING:
        (1) DAILY AP
        (2) 3 HR AP INDEX FOR CURRENT TIME
        (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
        (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
        (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
        (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
               TO CURRENT TIME
        (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
               TO CURRENT TIME
   MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
            CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
            MASS 17 IS Anomalous O ONLY.)

NOTES ON INPUT VARIABLES: 
   UT, Local Time, and Longitude are used independently in the
   model and are not of equal importance for every situation.  
   For the most physically realistic calculation these three
   variables should be consistent (STL=SEC/3600+GLONG/15).
   The Equation of Time departures from the above formula
   for apparent local time can be included if available but
   are of minor importance.

   F107 and F107A values used to generate the model correspond
   to the 10.7 cm radio flux at the actual distance of the Earth
   from the Sun rather than the radio flux at 1 AU. The following
   site provides both classes of values:
   ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/

   F107, F107A, and AP effects are neither large nor well
   established below 80 km and these parameters should be set to
   150., 150., and 4. respectively.

OUTPUT VARIABLES:
   D(1) - HE NUMBER DENSITY(CM-3)
   D(2) - O NUMBER DENSITY(CM-3)
   D(3) - N2 NUMBER DENSITY(CM-3)
   D(4) - O2 NUMBER DENSITY(CM-3)
   D(5) - AR NUMBER DENSITY(CM-3)                       
   D(6) - TOTAL MASS DENSITY(GM/CM3)
   D(7) - H NUMBER DENSITY(CM-3)
   D(8) - N NUMBER DENSITY(CM-3)
   D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
   T(1) - EXOSPHERIC TEMPERATURE
   T(2) - TEMPERATURE AT ALT

NOTES ON OUTPUT VARIABLES:
   TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.) 

   O, H, and N are set to zero below 72.5 km

   T(1), Exospheric temperature, is set to global average for
   altitudes below 120 km. The 120 km gradient is left at global
   average value for altitudes below 72 km.

   D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
   and GTD7D

     SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
     species labeled by indices 1-5 and 7-8 in output variable D.
     This includes He, O, N2, O2, Ar, H, and N but does NOT include
     anomalous oxygen (species index 9).

     SUBROUTINE GTD7D -- D(6) is the "effective total mass density
     for drag" and is the sum of the mass densities of all species
     in this model, INCLUDING anomalous oxygen.
   

   To get current values of SW: CALL TRETRV(SW)
=#
export gtd7!
"""
    gtd7!(output, input)

      SUBROUTINE GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)

# Inputs
- `input::OplMsis_Input`: Input struct, best constructed with `getMsisInput`

# Outputs (in place)
- `output::OplMsis_Output`: Output struct, best constructed with `getMsisOutput`
"""
function gtd7!(output::OplMsis_Output, input::OplMsis_Input)
    MN3 = 5
    ZN3 = SA[32.5, 20.0, 15.0, 10.0, 0.0]
    MN2 = 4
    ZN2 = SA[72.5, 55.0, 45.0, 32.5]
    ZMIX = 62.5

    gts3c = GTS3C(
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )
    meso7 = MESO7(
        MVector{5}(zeros(5)), MVector{4}(zeros(4)),
        MVector{5}(zeros(5)), MVector{2}(zeros(2)),
        MVector{2}(zeros(2)), MVector{2}(zeros(2)),
    )
    dmix = DMIX(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    lpoly = LPOLY(
        zeros(9, 4), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, zeros(4), 0.0
    )
    parmb = PARMB(0.0, 0.0)

    # Latitude variation of gravity (none for SW[2] = 0)
    xlat = input.lat
    if input.switches.SW[2] == 0
        xlat = 45.0
    end
    _glatf!(parmb, xlat)
    xmm = PDM[5, 3]

    # Thermosphere/mesosphere (above zn2[1])
    altt = max(input.alt, ZN2[1])
    mss = input.mass

    # Only calculate N2 in thermosphere if alt in mixed region
    if input.alt < ZMIX && input.mass > 0
        mss = 28
    end

    # Only calculate thermosphere if altitude above zn2[1] in mesosphere
    if input.alt > ZN2[1]
        ds, ts = gts7!(input, gts3c, meso7, dmix, lpoly, parmb)
        DM28M = DM28
        if input.switches.IMR == 1
            DM28M = DM28 * 1.0e6
        end
        mssl = mss
        output.T[1] = ts[1]
        output.T[2] = ts[2]
        for i in 1:9
            output.D[i] = ds[i]
        end
        return nothing
    end

    # Lower mesosphere / upper stratosphere
    #  Temperature at nodes and gradients at end nodes
    #  Inverse temperature a linear function of spherical harmonics

    meso7.TGN2[1] = meso7.TGN1[2]
    meso7.TN2[1] = meso7.TN1[5]
    meso7.TN2[2] = PMA[1, 1] * PAVGM[1] /
        (1.0 - input.switches.SW[20] * _glob7s(view(PMA, :, 1), input, lpoly))
    meso7.TN2[3] = PMA[1, 2] * PAVGM[2] /
        (1.0 - input.switches.SW[20] * _glob7s(view(PMA, :, 2), input, lpoly))
    meso7.TN2[4] = PMA[1, 3] * PAVGM[3] /
        (
        1.0 - input.switches.SW[20] * input.switches.SW[22] *
            _glob7s(view(PMA, :, 2), input, lpoly)
    )
    meso7.TGN2[2] = PAVGM[9] * PMA[1, 10] * (
        1.0 + input.switches.SW[20] * input.switches.SW[22] *
            _glob7s(view(PMA, :, 10), input, lpoly)
    ) * meso7.TN2[4]^2 / (PMA[1, 3] * PAVGM[3])^2
    meso7.TN3[1] = meso7.TN2[4]

    # Lower stratospher and troposphere
    #     Temperature at nodes and gradients at end nodes
    #     Inverse temperature a linear function of spherical harmonics
    #     Only calculate nodes if input changed
    if !(input.alt > ZN3[1])
        meso7.TGN3[1] = meso7.TGN2[2]
        meso7.TN3[2] = PMA[1, 4] * PAVGM[4] /
            (1.0 - input.switches.SW[22] * _glob7s(view(PMA, :, 4), input, lpoly))
        meso7.TN3[3] = PMA[1, 5] * PAVGM[5] /
            (1.0 - input.switches.SW[22] * _glob7s(view(PMA, :, 5), input, lpoly))
        meso7.TN3[4] = PMA[1, 6] * PAVGM[6] /
            (1.0 - input.switches.SW[22] * _glob7s(view(PMA, :, 6), input, lpoly))
        meso7.TN3[5] = PMA[1, 7] * PAVGM[7] /
            (1.0 - input.switches.SW[22] * _glob7s(view(PMA, :, 7), input, lpoly))
        meso7.TGN3[2] = PMA[1, 8] * PAVGM[8] *
            (1.0 - input.switches.SW[22] * _glob7s(view(PMA, :, 8), input, lpoly)) *
            meso7.TN3[5]^2 / (PMA[1, 7] * PAVGM[7])^2
    end

    if input.mass == 0
        tz, dd = _densm(
            input.alt, 1.0, 0.0, 0.0, MN3, ZN3, meso7.TN3,
            meso7.TGN3, MN2, ZN2, meso7.TN2, meso7.TGN2
        )
        output.T[2] = tz
        return
    end

    # Linear transition to full mixing below ZN2[1]
    DMC = 0.0
    if input.alt > ZMIX
        DMC = 1.0 - (ZN2[1] - input.alt) / (ZN2[1] - ZMIX)
    end
    DZ28 = ds[3]

    # N2 density
    DMR = ds[3] / DM28M - 1.0
    tz, output.D[3] = _densm(
        input.alt, DM28M, xmm, 0.0, MN3, ZN3, meso7.TN3,
        meso7.TGN3, MN2, ZN2, meso7.TN2, meso7.TGN2
    )
    output.D[3] = output.D[3] * (1.0 + DMR * DMC)

    # HE Density
    output.D[1] = 0.0
    if !(input.mass != 4 && input.mass != 48)
        DMR = ds[1] / (DZ28 * PDM[2, 1]) - 1.0
        output.D[1] = output.D[3] * PDM[2, 1] * (1.0 + DMR * DMC)
    end

    # O Density
    output.D[2] = 0.0
    output.D[9] = 0.0

    # O2 Density
    output.D[4] = 0.0
    if !(input.mass != 32 && input.mass != 48)
        DMR = ds[4] / (DZ28 * PDM[2, 4]) - 1.0
        output.D[4] = output.D[3] * PDM[2, 4] * (1.0 + DMR * DMC)
    end

    # AR density
    output.D[5] = 0.0
    if !(input.mass != 40 && input.mass != 48)
        DMR = ds[5] / (DZ28 * PDM[2, 5]) - 1.0
        output.D[5] = output.D[3] * PDM[2, 5] * (1.0 + DMR * DMC)
    end

    # Hydrogen density
    output.D[7] = 0.0

    # Atomic nitrogen density
    output.D[8] = 0.0

    # Total mass density
    if input.mass == 48
        output.D[6] = 1.66e-24 * (
            4.0 * output.D[1] +
                16.0 * output.D[2] +
                28.0 * output.D[3] +
                32.0 * output.D[4] +
                40.0 * output.D[5] +
                output.D[7] + 14.0 * output.D[8]
        )
        if input.switches.IMR == 1
            output.D[6] /= 1000.0
        end
    end

    return nothing
end
#=
=#

"""
    gts7!()

      SUBROUTINE GTS7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)

Thermospheric portion of NRLMSISE-00

"""
function gts7!(
        input::OplMsis_Input,
        gts3c::GTS3C,
        meso7::MESO7,
        dmix::DMIX,
        lpoly::LPOLY,
        parmb::PARMB,
    )
    MT = SA[48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17]
    ALTL = SA[200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0]
    MN1 = 5
    ZN1 = [120.0, 110, 100, 90, 72.5]
    DR = 1.72142e-2
    ALPHA = SA[-0.38, 0, 0, 0, 0.17, 0, -0.38, 0, 0]
    D = MVector{9}(zeros(9))
    T = MVector{2}(zeros(2))

    yrd = input.iyd
    za = PDL[16, 2]
    ZN1[1] = za

    # Tinf variations not important below za or zn1[1]
    Tinf = 0.0
    if input.alt > ZN1[1]
        Tinf = PTM[1] * PT[1] * (
            1.0 + input.switches.SW[16] *
                _globe7!(input, PT, lpoly)
        )
    else
        Tinf = PTM[1] * PT[1]
    end
    T[1] = Tinf

    # Gradient variations not important below ZN1[5]
    if input.alt > ZN1[5]
        g0 = PTM[4] * PS[1] * (
            1.0 + input.switches.SW[19] *
                _globe7!(input, PS, lpoly)
        )
    else
        g0 = PTM[4] * PS[1]
    end

    TLB = PTM[2] * (
        1.0 + input.switches.SW[17] *
            _globe7!(input, view(PD, :, 4), lpoly)
    ) * PD[1, 4]
    S = g0 / (Tinf - TLB)

    # Lower thermosphere temp variations not significant for density above
    # 300 km

    if input.alt < 300
        meso7.TN1[2] = PTM[7] * PTL[1, 1] /
            (1.0 - input.switches.SW[18] * _glob7s(view(PTL, :, 1), input, lpoly))
        meso7.TN1[3] = PTM[3] * PTL[1, 2] /
            (1.0 - input.switches.SW[18] * _glob7s(view(PTL, :, 2), input, lpoly))
        meso7.TN1[4] = PTM[8] * PTL[1, 3] /
            (1.0 - input.switches.SW[18] * _glob7s(view(PTL, :, 3), input, lpoly))
        meso7.TN1[5] = PTM[5] * PTL[1, 4] / (
            1.0 - input.switches.SW[18] *
                input.switches.SW[20] * _glob7s(view(PTL, :, 4), input, lpoly)
        )
        meso7.TGN1[2] = PTM[9] * PMA[1, 9] * (
            1.0 + input.switches.SW[18] *
                input.switches.SW[20] * _glob7s(view(PMA, :, 9), input, lpoly)
        ) * meso7.TN1[5] * meso7.TN1[5] / (PTM[5] * PTL[1, 4])^2
    else
        meso7.TN1[2] = PTM[7] * PTL[1, 1]
        meso7.TN1[3] = PTM[3] * PTL[1, 2]
        meso7.TN1[4] = PTM[8] * PTL[1, 3]
        meso7.TN1[5] = PTM[5] * PTL[1, 4]
        meso7.TGN1[2] = PTM[9] * PMA[1, 9] * meso7.TN1[5] * meso7.TN1[5] /
            (PTM[5] * PTL[1, 4])^2
    end

    if input.mass == 0
        if input.switches.IMR
            for i in 1:9
                D[i] *= 1.0e6
            end
            D[6] /= 1000
        end
        return (D, T)
    end

    #N2 variation factor at Zlb
    G28 = input.switches.SW[21] * _globe7!(input, view(PD, :, 3), lpoly)

    # Variation of turbopause height
    zhf = PDL[25, 2]
    T[1] = Tinf
    xmm = PDM[5, 3]
    z = input.alt

    skip = 0
    for j in 1:11
        if input.mass == MT[j]
            skip = j
            break
        end
    end
    if skip == 0
        if input.switches.IMR
            for i in 1:9
                D[i] *= 1.0e6
            end
            D[6] /= 1000
        end
        return (D, T)
    end
    # NOTE: Skip here is `J` in the original, for the GOTO between 17 & 20.

    # N2 Density
    if !(z > ALTL[6] && input.mass != 28 && input.mass != 48)
        # Diffusive dnesity at Zlb
        DB28 = PDM[1, 3] * exp(G28) * PD[1, 3]
        # Diffusive density at Alt
        T[2], D[3] = _densu!(
            z, DB28, Tinf, TLB, 28.0, ALPHA[3], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        DD = D[3]
        # Turbopause
        ZH28 = PDM[3, 3] * zhf
        ZHM28 = PDM[4, 3] * PDL[6, 2]
        XMD = 28.0 - xmm
        # Mixed density at Zlb
        tz, B28 = _densu!(
            ZH28, DB28, Tinf, TLB, XMD, ALPHA[3] - 1, 0.0, PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )

        if !(z > ALTL[3] || input.switches.SW[15] == 0)
            # Mixed density at altitude
            tz, DM28 = _densu!(
                z, B28, Tinf, TLB, xmm, ALPHA[3], tz, PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            # Net density at altitude
            D[3] = _dnet(D[3], DM28, ZHM28, xmm, 28.0)
        end
    end

    if skip == 1 || skip == 3
        # He Density

        # Density variation factor at Zlb
        G4 = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 1), lpoly)

        # Diffusive density at Zlb
        DB04 = PDM[1, 1] * exp(G4) * PD[1, 1]

        # Diffusive density at altitude
        T[2], D[1] = _densu!(
            z, DB04, Tinf, TLB, 4.0, ALPHA[1], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        DD = D[1]

        if z > ALTL[1] || input.switches.SW[15] == 0
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end
        # Turbopause
        ZH04 = PDM[3, 1]
        # Mixed density at Zlb
        T[2], B04 = _densu!(
            ZH04, DB04, Tinf, TLB, 4.0 - xmm, ALPHA[1] - 1.0, T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        # Mixed density at Alt
        T[2], DM04 = _densu!(
            z, B04, Tinf, TLB, xmm, 0.0, T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        ZHM04 = ZHM28

        # Net density at Alt
        D[1] = _dnet(D[1], DM04, ZHM04, xmm, 4.0)

        # Correction to specified mixing ratio at ground
        RL = log(B28 * PDM[2, 1] / B04)
        ZC04 = PDM[5, 1] * PDL[1, 2]
        HC04 = PDM[6, 1] * PDL[2, 2]

        # Net density corrected at alt
        D[1] = D[1] * _ccor(z, RL, HC04, ZC04)

        if input.mass != 48
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end
    end
    #  GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
    o_included = [1, 3, 4, 9]
    if skip in o_included
        # O Density

        # Density variation factor at Zlb
        G16 = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 2), lpoly)

        # Diffusive density at Zlb
        DB16 = PDM[1, 2] * exp(G16) * PD[1, 2]
        # Diffusive density at Alt
        T[2], D[2] = _densu!(
            z, DB16, Tinf, TLB, 16.0, ALPHA[2], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        DD = D[2]

        if !(z > ALTL[2] || input.switches.SW[15] == 0)
            # Turbopause
            ZH16 = PDM[3, 2]

            # Mixed density at Zlb
            T[2], B16 = _densu!(
                ZH16, DB16, Tinf, TLB, 16 - xmm, ALPHA[2] - 1.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            # Mixed density at Alt
            T[2], DM16 = _densu!(
                z, B16, Tinf, TLB, xmm, 0.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            ZHM16 = ZHM28
            # Net density at Alt
            D[2] = _dnet(D[2], DM16, ZHM16, xmm, 16.0)

            RL = PDM[2, 2] * PDL[17, 2] *
                (1.0 + input.switches.SW[1] * PDL[24, 1] * (input.f107a - 150.0))
            HC16 = PDM[6, 2] * PDL[4, 2]
            ZC16 = PDM[5, 2] * PDL[3, 2]
            HC216 = PDM[6, 2] * PDL[5, 2]
            D[2] = D[2] * _ccor2(z, RL, HC16, ZC16, HC216)

            # Chemistry correction
            HCC16 = PDM[8, 2] * PDL[14, 2]
            ZCC16 = PDM[7, 2] * PDL[13, 2]
            RC16 = PDM[4, 2] * PDL[15, 2]

            # Net density corrected at Alt
            D[2] = D[2] * _ccor(z, RC16, HCC16, ZCC16)
        end
        if input.mass != 48 && input.mass != 49
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end
    end
    #  GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
    # O2
    o2_included = [1, 3, 4, 6, 9]
    if skip in o2_included
        # Density variation factor at Zlb
        G32 = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 5), lpoly)

        # Diffusive density at Zlb
        DB32 = PDM[1, 4] * exp(G32) * PD[1, 5]

        # Diffusive density at Alt
        T[2], D[4] = _densu!(
            z, DB32, Tinf, TLB, 32.0, ALPHA[4], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        if input.mass == 49
            DD += 2.0 * D[4]
        else
            DD = D[4]
        end

        if input.switches.SW[15] == 0
            if input.mass != 48
                if input.switches.IMR
                    for i in 1:9
                        D[i] *= 1.0e6
                    end
                    D[6] /= 1000
                end
                return (D, T)
            end
        end

        if z <= ALTL[4]
            # Turbopause
            ZH32 = PDM[3, 4]

            # Mixed density at Zlb
            T[2], B32 = _densu!(
                ZH32, DB32, Tinf, TLB, 32.0 - xmm, ALPHA[4] - 1.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )

            # Mixed density at Alt
            T[2], DM32 = _densu!(
                z, B32, Tinf, TLB, xmm, 0.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            ZHM32 = ZHM28

            # Net density at Alt
            D[4] = _dnet(D[4], DM32, ZHM32, xmm, 32.0)

            # Correction to specified mixing ratio at ground
            RL = log(B28 * PDM[2, 4] / B32)
            HC32 = PDM[6, 4] * PDL[8, 2]
            ZC32 = PDM[5, 4] * PDL[7, 2]
            D[4] = D[4] * _ccor(z, RL, HC32, ZC32)
        end
        # Correction for general departure from diffusive equilibrium above Zlb
        HCC32 = PDM[8, 4] * PDL[23, 2]
        HCC232 = PDM[8, 4] * PDL[23, 1]
        ZCC32 = PDM[7, 4] * PDL[22, 2]
        RC32 = PDM[4, 4] * PDL[24, 2] *
            (1.0 + input.switches.SW[1] * PDL[24, 2] * (input.f107a - 150.0))

        # Net density corrected at Alt
        D[4] *= _ccor2(z, RC32, HCC32, ZCC32, HCC232)

        if input.mass != 48
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end
    end
    #  GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
    # AR
    AR_included = [1, 3, 4, 6, 7, 9]
    if skip in AR_included
        # Density variation factor at Zlb
        G40 = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 6), lpoly)

        # Diffusive density at Zlb
        DB40 = PDM[1, 5] * exp(G40) * PD[1, 6]

        # Diffusive density at Alt
        T[2], D[5] = _densu!(
            z, DB40, Tinf, TLB, 40.0, ALPHA[5], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        DD = D[5]

        if !(z > ALTL[5] || input.switches.SW[15] == 0)
            # Turbopause
            ZH40 = PDM[3, 5]

            # Mixed density at Zlb
            T[2], B40 = _densu!(
                ZH40, DB40, Tinf, TLB, 40.0 - xmm, ALPHA[5] - 1.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )

            # Mixed density at Alt
            T[2], DM40 = _densu!(
                ZH40, B40, Tinf, TLB, xmm, 0.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            ZHM40 = ZHM28

            # Net density at altitude
            D[5] = _dnet(D[5], DM40, ZHM40, xmm, 40.0)

            # Correction to specified mixing ratio at ground
            RL = log(B28 * PDM[2, 5] / B40)
            HC40 = PDM[6, 5] * PDL[10, 2]
            ZC40 = PDM[5, 5] * PDL[9, 2]

            # Net density corrected at alt
            D[5] *= _ccor(z, RL, HC40, ZC40)

        end
        if input.mass != 48
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end

    end

    #  GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
    # Hydrogen
    H_included = [1, 3, 4, 6, 7, 8, 9]
    if skip in H_included
        # Density variation factor at Zlb
        G1 = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 7), lpoly)

        # Diffusive density at Zlb
        DB01 = PDM[1, 6] * exp(G1) * PD[1, 7]

        # Diffusive density at Alt
        T[2], D[7] = _densu!(
            z, DB01, Tinf, TLB, 1.0, ALPHA[7], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        DD = D[7]

        if !(z > ALTL[7] || input.switches.SW[15] == 0)
            # Turbopuase
            ZH01 = PDM[3, 6]
            # Mixed density at Zlb
            T[2], B01 = _densu!(
                ZH01, DB01, Tinf, TLB, 1.0 - xmm, ALPHA[7] - 1.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )

            # Mixed density at Alt
            T[2], DM01 = _densu!(
                z, B01, Tinf, TLB, xmm, 0.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            ZHM01 = ZHM28

            # Net density at alt
            D[7] = _dnet(D[7], DM01, ZHM01, xmm, 1.0)

            # Correction to specified mixing ratio at ground
            RL = log(B28 * PDM[2, 6] * abs(PDL[18, 2]) / B01)
            HC01 = PDM[6, 6] * PDL[12, 2]
            ZC01 = PDM[5, 6] * PDL[11, 2]
            D[7] = D[7] * _ccor(z, RL, HC01, ZC01)

            # Chemistry correction
            HCC01 = PDM[8, 6] * PDL[20, 2]
            ZCC01 = PDM[7, 6] * PDL[19, 2]
            RC01 = PDM[4, 6] * PDL[21, 2]

            # Net density corrected at Alt
            D[7] = D[7] * _ccor(z, RC01, HCC01, ZCC01)
        end
        if input.mass != 48
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end

    end
    #  GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
    # Atomic Nitrogen
    N_included = [1, 3, 4, 6, 7, 8, 9, 10]
    if skip in N_included
        # Density variation factor at Zlb
        G14 = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 8), lpoly)

        # Diffusive density at Zlb
        DB14 = PDM[1, 7] * exp(G14) * PD[1, 8]

        # Diffusive density at Alt
        T[2], D[8] = _densu!(
            z, DB14, Tinf, TLB, 14.0, ALPHA[8], T[2], PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        DD = D[8]

        if !(z > ALTL[8] || input.switches.SW[15] == 0)
            # Turbopause
            ZH14 = PDM[3, 7]

            # Mixed density at Zlb
            T[2], B14 = _densu!(
                ZH14, DB14, Tinf, TLB, 14.0 - xmm, ALPHA[8] - 1.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )

            # Mixed density at Alt
            T[2], DM14 = _densu!(
                z, B14, Tinf, TLB, xmm, 0.0, T[2], PTM[6], S,
                MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
            )
            ZHM14 = ZHM28

            # Net density at Alt
            D[8] = _dnet(D[8], DM14, ZHM14, xmm, 14.0)

            # Correction to specified mixing ratio at ground
            RL = log(B28 * PDM[2, 7] * abs(PDL[3, 1]) / B14)
            HC14 = PDM[6, 7] * PDL[2, 1]
            ZC14 = PDM[5, 7] * PDL[1, 1]
            D[8] = D[8] * _ccor(z, RL, HC14, ZC14)

            # Chemistry correction
            HCC14 = PDM[8, 7] * PDL[5, 1]
            ZCC14 = PDM[7, 7] * PDL[4, 1]
            RC14 = PDM[4, 7] * PDL[6, 1]

            # Net density corrected at Alt
            D[8] = D[8] * _ccor(z, RC14, HCC14, ZCC14)
        end
        if input.mass != 48
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end
    end

    #  GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
    # AnomO Nitrogen
    AnomO_included = [1, 3, 4, 6, 7, 8, 9, 10, 11]
    if skip in AnomO_included
        G16H = input.switches.SW[21] *
            _globe7!(input, view(PD, :, 9), lpoly)
        DB16H = PDM[1, 8] * exp(G16H) * PD[1, 9]
        THO = PDM[10, 8] * PDL[7, 1]
        T2, DD = _densu!(
            z, DB16H, THO, THO, 16.0, ALPHA[9], T2, PTM[6], S,
            MN1, ZN1, meso7.TN1, meso7.TGN1, parmb.RE, parmb.GSURF
        )
        ZSHT = PDM[6, 8]
        ZMHO = PDM[5, 8]
        ZSHO = _scalh(ZMHO, 16.0, THO, parmb.GSURF, parmb.RE)
        D[9] = DD * exp(-ZSHT / ZSHO * (exp(-(z - ZMHO) / ZSHT) - 1.0))
        if input.mass != 48
            if input.switches.IMR
                for i in 1:9
                    D[i] *= 1.0e6
                end
                D[6] /= 1000
            end
            return (D, T)
        end
    end

    # Total Mass Density
    D[6] = 1.66e-24 * (4.0 * D[1] + 16.0 * D[2] + 28.0 * D[3] + 32.0 * D[4] + 40.0 * D[5] + D[7] + 14.0 * D[8])
    DB48 = 1.66e-24 * (4.0 * DB04 + 16.0 * DB16 + 28.0 * DB28 + 32.0 * DB32 + 40.0 * DB40 + DB01 + 14.0 * DB14)

    if input.mass != 48
        if input.switches.IMR
            for i in 1:9
                D[i] *= 1.0e6
            end
            D[6] /= 1000
        end
    end
    return D, T
end
#=
=#
