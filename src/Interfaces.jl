#=
OUTPUT VARIABLES:
   D(1) - HE NUMBER DENSITY(CM-3)
   D(2) - O NUMBER DENSITY(CM-3)
   D(3) - N2 NUMBER DENSITY(CM-3)
   D(4) - O2 NUMBER DENSITY(CM-3)
   D(5) - AR NUMBER DENSITY(CM-3)                       
   D(6) - TOTAL MASS DENSITY(GM/CM3) [includes anomalous oxygen]
   D(7) - H NUMBER DENSITY(CM-3)
   D(8) - N NUMBER DENSITY(CM-3)
   D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
   T(1) - EXOSPHERIC TEMPERATURE
   T(2) - TEMPERATURE AT ALT
=#
mutable struct OplMsis_Output
    D::MVector{9, Float64}
    T::MVector{2, Float64}
end

function getMsisOutput()
    return OplMsis_Output(MVector{9, Float64}(undef), MVector{2, Float64}(undef))
end
export getDensity
"""
    density = getDensity(jd, lat, lon, alt)

Determine the total atmospheric density at a given epoch and location.

Follows the NRLMSISE-00 Model, specifically running `GTD7D` to determine total 
density with anomalous oxygen.

# Inputs
- `jd::JulianDate`: Epoch, in OPLSOFA's JulianDate format.
- `lat::Float64`: Geodetic latitude, in degrees
- `lon::Float64`: Longitude, in degrees
- `alt::Float64`: Altitude, in km

# Output
- `ρ::Float64`: Total, local atmospheric density, in kg/m^3
"""
function getDensity(jd::JulianDate, lat::Float64, lon::Float64, alt::Float64)
    if lat > 90 || lat < -90
        error("Invalid latitude value, must be between -90 and 90 degrees")
    end
    lon = lon - 360.0 * round(lon / 360, RoundNearest)
    if alt < 0
        error("Positive altitudes are required")
    end
    if alt > 1000.0
        # TODO: determine if there's a better cutoff
        return 0.0
    end

    switches = getMsisSwitches()
    # switches.SW[9] = -1

    input = getMsisInput(jd, alt, lat, lon, switches = switches)
    #
    # input.f107 = 72.2
    # input.f107a = 71.1
    # input.ap = 4.2
    #
    output = getMsisOutput()

    gtd7d!(output, input)

    return output.D[6]
end
