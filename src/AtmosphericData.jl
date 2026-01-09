export getF107
function getF107(jd::JulianDate)
    mjd = _get_base_mjd(jd)
    if mjd < MIN_SW
        error("Invalid date, no data prior to 1957")
    end
    if mjd in SW_KEYS
        return SW[mjd].f107
    end
    #TODO implement predictive data
    if mjd > MAX_SW
        return SW[MAX_SW].f107
    end
    return 0
end

export getF107a
function getF107a(jd::JulianDate)
    mjd = _get_base_mjd(jd)
    if mjd < MIN_SW
        error("Invalid date, no data prior to 1957")
    end
    if mjd in SW_KEYS
        return SW[mjd].f107a
    end
    #TODO implement predictive data
    if mjd > MAX_SW
        return SW[MAX_SW].f107a
    end
    return 0
end

export getApAvg
function getApAvg(jd::JulianDate)
    mjd = _get_base_mjd(jd)
    if mjd < MIN_SW
        error("Invalid date, no data prior to 1957")
    end
    if mjd in SW_KEYS
        return SW[mjd].ap_avg
    end
    #TODO implement predictive data
    if mjd > MAX_SW
        return SW[MAX_SW].ap_avg
    end
    return 0
end
# getAp

export getAp
function getAp(jd::JulianDate)
    usejd = jd isa JDate ?
        sum(jdate_to_mjdate(jd).epoch) :
        sum(deepcopy(jd).epoch)

    basemjd = floor(Int, usejd)
    frac = usejd - basemjd
    rfrac = floor(Int, frac * 8) + 1

    if basemjd < MIN_SW
        error("Invalid date, no data prior to 1957")
    end
    if basemjd in SW_KEYS
        return SW[basemjd].ap[rfrac]
    end
    #TODO implement predictive data
    if basemjd > MAX_SW
        return SW[MAX_SW].ap[rfrac]
    end
end
# getApInput

export getApInput
function getApInput(jd::JulianDate)
    #=
    From the original Fortran, returns:
    C           - ARRAY CONTAINING:
    C             (1) DAILY AP
    C             (2) 3 HR AP INDEX FOR CURRENT TIME
    C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
    C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
    C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
    C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
    C                    TO CURRENT TIME
    C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
    C                    TO CURRENT TIME
    =#
    apVec = zeros(7)

    usejd = jd isa JDate ?
        sum(jdate_to_mjdate(jd).epoch) :
        sum(deepcopy(jd).epoch)

    basemjd = Int(floor(usejd))

    if basemjd < MIN_SW + 3
        error("Invalid date, no data prior to 1957")
    end
    #TODO implement predictive data
    if basemjd > MAX_SW
        basemjd = MAX_SW
    end

    frac = usejd - basemjd
    rfrac = floor(Int, frac * 8) + 1

    apVec[1] = SW[basemjd].ap_avg

    rfrac += 24
    fullList = [
        SW[basemjd - 3].ap;
        SW[basemjd - 2].ap;
        SW[basemjd - 1].ap;
        SW[basemjd].ap;
    ]
    apVec[2] = fullList[rfrac]
    apVec[3] = fullList[rfrac - 1]
    apVec[4] = fullList[rfrac - 2]
    apVec[5] = fullList[rfrac - 3]
    apVec[6] = mean(fullList[(rfrac - 11):(rfrac - 4)])
    apVec[7] = mean(fullList[(rfrac - 19):(rfrac - 12)])

    return apVec
end
