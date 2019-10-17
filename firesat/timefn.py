import datetime
import numpy as np
from numpy import datetime64


def julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day int."""
    janfeb = month < 3
    return (
        day
        + np.floor_divide(1461 * (year + 4800 - janfeb), 4)
        + np.floor_divide(367 * (month - 2 + janfeb * 12), 12)
        - np.floor_divide(3 * np.floor_divide(year + 4900 - janfeb, 100), 4)
        - 32075
    )


def julian_date(yr, mo=1, dy=1, hr=0, mn=0, sec=0.0):
    """Given a proleptic Gregorian calendar date, return a Julian date float."""
    if isinstance(yr, datetime.datetime) or isinstance(yr, datetime64):
        if isinstance(yr, datetime64):
            dt = yr.astype(datetime.datetime)
        else:
            dt = yr
        yr, mo, dy = dt.year, dt.month, dt.day
        hr, mn, sec = dt.hour, dt.minute, dt.second
        sec += dt.microsecond * (10 ** -6)
    return julian_day(yr, mo, dy) - 0.5 + (sec + mn * 60.0 + hr * 3600.0) / 86400.0


def jdt_tsince(tstart, tsince):
    """Return a vector of julian dates from tstart with points at tsince

    Args:
        tstart : float, julian date
        tsince: float (n), vector of minutes past tstart to calculate the julian date
    Output:
        jdt : float (n), vector of julian date ouputs. Can be inputted into solar functions
    References:
        Rhodes, python-sgp4/sgp4/ext.py
        Vallado, 'Revisiting Spacetrack Report #3'
    """
    return tstart + (tsince * 60.0) / 86400.0


def days2mdhms(year, days):
    """This procedure converts the day of the year, days, to the equivalent month, day, hour, minute and second.

    Args:
        year : int, year between 1900 - 2100
        days : float, julian day of the year between 0.0 - 366.0
    Outputs:
        mon : int, month, 1 .. 12
        day : int, day, 1 .. 28,29,30,31
        hr : int, hour, 0 .. 23
        min : int, minute 0 .. 59
        sec : float, second, 0.0 .. 59.999
    References:
        Vallado
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # # int array containing the number of days per month
    # lmonth = [(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
    #           (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)]
    # dayofyr = np.floor_divide(days, 1.0).astype(int)  # day of year
    # # dayofyr = int(days // 1.0)  # day of year
    # # find month and day of month
    # lpyr_idx = (year % 4) == 0
    # i = 1
    # inttemp = 0
    # while dayofyr > inttemp + lmonth[i - 1] and i < 12:
    #     inttemp = inttemp + lmonth[i - 1]
    #     i += 1
    # mon = i
    # day = dayofyr - inttemp
    # # find hours minutes and seconds
    # temp = (days - dayofyr) * 24.0
    # hr = np.floor_divide(temp, 1.0).astype(int)
    # # hr = int(temp // 1.0)
    # temp = (temp - hr) * 60.0
    # minute = np.floor_divide(temp, 1.0).astype(int)
    # # minute = int(temp // 1.0)
    # sec = (temp - minute) * 60.0
    # return mon, day, hr, minute, sec

    # non vectorized version...
    lmonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    dayofyr = int(days // 1.0)  # day of year
    # find month and day of month
    if (year % 4) == 0:
        # int array containing the number of days per month
        lmonth = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    i = 1
    inttemp = 0
    while dayofyr > inttemp + lmonth[i - 1] and i < 12:
        inttemp = inttemp + lmonth[i - 1]
        i += 1
    mon = i
    day = dayofyr - inttemp
    # find hours minutes and seconds
    temp = (days - dayofyr) * 24.0
    hr = int(temp // 1.0)
    temp = (temp - hr) * 60.0
    minute = int(temp // 1.0)
    sec = (temp - minute) * 60.0
    return mon, day, hr, minute, sec


def invjday(jd):
    """This procedure finds the year, month, day, hour, minute and second given the
    julian date. jd can be ut1, tdt, tdb, etc.

    Args:
        jd : float, julian date, days from 4713 BCE
    Outputs:
        year : int, year between 1900 - 2100
        mon : int, month between 1 - 12
        day : int, day between 1 - 31
        hr : int, hour between 0 - 23
        min : int, minute between 0 - 59
        sec : float, second between 0.0 - 59.999
    References:
        Vallado, 2007, 208, alg 22, ex 3-13
        Rhodes, python-sgp4/sgp4/ext.py
    """
    # find year and days of the year
    jd = np.atleast_1d(jd)
    temp = jd - 2415019.5
    tu = temp / 365.25  # julian centuries from 0 h jan 0, 1900
    year = 1900 + np.floor_divide(tu, 1.0).astype(int)
    # if type(year) == np.int64:
    #     # not hasattr(year, '__getitem__'):
    #     print('hello')
    #     year = np.asarray(year)
    #     print(type(year))
    # leapyrs = int(((year - 1901) * 0.25) // 1.0)  # number of leap years from 1900
    leapyrs = np.floor_divide(((year - 1901) * 0.25), 1.0).astype(int)  # number of leap years from 1900
    # optional nudge by 8.64x10-7 sec to get even outputs
    # day of year plus fractional portion of a day
    days = temp - ((year - 1900) * 365.0 + leapyrs) + 0.00000000001
    # check for case of beginning of a year
    day1_idx = days < 1.0
    year[day1_idx] = year[day1_idx] - 1
    leapyrs[day1_idx] = np.floor_divide((year[day1_idx] - 1901) * 0.25, 1.0).astype(int)
    days[day1_idx] = temp[day1_idx] - ((year[day1_idx] - 1900) * 365.0 + leapyrs[day1_idx])
    # find remaing data
    jd_size = len(jd)
    mon = np.empty(jd_size)
    day = np.empty(jd_size)
    hr = np.empty(jd_size)
    minute = np.empty(jd_size)
    sec = np.empty(jd_size)
    for i in range(jd_size):
        mon[i], day[i], hr[i], minute[i], sec[i] = days2mdhms(year[i], days[i])
    sec = sec - 0.00000086400
    return year, mon, day, hr, minute, sec
