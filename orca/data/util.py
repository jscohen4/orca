cfs_tafd = 2.29568411*10**-5 * 86400 / 1000
tafd_cfs = 1000 / 86400 * 43560
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
first_of_month = [1,31,59,90,119,150,180,211,242,272,303,334]
z_table_transform = [-1.645, -1.28, -1.035, -0.84, -0.675, -0.525, -0.385, -0.253, -0.125, 0, 0.125, 0.253, 0.385, 0.525, 0.675, 0.84, 1.035, 1.28, 1.645]

def water_day(d):
  return d - 274 if d >= 274 else d + 91
  
def water_year_day(d):  #obtain day of water year, which begins on October 1st
  if d.is_leap_year:
    if d.dayofyear >= 275:
      return d.dayofyear - 274
    elif d.dayofyear <= 274 and d.dayofyear >= 59:  
      return d.dayofyear + 92
    else:
      return d.dayofyear + 92
  elif not d.is_leap_year:
    if d.dayofyear >= 274:
      return d.dayofyear - 273
    else:
      return d.dayofyear + 92

