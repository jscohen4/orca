cfs_tafd = 2.29568411*10**-5 * 86400 / 1000
tafd_cfs = 1000 / 86400 * 43560
cfsd_mafd = 2.29568411*10**-5 * 86400 / 10 ** 6
first_of_month = [1,31,59,90,119,150,180,211,242,272,303,334]

def water_day(d):
  return d - 274 if d >= 274 else d + 91

def water_month(m):
  return m - 9 if m >= 9 else m + 3
for m in range(1,12):
	print water_month(m)
