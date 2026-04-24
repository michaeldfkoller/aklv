import markovlib as ml

wid = ml.WIDDOWLV()
wid.vSetStopTime(45)    # insured's entry age
wid.vSetStartTime(100)  # projection horizon

# Set constant mortality and improvement
for age in range(130):
    wid.dSetQx(age, 0.005 + age * 0.0003)   # insured (male)
    wid.dSetQy(age, 0.003 + age * 0.0002)   # survivor (female)
    wid.dSetFx(age, -0.015)
    wid.dSetFy(age, -0.020)
    wid.dSetHx(age, 0.70)                   # 70% chance of widow
    wid.dSetYx(age, float(age - 3))         # widow 3 yrs younger

# Benefit: unit widow's pension (constant)
# wid.vLeistReset()  # sets all benefits to 1.0

for t in range(2500):
    wid.dSetDisc(t, 1.0 / 1.03)

print(f"Reserve at entry age 45: {wid.dGetDK(45):.4f}")
print(f"Reserve at age 60:       {wid.dGetDK(60):.4f}")
