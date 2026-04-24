import markovlv as ml

a2 = ml.ANNUITYLV2()
a2.iSetTable1("CH-QX-ERM-1995")   # male
a2.iSetTable2("CH-QX-ERF-1995")   # female
a2.vSetStopTime(65)                # primary life entry age
a2.vSetStartTime(121)
a2.vSetSAge1(65)                   # payments start at 65 for x
a2.vSetSAge2(65)                   # same for y (relative to y's age)
a2.dSetY_Minus_X(62, 65)           # y is 3 years younger (y-x = -3)
a2.dSetBenefit(0, 1.0)             # both alive: full benefit
a2.dSetBenefit(1, 0.6)             # y alone: 60%
a2.dSetBenefit(2, 0.6)             # x alone: 60%

for t in range(2500):
    a2.dSetDisc(t, 1.0 / 1.03)

print(f"DK (both alive, t=65): {a2.dGetDK(65, 0):.4f}")
print(f"DK (y alone,    t=70): {a2.dGetDK(70, 1):.4f}")
print(f"DK (x alone,    t=70): {a2.dGetDK(70, 2):.4f}")
