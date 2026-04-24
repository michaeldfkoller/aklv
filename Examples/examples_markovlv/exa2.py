import markovlv as ml

ann = ml.ANNUITYLV()
ann.iSetTable("CH-QX-ERF-1995")  # Swiss female annuity table 1995
ann.vSetStopTime(65)              # valuation age
ann.vSetSAge(65)                  # payments start at 65
ann.vSetG(5)                     # 10-year guarantee period

# Fill discount at 3% across all positions (table sets its tech rate,
# but we can override)
for t in range(2500):
    ann.dSetDisc(t, 1.0 / 1.03)

dk = ann.dGetDK(65)
print(f"EPV of unit annuity at age 65: {dk:.4f}")

# Expected cash flows for first 5 years
for t in range(65, 90):
    print(f"  CF({t}): {ann.dGetCF(t):.6f}")
