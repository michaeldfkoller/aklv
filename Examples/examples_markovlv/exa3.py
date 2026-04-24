import markovlib as ml

cap = ml.CAPITALLV()
cap.iSetTable("CH-QX-EKM-1995")  # Swiss male capital table 1995
cap.vSetStopTime(35)              # entry age
# start_time is set automatically by iSetTable to omega+1

# Set contract term: maturity at 65
cap.vSetStartTime(65)
cap.vSetDeath(100000.0)
cap.vSetSurvival(65, 100000.0)   # survival benefit at maturity
cap.vSetPremium(2500.0)

# Override discount at 2.5%
for t in range(2500):
    cap.dSetDisc(t, 1.0 / 1.025)

V = cap.dGetDK(35)
print(f"Net reserve at age 35 (entry): {V:,.2f}")
print(f"Net reserve at age 50:         {cap.dGetDK(50):,.2f}")

