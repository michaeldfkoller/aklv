import markovlv as ml

q    = 0.01    # annual mortality
i    = 0.03    # technical interest
disc = 1.0 / (1.0 + i)
n    = 5       # term in years

m = ml.MARKOVLV(10,2,1)
m.vSetNrStates(2)
m.vSetStopTime(0)       # entry age (relative time)
m.vSetStartTime(n)      # contract end

for t in range(n):
    m.dSetPij(t, 0, 0, 1.0 - q)  # alive -> alive
    m.dSetPij(t, 0, 1, q)         # alive -> dead
    m.dSetPij(t, 1, 1, 1.0)       # dead absorbing
    m.dSetDisc(t, 0, 0, disc)
    m.dSetDisc(t, 1, 1, disc)
    m.dSetPre(t,  0, 0, -500.0)   # premium (negative = outflow)
    m.dSetPost(t, 0, 1, 100_000)  # death benefit

V0 = m.dGetDK(0, 0, 1)           # reserve at t=0, alive state
print(f"Reserve at issue: {V0:,.2f}")
# Expected: EPV(benefits) - EPV(premiums)
