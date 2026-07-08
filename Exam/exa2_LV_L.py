import markovlv as ml

license_str = open("./eth_summer_2027.lic").read().strip()
if ml.iSetLicense(license_str) != 1:
    raise RuntimeError("markovlv:␣invalid␣or␣expired␣license")

SEP = "-" * 56

# -- helper --------------------------------------------------------------------
def term_insurance(q, i, n, S, P):
    """Return net reserve at t=0 (alive) for a simple term policy."""
    disc = 1.0 / (1.0 + i)
    m = ml.MARKOVLV(n + 2, 2, 1)
    m.vSetNrStates(2)
    m.vSetStopTime(0)
    m.vSetStartTime(n)
    m.lSetFolgezustand(0, 0)
    m.lSetFolgezustand(1, 1)
    for t in range(n):
        m.dSetPij(t, 0, 0, 1.0 - q)
        m.dSetPij(t, 0, 1, q)
        m.dSetPij(t, 1, 1, 1.0)
        m.dSetDisc(t, 0, 0, disc)
        m.dSetDisc(t, 1, 1, disc)
        m.dSetPre (t, 0, 0, -P)       # premium paid at start of year (negative = inflow)
        m.dSetPost(t, 0, 1,  S)       # death benefit at end of year
    return m.dGetDK(0, 0, 1)

# -- Example 1: baseline -------------------------------------------------------
print("MARKOVLV  --  5-year term insurance, baseline")
print(SEP)
q0, i0, n0, S0, P0 = 0.01, 0.03, 5, 100_000, 500
V = term_insurance(q0, i0, n0, S0, P0)
print(f"  q={q0:.3f}  i={i0:.3f}  n={n0}  S={S0:,.0f}  P={P0}")
print(f"  Net reserve at issue  V_0 = {V:,.2f}")
print()

# -- Example 2: vary mortality q -----------------------------------------------
print("Sensitivity  --  vary annual mortality q  (i=3%, n=5, S=100k, P=500)")
print(SEP)
print(f"  {'q':>6}  {'V_0':>12}")
for q in [0.002, 0.005, 0.010, 0.020, 0.040]:
    print(f"  {q:6.3f}  {term_insurance(q, 0.03, 5, 100_000, 500):12.2f}")
print()

# -- Example 3: vary interest i ------------------------------------------------
print("Sensitivity  --  vary interest rate i  (q=1%, n=5, S=100k, P=500)")
print(SEP)
print(f"  {'i':>6}  {'V_0':>12}")
for i in [0.01, 0.02, 0.03, 0.04, 0.05]:
    print(f"  {i:6.3f}  {term_insurance(0.01, i, 5, 100_000, 500):12.2f}")
print()

# -- Example 4: net premium by term n -----------------------------------------
# Net premium P* satisfies EPV(benefits) = P* * EPV(unit-premium annuity)
# EPV of unit-premium annuity = term_insurance(q,i,n, S=0, P=1)
#   -> dSetPre = -1 per year -> V_0 = -EPV(annuity-due) -> multiply by -1
print("Net premium by term  (q=1%, i=3%, S=100k)")
print(SEP)
print(f"  {'n':>4}  {'EPV(ben)':>12}  {'EPV(ann-due)':>13}  {'P_net':>8}  {'V_0(check)':>12}")
for n in [5, 10, 15, 20, 30]:
    q = 0.01
    epv_ben  = term_insurance(q, 0.03, n, 100_000, 0)  # P=0: pure benefit EPV
    epv_ann  = -term_insurance(q, 0.03, n, 0, 1)       # P=1: premium income -> negative V -> negate
    P_net    = epv_ben / epv_ann if epv_ann > 0 else 0
    V_check  = term_insurance(q, 0.03, n, 100_000, P_net)  # should ~= 0
    print(f"  {n:4d}  {epv_ben:12.2f}  {epv_ann:13.5f}  {P_net:8.2f}  {V_check:12.6f}")
print()

# -- Example 5: reserve progression over time ---------------------------------
print("Reserve progression  (q=1%, i=3%, n=10, S=100k, P=600)")
print(SEP)
disc = 1.0 / 1.03
q, n, S, P = 0.01, 10, 100_000, 600
m = ml.MARKOVLV(n + 2, 2, 1)
m.vSetNrStates(2)
m.vSetStopTime(0)
m.vSetStartTime(n)
m.lSetFolgezustand(0, 0); m.lSetFolgezustand(1, 1)
for t in range(n):
    m.dSetPij(t, 0, 0, 1.0 - q); m.dSetPij(t, 0, 1, q); m.dSetPij(t, 1, 1, 1.0)
    m.dSetDisc(t, 0, 0, disc);    m.dSetDisc(t, 1, 1, disc)
    m.dSetPre (t, 0, 0, -P);      m.dSetPost(t, 0, 1,  S)
print(f"  {'t':>3}  {'V_t(alive)':>14}  {'RP_t':>10}  {'SP_t':>10}  {'CF_t':>10}")
for t in range(n + 1):
    V  = m.dGetDK(t, 0, 1)
    rp = m.dGetRP(t, 0)
    sp = m.dGetSP(t, 0)
    cf = m.dGetCF(t,0,0) + m.dGetCF(t,0,1)
    print(f"  {t:3d}  {V:14.2f}  {rp:10.2f}  {sp:10.2f} {cf:10.2f}")
