"""
driver.py
=========
Comprehensive driver for the PrivateCredit Markov-chain valuation model.

Sections
--------
  1. Replication of private_credit_model.py scenarios (Monte-Carlo bootstrap)
  2. Expected cash-flow schedules (dGetCF) for all four toggle combinations
  3. Fair-value comparison table (exact values)
  4. Notch-downgrade sensitivity table (0, 1, 2 notches)
  5. Cash-flow profile figures (matplotlib, one figure per PIK×extension)
  6. LaTeX table data export for documentation

Instrument (base case)
-----------------------
  Nominal    : 100,000
  Coupon     : see per-rating table
  Frequency  : annual (freq=1)
  Maturity   : 5 years
  LGD        : 40%
  RF         : 4.5%
  Extension  : 40% at T, 40% at T+1, 20% at T+2

Run
---
  cd MarkovLib/PrivateCredit
  python driver.py
"""

import sys, os, math, random, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from PrivateCredit import (
    PrivateCredit, RATING_LABELS, M_STATES, D_IDX,
    SP_BOOTSTRAP, BOOTSTRAP_YEARS, _normalise,
)

# ============================================================
# Instrument parameters
# ============================================================

N_YEARS   = 5
LGD       = 0.40
RF        = 0.045
FREQ      = 1           # annual
N_SIM     = 20_000
SEED      = 42

# Rating → (init_state_index, annual_coupon_rate, label)
RATINGS = [
    ("BBB", 3, 0.080, "BBB"),
    ("BB",  4, 0.090, "BB"),
    ("B",   5, 0.105, "B"),
    ("CCC", 6, 0.130, "CCC"),
]

# 3-tranche extension schedule: 40% at T, 40% at T+1, 20% at T+2
EXT_SCHEDULE = (0.40, 0.40, 0.20)

# Four PIK × extension toggle combinations
TOGGLES = [
    ("Cash-pay / No ext",  dict(pik=False, extension=False)),
    ("Cash-pay / Ext",     dict(pik=False, extension=True,  ext_schedule=EXT_SCHEDULE)),
    ("PIK / No ext",       dict(pik=True,  extension=False)),
    ("PIK / Ext",          dict(pik=True,  extension=True,  ext_schedule=EXT_SCHEDULE)),
]

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


# ============================================================
# Helpers
# ============================================================

def make_instrument(rating_key, toggle_kw, notch=0, nominal=100_000.0):
    """Build a PrivateCredit object for a given rating, toggle set and notch."""
    _, i0, coup, _ = next(r for r in RATINGS if r[0] == rating_key)
    eff_state = min(i0 + notch, M_STATES - 1)
    return PrivateCredit(
        nominal     = nominal,
        coupon_rate = coup,
        n_years     = N_YEARS,
        freq        = FREQ,
        init_rating = eff_state,
        lgd         = LGD,
        rf          = RF,
        **toggle_kw,
    )


# ============================================================
# Section 1: Bootstrap replication of private_credit_model.py
# ============================================================

def _normalise_8(mat):
    """Normalise 8-state (no C-state) bootstrap matrix."""
    out = []
    for row in mat:
        s = sum(row)
        out.append([v / s for v in row])
    return out


def bootstrap_loss_distribution(
    init_state,
    nominal,
    coupon,
    recovery,
    n_years,
    rf,
    n_sim      = N_SIM,
    seed       = SEED,
    pik        = False,
    extension  = False,
    ext_sched  = EXT_SCHEDULE,
    notch      = 0,
):
    """
    Bootstrap loss distribution mimicking private_credit_model.py.

    Parameters
    ----------
    init_state : 0-based index into 8-state bootstrap matrix (AAA..CCC)
    nominal    : face value
    coupon     : annual coupon rate
    recovery   : 1 - LGD
    n_years    : maturity in years (annual periods)
    rf         : annual risk-free rate
    n_sim      : number of Monte Carlo paths
    seed       : RNG seed
    pik        : PIK mode
    extension  : maturity extension active
    ext_sched  : extension fractions tuple (f0, f1, ...) summing to 1
    notch      : initial-state downgrade notches applied

    Returns
    -------
    sorted list of losses (Nominal - PV received)
    """
    rng  = random.Random(seed)
    disc = 1.0 / (1.0 + rf)
    N    = nominal
    T    = n_years
    n_ext = len(ext_sched) - 1 if extension else 0
    total_T = T + n_ext

    # Map 7-state (AAA..CCC) to 8-state with D at index 7
    D_bs = 7  # default index in 8-state bootstrap matrix

    losses = []
    for _ in range(n_sim):
        # Sample annual matrices for each period
        years_drawn = [rng.choice(BOOTSTRAP_YEARS) for _ in range(total_T)]

        cur_state   = min(init_state + notch, D_bs - 1)
        pv_received = 0.0
        defaulted   = False
        pik_accum   = N  # outstanding principal under PIK

        for t in range(total_T):
            if defaulted:
                break
            yr  = years_drawn[t]
            mat = _normalise_8(SP_BOOTSTRAP[yr])
            row = mat[cur_state]

            u   = rng.random()
            cum = 0.0
            next_state = cur_state
            for j, p in enumerate(row):
                cum += p
                if u <= cum:
                    next_state = j
                    break

            df = disc ** (t + 1)

            if t < T:
                if cur_state < D_bs:
                    if not pik:
                        pv_received += N * coupon * df  # cash coupon
                    else:
                        pik_accum = N * (1 + coupon) ** (t + 1)  # track accrual
                if next_state == D_bs and not defaulted:
                    outstanding = pik_accum if pik else N
                    pv_received += outstanding * recovery * df
                    defaulted = True

            else:
                # Extension period t = T, T+1, ...
                k    = t - T   # 0-based extension index
                frac = ext_sched[k + 1] if extension else 0.0

                if not defaulted and cur_state < D_bs:
                    mat_amt  = pik_accum if pik else N
                    if next_state < D_bs:
                        pv_received += mat_amt * frac * df
                    else:
                        pv_received += mat_amt * frac * recovery * df
                        defaulted = True

            # Maturity payment at T (first tranche)
            if t == T - 1 and not defaulted and cur_state < D_bs:
                mat_amt = pik_accum if pik else N
                if next_state < D_bs:
                    pv_received += mat_amt * ext_sched[0] * df
                elif not defaulted:
                    # Default at maturity
                    pv_received += mat_amt * ext_sched[0] * recovery * df
                    defaulted = True

            cur_state = next_state

        losses.append(N - pv_received)

    return sorted(losses)


def var_at(losses, alpha):
    idx = max(0, min(int(math.ceil(alpha * len(losses))) - 1, len(losses) - 1))
    return losses[idx]


def cvar_at(losses, alpha):
    idx  = int(math.ceil(alpha * len(losses)))
    tail = losses[idx:]
    return sum(tail) / len(tail) if tail else losses[-1]


def run_bootstrap():
    print("\n" + "="*72)
    print("  Section 1 — Bootstrap loss distribution (replica of private_credit_model.py)")
    print("="*72)

    rows = []
    for rtag, i0, coup, rlabel in RATINGS:
        for tlab, tkw in TOGGLES[:2]:   # cash-pay only to keep output concise
            pik  = tkw["pik"]
            ext  = tkw.get("extension", False)
            esched = tkw.get("ext_schedule", (1.0,)) if ext else (1.0,)
            losses = bootstrap_loss_distribution(
                init_state = i0,
                nominal    = 100_000.0,
                coupon     = coup,
                recovery   = 1.0 - LGD,
                n_years    = N_YEARS,
                rf         = RF,
                pik        = pik,
                extension  = ext,
                ext_sched  = esched if ext else (1.0,),
            )
            el   = sum(losses) / len(losses)
            v95  = var_at(losses, 0.95)
            v99  = var_at(losses, 0.99)
            cv99 = cvar_at(losses, 0.99)
            row  = dict(rating=rlabel, scenario=tlab, EL=el, VaR95=v95,
                        VaR99=v99, CVaR99=cv99)
            rows.append(row)
            print(f"  {rlabel:>3} | {tlab:<22} | EL={el:8.1f} | "
                  f"VaR95={v95:8.1f} | VaR99={v99:8.1f} | CVaR99={cv99:8.1f}")

    csv_path = os.path.join(OUT_DIR, "bootstrap_summary.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved to: {csv_path}")
    return rows


# ============================================================
# Section 2: Expected cash flows per period
# ============================================================

def run_cashflow_tables():
    print("\n" + "="*72)
    print("  Section 2 — Expected cash flows per period (dGetCF)")
    print("="*72)
    for rtag, i0, coup, rlabel in RATINGS:
        for tlab, tkw in TOGGLES:
            pc = make_instrument(rtag, tkw)
            pc.print_cash_flows(label=f"{rlabel} | {tlab}")


# ============================================================
# Section 3: Fair-value comparison table
# ============================================================

def compute_fair_value_table():
    """
    Returns dicts keyed [rating][toggle_label]: fair_value, maturity_amt,
    macaulay_duration, modified_duration.
    """
    fv_table  = {}
    mat_table = {}
    mac_table = {}
    mod_table = {}
    for rtag, i0, coup, rlabel in RATINGS:
        fv_table[rlabel]  = {}
        mat_table[rlabel] = {}
        mac_table[rlabel] = {}
        mod_table[rlabel] = {}
        for tlab, tkw in TOGGLES:
            pc = make_instrument(rtag, tkw)
            fv_table[rlabel][tlab]  = pc.dGetDK(0)
            mat_table[rlabel][tlab] = pc._maturity_amount()
            mac_table[rlabel][tlab] = pc.macaulay_duration()
            mod_table[rlabel][tlab] = pc.modified_duration()
    return fv_table, mat_table, mac_table, mod_table


def run_fair_value_table(fv_table, mat_table, mac_table, mod_table):
    print("\n" + "="*88)
    print("  Section 3 — Fair-value and duration comparison (exact values)")
    print("="*88)
    hdr = (f"  {'Rating':<5}  {'Scenario':<24}  {'Mat.Amt':>12}  "
           f"{'Fair Value':>12}  {'FV/Nom':>9}  {'D_mac':>7}  {'D_mod':>7}")
    print(hdr)
    print("  " + "-"*84)
    for rtag, i0, coup, rlabel in RATINGS:
        for tlab, tkw in TOGGLES:
            fv  = fv_table[rlabel][tlab]
            mat = mat_table[rlabel][tlab]
            mac = mac_table[rlabel][tlab]
            mod = mod_table[rlabel][tlab]
            print(f"  {rlabel:<5}  {tlab:<24}  {mat:>12,.2f}  "
                  f"{fv:>12,.2f}  {fv/100_000:>9.4%}  {mac:>7.4f}  {mod:>7.4f}")
    print()


# ============================================================
# Section 4: Notch-downgrade sensitivity
# ============================================================

def compute_notch_table():
    """
    Returns notch_table[rating][scenario][notch] = {fv, mac, mod}.
    """
    notch_table = {}
    for rtag, i0, coup, rlabel in RATINGS:
        notch_table[rlabel] = {}
        for tlab, tkw in TOGGLES:
            notch_table[rlabel][tlab] = {}
            for notch in range(3):
                pc = make_instrument(rtag, tkw, notch=notch)
                notch_table[rlabel][tlab][notch] = {
                    "fv":  pc.dGetDK(0),
                    "mac": pc.macaulay_duration(),
                    "mod": pc.modified_duration(),
                }
    return notch_table


def run_notch_table(notch_table):
    print("\n" + "="*100)
    print("  Section 4 — Notch-downgrade sensitivity (0 / 1 / 2 notches)")
    print("  Columns: Fair Value | Macaulay Dur | Modified Dur   (per notch level)")
    print("  (Downgrade applied to initial rating; coupon unchanged)")
    print("="*100)
    hdr = (f"  {'Rating':<5}  {'Scenario':<24}  "
           f"{'FV(0)':>10}  {'Mac':>6}  {'Mod':>6}  "
           f"{'FV(1)':>10}  {'Mac':>6}  {'Mod':>6}  "
           f"{'FV(2)':>10}  {'Mac':>6}  {'Mod':>6}")
    print(hdr)
    print("  " + "-"*96)
    for rtag, i0, coup, rlabel in RATINGS:
        for tlab, tkw in TOGGLES:
            d = notch_table[rlabel][tlab]
            row = f"  {rlabel:<5}  {tlab:<24}"
            for n in range(3):
                row += (f"  {d[n]['fv']:>10,.2f}  {d[n]['mac']:>6.3f}  {d[n]['mod']:>6.3f}")
            print(row)
    print()


# ============================================================
# Section 5: Cash-flow profile figures
# ============================================================

def plot_cf_profiles(fig_prefix="cf_profile"):
    """
    4 figures (one per PIK × extension combo), each showing expected CF
    per period for BBB / BB / B / CCC as grouped bars.
    Returns list of saved file paths.
    """
    colors = {
        "BBB": "#2166ac",
        "BB":  "#4dac26",
        "B":   "#d01c8b",
        "CCC": "#f1a340",
    }
    saved = []

    for fig_idx, (tlab, tkw) in enumerate(TOGGLES):
        # Build per-rating instruments and collect CFs + durations
        instruments = {}
        all_cfs     = {}
        for rtag, i0, coup, rlabel in RATINGS:
            pc = make_instrument(rtag, tkw)
            instruments[rlabel] = pc
            all_cfs[rlabel]     = pc.expected_cash_flows()

        # Common period range
        n_periods = max(len(v) for v in all_cfs.values())
        years     = [(t + 1) / FREQ for t in range(n_periods)]

        fig, (ax_cf, ax_dur) = plt.subplots(
            1, 2, figsize=(14, 5),
            gridspec_kw={"width_ratios": [3, 1]}
        )

        # ---- left: CF bars ----
        n_r     = len(RATINGS)
        width   = 0.18
        offsets = np.linspace(-(n_r - 1) * width / 2, (n_r - 1) * width / 2, n_r)

        for oi, (rtag, i0, coup, rlabel) in enumerate(RATINGS):
            cf_vals = [0.0] * n_periods
            for t, yr, cf in all_cfs[rlabel]:
                if t < n_periods:
                    cf_vals[t] = cf
            x = np.array(range(n_periods), dtype=float) + offsets[oi]
            ax_cf.bar(x, cf_vals, width=width, label=rlabel,
                      color=colors[rlabel], alpha=0.82, edgecolor="white")

        ax_cf.set_xticks(range(n_periods))
        ax_cf.set_xticklabels([f"t={t}\n({yr:.1f}y)" for t, yr in enumerate(years)],
                              fontsize=7.5)
        ax_cf.set_xlabel("Period (t) / Year")
        ax_cf.set_ylabel("Expected CF (EUR)")
        ax_cf.set_title(f"Expected Cash-Flow Profile — {tlab}\n"
                        f"(N=100,000  5y  LGD=40%  rf=4.5%)", fontsize=10)
        ax_cf.legend(title="Init. Rating", loc="upper left", fontsize=8)
        ax_cf.yaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f"{x/1000:.0f}K" if abs(x) >= 1000 else f"{x:.0f}"
        ))
        ax_cf.grid(axis="y", linestyle="--", alpha=0.4)

        # ---- right: duration bar chart ----
        rlabels_list = [r[3] for r in RATINGS]
        mac_vals = [instruments[rl].macaulay_duration() for rl in rlabels_list]
        mod_vals = [instruments[rl].modified_duration()  for rl in rlabels_list]

        x_d = np.arange(len(rlabels_list))
        w_d = 0.35
        ax_dur.bar(x_d - w_d/2, mac_vals, w_d, label="Macaulay",
                   color=[colors[r] for r in rlabels_list], alpha=0.85, edgecolor="white")
        ax_dur.bar(x_d + w_d/2, mod_vals, w_d, label="Modified",
                   color=[colors[r] for r in rlabels_list], alpha=0.45,
                   edgecolor=[colors[r] for r in rlabels_list], linewidth=1.2)
        ax_dur.set_xticks(x_d)
        ax_dur.set_xticklabels(rlabels_list, fontsize=9)
        ax_dur.set_ylabel("Duration (years)")
        ax_dur.set_title("Macaulay / Modified\nDuration", fontsize=10)
        ax_dur.grid(axis="y", linestyle="--", alpha=0.4)
        # Annotate bars with values
        for xi, (mac, mod) in enumerate(zip(mac_vals, mod_vals)):
            ax_dur.text(xi - w_d/2, mac + 0.02, f"{mac:.2f}", ha="center",
                        va="bottom", fontsize=7, color="black")
            ax_dur.text(xi + w_d/2, mod + 0.02, f"{mod:.2f}", ha="center",
                        va="bottom", fontsize=7, color="dimgray")

        # Custom legend for duration plot
        from matplotlib.patches import Patch
        ax_dur.legend(handles=[
            Patch(facecolor="gray", alpha=0.85, label="Macaulay"),
            Patch(facecolor="gray", alpha=0.45, label="Modified"),
        ], fontsize=8, loc="upper right")

        plt.tight_layout()
        fname = os.path.join(
            OUT_DIR,
            f"{fig_prefix}_{fig_idx+1}_{tlab.replace('/', '_').replace(' ','')}.pdf"
        )
        plt.savefig(fname, bbox_inches="tight")
        plt.close(fig)
        saved.append(fname)
        print(f"  Saved figure: {fname}")

    return saved


# ============================================================
# Section 6: Export LaTeX table data to JSON for documentation
# ============================================================

def export_table_data(fv_table, mat_table, mac_table, mod_table, notch_table):
    """Write all computed values to a JSON file for use in the LaTeX source."""
    data = {
        "fair_values":        fv_table,
        "maturity_amt":       {r: {t: mat_table[r][t] for t in mat_table[r]}
                               for r in mat_table},
        "macaulay_duration":  {r: {t: mac_table[r][t] for t in mac_table[r]}
                               for r in mac_table},
        "modified_duration":  {r: {t: mod_table[r][t] for t in mod_table[r]}
                               for r in mod_table},
        "notch":              notch_table,
        "parameters":   {
            "n_years": N_YEARS, "freq": FREQ, "lgd": LGD, "rf": RF,
            "ext_schedule": list(EXT_SCHEDULE),
            "ratings": {r[0]: {"state": r[1], "coupon": r[2]} for r in RATINGS},
        },
    }
    path = os.path.join(OUT_DIR, "table_data.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  Saved table data: {path}")
    return data


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("\nPrivateCredit Driver — MarkovLib  (refined, June 2026)")
    print("="*72)
    print(f"  N=100,000 | 5y annual | LGD={LGD:.0%} | rf={RF:.1%}")
    print(f"  Extension schedule: {EXT_SCHEDULE}  (40% at T, 40% at T+1, 20% at T+2)")
    print(f"  Ratings: BBB ({RATINGS[0][2]:.1%}), BB ({RATINGS[1][2]:.1%}), "
          f"B ({RATINGS[2][2]:.1%}), CCC ({RATINGS[3][2]:.1%})")
    print("="*72)

    # Section 1: Bootstrap (mimics private_credit_model.py)
    boot_rows = run_bootstrap()

    # Section 2: CF schedules
    run_cashflow_tables()

    # Section 3: Fair values + durations
    fv_table, mat_table, mac_table, mod_table = compute_fair_value_table()
    run_fair_value_table(fv_table, mat_table, mac_table, mod_table)

    # Section 4: Notch downgrade
    notch_table = compute_notch_table()
    run_notch_table(notch_table)

    # Section 5: Figures
    print("\n" + "="*72)
    print("  Section 5 — CF profile figures")
    print("="*72)
    saved_figs = plot_cf_profiles()

    # Section 6: Export data
    print("\n" + "="*72)
    print("  Section 6 — Export table data for LaTeX")
    print("="*72)
    tdata = export_table_data(fv_table, mat_table, mac_table, mod_table, notch_table)

    print("\nDone.  All outputs saved to:", OUT_DIR)
