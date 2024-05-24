"""
stores all the needed variable in a hopefully organized way
"""

from symbolica import Expression, Transformer, set_license_key

x = Expression.symbol("x")
e = (1 + x) ** 2
r = e.expand()
print(r)
ii = Expression.I
pi = Expression.PI
set_license_key("2a14d489-bbd8-5948-8de8-c6fcb28d2a49")
# generic variables
(a, b, c, d, e, g, h, i, j, k, l, m, mm, n, o, p, q, r, s, t) = Expression.symbols(
    "a",
    "b",
    "c",
    "d",
    "e",
    "g",
    "h",
    "i",
    "j",
    "k",
    "l",
    "m",
    "mm",
    "n",
    "o",
    "p",
    "q",
    "r",
    "s",
    "t",
)

X_, Y_, Z_, P_, Q_ = Expression.symbols("X_", "Y_", "Z_", "P_", "Q_")
# f = Expression.symbol("f")
alpha, xi = Expression.symbols("alpha", "xi")
(u, w, v, x, y, z) = Expression.symbols("u", "w", "v", "x", "y", "z")
(delt, k1, D, eps, rep, g0) = Expression.symbols("delt", "k1", "D", "eps", "rep", "g0")

(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, q1, q2, mom) = Expression.symbols(
    "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "q1", "q2", "mom"
)  # momentums
(pd, ph, pb) = Expression.symbols("pd", "ph", "pb")  # momentums
# matching variables
l___, r___, l1___, r1___, l2___, r2___, l3___, r3___, l4___, r4___ = Expression.symbols(
    "l___",
    "r___",
    "l1___",
    "r1___",
    "l2___",
    "r2___",
    "l3___",
    "r3___",
    "l4___",
    "r4___",
)
l__ = Expression.symbol("l__")
(x_, y_, z_, p_, k_, q_, p1_, p2_, p3_, a_, s_) = Expression.symbols(
    "x_", "y_", "z_", "p_", "k_", "q_", "p1_", "p2_", "p3_", "a_", "s_"
)
(
    x1_,
    x2_,
    x3_,
    x4_,
    x5_,
    x6_,
    x7_,
    x8_,
    x9_,
    x10_,
    x11_,
    x12_,
    x13_,
    x14_,
    x15_,
    x16_,
    x17_,
) = Expression.symbols(
    "x1_",
    "x2_",
    "x3_",
    "x4_",
    "x5_",
    "x6_",
    "x7_",
    "x8_",
    "x9_",
    "x10_",
    "x11_",
    "x12_",
    "x13_",
    "x14_",
    "x15_",
    "x16_",
    "x17_",
)
# SUN
(ca, cf, sunn, color, dirac) = Expression.symbols(
    "ca", "cf", "sunn", "color", "dirac"
)  # variables
sd = Expression.symbol("sd", is_symmetric=True)
sdf = Expression.symbol("sdf", is_symmetric=True)
(sunf, sund, sunt, suntf, suntr, suntrrep, adj) = Expression.symbols(
    "sunf", "sund", "sunt", "suntf", "suntr", "suntrrep", "adj"
)  # functions

(mul_, mur_, mul1_, mur1_, mul2_, mur2_, mul3_, mur3_, mul4_, mur4_) = (
    Expression.symbols(
        "mul_",
        "mur_",
        "mul1_",
        "mur1_",
        "mul2_",
        "mur2_",
        "mul3_",
        "mur3_",
        "mul4_",
        "mur4_",
    )
)

(mu, muh, mub) = Expression.symbols("mu", "muh", "mub")  # lorentz indices

(d1, d2, d3, d4, d5, d1ext, d2ext, d3ext, d4ext) = Expression.symbols(
    "d1", "d2", "d3", "d4", "d5", "d1ext", "d2ext", "d3ext", "d4ext"
)  # dirac/spinor indices
(f1, f2, f3, f4, f5, f1ext, f2ext, f3ext, f4ext) = Expression.symbols(
    "f1", "f2", "f3", "f4", "f5", "f1ext", "f2ext", "f3ext", "f4ext"
)  # dirac/spinor indices

# Lorentz and dirac structure
(mom) = Expression.symbol("mom")  # variables
sd = Expression.symbol("sd", is_symmetric=True)
mt = Expression.symbol("mt", is_symmetric=True)
sp = Expression.symbol("sp", is_symmetric=True)
(fv, fvd, fve) = Expression.symbols("fv", "fvd", "fve")  #  lorentz functions

(spd, spe, spp, sppd, sppe) = Expression.symbols(
    "spd", "spe", "spp", "sppd", "sppe"
)  # lorentz functions

(gs, gsd, gse, gsr, gsdr, gser, g5, g55, didelta, dchn) = Expression.symbols(
    "gs", "gsd", "gse", "gsr", "gsdr", "gser", "g5", "g55", "didelta", "dchn"
)  #  dirac functions

(fvr, fvdr, fver, spr, sprr, spdr, spdrr) = Expression.symbols(
    "fvr", "fvdr", "fver", "spr", "sprr", "spdr", "spdrr"
)  # lorentz functions

(sper, sperr, sca, tens, ga, gaa, gae, gad) = Expression.symbols(
    "sper", "sperr", "sca", "tens", "ga", "gaa", "gae", "gad"
)  # lorentz functions
spn, fun = Expression.symbols("spn", "fun")
# weinber operator functions
# Flow time
(tinsert, mm, ss, nn, cbb, caa, cc, m_, ca_, c_, cb_, ta, tb) = Expression.symbols(
    "tinsert",
    "mm",
    "ss",
    "nn",
    "cbb",
    "caa",
    "cc",
    "m_",
    "ca_",
    "c_",
    "cb_",
    "ta",
    "tb",
)

# Generic functions
(exp, lc, masterint, gamma, ldot, mdot, log) = Expression.symbols(
    "exp", "lc", "masterint", "gamma", "ldot", "mdot", "log"
)
dtr, dtrr, dtr5 = Expression.symbols("dtr", "dtrr", "dtr5")
