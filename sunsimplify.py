import sys
import time

import numpy as np
from symbolica import Expression, Transformer, get_version, set_license_key

set_license_key("2a14d489-bbd8-5948-8de8-c6fcb28d2a49")
sys.path.append("/home/oscar/Desktop/PhD/BFM-GF/LEFT/Symbolica/")
from myfunctions import freeq
from variables import *

# x = Expression.symbol("x")
# e = (1 + x) ** 2
# r = e.expand()
# print(r)

T = Transformer()


def suntwo() -> Transformer:
    """
    simplifies sunf and sund  structure constants when two indices are contracted
    """
    trans = T.repeat(
        T.repeat(
            T.replace_all(  # order sunf so next formula always triggers, otherwise we could get wrong ordering
                sunf(l___, x_, y_, r___),
                -sunf(l___, y_, x_, r___),
                x_.req_cmp_gt(y_, True),
            )
        ),
        T.repeat(
            T.replace_all(  # order sund so next formula always triggers, otherwise we could get wrong ordering
                sund(l___, x_, y_, r___),
                sund(l___, y_, x_, r___),
                x_.req_cmp_gt(y_, True),
            )
        ),
        T.replace_all(
            sunf(l___, x_, l2___, y_, r___) * sund(l3___, x_, r2___, y_, r3___), 0
        ),
        T.replace_all(  # use formula for CA
            sunf(l___, x_, l2___, y_, r___) * sunf(l1___, x_, r2___, y_, r1___),
            ca
            * sd(l___, r___, l1___, r1___, l2___, r2___)
            * (-1) ** (l2___.transform().nargs())
            * (-1) ** (r2___.transform().nargs()),
        ),
        T.replace_all(  # use formula for sund
            sund(l___, x_, l2___, y_, r___) * sund(l1___, x_, r2___, y_, r1___),
            (ca**2 - 4) / ca * sd(l___, r___, l1___, r1___, l2___, r2___),
        ),
        T.replace_all(  # sunf^2
            sunf(x_, y_, z_) ** 2,
            ca * (ca**2 - 1),
        ).expand(),
        T.replace_all(  # sund^2
            sunf(x_, y_, z_) ** 2,
            ca * (ca**2 - 1),
        ),
        T.repeat(
            T.replace_all(  # sd*sd
                sd(l___, x_, r___) * sd(l2___, x_, r2___),
                sd(l2___, l___, r___, r2___),
            )
        ),
        T.repeat(T.replace_all(sd(x_, x_), ca**2 - 1)),
        T.repeat(T.replace_all(sd(x_, y_) ** 2, ca**2 - 1)),  # sd(a,a)
        T.repeat(
            T.replace_all(  # sd*sunf
                sd(l___, x_, r___) * sunf(l2___, x_, r2___),
                sunf(l2___, l___, r___, r2___),
            )
        ),
        T.repeat(
            T.replace_all(  # sd*sund
                sd(l___, x_, r___) * sund(l2___, x_, r2___),
                sund(l2___, l___, r___, r2___),
            )
        ),
    )

    return trans


def sunone() -> Transformer:
    """
    simplifies sunf and sund structure constants when only one index is contracted
    """
    trans = (
        T.repeat(  # sunf contractions
            T.replace_all(
                sunf(x1_, x2_, x3_) * sunf(x1_, x4_, x5_),
                2
                * (
                    -suntr(x2_, x3_, x4_, x5_)
                    + suntr(x2_, x3_, x5_, x4_)
                    + suntr(x3_, x2_, x4_, x5_)
                    - suntr(x3_, x2_, x5_, x4_)
                ),
            ),
            T.replace_all(
                sunf(l___, x_, r___) * sunf(l2___, x_, r2___),
                sunf(x_, l___, r___)
                * sunf(x_, l2___, r2___)
                * (-1) ** (l___.transform().nargs() + l2___.transform().nargs()),
            ),
        )
        .repeat(  # sund contractions
            T.replace_all(
                sund(x1_, x2_, x3_) * sund(x1_, x4_, x5_),
                -4
                * (
                    sd(x2_, x3_) * sd(x4_, x5_) / (2 * ca)
                    + (
                        suntr(x2_, x3_, x4_, x5_)
                        + suntr(x2_, x3_, x5_, x4_)
                        + suntr(x3_, x2_, x4_, x5_)
                        + suntr(x3_, x2_, x5_, x4_)
                    )
                    / 2
                ),
            ),
            T.replace_all(
                sund(l___, x_, r___) * sund(l2___, x_, r2___),
                sund(x_, l___, r___) * sund(x_, l2___, r2___),
            ),
        )
        .chain(sunttraces())  # we simplify the traces now to kill terms earlier
    )

    return trans


def sunttraces(conv: Expression = 1) -> Transformer:
    """
    simplifies products of SU(N) generators with contracted indices
    """
    trans = T.repeat(
        sunsimple(conv),
        T.replace_all(  # equal trace contraction
            suntr(l___, x_, l2___, x_, r___),
            conv**2 * suntr(l___, l2___, r___) / (2 * ca)
            - conv**2 * suntr(l___, r___) * suntr(l2___) / 2,
        ),
        sunsimple(conv),
        T.replace_all(  # different trace contraction
            suntr(l___, x_, r___) * suntr(l2___, x_, r2___),
            conv**2 * suntr(l___, r___) * suntr(l2___, r2___) / (2 * ca)
            - conv**2 * suntr(l___, r2___, l2___, r___) / 2,
        ),
        sunsimple(conv),
        T.replace_all(  # different trace contraction with a square
            suntr(x_, r___) ** 2,
            conv**2 * suntr(r___) * suntr(r___) / (2 * ca)
            - conv**2 * suntr(r___, r___) / 2,
        ),
        T.replace_all(suntr(), ca),  # trace I = N,
    )
    return trans


def sunsimple(conv: Expression = 1) -> Transformer:
    """
    performs SU(N) simplifications of two or three generators
    """
    trans = T.repeat(T.replace_all(suntr(x_), 0), twotraces(conv), threetraces(conv))
    return trans


def tonumberofcolors() -> Transformer:
    """
    replaces cf->(N^2-1)/(2N) and ca->N
    """
    trans = T.chain(
        T.repeat(T.replace_all(cf, (n**2 - 1) / (2 * n))),
        T.repeat(T.replace_all(ca, n)),
    )
    return trans


def twotraces(conv: Expression = 1) -> Transformer:
    """
    simplifies traces of two SU(N) generators
    """
    trans = T.repeat(
        T.replace_all(suntr(x_, x_), -cf * ca * conv**2),
        T.replace_all(suntr(l___, x_, x_, r___), -cf * conv**2 * suntr(l___, r___)),
        T.replace_all(
            suntr(x_, l___, x_),
            -cf * conv**2 * suntr(l___),  # trace cyclicity of the above formula
        ),  # cyclicity
        T.replace_all(suntr(x_, y_), -(conv**2) * sd(x_, y_) / 2),
        T.replace_all(
            sd(l___, x_, r___) * suntr(l2___, x_, r2___),
            suntr(l2___, l___, r___, r2___),
        ),
        T.replace_all(
            suntr(x_), 0
        ),  # check again in case it triggers due to the previous replacement
        T.replace_all(suntr(), ca),  # trace I = N
        T.replace_all(  # sd*sd
            sd(l___, x_, r___) * sd(l2___, x_, r2___),
            sd(l2___, l___, r___, r2___),
        ),
        T.replace_all(sd(x_, x_), ca**2 - 1),
        T.replace_all(sd(x_, y_) ** 2, ca**2 - 1),
    )

    return trans


def threetraces(conv: Expression = 1) -> Transformer:
    """
    simplifies traces of three SU(N) generators
    """
    trans = T.repeat(
        T.replace_all(
            suntr(x_, y_, l___, y_), suntr(y_, x_, y_, l___)
        ),  # cyclicity at the start
        T.replace_all(
            suntr(x_, l___, x_, y_), suntr(l___, x_, y_, x_)
        ),  # cyclicity at the end
        T.replace_all(
            suntr(l___, x_, y_, x_, r___), conv**2 * suntr(l___, y_, r___) / (2 * ca)
        ),  # replacement
    )

    return trans


def sunfdt(conv: Expression = 1) -> Transformer:
    """
    simplifies products of SU(N) generators with sunf and/or sund
    """
    trans = (
        T.repeat(
            T.replace_all(  # sunf*sund
                sund(l___, x_, r___) * sunf(l2___, x_, r2___),
                (-1) ** (l2___.transform().nargs())
                * 2
                * ii
                * (
                    -suntr(l2___, r2___, l___, r___)
                    + suntr(r2___, l2___, l___, r___)
                    - suntr(l2___, r2___, r___, l___)
                    + suntr(l2___, r2___, l___, r___)
                ),
            )
        )
        .replace_all(
            sunf(l___, x_, r___) * suntr(l2___, x_, r2___),
            (-1) ** (l___.transform().nargs())
            * sunf(x_, l___, r___)
            * suntr(l2___, x_, r2___),
        )
        .replace_all(
            sunf(x_, y_, z_)
            * suntr(l2___, x_, r2___),  # the extra factor of conv is because of Fierz
            conv**3 * suntr(r2___, l2___, y_, z_)
            - conv**3
            * suntr(
                r2___, l2___, z_, y_
            ),  # conv instead of conv**2 is not a typo, there's just a factor of i here
        )
        .expand()
        .replace_all(
            sund(l___, x_, r___) * suntr(l2___, x_, r2___),
            sund(x_, l___, r___) * suntr(l2___, x_, r2___),
        )
        .replace_all(
            sund(x_, y_, z_) * suntr(l2___, x_, r2___),
            (conv**3) * ii * sd(y_, z_) * suntr(l2___, r2___) / (ca)
            + ii * conv**3 * suntr(r2___, l2___, y_, z_)
            + ii * conv**3 * suntr(r2___, l2___, z_, y_),
        )
        .chain(sunttraces())
    )

    return trans


def colorsimplify(conv: Expression = 1) -> Transformer:
    """
    performs the necessary SU(N) simplifications after having projected the color structures
    """
    trans = T.chain(
        adjointcontract(),  # contracts all adjoint indices
        suntwo(),  # simplifies sunf and sund with two or three contracted indices
        sunttraces(conv),  # simplifies products of SU(N) traces with contracted indices
        sunone(),  # simplifies sunf and sund with one contracted index
        sunfdt(
            conv
        ),  # we simplify the traces now to kill terms earlier  # simplifies products of SU(N) generators with sunf and/or sund
        sunsimple(conv),
    )
    return trans


def factorcolor() -> Transformer:
    """
    multiplies and expands factors that contain SUN structures, but leaves the rest unexpanded
    """
    trans = T.chain(
        T.replace_all(x_, color(x_)),
        T.repeat(
            T.replace_all(
                color(x_ * r___),
                x_ * color(r___),
                x_.req(
                    lambda m: freeq(
                        m,
                        [
                            sunf(l___),
                            adj(l___),
                            sd(l___),
                            sunf(l___),
                            sund(l___),
                            sunt(l___),
                            suntr(l___),
                            sdf(l___),
                        ],
                    )
                ),
            )
        ),
        T.replace_all(color(x_), color(x_.transform().expand())),
    )

    return trans


def fcjoin() -> Transformer:
    """
    performs multiplications of SU(N) structures with fundamental SU(N) indices
    """
    trans = T.chain(
        T.repeat(
            T.replace_all(
                sdf(l___, x_, r___) * sdf(l1___, x_, r1___),
                sdf(l___, l1___, r1___, r___),
            )
        ),
        T.replace_all(
            suntf(l__, x_, y_),
            suntf(mdot(l__), x_, y_),
        ),
        T.repeat(
            T.replace_all(
                sdf(l___, x1_, r___) * suntf(x2_, l2___, x1_, r2___),
                suntf(x2_, l2___, l___, r___, r2___),
            )
        ),
        T.repeat(
            T.replace_all(
                suntf(mdot(x1_), x2_, x3_) * suntf(mdot(x4_), x3_, x5_),
                suntf(mdot(x1_, x4_), x2_, x5_),
            )
        ),
        T.repeat(
            T.replace_all(
                suntf(mdot(l___, mdot(l1___), r___), x2_, x3_),
                suntf(mdot(l___, l1___, r___), x2_, x3_),
            )
        ),
    )

    return trans


def tocolortraces() -> Transformer:
    trans = T.chain(
        T.repeat(T.replace_all(suntf(mdot(l___), fun(x_), fun(x_)), suntr(l___)))
    )
    return trans


def adjointcontract() -> Transformer:
    trans = T.chain(
        T.repeat(
            T.replace_all(
                sd(l___, x_, r___) * sd(l1___, x_, r1___), sd(l___, l1___, r1___, r___)
            )
        ),
        T.repeat(
            T.replace_all(
                sd(l___, x_, r___) * suntr(l1___, x_, r1___),
                suntr(l1___, l___, r___, r1___),
            )
        ),
    )
    return trans


def fullcoloralgorithm(conv: Expression = 1) -> Transformer:
    """
    takes care of all the SU(N) simplifications: contraction of fundamental indices, writing in terms of traces, and performs the color algorithm.
    conv=ii => hermitian generators, conv=1=> antihermitian generators (lattice convention)
    """
    trans = T.chain(
        fcjoin(),
        factorcolor(),
        fcjoin(),
        tocolortraces(),
        T.replace_all(color(y_), y_.transform().chain(colorsimplify(conv), T.expand())),
    )
    return trans


def checksunsimplifynonhermitian():
    """
    # checks that we get the right results for some chosen cases
    """
    checks = [
        sunf(adj(1), adj(2), adj(3)) * sunf(adj(1), adj(2), adj(4)),
        sunf(adj(1), adj(2), adj(3))
        * sunf(adj(1), adj(4), adj(5))
        * sunf(adj(2), adj(3), adj(6))
        * sunf(adj(4), adj(5), adj(6)),
        sunf(adj(1), adj(2), adj(3))
        * sunf(adj(1), adj(4), adj(5))
        * sunf(adj(2), adj(4), adj(6))
        * sunf(adj(3), adj(5), adj(6)),
        suntr(adj(1), adj(2), adj(1), adj(3)),
        sunf(adj(1), adj(2), adj(4)) * suntr(adj(1), adj(2), adj(3)),
        sund(adj(1), adj(2), adj(4)) * suntr(adj(1), adj(2), adj(3)),
        sund(adj(1), adj(2), adj(3)) * sund(adj(1), adj(2), adj(4)),
    ]
    results = [
        ca * sd(adj(3), adj(4)),
        ca**4 - ca**2,
        ca**4 / 2 - ca**2 / 2,
        -sd(adj(2), adj(3)) / (4 * ca),
        -ca * sd(adj(3), adj(4)) / 4,
        ii * ((ca**2 - 4) / (4 * ca)) * sd(adj(3), adj(4)),
        ((ca**2 - 4) / ca) * sd(adj(3), adj(4)),
    ]
    for i in range(len(checks)):
        start = time.time()
        hey = checks[i].transform().chain(fullcoloralgorithm(1)).execute()
        hey = hey - results[i]
        hey = hey.expand()
        hey = hey.transform().chain(tonumberofcolours()).execute()
        hey = hey.expand()
        end = time.time()
        dif = end - start
        print(i, hey, dif)
    return None


def checksunsimplifyhermitian():
    """
    # checks that we get the right results for some chosen cases, hermitian generator conventions
    """
    checks = [
        suntr(adj(1), adj(2)),
        suntr(adj(1), adj(2), adj(1), adj(3)),
        sunf(adj(1), adj(2), adj(4)) * suntr(adj(1), adj(2), adj(3)),
        sund(adj(1), adj(2), adj(4)) * suntr(adj(1), adj(2), adj(3)),
        sunf(adj(1), adj(2), adj(3)) * sunf(adj(1), adj(2), adj(4)),
        sund(adj(1), adj(2), adj(3)) * sund(adj(1), adj(2), adj(4)),
    ]
    results = [
        sd(adj(1), adj(2)) / 2,
        -sd(adj(2), adj(3)) / (4 * ca),
        ii * ca * sd(adj(3), adj(4)) / 4,
        ((ca**2 - 4) / (4 * ca)) * sd(adj(3), adj(4)),
        ca * sd(adj(3), adj(4)),
        ((ca**2 - 4) / ca) * sd(adj(3), adj(4)),
    ]
    for i in range(len(checks)):
        start = time.time()
        hey = checks[i].transform().chain(fullcoloralgorithm(ii)).execute()
        hey = hey - results[i]
        hey = hey.expand()
        hey = hey.transform().chain(tonumberofcolours()).execute()
        hey = hey.expand()
        end = time.time()
        dif = end - start
        print(i, hey, dif)
    return None


# checksunsimplifynonhermitian()
# checksunsimplifyhermitian()
