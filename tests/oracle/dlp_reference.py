"""Independent reference implementation of P^2 non-emptiness (E12-M0).

CHARTER
=======
This module is a slow, naive transcription of published THEOREM STATEMENTS. It
must remain independent of the library it is used to test:

  * it imports nothing from the package under test (an integrity test asserts the
    literal package name never appears in this file's source);
  * it uses exact ``fractions.Fraction`` arithmetic only -- there is deliberately
    no floating-point value and no standard-library maths import anywhere,
    including the rank cutoff, which is applied as an INTEGER denominator bound,
    never via a square root;
  * it is written to be obviously-correct, not fast. Naive loops, no memoisation.
    A second per class is fine; a reader must be able to check it against the
    theorem statement line by line.

PRIMARY SOURCES (quoted, not recalled)
======================================
[CHW]  Coskun-Huizenga-Woolf, "The effective cone of the moduli space of sheaves
       on the plane", arXiv:1401.1613, Section 2, Theorem 2.2:

         "There exists a positive dimensional moduli space of semistable sheaves
          M(xi) with Chern character xi if and only if c_1 = r mu in Z,
          chi = r(P(mu) - Delta) in Z and Delta >= delta(mu)."

       Same section: "exceptional bundles are precisely the stable bundles E on
       P^2 with Delta(E) < 1/2"; they "are rigid and their moduli spaces consist
       of a single reduced point"; and  Delta = (1/2) mu^2 - ch_2 / r,  mu = ch_1/r.

[CH]   Coskun-Huizenga, "Existence of semistable sheaves on Hirzebruch surfaces",
       arXiv:1907.06739, Example 1.9: "Non-exceptional mu-stable sheaves of
       character (r, mu, Delta) exist if and only if Delta >= delta(mu)."
       Example 1.14 calls a direct sum of copies of an exceptional bundle a
       SEMIEXCEPTIONAL bundle.

[DLP]  The Drezet-Le Potier delta curve, as tabulated in the Coskun-Huizenga
       survey and in Le Potier's lectures:

         P(m)          = (m^2 + 3 m + 2) / 2                (Euler char of O(m) on P^2)
         Delta_alpha   = (1 - 1/rho^2) / 2                  (rho = rank of E_alpha)
         alpha . beta  = (alpha + beta)/2
                          + (Delta_beta - Delta_alpha) / (3 + alpha - beta)
         delta(mu)     = max over the exceptional slopes alpha of the value the
                         bundle E_alpha forces at mu, clamped below by 1/2.

       The exceptional slopes are exactly the epsilon-recursion image seeded by
       the integers (rank-1 line bundles) under the binary operation ``alpha.beta``.

CONSEQUENCE -- the statement implemented here (P^2 only):

  M(xi) != empty
     <=>  c_1 in Z  AND  chi in Z  AND
          ( Delta >= delta(mu)  OR  xi = m * ch(E), E exceptional, m >= 1 ).

The certified rank cutoff (re-derived in the E12 anchor): an exceptional slope of
rank rho lifts delta above 1/2 at mu only inside its control interval of half-width
x_rho = (3 - S)/2, where S is the (irrational) square root of (9 - 4/rho^2).  If
alpha != mu and q = denominator(mu), then |mu - alpha| >= 1/(q rho); together with
x_rho <= 0.382/rho^2 this forces rho < 0.382 q, and the only remaining case alpha = mu
needs rho = q.  Hence enumerating exceptional slopes of denominator <= denominator(mu)
is EXACT.  We apply that purely as an integer denominator bound; the square root above
is never evaluated -- this file contains no float and no irrational arithmetic at all.
"""

from enum import Enum
from fractions import Fraction
from typing import Union

Number = Union[int, Fraction]


class Status(Enum):
    """Verdict of the P^2 non-emptiness theorem for a Chern character."""

    NONEMPTY = "nonempty"   # a positive-dimensional M(xi), or a single reduced point
    EMPTY = "empty"         # a valid character, but no semistable sheaf exists
    INVALID = "invalid"     # not the Chern character of any sheaf (chi not in Z)


# --------------------------------------------------------------------------- #
# Exact floor / ceil on a Fraction, using only integer and Fraction ops.       #
# --------------------------------------------------------------------------- #
def _floor(x: Number) -> int:
    x = Fraction(x)
    # Python's // is floor division and is exact on ints.
    return x.numerator // x.denominator


def _ceil(x: Number) -> int:
    x = Fraction(x)
    return -((-x.numerator) // x.denominator)


# --------------------------------------------------------------------------- #
# Theorem constants (verbatim).                                                #
# --------------------------------------------------------------------------- #
def _P(m: Number) -> Fraction:
    """Euler characteristic polynomial  P(m) = (m^2 + 3 m + 2) / 2  ([CHW] Sec.2)."""
    m = Fraction(m)
    return (m * m + 3 * m + 2) / 2


def _delta_of_slope(alpha: Number) -> Fraction:
    """Discriminant of the exceptional bundle of slope ``alpha``.

    An exceptional bundle of slope alpha = p/q in lowest terms has rank q, and
    Delta_alpha = (1 - 1/q^2)/2  ([CHW] Sec.2: Delta_alpha = (1 - 1/r_alpha^2)/2).
    """
    q = Fraction(alpha).denominator
    return (Fraction(1) - Fraction(1, q * q)) / 2


def _mediant(a: Number, b: Number) -> Fraction:
    """The epsilon-recursion binary operation  alpha . beta  ([DLP]).

    alpha . beta = (alpha + beta)/2 + (Delta_beta - Delta_alpha)/(3 + alpha - beta).

    For adjacent slopes in a unit interval |a - b| < 1, so 3 + a - b lies in (2, 4)
    and the denominator never vanishes.
    """
    a = Fraction(a)
    b = Fraction(b)
    return (a + b) / 2 + (_delta_of_slope(b) - _delta_of_slope(a)) / (3 + a - b)


# --------------------------------------------------------------------------- #
# The exceptional slopes as the epsilon-recursion image.                       #
# --------------------------------------------------------------------------- #
def _exceptional_slopes(lo: Number, hi: Number, max_denominator: int) -> set:
    """All exceptional slopes in ``[lo, hi]`` with denominator <= ``max_denominator``.

    Seed: every integer is an exceptional (rank-1 line bundle) slope.  Between two
    adjacent slopes a < b we insert g = a . b (the epsilon-recursion child) and
    recurse into (a, g) and (g, b).  The child's denominator is strictly larger
    than both parents' (the DLP/Markov tree property), so the denominator bound
    both prunes exactly and guarantees termination.  Written as an explicit stack;
    no memoisation, no cleverness.
    """
    lo = Fraction(lo)
    hi = Fraction(hi)
    generated = set()

    integers = list(range(_floor(lo) - 1, _ceil(hi) + 2))
    for n in integers:
        generated.add(Fraction(n))

    stack = []
    for i in range(len(integers) - 1):
        stack.append((Fraction(integers[i]), Fraction(integers[i + 1])))

    while stack:
        a, b = stack.pop()
        g = _mediant(a, b)
        if g.denominator > max_denominator:
            continue                # denominators strictly increase: exact prune
        if g in generated:
            continue                # safety net; the DLP tree has no repeats
        generated.add(g)
        stack.append((a, g))
        stack.append((g, b))

    return {s for s in generated if lo <= s <= hi}


# --------------------------------------------------------------------------- #
# Public reference API.                                                        #
# --------------------------------------------------------------------------- #
def reference_discriminant(r: int, c1: int, ch2: Number) -> Fraction:
    """CH discriminant  Delta = (1/2) mu^2 - ch_2 / r,  mu = c_1 / r  ([CHW] Sec.2).

    On P^2 (H^2 = 1) this is *the* discriminant the non-emptiness theorem compares
    against; it is polarization-independent.
    """
    mu = Fraction(c1, r)
    return Fraction(1, 2) * mu * mu - Fraction(ch2) / r


def reference_delta(mu: Number) -> Fraction:
    """The Drezet-Le Potier delta curve at slope ``mu`` ([DLP]).

    delta(mu) = max( 1/2 ,  max over exceptional alpha of  P(-|mu - alpha|) - Delta_alpha ),
    where -- crucially -- alpha contributes at mu only inside its control interval.

    The control intervals of distinct exceptional bundles are DISJOINT: mu lies in
    at most one, whose bundle gives the peak value; for every other exceptional alpha
    the value P(-|mu-alpha|) - Delta_alpha is < 1/2.  The near edge of the control
    interval is exactly where that value equals 1/2, and it always sits at
    |mu - alpha| < 3/2 (the vertex of the parabola m |-> P(-m) is at m = 3/2).  The
    parabola is symmetric about m = 3/2, so restricting to |mu - alpha| <= 3/2 keeps
    every genuine (near-side) contribution and discards the spurious far-side branch
    without ever evaluating a square root.

    Enumerating exceptional slopes of denominator <= denominator(mu) is exact (see
    the module docstring's certified rank cutoff).
    """
    mu = Fraction(mu)
    half = Fraction(3, 2)
    best = Fraction(1, 2)
    for alpha in _exceptional_slopes(mu - half, mu + half, mu.denominator):
        d = abs(mu - alpha)
        if d <= half:                      # near side only: exact, float-free
            v = _P(-d) - _delta_of_slope(alpha)
            if v > best:
                best = v
    return best


def _exceptional_ch2(r: int, c1: int) -> Fraction:
    """The unique ch_2 an exceptional bundle of rank r and c_1 must carry.

    From Delta_alpha = (1 - 1/r^2)/2 and Delta = (1/2)mu^2 - ch_2/r one gets
    ch_2 = (c_1^2 - r^2 + 1) / (2 r).
    """
    return Fraction(c1 * c1 - r * r + 1, 2 * r)


def reference_is_exceptional(r: int, c1: int, ch2: Number) -> bool:
    """True iff (r, c_1, ch_2) is the Chern character of an exceptional bundle on P^2.

    Necessary and sufficient conditions, all checked exactly:
      1. rank equals the slope denominator (mu = c_1/r must already be in lowest terms);
      2. ch_2 equals the unique exceptional value (c_1^2 - r^2 + 1)/(2 r);
      3. mu is a member of the epsilon-recursion image (denominator == r).

    Conditions (1)+(2) are only NECESSARY (they encode chi(E, E) = 1); membership
    in the epsilon-image is what makes it sufficient.  So this REJECTS the impostors
    whose rank is a Markov number but whose slope is not an epsilon-slope, e.g.
    3/10 and 133/610.
    """
    if r <= 0:
        return False
    mu = Fraction(c1, r)
    if mu.denominator != r:
        return False
    if Fraction(ch2) != _exceptional_ch2(r, c1):
        return False
    return mu in _exceptional_slopes(_floor(mu), _floor(mu) + 1, r)


def reference_is_semiexceptional(r: int, c1: int, ch2: Number) -> bool:
    """True iff xi = m * ch(E) for an exceptional bundle E and an integer m >= 1.

    [CH] Ex.1.14: a semiexceptional bundle is a direct sum of copies of an
    exceptional bundle.  m = 1 makes every exceptional bundle semiexceptional too.
    """
    ch2 = Fraction(ch2)
    for m in range(1, r + 1):
        if r % m != 0 or c1 % m != 0:
            continue
        if reference_is_exceptional(r // m, c1 // m, ch2 / m):
            return True
    return False


def _chi(r: int, c1: int, ch2: Number) -> Fraction:
    """Euler characteristic  chi(E) = ch_2 + (3/2) c_1 + r  on P^2 (Riemann-Roch).

    (td(P^2) = 1 + (3/2)H + [pt], H^2 = 1.)  Note the identity
    chi + c_2 = c_1(c_1 + 3)/2 + r, which is an integer for all integer (r, c_1);
    hence chi in Z  <=>  c_2 in Z.  We test Theorem 2.2's literal clause chi in Z.
    """
    return Fraction(ch2) + Fraction(3 * c1, 2) + r


def reference_nonempty(r: int, c1: int, ch2: Number) -> Status:
    """Verdict of [CHW] Thm 2.2 (+ [CH] Ex.1.9/1.14) for the character (r, c_1, ch_2).

        chi not in Z                      -> INVALID (not a Chern character)
        Delta >= delta(mu)                -> NONEMPTY
        xi = m * ch(E), E exceptional     -> NONEMPTY (semiexceptional, incl. m = 1)
        otherwise                         -> EMPTY
    """
    if r <= 0:
        return Status.INVALID
    if _chi(r, c1, ch2).denominator != 1:
        return Status.INVALID
    disc = reference_discriminant(r, c1, ch2)
    mu = Fraction(c1, r)
    if disc >= reference_delta(mu):
        return Status.NONEMPTY
    if reference_is_semiexceptional(r, c1, ch2):
        return Status.NONEMPTY
    return Status.EMPTY
