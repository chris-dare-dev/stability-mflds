"""Second independent P^2 reference (E13-M4): exceptional bundles by K-THEORETIC MUTATION.

CHARTER -- why a second oracle exists
=====================================
``tests/oracle/dlp_reference.py`` is import-independent of the package, but it
generates the exceptional slopes by the SAME algorithm production uses: the
Drezet-Le Potier epsilon-mediant interval subdivision, seeded on the integers.
The E13 adversarial re-audit named that a COMMON-MODE risk: a shared
misunderstanding of the epsilon-recursion (the formula, the seeding, the tree
walk, the cutoff) would reproduce itself identically on both sides and the
differential gate would be blind to it.

This module is the hardening.  It generates the same objects from a DIFFERENT
THEOREM, in different coordinates, with different arithmetic:

  * classes live in the numerical K-group K(P^2) = Z^3, recorded as INTEGER
    triples ``(r, c1, chi)`` with ``chi`` the Euler characteristic -- no slope,
    no Fraction, appears anywhere in the generator;
  * generation is by MUTATION of full exceptional collections -- pure integer
    linear algebra:  ``L_A B = chi(A,B).A - B``,  ``R_B A = chi(A,B).B - A``
    ([GR]) -- starting from the line-bundle collection
    ``(O, O(1), O(2)) = ((1,0,1), (1,1,3), (1,2,6))``;
  * completeness is [GR]'s constructibility theorem (every exceptional bundle
    on P^2 arises from that collection by iterated mutation), not the
    epsilon-image description;
  * NO mediant, NO epsilon-recursion, NO dyadic tree, and NO ``Delta_alpha``
    rank formula: the discriminant of a generated class is COMPUTED from its
    own ``(r, c1, chi)``, and ``chi(X, X) = 1`` is asserted, never assumed.

The two generators can therefore agree only if both theorems are transcribed
correctly.  ``tests/test_mutation_oracle.py`` additionally pins a THIRD,
purely number-theoretic recursion against both -- Springborn's Markov
fractions ([V]) -- and sweeps the impostor family (classes with ``chi = 1``,
integral ``c2``, Markov rank, and even ``p^2 = -1 (mod q)``, that are still
NOT exceptional, e.g. ``(610, 133, -581/2)``).

Like ``dlp_reference`` it imports nothing from the package under test, uses
exact int/Fraction arithmetic only (no float, no square root, no ``math``
import), and is written to be obviously-correct, not fast.

PRIMARY SOURCES (quoted, not recalled)
======================================
[GR]  Gorodentsev-Rudakov, "Exceptional vector bundles on projective spaces",
      Duke Math. J. 54 (1987), 115-130.  Mutations of exceptional pairs and
      collections; on P^2 every exceptional bundle is obtained from the
      line-bundle collection by iterated mutations (constructibility).  In the
      numerical K-group a mutation acts by
          [L_A B] = chi(A, B) [A] - [B],     [R_B A] = chi(A, B) [B] - [A]
      (the class of the mutated object, up to the shift sign, which is
      normalized here to positive rank).

[R]   Rudakov, "Exceptional vector bundles on P^2 and Markov numbers",
      Izv. Akad. Nauk SSSR Ser. Mat. 52 (1988); Engl. transl. Math. USSR-Izv.
      32 (1989), 99-112.  The ranks (a, b, c) of a full exceptional collection
      on P^2 satisfy the Markov equation  a^2 + b^2 + c^2 = 3abc,  and an
      exceptional bundle is determined by (rank, slope).  The Markov equation
      is asserted on EVERY collection this generator visits.

[CHW] Coskun-Huizenga-Woolf, arXiv:1401.1613, Thm 2.2 (the non-emptiness
      verdict), and [CH] Coskun-Huizenga, arXiv:1907.06739, Ex. 1.9/1.14
      (the (semi)exceptional disjunct) -- exactly as in dlp_reference.

[V]   Veselov, "Markov fractions and the slopes of the exceptional bundles on
      P^2", arXiv:2501.06779 (Markov fractions after B. Springborn; used as
      the third-recursion pin in the tests, not imported here).

DERIVED CLOSED FORMS (each re-derived and pinned in tests/test_mutation_oracle.py)
==================================================================================
Riemann-Roch on P^2 gives  chi(E) = ch2 + (3/2) c1 + r,  so the coordinates
convert by  ch2 = chi - (3/2) c1 - r.  Substituting into the RR Euler pairing
chi(E,F) = r_E ch2_F + r_F ch2_E - c_E c_F + (3/2)(r_E c_F - r_F c_E) + r_E r_F
collapses to the ALL-INTEGER form used here:

    chi(E, F) = r_E chi_F + r_F chi_E - r_E r_F - 3 r_F c_E - c_E c_F.

Twisting by O(n) (an equivalence, so it preserves exceptionality):

    (r, c, chi) . O(n) = (r, c + n r, chi + n c + (n(n+3)/2) r),

and in (r, c1, ch2) coordinates  ch2 -> ch2 + n c1 + (n^2/2) r.
Pinned against  chi(O(a), O(b)) = (b-a+2)(b-a+1)/2  (cohomology of line
bundles on P^2) over an integer grid.
"""

from fractions import Fraction
from typing import Union

Number = Union[int, Fraction]

#: The standard line-bundle collection (O, O(1), O(2)) in (r, c1, chi) coordinates:
#: chi(O(n)) = (n+2)(n+1)/2 gives 1, 3, 6.
_ROOT = ((1, 0, 1), (1, 1, 3), (1, 2, 6))

#: Hard iteration bound for the collection walk -- a tripwire, far above any
#: real traversal (64 canonical collections suffice for rank <= 610).
_MAX_STEPS = 500000


# --------------------------------------------------------------------------- #
# Exact floor / ceil on a Fraction (integer ops only; same idiom as             #
# dlp_reference -- no math import anywhere in this file).                       #
# --------------------------------------------------------------------------- #
def _floor(x: Number) -> int:
    x = Fraction(x)
    return x.numerator // x.denominator


def _ceil(x: Number) -> int:
    x = Fraction(x)
    return -((-x.numerator) // x.denominator)


# --------------------------------------------------------------------------- #
# The integer Euler form and the twist, in (r, c1, chi) coordinates.            #
# --------------------------------------------------------------------------- #
def mutation_euler(E, F) -> int:
    """``chi(E, F)`` on K(P^2) in (r, c1, chi) coordinates (module docstring).

    All-integer: derived from the RR pairing by eliminating ch2 via
    ``ch2 = chi - (3/2) c1 - r``; the halves cancel exactly.
    """
    rE, cE, xE = E
    rF, cF, xF = F
    return rE * xF + rF * xE - rE * rF - 3 * rF * cE - cE * cF


def _twist(X, n: int):
    """``X (x) O(n)``: (r, c, chi) -> (r, c + n r, chi + n c + (n(n+3)/2) r).

    ``n(n+3)`` is always even, so this is integer arithmetic.  Twisting is an
    equivalence of categories: it preserves exceptionality and collections.
    """
    r, c, x = X
    return (r, c + n * r, x + n * c + ((n * (n + 3)) // 2) * r)


def _normalize_class(X):
    """Twist ``X`` so its slope ``c/r`` lies in [0, 1) (r > 0 assumed)."""
    r, c, _ = X
    return _twist(X, -(c // r))


def _canonical_collection(col):
    """Twist the ordered collection so the FIRST member's slope lies in [0, 1).

    Twisting acts diagonally on a collection and preserves it; this fixes the
    twist ambiguity so the visited-set walk terminates.
    """
    r0, c0, _ = col[0]
    n = -(c0 // r0)
    return tuple(_twist(X, n) for X in col)


def _norm_sign(X):
    """Normalize the mutated class to positive rank (the shift ambiguity)."""
    r, c, x = X
    if r < 0:
        return (-r, -c, -x)
    if r == 0:
        raise AssertionError(f"mutation produced a rank-0 class: {X!r}")
    return X


def _lmut(A, B):
    """[L_A B] = chi(A,B) [A] - [B]  ([GR]), sign-normalized."""
    k = mutation_euler(A, B)
    return _norm_sign((k * A[0] - B[0], k * A[1] - B[1], k * A[2] - B[2]))


def _rmut(A, B):
    """[R_B A] = chi(A,B) [B] - [A]  ([GR]), sign-normalized."""
    k = mutation_euler(A, B)
    return _norm_sign((k * B[0] - A[0], k * B[1] - A[1], k * B[2] - A[2]))


# --------------------------------------------------------------------------- #
# The mutation walk over full exceptional collections.                          #
# --------------------------------------------------------------------------- #
def mutation_collections(rank_max: int):
    """All canonical full exceptional collections reachable by mutation, pruned
    at ``rank_max``.

    Walks the four one-step mutations of ``(A, B, C)``:

        (L_A B, A, C),  (A, L_B C, B),  (B, R_B A, C),  (A, C, R_C B),

    deduplicating collections modulo twist (:func:`_canonical_collection`) and
    pruning a move when the NEWLY created class has rank > ``rank_max``.  By
    [R] the rank triples are Markov triples and the tree's maxima strictly
    increase away from the root (1, 1, 1), so the prune is exact: every
    exceptional bundle of rank <= rank_max (mod twist) appears in a retained
    collection.  Two tripwires run on every visited collection:

      * the Markov equation  a^2 + b^2 + c^2 = 3abc  on the ranks ([R]);
      * ``chi(X, X) = 1`` for each member (exceptionality is unit self-pairing).

    Returns the set of canonical collections (tuples of (r, c1, chi) triples).
    """
    if rank_max < 1:
        return set()
    root = _canonical_collection(_ROOT)
    seen = {root}
    stack = [root]
    steps = 0
    while stack:
        steps += 1
        if steps > _MAX_STEPS:
            raise AssertionError("mutation walk exceeded the iteration tripwire")
        A, B, C = stack.pop()
        a, b, c = A[0], B[0], C[0]
        if a * a + b * b + c * c != 3 * a * b * c:
            raise AssertionError(
                f"Markov equation failed on collection ranks {(a, b, c)!r} -- "
                "the mutation transcription is broken ([R])")
        for X in (A, B, C):
            if mutation_euler(X, X) != 1:
                raise AssertionError(f"chi(X,X) != 1 for generated class {X!r}")
        for col in (
            (_lmut(A, B), A, C),
            (A, _lmut(B, C), B),
            (B, _rmut(A, B), C),
            (A, C, _rmut(B, C)),
        ):
            new = [X for X in col if X not in (A, B, C)]
            if new and new[0][0] > rank_max:
                continue                      # exact prune: maxima increase
            colc = _canonical_collection(col)
            if colc not in seen:
                seen.add(colc)
                stack.append(colc)
    return seen


def _to_rch2(X):
    """(r, c1, chi) -> (r, c1, ch2) with ch2 = chi - (3/2) c1 - r (exact Fraction)."""
    r, c, x = X
    return (r, c, Fraction(x) - Fraction(3 * c, 2) - r)


def mutation_exceptional_classes(rank_max: int):
    """Every exceptional class of rank <= ``rank_max``, slope normalized to [0, 1).

    Returns a set of ``(r, c1, ch2)`` triples (``ch2`` an exact ``Fraction``).
    Completeness is [GR] constructibility + the Markov-tree prune ([R]).
    """
    classes = set()
    for col in mutation_collections(rank_max):
        for X in col:
            Xn = _normalize_class(X)
            if Xn[0] <= rank_max:
                classes.add(_to_rch2(Xn))
    return classes


def mutation_exceptional_in_window(lo: Number, hi: Number, rank_max: int):
    """Exceptional classes of rank <= ``rank_max`` with slope in ``[lo, hi]``.

    Twists of the normalized base set: twisting by O(n) shifts the slope by n
    and is an equivalence, so this loses (and invents) nothing.
    """
    lo = Fraction(lo)
    hi = Fraction(hi)
    out = set()
    for (r, c1, ch2) in mutation_exceptional_classes(rank_max):
        alpha0 = Fraction(c1, r)
        for n in range(_floor(lo - alpha0), _ceil(hi - alpha0) + 1):
            alpha = alpha0 + n
            if lo <= alpha <= hi:
                out.add((r, c1 + n * r,
                         ch2 + n * c1 + Fraction(n * n, 2) * r))
    return out


# --------------------------------------------------------------------------- #
# Membership, the delta curve, and the non-emptiness verdict.                   #
# --------------------------------------------------------------------------- #
def mutation_is_exceptional(r: int, c1: int, ch2: Number) -> bool:
    """True iff ``(r, c1, ch2)`` is the Chern character of an exceptional bundle.

    Twist-normalizes the queried class into slope [0, 1) and tests FULL-TRIPLE
    membership in the mutation-generated set at ``rank_max = r``.  No slope
    denominator test and no exceptional-ch2 formula: membership of the whole
    triple is the theorem-level criterion ([GR] constructibility).
    """
    if r < 1:
        return False
    n = -(c1 // r)
    key = (r, c1 + n * r, Fraction(ch2) + n * c1 + Fraction(n * n, 2) * r)
    return key in mutation_exceptional_classes(r)


def mutation_is_semiexceptional(r: int, c1: int, ch2: Number) -> bool:
    """True iff ``(r, c1, ch2) = m . ch(E)`` for exceptional ``E``, ``m >= 1``
    ([CH] Ex. 1.14; ``m = 1`` includes the exceptional bundles themselves)."""
    ch2 = Fraction(ch2)
    for m in range(1, r + 1):
        if r % m != 0 or c1 % m != 0:
            continue
        if mutation_is_exceptional(r // m, c1 // m, ch2 / m):
            return True
    return False


def _P(m: Number) -> Fraction:
    """Euler polynomial  P(m) = (m^2 + 3m + 2)/2  ([CHW] Sec. 2; = chi(O(m)) at
    integer m, which the tests pin against :func:`mutation_euler`)."""
    m = Fraction(m)
    return (m * m + 3 * m + 2) / 2


def mutation_delta(mu: Number, rank_max: "int | None" = None) -> Fraction:
    """The Drezet-Le Potier curve ``delta(mu)`` over the MUTATION-generated set.

    Same theorem statement as ``dlp_reference.reference_delta`` (near-side
    branch on ``|mu - alpha| <= 3/2``, clamped below by 1/2), but with two
    deliberate independence choices:

      * the exceptional set comes from the mutation walk, not the
        epsilon-recursion;
      * ``Delta_alpha`` is computed from each generated class's OWN
        ``(r, c1, ch2)`` -- the rank formula ``(1 - 1/rho^2)/2`` is never
        transcribed here (it is forced by ``chi(X,X) = 1``, which the walk
        asserts);
      * the default rank cutoff is ``4 q + 64``, strictly ABOVE the certified
        sharp cutoff ``q = denominator(mu)`` that production and dlp_reference
        both use.  If the sharp-cutoff theorem (docs/CORRECTIONS.md Sec. 8)
        were wrong -- if some bundle of rank in ``(q, 4q + 64]`` actually
        bound -- the differential against production would now FAIL instead of
        agreeing blindly.  That retires the cutoff as a common mode too.
    """
    mu = Fraction(mu)
    if rank_max is None:
        rank_max = 4 * mu.denominator + 64
    half = Fraction(3, 2)
    best = Fraction(1, 2)
    for (r, c1, ch2) in mutation_exceptional_classes(rank_max):
        alpha0 = Fraction(c1, r)
        # Delta of the class itself -- never the (1 - 1/rho^2)/2 transcription.
        delta_alpha = alpha0 * alpha0 / 2 - Fraction(ch2) / r
        for n in range(_floor(mu - half - alpha0), _ceil(mu + half - alpha0) + 1):
            d = abs(mu - (alpha0 + n))
            if d <= half:                      # near-side branch only (exact)
                v = _P(-d) - delta_alpha
                if v > best:
                    best = v
    return best


def mutation_semistable_exists(r: int, c1: int, ch2: Number) -> bool:
    """[CHW] Thm 2.2 + [CH] Ex. 1.9/1.14 over the mutation-generated machinery.

        chi = ch2 + (3/2) c1 + r  not in Z   ->  False (not a Chern character)
        Delta >= delta(mu)                   ->  True
        (semi)exceptional                    ->  True
        otherwise                            ->  False
    """
    if r < 1:
        return False
    ch2 = Fraction(ch2)
    if (ch2 + Fraction(3 * c1, 2) + r).denominator != 1:
        return False
    mu = Fraction(c1, r)
    disc = mu * mu / 2 - ch2 / r
    if disc >= mutation_delta(mu):
        return True
    return mutation_is_semiexceptional(r, c1, ch2)
