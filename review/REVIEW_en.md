# Referee report on Zenodo v5

**Paper:** R. Chen, *The Geometry of Prime Vacuums: Legendre's Conjecture in Function Fields via Monodromy and Chebotarev Density*, Zenodo, 20 Feb 2026, [10.5281/zenodo.18705744](https://doi.org/10.5281/zenodo.18705744). The file reviewed is `paper/archive/Legendre_v5_zenodo.tex`.

**Recommendation:** major revision. The asymptotic statement is true, but it is not new. Three proofs have gaps, and the explicit threshold is unproved. All of these issues are fixed in the revised manuscript `paper/Legendre_FF_v6.tex`.

## 1. Summary of the submission

The paper considers the Legendre interval I_f = {f² + s : deg s ≤ d} for monic f ∈ F_q[t] of degree d, with f squarefree and char F_q > 2d. It claims:

- N_irr = q^{d+1}/(2d) + O_d(q^{d+1/2}).
- N_irr ≥ 1 for q > (8d)^{2d+6}.

The proof goes through full S_{2d} monodromy, an effective Chebotarev theorem, and a Katz Betti-number bound.

## 2. Major issues

**M1. The main theorem is a known special case.** Bank, Bary-Soroker and Rosenzweig ([Duke Math. J. 164 (2015)](https://arxiv.org/abs/1302.0625), Thm 2.3 and Cor 2.4) prove the factorization-type count in every short interval f + F_q[t]_{≤m} with 1 ≤ m < k. The only extra conditions are for p | k(k−1) and for p = 2. Taking k = 2d, m = d and f ↦ f² gives the main theorem for every odd q. That version needs no squarefree hypothesis and no condition p > 2d. Their abstract even mentions breaking the ε = 1/2 barrier, which is exactly the Legendre case. The submission cites this paper only as a possible route to "sharper bounds". The introduction and the "Relation to prior work" section must say plainly that the asymptotic is known.

**M2. The primitivity proof (Prop 2.3) has a gap.** The step "imprimitive ⇔ F = g∘h" is quoted from Müller and Turnwald. That equivalence is about Gal(G(t) − x / K(x)) with x transcendental and occurring only in the constant term. It does not hold for the Galois group of an arbitrary polynomial over a field. The proof never isolates a₀ to put itself in that setting. The closing argument ("a map A^n → A^{d+1} with n < d+1 is not surjective") treats the coefficients of g as free parameters in A^n. In fact they lie in an algebraic extension of F̄_q(a). The conclusion is true, but the proof as written does not establish it.

**M3. The transposition proof (Prop 2.4) has errors.**

- The local equation is written F(α_j + u) = f′(α_j)²u² + P(α_j) + O(u·a) + O(u³). It then asserts that the branches of the discriminant near a = 0 are {P(α_j) = 0}. When P(α_j) = 0 but P′(α_j) ≠ 0, the double root splits into two roots, so the branch is mis-identified.
- The variable a_d is dropped.
- "Picard–Lefschetz" is invoked in characteristic p without the hypotheses that make it applicable.
- The claim that p > 2d is needed for f′(α_j) ≠ 0 is wrong. Squarefreeness alone gives it.

**M4. The effective Chebotarev theorem (Thm 3.2) is mis-stated.** The bound is expressed through the Betti numbers of the permutation sheaf π_*Q̄_ℓ. Frobenius traces on that sheaf count fixed points, which correspond to linear factors. They do not detect 2d-cycles. Counting 2d-cycles needs every irreducible character of S_{2d} that is nonzero on the class; these are the 2d hook characters. The citation "Katz 1988, Thm 9.2.6" does not point to an effective Chebotarev theorem.

**M5. The explicit threshold is unproved.** The bound B(F) ≤ N(2δ+1)^n is attributed to Katz, *Sums of Betti numbers in arbitrary characteristic*, FFA 7 (2001), Thm 12 ([author's copy](https://web.math.princeton.edu/~nmk/BettiSum14.pdf)). Theorem 12 there bounds Betti numbers of L_ψ(f) ⊗ ⊗_j L_{ρ_j}(G_j) on open subsets of affine varieties. It says nothing about push-forwards of finite covers. Katz's general results for compatible systems (Thms 6 and 9) are not explicit. So the threshold q > (8d)^{2d+6} has no proof.

**M6. The error terms are conflated in Step 4.** The unspecified O_{N,n}(q^{n−1}) term of Theorem 3.2 is replaced by the Lang–Weil count of D(F_q). Those are different quantities. The D(F_q) count is not needed at all, because Step 3 already shows N_irr = N_irr°.

## 3. Minor issues

1. Equation (4.1) writes χ(P¹, j_!F) = 2r − Σ(drop_v + Sw_v). That is the formula for j_*. Section 4 is not used anywhere in the proof.
2. The citation [barysoroker2012] is wrong. The paper is *Proc. AMS* **137** (2009), 73–83, and it works over PAC fields, not through Chebotarev over finite fields.
3. The Arnold–Gusein-Zade–Varchenko entry has key `AGLV1993` but year 1988.
4. Seven bibliography entries are never cited: deligne1977, ritt1922, rosen2002, serre1992, weil1948, zannier2008, and grothendieck1977 in substance.
5. The body contains version-history remarks such as "erroneously used in earlier versions" and "incorrectly attempted in an earlier version". These belong in a changelog.
6. Remark 4.1 infers "N_irr ≥ 1 for all q and d" from three examples with q ≤ 11 and d ≤ 2. That is far too little evidence.
7. Section 1.2, on the "Titan Project" context, overstates the link to the conductor-rigidity papers.
8. The Jordan theorem citation, Dixon–Mortimer Thm 3.3E, appears to point at a different statement (the prime-cycle version). The "primitive + transposition ⇒ S_n" statement is Thm 3.3A.

## 4. Verification of the numerical claims

All numerical claims in v5 were reproduced by brute force (`code/test_fflib.py`):

| Claim in v5 | Recomputed |
|---|---|
| d=1, q=3,5,7,11: N_irr = 3, 10, 21, 55 | 3, 10, 21, 55 ✓ |
| d=2, f=t²+1, q=11: N_irr = 330 | 330 ✓ (= (q³−q)/4) |

## 5. What survives and what was added in the revision

- **Theorem A.** G_geom = S_{2d} for p > 2d, without the squarefree hypothesis. The proof is new and self-contained: it applies Lüroth with a₀ isolated, and uses a DVR inertia lemma at an explicitly constructed fibre with exactly one double root.
- **Theorem B.** The exact identity 2d·N_irr = #W(F_q) − E. Here W is an explicit twisted root variety with deg W ≤ (d−1)!, and 0 ≤ E < q^{d+1}/(q−1). This replaces the cohomological Chebotarev step entirely.
- **Theorem C.** An explicit, proved bound. It follows from the Cafure–Matera explicit Lang–Weil theorem, [FFA 12 (2006), Thm 7.1](https://arxiv.org/abs/math/0405302). For example the threshold is q > 1.65·10⁴ for d=4 and q > 7.3·10⁶ for d=5. A second route through Katz's Betti bound and Deligne's theorem is also proved.
- **Theorem D.** Closed formulas for d ≤ 3, with N_irr ≥ 1 for every admissible q. For d=3 the formula is a special case of Kuz'min's two-coefficient formula.
- **Computations.** 186 exact counts for d ≤ 6, all of I_f factored. There are 0 mismatches with Theorem D. The twisted identity was verified by enumerating F_{q^{2d}}, and an exhaustive search over all f was run for small q. No interval without an irreducible was found.
