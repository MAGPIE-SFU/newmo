# Notes

Here are some notes (from CC) about the code that builds the NFDS model including $x(t)$ and $G$ 

The overall flow is: 

`cleanmetadata.R`  which makes the file `meta_clean.csv`. This file
has the metadata, but with some new serotype names to match Nick's
invasiveness, the groups defined (now there are 87 of them) and the
ages fixed. 

`make_invasiveness_lookup.R` makes the file `invasiveness_lookup.csv`. It adjusts the
invasiveness weights, pulling towards the mean. This avoids having just one IPD case end
up contributing a large fraction of the synthetic carriage population in a year, when we
come to that.

The reweighting script, `reweight.R` creates the variable $x(t)$ (`xt` in the code), in
`xt_grouped.csv`. This is the synthetic carriage population. 

It, in turn, is needed to create the locus frequencies, because we can't know how frequent
a locus was in year $t$ without knowing the population composition. And locus frequencies
are needed for two reasons: first, we need to select which loci will go in the model (and
that's usually between 5 and 95% frequency). Also, we like to order the columns of G by
locus stability, and we can't know the stability before we know the frequencies.

Therefore, the $G$ script, `groupG.R`, reads the raw $G$ matrix (now in
`pa_matrix_final.csv`) and the $x(t)$ information. It groups $G$ by the ST-SC groups, sets
$G_{ij}$ to the fraction of isolates in group $i$ with gene $j$, computes locus
frequencies, saves them in `locus_freqs_deltas.csv`, then reorders $G$'s columns and saves
`Gnfds.csv`.

Finally, `create-SA-nfdsmodel.R` pulls the information together to make an NFDS model with
the initial condition $x_0$, $G$, the serotype (for record-keeping), the group, the
equilibrium locus frequencies and for context the logit frequency difference (see below).

| Script | Inputs | Outputs |
|---|---|---|
| `cleanmetadata.R` | `metadata.csv` | `meta_clean.csv` |
| `make_invasiveness_lookup.R` | `meta_clean.csv`, Nick's two xlsx tables | `invasiveness_lookup.csv` |
| `reweight.R` | `meta_clean.csv`, `invasiveness_lookup.csv` | $x(t)$ (`xt_grouped.csv`) |
| `groupG` | $x(t)$, `pa_matrix_final.csv`, `meta_clean.csv` | `locus_freqs_deltas.csv`, `Gnfds` |
| `create_SA_nfdsmodel.R` | Gnfds, metadata, locus freqs, xt | `samodel.rds`|

---
## 2026-10-20: serotype distribution questions
I have been experimenting with the model and looking at the serotype distributions in the extrapolated carriage population. Why is 35B supposedly giving over 20% of the carriage? This doesn't seem right and isn't consistent with other observations. Skosana et al https://www.scielo.org.za/scielo.php?script=sci_arttext&pid=S0256-95742021000200014 don't report such a skewed distribution. The 35B portion happens because 35B has a low invasiveness, so the IPD cases we have get pushed to large numbers in the extrapolated carriage population. We may want to decrease the parameter $k$ in the reweighting script. Claude gave me a bit of background about what the ESS is for our IPD samples. 

Each isolate $i$ carries weight $w_i = 1/v_{s(i)}$. We estimate a group’s share as $\hat p = \sum_{i \in g} w_i /\sum_i w_i$. Writing $W = \sum_i w_i$ and treating each isolate as independently in group $g$ with probability $p$:

$$\text{Var}(\hat p) = \frac{1}{W^2}\sum_i w_i^2 \text{Var}\bigl(\mathbb{1}[i \in g]\bigr) = p(1-p)\frac{\sum_i w_i^2}{\left(\sum_i w_i\right)^2} = \frac{p(1-p)}{\text{ESS}}$$

with $\mathrm{ESS} = (\sum w_i)^2 / \sum w_i^2$ (because the ESS is the number of independent samples you would need to get the observed variance). These ESS values turn out to be around 100 (183 in 2009 for the IPD samples and 90 by 2013). This suggests $k$ of more like 100 than my earlier 780 (which was based on getting the carriage data to be about 0.27 of the total - but in carriage years I guess it should be more). 
Anyway that won't help with 35B being 20% of the extrapolated carriage population in 2014, because there's no carriage data in 2014. For that, we might want to think about resampling $x(t)$ using resampled invasiveness values, fiddling some of the invasiveness values (e.g. by digging through Nick's paper some more), or smoothing the invasiveness values more such that we match the effective number of serotypes in typical post-vaccine carriage studies. 

---
## 2026-10-02: rare NVT groups
See the model testing notes. We need rare NVT to be able to rise. Need to change the way I make the groups. 
Also, should I not be making the groups based at least a little on how many carriage cases I might have expected from the IPD counts? Not just on the raw number of isolates? Sigh. 

Going to at least use only pre-vaccine counts to create the rare NVT groups. 

## 2026-09-25: IPD vs carriage in G

$G_{j\ell}$ is the fraction of all isolates in group $j$ that have locus $\ell$, pooling IPD and
carriage.  $x(t)$ is a
synthetic *carriage* population, but about 73% of the isolates behind G are IPD (4952 of 6777). If
the accessory genome within a group differs between invasive and carried isolates (which
they can) then the $G$ entries might not reflect the
circulating population. I decided not to worry about this. See code
`whatabout-Gentries-ipdmix.R`. 

43 of the 87 lineages have at least 10 isolates of each type (`both` in the code). If I substitute
carriage-only prevalences for those, the population locus frequencies move by a median
of 0.043. That's not great because the locus frequencies are often small, but it's not
huge. Of course, most genes aren't specifically associated with IPD or carriage, so maybe
that's fine.

Permuting the IPD/carriage labels within a lineage (keeping the sample sizes, but dropping
any association with IPD for a locus) gives a median change in frequency of 0.03. So most of the apparent shift
is noise in the carriage estimates, which rest on about 21 isolates per group and only from
2009-2013. I decided not to change this. The correction would come from the same  carriage data so
its own error would probably be on this scale anyway. 

Idea: However, we could consider looking in these data for IPD-associated genes, where they are
in the decotur network, and so on. I saved some code to create the locus log odds ratios
and sample numbers, for the 43 groups that have samples in both IPD and carriage. See
`whatabout-IPDgenes.R`. I left a couple of thoughts at the bottom of that script. Would
need a comparison to known IPD-associated (or non-associated, i guess) genes. 

Still to do: check whether dropping the top logOR_ipd loci (or bottom?) changes anything in the model
(since these are the ones for which the Gnfds entries reflect the circulating population
least well). 


---

## 2026-09-24: ranking loci by |logit(post) - logit(pre)|

I replaced $\delta = (f-e)^2/(1-e(1-e))$ (which previous authors used to order genes for NFDS) with logitdiff = |logit(post2011) - logit(pre2009)|,
sorted in ascending order to list the most stable loci first, restricted to loci with pre-vaccine frequency
between 0.05 and 0.95.

The $\delta$ statistic, $\delta = \frac{ (e-f)^2}{1 - e(1-e)}$, has a scale that depends  on
frequency. For example, two loci that both halve, say A from 0.06 to 0.03 and B from 0.50
to 0.25, differ 87-fold. If two loci both drop by the same absolute amount, $\delta$ is
nearly the same (0.5 to 0.47 and 0.06 to 0.03) which seems weird. The $1-e(1-e)$ denominator
only ranges from 0.75 to 1 so it does not really boost intermediate frequency genes much
(I think that is what they intended). 

I checked this by regressing log(measure) and log(pre-vaccine frequency) over the
loci in the model for 4 measures: $\delta$, $z_\ell = \frac{f_\ell - e_\ell}{\sqrt{\dfrac{e_\ell(1-e_\ell)}{n_{\mathrm{pre}}}
         + \dfrac{f_\ell(1-f_\ell)}{n_{\mathrm{post}}}}} $, logitdiff is |logit(e)  -
		 logit(f)| and the last is the log ratio. $e$ and $f$ are pre- and post- vaccine
		 frequency. 
		 Here are the slopes: 

```
    delta       +0.65
    |z|         +0.11
    logitdiff   -0.08
    |log(f/e)|  -0.72
```

$\delta$ and the log
ratio have opposite "preferences" for the baseline or starting frequency. $\delta$ calls low-frequency genes stable because they have
no absolute room to move, but log ratio tends to consider high-frequency genes stable because they
have little proportional room to move (at least, only in one direction). The logitdiff is
least dependent on frequency (as is $z$ which has the bonus of having a statistical
interpretation and depending on some counts, so it's a possibility). 

The logitdiff measure is also invariant if you swap presence for absence, so a change from 0.9 to 0.95 scores the same as
one from 0.1 to 0.05. I'm not that committed to the symmetry point (the idea that not having a gene is just as
much a "trait that could be under NFDS" as having one). But it's a thought. $\delta$ and
$z$ are also invariant to this switch but the log ratio one isn't. 

If I choose the most stable 200 genes using $\delta$ and the logitdiff, 160 of the top 200
agree.   The other 40 are more intermediate in frequency than your typical locus, which
seems good. 

Locus selection for the model: we remove genes with frequencies outsize the range $[0.05,
0.95]$  in the pre-vaccine period, because the model treats pre-vaccine
frequencies as the equilibrium. Filtering on post-vaccine
frequency could drop  loci that moved out of that window, which seems unfair (could drop
ones that move a lot, making them more likely to get out of the range, and then use a model that assumes they won't move too much). 
13 loci have frequency 1 in every lineage — rlmN, rnhB, prmA, pepA, thiN and friends,
i.e. housekeeping genes. They appear core now even though the pangenome's
core/accessory split called them accessory at some point. They are dropped from the NFDS model with 
`in_model = FALSE`. 

Model results are likely to depend on how the NFDS ranking is done. Consider removing DR
genes, IPD-associated genes and similar from the most stable set used in the NFDS term,
perhaps. 

Thought: compare the frequency dynamics to what we would expect without NFDS (started down
this path but didn't continue much; choosing the right null is a bit fraught), and to what
we might expect from drift (but this is hard because the frequencies are not changing by
inheritance with genes moving in / out like alleles getting inherited or not, they have
this lineage structure. Still could be interesting - what are the frequency dynamics, and
are they different from what we'd expect under .. (some alternative). 


---

## 2026-09-24 — uncertainty in x(t)


I looked at this but didn't use it (yet) to do anything. 
We can  resample the isolate counts within each year, age group, sample type combination
(multinomial, conditioning on the cell total, since the number sequenced was fixed by
design and what's uncertain is which groups they were), then re-run the whole reweighting.
There is code there that does 200 draws (not yet pushed). 

Group weights turn out to be quite uncertain — median CV 0.51 across group-years holding
more than 0.5% of the population, above 1 for the worst, which are groups resting on one
IPD case of a low-invasiveness serotype and then amplified by 1/v. Locus frequencies are
much better determined, since a locus frequency sums over many groups and the errors
mostly cancel. Comparing the count-based spread to the simple $\sqrt{f(1-f)/n}$
approximation, they agree within about a factor of two for most loci and only diverge at
extreme frequencies. 

The group-level uncertainties should be taking into account  for model fitting - let's not
try to fit $x(t)$ super precisely. Have put some standard errors in. 

This does not use uncertainty in the invasiveness values (held at their shrunk point
estimates, and probably the bigger uncertainty, but at least we did take that into account
in the reweighting), in k, tauA, tauU, and the fact that a
group absent from a year stays absent. 

Also — Nick has post-vaccine ("gamma") invasiveness values that I'm not using, and since
this whole thing is a pre/post comparison, invasiveness changing between the periods would
bias it directionally rather than just adding noise. Worth a sensitivity run.

---

## 2026-09-23 — Dirichlet blend for carriage + IPD

This all happens in the `reweight.R` file. 

I replaced the fixed $\omega = 0.27$ weight for carriage compared to IPD, and I dropped
the two pseudocounts, upon reflection. Instead, I used
`p_final = (N_car + k * p_C_aug) / (N_car + k)`, with `k = 780`. 

The old version was `(1-omega) p_C_aug + omega q_car_child`, which rearranges into exactly the
posterior mean of a Dirichlet-multinomial thing, with the IPD-derived estimate as the prior mean. 
Writing it in the way I had before with an omega value makes the
carriage weight `N_car/(N_car + k)`, in each year. This means, first, that this is not a
huge change, and also it lets the weight vary with the
sample size (which doesn't vary much here but this Dirichlet thing is statistically
better because its conjugate to the multinomial etc).

The interpretation of $k$ is the
number of carriage swabs an IPD-derived estimate counts for, basically. I don't know that either but 
k = 780 reproduces my earlier `omega = 0.27` for $N$ in the 300 range (290), and about
that fraction of our data is from carriage, so... 

I also dropped the pseudocounts. alpha was giving carriage-only groups a share of the
IPD-derived estimate (weird). beta was the giving IPD-only groups a
made-up carriage share, about 0.5% of `p_final` in carriage years. Both weird. Dropped. 

Note: we should consider looking at dependence on $k$. The carriage data are from two sites, different methods, small
fraction of all carriage in the country. These are all arguments for a big k, meaning less weight on carriage.
Eventually we could run the model at with k = 400, 780, 1500 and see if the fit is better,
stable, or even perhaps possible at some values of k but not others. 

Note: we don't use adult data. Set the relevant fraction to 0.01 just so as not to break
anything. Previous papers, I believe, only had kids data. 

---

## 2026-09-23 — grouping thresholds and rare_NVT

In `cleanmetadata.R`: a serotype-GPSC pair keeps its own group if it has at least 20 isolates
(min_gpsc_n), otherwise it goes to `<serotype>_other`. NVT serotypes with less than 40 isolates total
(`min_serotype_n`) are grouped into one `rare_NVT` group. Vaccine serotypes never pool. 
This means that we have 87 groups, where previously we had 141 or 142 or something. 

Obviously the smaller number of groups is an advantage. Also,  44% of group-years were
zero. Most of those zeros are sampling zeros, and probablby not real absences in
carriage. With
100-400 IPD cases a year spread over ~140 groups, a group that's carried but not very
invasive easily produces none in a given year. This new grouping doesn't eliminate that
but reduces it a big and means that a single IPD case can no longer dominate too much. 


Vaccine types never pool and every group belongs to exactly one of PCV7/PCV13/NVT so
vaccine status can always be recovered. Only 7F (with 29 isolates) would otherwise have
been pooled. 

Also changed: gpsc keeps its real value now. The old code overwrote rare GPSCs with
"other", which threw away the lineage identity for any later analysis. Only `group` gets
coarsened. 

A sanity check: VT overall portion drops after 2009 and again after 2011 as expected. NVTs go up in
frequency from about 0.3 to 0.85 pre- to post-vaccine. 

Open: about 25% of group-years are still zero, mostly moderately common groups missing in
single early years. No grouping seems to fix that, really. I think we should just set the
model up with pre-vaccine equilibrium frequencies and pre-vaccine overall population
portions, and then try to model only what went up vs down (a lot or a little,
maybe). Might not be worth trying to fit the whole time course given longer IPD sampling
pre-vaccine which we have no driveres for. 

Idea: Longer term, better grouping might
come from dimension reduction on the gene matrix, defining groups and locus sets that best
describe the dynamics rather than serotype - GPSC combinations. Would be a whole new
thing. 

---

## 2026-09-23 — invasiveness: coverage and shrinkage

`make_invasiveness_lookup.R` makes one invasiveness value per serotype per age group,
covering every serotype in meta_clean, with the estimates shrunk toward the overall mean
on the log scale.

Previously, a line `filter(!is.na(v))` was silently dropping every serotype that isn't in
Nick's tables. The big one was 15B/15C with 131 IPD cases and 142 carriage samples; 
 Nick has 15B/C. Also 11BC, 20A, 35D, 33D, 6G, 23, 24, 25F, and the
adult table was additionally missing 10B, 18F, 19B, 39, 46. The name fixes are  in
`cleanmetadata.R`, 
Invasiveness values for serotypes that are still  absent from Nick's tables fall back to
the geometric mean of the serogroup, then the other age table, then the overall mean,
and `v_source` notes which. Non-typeables get NA and only come in through carriage
anyway. 

Invasiveness estimates: Nick's child estimates span 3745-fold, and the lowest ones (19B at 1.2e-5,
35A, 33A, 10F) all have lower CI bounds around 1e-6, so there is  almost no IPD from
them. When we divide by them to get the implied carriage counts we can end up with very
noisy results. In 6 of 15 years the top group in the child IPD-derived estimate was due to
single IPD case, which seems a little wrong.. 

I used empirical Bayes on the log scale (DerSimonian-Laird random effects, SEs from the reported
95% intervals in Nick's tables) to move each estimate toward the overall mean in proportion to how uncertain
it is. (That sounds hard but it's just $\log(v_{adj})  = \mu + B_i (\log(v_{orig}) - \mu)$. When $B = 0$
the new invasiveness is the mean (here $\mu$ is the log of the raw means). When $B=1$ the
original (raw) invasiveness estimate is kept. $B$ comes from the uncertainty in the
invasiveness estimates. Where there is higher certainty, invasiveness estimates are not much affected  (14, 6A, 19A change by less than 10%), but the near-zero ones go up about
10-fold, and the range is much more reasonable (though still a factor of 180). After this, no year's top group
comes from a single case. 

Open: if Nick's estimates came out of a hierarchical model themselves then this shrinks
them twice and over-smooths..  should check his methods but .. probably even better to include his post-vaccine invasiveness numbers?
Invasiveness is per serotype so every GPSC within a serotype gets the same inv value. I believe Nick also has post-vaccine
invasiveness estimates, which we could apply to the post-vaccine years (would that make
the locus frequencies stabilize more? Perhaps! Maybe this is a priority). 

South Africa's HIV prevalence over this period could affect the invasiveness. 
And uncertainty in $v$ isn't propagated anywhere, and it is probably the largest contributor to uncertainty
in $x(t)$. We could draw new vs from Nick's intervals, thereby resampling $x(t)$, and
then.. I guess we would potentially have a new NFDS model each time? Anyway I did put
standard errors in the new reweighted $x(t)$ file (but they aren't based on the unc in $v$) so that's something we could think about
(don't penalize the model for a mismatch if the standard error is large, for example) . 
