# Speaking notes — `regression_showcase.pdf`

Purpose of the meeting: **come out with decisions**, not with approval. Every number below is in the deck
or in `Amoc/results/`; nothing here is an estimate from memory.

If the meeting is cut to ten minutes, show slides **2, 6, 7, 12, 13** and skip the rest.

---

## 1 — Title (20 seconds)

> Two response functions: how crop yields and household energy use react to weather. Today is only the
> econometrics — no climate scenario, no impact figures. I have four or five choices I would like your
> view on, and I have put them at the end.

---

## 2 — What is estimated (1 minute)

> The estimand is a weather-response function per crop and per fuel: the proportional change in yield, or
> in consumption, when a season comes out different from what that region normally gets. There is no
> treatment and no control group — the identifying assumption is that weather is assigned by nature, so
> a region's yield cannot cause its own rainfall deviation.
>
> Crops: NUTS3 region by crop by harvest year, seven crops, 1989 to 2023, between 7,000 and 18,500
> observations per crop. Energy: country by year, two fuels, 1990 to 2024.

Do not linger. The professor will want to know the identification, which is the last bullet: within
region, year-to-year deviations from that region's own path, net of anything common to Europe that year.

---

## 3 — Roadmap and what I need (30 seconds)

> The path is eleven more slides. What I actually need from you is on the last two: six open questions.
> The first is the one I would spend the meeting on — whether the region-specific trends are the right
> price to pay, given some regions have very short series.

Say this out loud. It converts the meeting from a presentation into a consultation, and it tells him
where to interrupt.

---

## 4 — Every choice, and why (1–2 minutes)

Do **not** read eight rows. Read three and say the rest are there for reference:

> Three examples. Region rather than country, because country averaging destroys the within-country
> weather variation — it is 45 times the observations. Fixed effects rather than random, because soil
> and irrigation are correlated with climate, and Hausman rejects at 10 to the minus 22. Country
> clustering rather than region, because the spatial correlation survives a year effect.

---

## 5 — Specification F (1 minute)

> One regression per crop. The year is split in two blocks: overwinter, from sowing to February, and the
> season, from March to harvest — November for maize, September for beet and sunflower, August for the
> cereals. Spring-sown crops have no overwinter block: bare soil, nothing to damage.
>
> Region fixed effect, a region-specific quadratic trend, a year fixed effect. So the coefficient is
> identified from a region's deviation from its own smooth path, net of the European year.

The reason for two blocks, if asked: half the sample is autumn-sown and 78% of frost days fall outside
March–July, so a fixed March–July window throws away the sowing and overwinter period.

---

## 6 — The results (2–3 minutes, the core of the talk)

Read **one cell** out loud, then the pattern:

> Winter barley, season GDD: minus 6.5 times ten to the minus four, with a clustered standard error of
> 2.0. In round units, on the row below: a hundred extra growing degree-days over the season costs
> 6.5% of the yield.
>
> Now read across. The four cool-season cereals are negative on season warmth, minus 3 to minus 6.5 per
> hundred degree-days. The three summer crops are not — for them the damage is extreme heat instead:
> maize minus 1.8% per ten degree-days above 28.

Then the overwinter block:

> What survives overwinter is frost: winter barley at 1.5 times ten to the minus three, bootstrap p
> 0.000, soft wheat at p 0.075. That is winterkill, and a fixed March–July window cannot see it.

**Volunteer the weak point before he finds it:**

> Durum's overwinter heat coefficient is minus 0.12 with a standard error of 0.11. That block is almost
> always zero for durum, so it is fitted on nearly nothing. I do not report it as a finding.

---

## 7 — The two crop families (1 minute)

> This split is not imposed anywhere: the same nine regressors and the same fixed effects, fitted crop by
> crop. The agronomy comes out of the coefficients. That is my main evidence that the design is
> measuring something real rather than fitting noise.

If he asks about maize's positive GDD: say it is small, +1.2% per hundred degree-days, and that it does
**not** survive a change of trend specification — it is on slide 13. Do not defend it.

---

## 8 — Energy specification (45 seconds)

> Same design, country level. Heating degree-days below 15, cooling degree-days above 24, population
> weighted. Log household price as a control. Country fixed effect and trend, year effect, clustered by
> country.

---

## 9 — Energy results (1–2 minutes)

> Columns one to three are the build-up; column three is the benchmark. Columns four and five drop the
> price.
>
> The point of those two: price is the one regressor not assigned by nature. On the same rows, dropping
> it moves heating degree-days from 6.52 to 6.56 times ten to the minus five — nothing. So the weather
> response does not rest on the one endogenous variable.
>
> The last column adds back the early 1990s and Finland and Norway, which have no price series. There the
> response is weaker. That is a statement about the early data, not about the price.

---

## 10 — Inference (1 minute)

> Standard errors clustered by country, because regions of one country share weather and policy. But
> that leaves 7 to 15 clusters, where the clustered t over-rejects. So the p-values come from a wild
> cluster bootstrap, 99,999 replications, Webb weights when there are fewer than 12 clusters.

---

## 11 — Which standard error (1–2 minutes)

> Same coefficient in every row, only the variance changes. i.i.d. is six times too small. Clustering by
> region is barely better — it assumes a drought stops at a NUTS3 border. Clustering by country is two to
> three times the region error, and that gap is the spatial correlation.
>
> Look at the bottom two rows: for soft wheat the clustered t gives 0.063 and the bootstrap 0.027; for
> durum, 0.078 against 0.022. With this few clusters they disagree, so which one you quote is a choice —
> I quote the bootstrap.

---

## 12 — The three things that could break it (2 minutes)

Say these before he does. It is the slide that buys credibility.

> One: RESET rejects linearity for five of seven crops. That is the only violation that touches the
> coefficients rather than the standard errors. Heteroskedasticity, serial correlation and
> cross-sectional dependence also fail, but those only invalidate the classical standard error.
>
> Two, and this is the real limit: only 2.4 to 3.8% of season-GDD variance survives the fixed effects,
> and the counterfactual I eventually want to apply moves that regressor by up to ten standard
> deviations — far outside anything observed. The between-region slope, the only independent check,
> disagrees with the within slope.
>
> Three: the panel is unbalanced, 52 to 71% of cells. Most of the missing ones are before a region enters,
> not gaps. Slovak soft wheat contributes two years and is absorbed entirely by its own trend, so 15
> clusters are reported and 14 identify.

---

## 13 — Question 1: the trends (3–4 minutes, the heart of the meeting)

Present the evidence, then stop talking.

> The design spends three parameters per region — fixed effect, t, t squared — on a median of 24 years,
> and 341 of 3,769 region-crop units have fewer than ten years. So I tested it.
>
> The cereals do not care: quadratic, linear or no trend all sit inside one standard error. Dropping
> every region with fewer than ten years changes nothing.
>
> But country-times-year effects are a different story: soft wheat halves, durum goes to zero, sugar beet
> flips sign. And the small positive GDD of the summer crops survives none of it.
>
> **So my question: would you insist on country-by-year fixed effects — identifying only from differences
> between regions within the same country and year? It is the most demanding version, and under it two of
> my seven results largely disappear.**

Then be quiet and let him answer. Write down what he says.

---

## 14 — The other five questions (3–5 minutes)

Ask them one at a time, in order. Do not rush to fill silence.

1. **Extrapolation.** Estimated on ±1 SD, applied to −10 SD. Is there any defensible bound, or should the
   exercise be restricted to bins inside support and reported as a sign rather than a magnitude?
2. **Functional form.** RESET rejects for five of seven. Bins, spline, or note it and move on? Bins cost
   precision I may not have with 7–15 clusters.
3. **Seven regressions or one?** Durum and sunflower have seven country clusters each. Pooling with
   crop-interacted weather gives more clusters and a formal test of the cool/warm split, at the cost of
   common fixed effects.
4. **Only the intensive margin.** Yield per hectare, area weights held fixed. Does a credible impact
   number need the extensive margin?
5. **Energy dynamics.** ρ = 0.60 for electricity and 0.44 for gas survive trends and year effects. Partial
   adjustment, or is a static model with clustered errors answering my question?

**Closing line — the most valuable question of the meeting:**

> Of those six, which would you fix first, and which would you not bother with?

---

## Cheat sheet for likely challenges

| He says | You answer |
|---|---|
| "Is weather really exogenous?" | Year-to-year deviations of a region's weather from its own mean cannot be caused by that region's yield. The region fixed effect removes permanent confounders — soil, irrigation — and Hausman confirms they are correlated with climate, which is why FE and not RE. What remains untestable is same-year confounding, a pest outbreak caused by a warm spring. |
| "Why not random effects?" | Hausman rejects at p ≤ 4e-22 for every crop. For energy it does not reject, but θ = 0.96 there — RE has collapsed onto FE, so the test has nothing to discriminate; FE is kept on identification grounds. |
| "Your within R² is tiny." | 0.04–0.10, and it should be: most of the variance is fixed effects and trends. The question is whether the weather coefficient is identified, not whether weather explains yields. |
| "Only 3% of the variation survives." | Correct, and that is the extrapolation problem on slide 12, not a flaw in the estimate. |
| "What about adaptation?" | Burke–Emerick long difference is in the appendix. Every long-run coefficient is insignificant — 640 and 464 cells on 12 and 9 clusters. It cannot reject that long-run equals short-run, and cannot confirm it. I report it as uninformative rather than as evidence of no adaptation. |
| "Why NUTS3 and not country?" | 45× the observations and within-country identification. Country-level FAOSTAT was the original proposal and would have destroyed exactly the variation the weather pipeline was built to capture. |
| "Spec A versus F?" | A is a fixed March–July window, shown as a column. E is one window per crop and fits worse for all seven. F is reported because the overwinter frost result exists and A cannot see it. |
| "Sugar beet has a strange sample." | Ends 2022, and it is the least balanced panel, 52% of cells, 115 of 581 regions with fewer than ten years. Its coefficients are the ones I would trust least. |
