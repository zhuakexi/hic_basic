# HMM for Sister-Chromatid Composition Inference from SNP Read Counts

## Overview

This document describes a simple Hidden Markov Model (HMM) for inferring the local composition of two sister chromatids along a chromosome in secondary spermatocytes.

At each SNP position, the hidden state represents the combination of the two sister chromatids:

- `PP`: both chromatids carry haplotype `P`
- `PM`: one chromatid carries `P`, the other carries `M`
- `MM`: both chromatids carry `M`
- `U`: uncertain

Here `P` and `M` are **internal symmetric labels**. They do **not** need to be interpreted directly as paternal/maternal origin inside the model. The user is responsible for defining `ref` and `alt` consistently in the input data.

The observed data at each SNP are the read counts:

- `ref`: number of reads supporting the reference allele
- `alt`: number of reads supporting the alternative allele
- `other`: number of reads supporting neither ref nor alt

The goal is to infer the most likely hidden state path along the chromosome using an HMM and the Viterbi algorithm.

---

## 1. Hidden states

Let the hidden state at SNP position `t` be

$$s_t \in \{PP, PM, MM, U\}$$

with the following interpretation:

- `PP`: both sister chromatids are of haplotype `P`
- `PM`: one sister chromatid is `P`, the other is `M`
- `MM`: both sister chromatids are of haplotype `M`
- `U`: the local state cannot be reliably determined

The states `PP`, `PM`, and `MM` are the biological states of interest.  
The state `U` is introduced to absorb low-quality, conflicting, or sparse observations.

---

## 2. Observations

At SNP position `t`, the observation is the read-count vector

$$o_t = (r_t, a_t, x_t)$$

where:

- `r_t`: ref read count
- `a_t`: alt read count
- `x_t`: other read count

and the total depth is

$$n_t = r_t + a_t + x_t$$

These counts are treated as direct quantitative observations in the emission model.

---

## 3. HMM structure

The model consists of:

1. an initial state distribution \(\pi\)
2. a transition matrix \(A\)
3. an emission model \(B\)

The full parameter set is

$$\lambda = (\pi, A, B)$$

The joint probability of a hidden-state path $S = (s_1, \dots, s_T)$ and observations $O = (o_1, \dots, o_T)$ is

$$P(S,O \mid \lambda) = \pi_{s_1} \, b_{s_1}(o_1) \prod_{t=2}^{T} a_{s_{t-1},s_t} \, b_{s_t}(o_t)$$

where:

- $\pi_i = P(s_1 = i)$
- $a_{ij} = P(s_t = j \mid s_{t-1} = i)$
- $b_j(o_t) = P(o_t \mid s_t = j)$

---

## 4. Initial state distribution

The initial state distribution is

$$\pi = (\pi_{PP}, \pi_{PM}, \pi_{MM}, \pi_U)$$

with

$$\pi_{PP} + \pi_{PM} + \pi_{MM} + \pi_U = 1$$

A simple default choice is to assign a very small probability to `U` and distribute most mass among the biological states.

Example:

\[
\pi = (0.3, 0.4, 0.3, 0.0)
\]

If prior biological knowledge suggests that one state is more common at chromosome starts, the initial distribution can be adjusted accordingly.

---

## 5. Transition matrix

The transition matrix defines how hidden states change between adjacent SNPs:

$$A = (a_{ij}), \quad a_{ij} = P(s_t=j \mid s_{t-1}=i)$$

A simple and biologically intuitive transition structure is:

$$
A = \begin{bmatrix}
1-\alpha-\gamma & \alpha & 0 & \gamma \\
\beta/2 & 1-\beta-\gamma & \beta/2 & \gamma \\
0 & \alpha & 1-\alpha-\gamma & \gamma \\
\rho & \rho & \rho & 1-3\rho
\end{bmatrix}
$$

Interpretation:

- `PP -> PM` and `MM -> PM` are allowed with probability controlled by $\alpha$
- `PM -> PP` and `PM -> MM` are allowed with probability controlled by $\beta$
- entry into `U` is controlled by $\gamma$
- recovery from `U` is controlled by $\rho$
- direct `PP -> MM` and `MM -> PP` transitions are set to 0 in the simplest model

This structure encodes the idea that:

- states usually persist over neighboring SNPs
- changes are relatively rare
- `U` should appear only occasionally
- a direct jump from `PP` to `MM` is less plausible than passing through `PM`

If needed, direct `PP <-> MM` transitions can later be assigned a small nonzero probability.

---

## 6. Emission model

### 6.1 Principle

The emission model describes the probability of observing the read-count vector

$$(r_t,a_t,x_t)$$

given a hidden state.

Because the observations are counts in three categories (`ref`, `alt`, `other`), the most natural first model is a multinomial distribution.

---

### 6.2 Multinomial emission

For each hidden state `s`, define a probability vector

$$\theta_s = (p^{(s)}_{\text{ref}}, p^{(s)}_{\text{alt}}, p^{(s)}_{\text{other}})$$

with

$$p^{(s)}_{\text{ref}} + p^{(s)}_{\text{alt}} + p^{(s)}_{\text{other}} = 1$$

Then the observation at position `t` is modeled as

$$(r_t,a_t,x_t) \mid s_t=s \sim \text{Multinomial}(n_t, \theta_s)$$

and therefore

$$b_s(o_t) = P(o_t \mid s_t=s) = \frac{n_t!}{r_t!a_t!x_t!} (p^{(s)}_{\text{ref}})^{r_t} (p^{(s)}_{\text{alt}})^{a_t} (p^{(s)}_{\text{other}})^{x_t}$$

---

### 6.3 Intuitive emission patterns

A simple intuitive choice is:

- `PP`: mostly ref reads
- `PM`: roughly balanced ref and alt
- `MM`: mostly alt reads
- `U`: enriched in `other`, or broadly uninformative

Example parameterization:

$$\theta_{PP} = (0.95, 0.03, 0.02)$$

$$\theta_{PM} = (0.49, 0.49, 0.02)$$

$$\theta_{MM} = (0.03, 0.95, 0.02)$$

$$\theta_U = (0.25, 0.25, 0.50)$$

Interpretation:

- in `PP`, ref dominates
- in `MM`, alt dominates
- in `PM`, ref and alt are both expected
- in `U`, the signal is weak, mixed, or enriched for `other`

This parameterization assumes the user has defined `ref` and `alt` consistently with the intended orientation of `P` and `M`.

---

### 6.4 Why use count emissions directly

Using count emissions is preferable to first converting each SNP into a discrete label like `PP-like` or `PM-like`, because the count model retains:

- total depth information
- allele balance information
- the extent of `other` signal
- confidence differences between SNPs

For example, the two observations below should not be treated equally:

- `(10, 9, 0)`
- `(2, 1, 9)`

The first strongly supports `PM`, while the second is much more compatible with `U`.

---

### 6.5 Optional refinement: Dirichlet-multinomial

If the data show more variability than a multinomial model can explain, the emission model can be upgraded to a Dirichlet-multinomial:

$$(r_t,a_t,x_t)\mid s_t=s \sim \text{Dirichlet-Multinomial}(n_t,\alpha_s)$$

This allows extra dispersion and is often more realistic for sequencing data, but the multinomial version is the simplest starting point.

---

## 7. What probabilities need to be estimated

The model requires the following probabilities or parameters.

### 7.1 Initial state probabilities

$$\pi_{PP}, \pi_{PM}, \pi_{MM}, \pi_U$$

These describe the prior probabilities at the first SNP.

---

### 7.2 Transition probabilities

The transition matrix requires estimating:

- $\alpha$: transition probability between a pure state and the mixed state
- $\beta$: probability of leaving `PM`
- $\gamma$: probability of entering `U`
- $\rho$: probability of leaving `U`

Equivalent parameterizations are also possible; the key is to estimate the relative tendency of:

- state persistence
- biological state switching
- entering and leaving uncertain regions

---

### 7.3 Emission probabilities

For the multinomial model, each state requires estimating a 3-component probability vector:

$$\theta_{PP}, \theta_{PM}, \theta_{MM}, \theta_U$$

For example:

$$\theta_{PP} = (p^{(PP)}_{\text{ref}}, p^{(PP)}_{\text{alt}}, p^{(PP)}_{\text{other}})$$

and similarly for the other states.

These parameters describe the expected read-count composition under each hidden state.

---

## 8. How to estimate the probabilities

This section describes simple and practical ways to estimate the required parameters.

### 8.1 Estimating initial probabilities

The initial probabilities can be set in one of three ways:

#### Option 1: fixed defaults

Use biologically reasonable defaults, for example:

$$\pi = (0.3, 0.4, 0.3, 0.0)$$

This is the simplest option.

#### Option 2: empirical frequency from rough classification

Perform a rough preliminary SNP classification and use the state frequencies near chromosome starts.

#### Option 3: treat as non-critical

Because the initial position usually has little influence on a long chromosome, \(\pi\) can often be fixed to simple values and left unoptimized.

---

### 8.2 Estimating transition probabilities

A practical approach is to estimate transitions from a rough first-pass segmentation.

#### Step 1: local rough classification

For each SNP or local window, compute summary statistics such as:

- ref fraction:
  $$f_t = \frac{r_t}{r_t+a_t+x_t}$$
- alt fraction
- other fraction:
  $$o_t = \frac{x_t}{r_t+a_t+x_t}$$

Then assign rough labels using thresholds, for example:

- `PP-like` if ref fraction is very high
- `MM-like` if alt fraction is very high
- `PM-like` if ref and alt are both substantial and balanced
- `U-like` if depth is too low or `other` is too high

This preliminary labeling is only used to estimate parameters.

#### Step 2: count empirical transitions

From the rough labels, count how often each transition occurs between adjacent SNPs or between adjacent windows.

For example:

- number of `PP -> PM`
- number of `PM -> MM`
- number of `PP -> U`
- number of `U -> PM`

Then normalize row-wise to obtain approximate transition probabilities.

#### Step 3: smooth the estimates

To avoid zeros caused by small sample size, add pseudocounts before normalization.

This gives a stable initial transition matrix for the HMM.

---

### 8.3 Estimating emission probabilities

The emission probabilities can be estimated from curated SNP subsets or from preliminary state assignments.

#### Option 1: estimate from high-confidence subsets

Construct high-confidence SNP sets:

- `PP`-like: sites with very high ref fraction and very low `other`
- `MM`-like: sites with very high alt fraction and very low `other`
- `PM`-like: sites with near-balanced ref and alt and low `other`
- `U`-like: sites with low depth or high `other`

Then pool counts within each set and estimate:

$$\hat{\theta}_s = \left( \frac{\sum r}{\sum n}, \frac{\sum a}{\sum n}, \frac{\sum x}{\sum n} \right)$$

for each state `s`.

This is intuitive and easy to implement.

#### Option 2: initialize manually

Use fixed values such as:

$$\theta_{PP} = (0.95,0.03,0.02)$$

$$\theta_{PM} = (0.49,0.49,0.02)$$

$$\theta_{MM} = (0.03,0.95,0.02)$$

$$\theta_U = (0.25,0.25,0.50)$$

This is suitable for an initial proof of concept.

#### Option 3: iterative refinement

After one round of Viterbi decoding, collect SNPs assigned to each state and re-estimate the emission parameters from those assignments. Repeat if necessary.

This is a simple practical approximation to iterative training.

---

### 8.4 Estimating `U`-related probabilities

The uncertain state is special and should usually be sparse.

A practical way to estimate `U` parameters is:

- define high-confidence `U` candidates using low depth, high `other`, or highly conflicting counts
- estimate `\theta_U` from those candidates
- estimate `\gamma` and `\rho` from the lengths and frequencies of these rough uncertain segments

Usually:

- `\gamma` should be small
- `\rho` should be larger than `\gamma` unless `U` regions are expected to be long

---

### 8.5 Recommended estimation workflow

A simple and robust workflow is:

1. compute rough SNP- or window-level labels using thresholds
2. estimate emission parameters from high-confidence subsets
3. estimate transition probabilities from rough label transitions
4. construct the HMM
5. run Viterbi decoding
6. optionally re-estimate parameters from the decoded path and rerun

This gives an interpretable and practical estimation procedure without needing full EM training at the start.

---

## 9. Dynamic programming and Viterbi algorithm

The goal is to infer the most likely hidden-state sequence

$$\hat{S} = (\hat{s}_1,\hat{s}_2,\dots,\hat{s}_T)$$

given the observations.

The Viterbi algorithm computes the single most likely hidden path.

---

### 9.1 Viterbi variable

Define

$$\delta_t(j) = \max_{s_1,\dots,s_{t-1}} P(s_1,\dots,s_{t-1}, s_t=j, o_1,\dots,o_t \mid \lambda)$$

This is the probability of the most likely path ending in state `j` at position `t`, together with the first `t` observations.

---

### 9.2 Initialization

At the first SNP:

$$\delta_1(j) = \pi_j \, b_j(o_1)$$

$$\psi_1(j) = 0$$

where $\psi_t(j)$ stores the best previous state for traceback.

---

### 9.3 Recursion

For `t = 2, 3, ..., T`:

$$\delta_t(j) = \max_i \left[ \delta_{t-1}(i)\, a_{ij} \right] \, b_j(o_t)$$

$$\psi_t(j) = \arg\max_i \left[ \delta_{t-1}(i)\, a_{ij} \right]$$

This step compares all possible previous states and keeps only the best one for each current state.

---

### 9.4 Termination

At the last SNP:

$$P^* = \max_j \delta_T(j)$$

$$\hat{s}_T = \arg\max_j \delta_T(j)$$

---

### 9.5 Traceback

Recover the full path by following the back-pointers:

$$\hat{s}_{t-1} = \psi_t(\hat{s}_t) \quad \text{for } t=T,T-1,\dots,2$$

This produces the most likely hidden-state sequence.

---

## 10. What the algorithm is doing intuitively

At each SNP, the algorithm evaluates all hidden states and asks:

- how compatible is this state with the local read counts?
- how compatible is this state with the previously inferred state?
- among all ways of reaching this state, which previous path is best?

In this way, the HMM balances two factors:

1. **local evidence** from SNP read counts
2. **global smoothness** from state transitions

This avoids overreacting to isolated noisy SNPs.

For example:

- a single SNP with counts `(12,11,0)` in the middle of a long `PP` region may not be enough to force a real `PP -> PM -> PP` transition
- several consecutive SNPs with balanced counts may together support a true `PM` segment

---

## 11. Why this HMM improves over per-SNP independent calling

Independent per-SNP calling ignores continuity along the chromosome.  
The HMM improves inference by enforcing path consistency.

Advantages include:

- suppression of isolated noisy SNP calls
- smoother and more biologically plausible segments
- explicit handling of uncertain regions
- principled integration of read depth, allele balance, and continuity

For example:

- `(20,1,0)` strongly supports `PP`
- `(10,10,0)` strongly supports `PM`
- `(2,1,9)` is much better explained by `U`

The HMM combines these local signals with state persistence, producing a more robust chromosome-wide inference.

---

## 12. Practical notes

### 12.1 Use log probabilities

In implementation, probabilities should be converted to log scale to avoid numerical underflow:

$$\log \delta_t(j) = \max_i \left[ \log \delta_{t-1}(i) + \log a_{ij} \right] + \log b_j(o_t)$$

---

### 12.2 SNP spacing

In the simplest model, transition probabilities depend only on adjacent SNP order, not physical distance.  
A later refinement can make transitions depend on genomic distance, but this is not required in the first version.

---

### 12.3 State `U`

The `U` state is useful, but it should not dominate the path.  
If `U` is too easy to enter, the result becomes fragmented and uninformative.  
Therefore, entry into `U` should usually be penalized.

---

## 13. Summary

This model provides a simple and interpretable HMM for inferring sister-chromatid composition from SNP read counts.

### Hidden states

$$\{PP, PM, MM, U\}$$

### Observations

$$(r_t, a_t, x_t)$$

### Emission model

$$(r_t,a_t,x_t) \mid s_t=s \sim \text{Multinomial}(n_t,\theta_s)$$

### Inference

Use the Viterbi algorithm to infer the most likely chromosome-wide hidden-state path.

### Parameters to estimate

- initial probabilities $\pi$
- transition probabilities $A$
- emission probabilities $\theta_s$

### Recommended estimation strategy

1. rough preliminary labeling
2. estimate emissions from high-confidence subsets
3. estimate transitions from rough transition counts
4. run Viterbi decoding
5. optionally refine parameters and rerun

This model is intentionally simple and serves as a clean starting point for more advanced extensions.
