# Understanding scVI and Variational Autoencoders for scRNA-seq

## 1. What is scVI?

**scVI** (single-cell Variational Inference) is a deep generative model for single-cell RNA-seq data.  
It learns a low-dimensional latent space that denoises the data, corrects batch effects, and models gene expression counts using a probabilistic framework (VAE - Variational Autoencoder).

---

## 2. Challenges Addressed by scVI

- High technical noise and dropout in scRNA-seq  
- Batch effects across samples or conditions  
- Nonlinear gene expression relationships  
- Dimensionality reduction for clustering, DE, and visualization  

---

## 3. High-Level VAE Architecture

- **Encoder**: Maps gene expression to a latent representation `z`  
- **Latent space `z`**: Represents cell state, denoised and compressed  
- **Decoder**: Reconstructs gene expression from `z`  
- **KL Divergence**: Regularizes `z` to stay close to a prior  

---

## 4. Mathematical Core

### Posterior Approximation

Instead of directly computing the intractable posterior:

\[
p(z | x) = \frac{p(x | z) p(z)}{p(x)}
\]

We approximate it using:

\[
q(z | x, s) = \mathcal{N}(\mu(x, s), \text{diag}(\sigma^2(x, s)))
\]

### Encoder Output

- \(\mu(x, s)\): Mean vector of latent space  
- \(\sigma^2(x, s)\): Variance vector (diagonal covariance)  

### Decoder Model

\[
x \sim \text{NB}(\mu = f_\theta(z, s), \theta)
\]

Negative Binomial is used to capture overdispersion and dropout.

### ELBO (Loss Function)

\[
\text{ELBO} = \mathbb{E}_{q(z|x)}[\log p(x|z)] - D_{\text{KL}}(q(z|x) \parallel p(z))
\]

Encourages both reconstruction and a regularized latent space.

---

## 5. Key Tricks

### Reparameterization Trick

\[
z = \mu + \sigma \cdot \epsilon, \quad \epsilon \sim \mathcal{N}(0, I)
\]

Allows gradient-based optimization despite stochastic sampling.

**Where the stochastic sampling occurs**:  
During training, the encoder outputs parameters of a distribution (not a point).  
To train the model, we need to sample \(z \sim q(z|x)\) for each cell — this is the stochastic step.  
To keep gradients flowing during this sampling, we use the reparameterization trick.

### Diagonal Gaussian

- Simplifies the model  
- Assumes no correlation between latent dimensions  
- Enables fast KL computation  

---

## 6. Posterior Collapse

- **What**: Encoder learns \(q(z|x) \approx p(z)\), ignoring `x`  
- **Why**: Powerful decoder reconstructs `x` without needing latent `z`  
- **Fixes**: KL annealing, limited decoder capacity, dropout regularization  

---

## 7. KL Divergence Role

- Penalizes divergence of \(q(z|x)\) from \(p(z)\)  
- Keeps latent space compact and smooth  
- Chosen because it has an analytic form for Gaussians and fits naturally into the variational inference framework  

---

## 8. Visual Intuition

- **Prior**: Standard Gaussian (\(\mathcal{N}(0, I)\)) in latent space  
- **Posterior**: Encoder’s Gaussian centered at \(\mu(x)\)  
- **KL**: Distance between these distributions  

Two diagrams were created:

1. A normal posterior distribution offset from the prior  
2. A collapsed posterior nearly matching the prior (posterior collapse)  

---

## 9. Summary Table

| Concept                | Role                                               |
|------------------------|----------------------------------------------------|
| **Latent `z`**         | Denoised, compressed representation of cells       |
| **Encoder**            | Maps gene counts to a distribution over `z`        |
| **Decoder**            | Predicts gene expression from `z` and batch info   |
| **Reparameterization** | Enables sampling + gradient flow                   |
| **KL Divergence**      | Keeps latent space regularized                     |
| **Posterior collapse** | Model fails to learn useful encodings              |
