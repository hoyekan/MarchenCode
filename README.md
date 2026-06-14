# MarchenCode

This repo contains GPU-optimized Python code to compute the Green's functions from single-sided surface reflectivity data. **marchenko_gpu_optimized.ipynb** computes the **Marchenko Green's functions** and **focusing functions** from the **input reflectivity**, **direct arrival**, **filter (theta)**, and **ground-truth Green’s functions** (for comparison). **marchenko_imaging_gpu_optimized.ipynb** implements **Marchenko-based imaging** using the **decomposed Green’s functions**, specifically the **upgoing Green’s function** $g^-$ and the **downgoing direct arrival** (approximated as $f_0^+$). I embarked on this project in an effort to understand the Marchenko methods. For learners proficient in MATLAB, you can download the original version of the code in MATLAB [here](https://wiki.seg.org/wiki/Software:Marchenko_for_imaging) by [`Angus Lomas and Andrew Curtis`](https://www.geos.ed.ac.uk/~acurtis/assets/Lomas_Curtis_Geop_March_2019.pdf)

---

## Green's Functions Estimation Using Marchenko Methods - GPU-optimized CuPy code

### Marchenko Method Formulas

The Marchenko method reconstructs the **Green's functions** (both upgoing and downgoing) between a virtual source and surface receivers using single-sided reflection data and an estimate of the direct arrival.

The Marchenko method involves several key equations to estimate Green's functions from seismic reflection data. Here are the formulas corresponding to each section of the code below:

---

### **1. Initial Focusing Functions (Section 4)**

#### a. Downgoing focusing function $f_0^+$ :

$$
f_0^+(\mathbf{x}_f, \mathbf{x}_s, t) = T_d(\mathbf{x}_r, \mathbf{x}_f, -t) \tag{1}
$$

- Where:
  - $T_d$ = Direct arrival wavefield (one-way travel time) from the subsurface focusing point (virtual source) to the surface receivers
  - $T_d(\mathbf{x}_r, \mathbf{x}_f, -t)$ is read as the time-reversed direct arrival at the receiver location $\mathbf{x}_r$ from the virtual source $\mathbf{x}_f$ in the subsurface. Other notations in this note can be read in a similar way.
  - $\mathbf{x}_r$ = Receiver position
  - $\mathbf{x}_f$ = Focal point position
- Implementation: Time reversal of direct arrival (Time-reversed direct arrival)


#### b. Upgoing focusing function $f_0^-$ :

$$
f_0^-(\mathbf{x}_r, \mathbf{x}_f, t) = \Theta(t) \ast \int_S R(\mathbf{x}_r, \mathbf{x}_s, t) \ast f_0^+(\mathbf{x}_f, \mathbf{x}_s, t)  d\mathbf{x}_s
$$

- Where:
  - $R$ = Reflection response
  - $\Theta$ = Time-window function. It is a muting filter that removes acausal or undesired energy (e.g., above the direct arrival time).
  - $\ast$ = Convolution operator
- Implementation: Convolution of reflection response with $f_0^+$. In frequency domain, this becomes a matrix product.

---

### **2. Iterative Focusing Functions (Section 5)**

For k = 1 to N iterations:
- Downgoing component:
$$
m_{k}^{+}(\mathbf{x}_f, \mathbf{x}_s, t) = \Theta(t) \ast \int_S R(\mathbf{x}_r, \mathbf{x}_s, t) \ast f_{k-1}^{-}(\mathbf{x}_r, \mathbf{x}_f, -t) d\mathbf{x}_r
$$

For the first iteration, for instance,
$$
m_{1}^{+}(\mathbf{x}_f, \mathbf{x}_s, t) = \Theta(t) \ast \int_S R(\mathbf{x}_r, \mathbf{x}_s, t) \ast f_{0}^{-}(\mathbf{x}_r, \mathbf{x}_f, -t) d\mathbf{x}_r
$$

* Where: $\mathbf{x}_f$: focusing point , $\mathbf{x}_s$: source location (or source side in reflection response), and $\mathbf{x}_r$: receiver location (or receiver side in reflection response)

Implementation: Convolve reflectivity and (time-reversed) fk_minus

---

- Updating upgoing function:
  $$
  f_{k}^{-}(\mathbf{x}_r, \mathbf{x}_f, t) = f_{0}^{-}(\mathbf{x}_r, \mathbf{x}_f, t) + \Theta(t) \ast \int_S R(\mathbf{x}_r, \mathbf{x}_s, t) \ast m_{k}^{+}(\mathbf{x}_f, \mathbf{x}_s, -t) d\mathbf{x}_s
  $$

    $$
    f_k^-(\mathbf{x}_r, \mathbf{x}_f, t) = \Theta(t) \ast \left[f_0^-(\mathbf{x}_r, \mathbf{x}_f, t) + \int R(\mathbf{x}_r, \mathbf{x}_s, t) \ast m_k^+(\mathbf{x}_f, \mathbf{x}_s, -t) \, d\mathbf{x}_s \right]  \tag{4}
    $$
    
- Downgoing function update:
  
$$
f_{k}^{+}(\mathbf{x}_f, \mathbf{x}_s, t) = f_{0}^{+}(\mathbf{x}_f, \mathbf{x}_s, t) + m_{k}^{+}(\mathbf{x}_f, \mathbf{x}_s, t)
$$

$$
f_k^+( \mathbf{x}_f, \mathbf{x}_s, t) = f_0^+( \mathbf{x}_f, \mathbf{x}_s, t) + \Theta(t) \ast \int_S R( \mathbf{x}_r, \mathbf{x}_s, t) \ast f_{k-1}^-( \mathbf{x}_r, \mathbf{x}_f, -t)  d\mathbf{x}_r
$$

$$
f_k^+( \mathbf{x}_f, \mathbf{x}_s, t) = f_0^+( \mathbf{x}_f, \mathbf{x}_s, t) + \Theta(t) \ast \int_S R( \mathbf{x}_r, \mathbf{x}_s, t) \ast f_{k-1}^-(\mathbf{x}_r, \mathbf{x}_f, -t)  d\mathbf{x}_r
$$

- Where $k$ = Iteration index

---

### **3. Green's Functions (Section 6)**

#### a. Downgoing Green's function ($G^+$):

$$
G^+( \mathbf{x}_r, \mathbf{x}_f, t) = f_k^+( \mathbf{x}_r, \mathbf{x}_f, -t) - \int_S R( \mathbf{x}_r, \mathbf{x}_s, t) \ast f_k^-( \mathbf{x}_s, \mathbf{x}_f, -t)  d\mathbf{x}_s
$$

#### b. Upgoing Green's function ($G^-$):

$$
G^-( \mathbf{x}_r, \mathbf{x}_f, t) = -f_k^-( \mathbf{x}_r, \mathbf{x}_f, t) + \int_S R( \mathbf{x}_r, \mathbf{x}_s, t) \ast f_k^+( \mathbf{x}_s, \mathbf{x}_f, t)  d\mathbf{x}_s
$$

Note: In numerical implementations, the convolution is often performed in the frequency domain. which is a simple multiplication

#### c. Total Green's function:

$$
G_{\text{total}}(\mathbf{x}_r, \mathbf{x}_s, t) =  G^+(\mathbf{x}_f, \mathbf{x}_s, t)  +  G^-(\mathbf{x}_r, \mathbf{x}_f, t) \tag{8}
$$

$$
G_{\text{total}} = G^{-} + G^{+}
$$

The upgoing Green's Function can also be extracted by:

$$
G^-(t) = G(t) - G^+(t)
$$

Where:

* $G(t)$: Total Green’s function
* $G^+(t)$: Previously computed downgoing part

This is derived from the assumption that the total Green's function is a sum of its downgoing and upgoing parts.

---

### **4. Key Processing Steps**
For whatever reason, you may choose to design your own taper function using the function below (or a modified function of it) instead of the built-in `tukey` in the SciPy library.

1. **Tapering Window (Tukey)**:

    $$
   w_{\text{tukey}}(i) = \begin{cases}
   \frac{1}{2} \left[1 + \cos\left(\pi \frac{2i}{\alpha(N-1)} - 1\right)\right] & 0 \leq i \leq \frac{\alpha(N-1)}{2} \\
   1 & \frac{\alpha(N-1)}{2} < i < (N-1) - \frac{\alpha(N-1)}{2} \\
   \frac{1}{2} \left[1 + \cos\left(\pi \frac{2i}{\alpha(N-1)} - \frac{2}{\alpha} + 1\right)\right] & (N-1) - \frac{\alpha(N-1)}{2} \leq i \leq (N-1)
   \end{cases}  \tag{9}
   $$

   more compact

   $$
   w_{\text{tukey}}(i) = \begin{cases}
   \frac{1}{2} \left[1 + \cos\left(\pi \frac{2i}{\alpha(N-1)} - 1\right)\right] & 0 \leq i \leq \frac{\alpha(N-1)}{2} \\
   1 & \frac{\alpha(N-1)}{2} < i < (N-1)(1-\frac{\alpha}{2}) \\
   \frac{1}{2} \left[1 + \cos\left(\pi \frac{2i}{\alpha(N-1)} - \frac{2}{\alpha} + 1\right)\right] & (N-1)(1-\frac{\alpha}{2}) \leq i \leq (N-1)
   \end{cases}  \tag{10}
   $$

   - Where $\alpha$ = taper fraction (`tp`)
      
#### Comparison with Other Windows
| Window Type | Best For | Spectral Leakage |
|-------------|----------|------------------|
| Tukey       | General purpose with controllable tapering | Moderate |
| Hann        | General purpose | Good |
| Hamming     | Narrowband analysis | Fair |
| Rectangular | Transient detection | Poor |

The Tukey window provides a good compromise between frequency resolution and spectral leakage control, especially when you need to adjust the amount of tapering for your specific application. The Tukey function above is controlled by the `alpha` parameter. The `alpha` parameter controls the transition:
- `alpha=0`: Equivalent to a rectangular window (no tapering)
- `alpha=1`: Equivalent to a Hann window (full cosine taper)
       
2. **Pre-direct Arrival Mute**:

 $$
  G_{\text{final}} = G_{\text{total}} \cdot (1 - \Theta(t))
 $$
   - Removes artifacts before first arrival

---

### **5. Normalization**
$$
G_{\text{norm}} = \frac{G}{\max(|G_{\text{total}}|)}
$$

$$
G_{\text{norm}} = \frac{G}{\max_{t,\mathbf{x}_r}( |G_{\text{total}}|)}
$$

- Ensures consistent amplitude scaling
  
---

### **6. Reference True Green's Function (`true`)**

This is typically computed via direct **forward modeling** (e.g., finite difference, finite element, or Kirchhoff modeling) using the known subsurface model $c(x, z)$ (velocity field):

  $$
  \left( \frac{1}{c^2(x, z)} \frac{\partial^2}{\partial t^2} - \nabla^2 \right) G(x, z, t) = \delta(x - x_s)\delta(z - z_s)\delta(t)  \tag{11}
  $$

Where:

* $c(x, z)$: Velocity model
* $\delta$: Dirac delta function at the source location $(x_s, z_s)$

This equation ensures that `true` serves as the ground truth against which the retrieved Green’s functions are compared.

---

### **Variable Correspondence Table**

| Math Symbol | Code Variable | Description |
|------------|---------------|-------------|
| $R$        | `sg`          | Reflection response |
| $T_d$      | `direct`      | Direct arrival |
| $\Theta$   | `theta`       | Time window function (causality constraint) |
| $f_k^+$    | `fk_plus`     | Downgoing focusing function |
| $f_k^-$    | `fk_minus`    | Upgoing focusing function |
| $G^+$      | `g_plus`      | Downgoing Green's function |
| $G^-$      | `g_minus`     | Upgoing Green's function |
| $w$        | `tap`         | Taper window |
| $\alpha$   | `tp`          | Taper fraction |
| $N$        | `ns`          | Number of receivers |
| $\mathbf{x}_r$          | `-`     | Receiver position vector |
| $\mathbf{x}_f$        | `-`         | Focal point position vector (focusing point) |
| $\mathbf{x}_s$         | `-`          | Source position vector |
| $S$                   | `-`          | Measurement Surface |
| $\ast$        | `-`         | Convolution Operator |
| $k$         | `itr `          | Iteration index |
| $t$                   | `-`          | time |

These equations implement the Marchenko method to retrieve Green's functions from single-sided reflection measurements, enabling redatuming without detailed velocity model information. The iterative scheme progressively removes artifacts caused by internal multiples.

### **Summary**

* The **direct arrival** $T_d$ gives $f_0^+$.
* Using $f_0^+$ and the **reflection response** $R$, we compute $f_0^-$.
* We iteratively update $f_k^+$ and $f_k^-$ using $R$ and a time-reversed convolution process.
* The final Green's functions $G^-$ and $G^+$ are computed from these focusing functions and $R$, from which **total Green's function** is constructed.

In summary, Marchenko methods allow for the reconstruction of $G(t, \mathbf{x}_r, \mathbf{x}_f)$ using only:
- Surface measurements ($\mathbf{x}_s$ and $\mathbf{x}_r$ at the surface)
- An estimate of the direct arrival
- A time window (causality constraint)




