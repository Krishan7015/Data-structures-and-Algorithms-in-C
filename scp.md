# SEMICONDUCTOR PHYSICS - B.TECH 1ST Y# UNIT 1: CRYSTAL STRUCTURE AND BONDING

## Q1: HCP and Diamond Structure - Packing Fractions

### Hexagonal Close-Packed (HCP) Structure

**Structure Description:**
- HCP is a hexagonal lattice with two layers of atoms stacked in an ABAB pattern
- Each atom touches 12 nearest neighbors (coordination number = 12)
- Lattice parameters: basal plane spacing = a, height between layers = c
- For ideal HCP: c/a = √(8/3) ≈ 1.633

**Packing Fraction Derivation:**

In HCP unit cell:
- Number of atoms per unit cell = 6 × (1/6) [corner] + 3 × (1/2) [face] + 2 [inside] = 2 atoms

Volume of one atom = (4/3)πr³

Total volume of atoms in unit cell = 2 × (4/3)πr³ = (8/3)πr³

For hexagonal structure with ideal c/a ratio:
- Base area of hexagonal unit cell = (√3/2)a²
- Height = c = a√(8/3)
- Unit cell volume = (√3/2)a² × a√(8/3) = (√3/2)a³ × √(8/3)

From geometry, atoms touch along the basal plane: a = 2r

Substituting and simplifying:

**Packing Fraction of HCP = (π√3)/(6√2) ≈ 0.74 or 74%**

This is the highest packing fraction for any structure.

### Diamond Structure

**Structure Description:**
- Diamond is a cubic structure (face-centered cubic - FCC basis)
- Each carbon atom is tetrahedrally bonded to four other carbon atoms
- Coordination number = 4
- Lattice constant = a

**Packing Fraction Derivation:**

In diamond cubic unit cell:
- FCC atoms at corners: 8 × (1/8) = 1
- FCC atoms at face centers: 6 × (1/2) = 3
- Atoms inside unit cell at tetrahedral positions: 4
- **Total atoms per unit cell = 8**

Volume of one atom = (4/3)πr³

Total volume of atoms = 8 × (4/3)πr³ = (32/3)πr³

From tetrahedral bonding geometry:
- Distance between nearest neighbors (C-C bond) = a√3/4
- This equals 2r (atom diameter)
- Therefore: r = a√3/8

Unit cell volume = a³

Packing Fraction = [(32/3)πr³]/a³
= [(32/3)π(a√3/8)³]/a³
= [(32/3)π × (a³ × 3√3/512)]/a³
= [(32/3)π × (3√3/512)]

**Packing Fraction of Diamond ≈ 0.34 or 34%**

**Comparison Table:**
| Parameter | HCP | Diamond |
|-----------|-----|---------|
| Coordination Number | 12 | 4 |
| Packing Fraction | 74% | 34% |
| Structure Type | Close-packed | Open structure |
| Bonding | Metallic | Covalent |

---

## Q2: Miller Indices and Interplanar Spacing

### What are Miller Indices?

Miller indices (hkl) are a notation system to identify crystal planes and directions:

**Rules for determining Miller Indices:**
1. Find the intercepts of the plane with x, y, z axes in terms of lattice constants (a, b, c)
2. Take reciprocals of these intercepts
3. Clear fractions by multiplying by LCM
4. Enclose in parentheses: (hkl)

**Special Cases:**
- If plane is parallel to an axis, intercept = ∞, reciprocal = 0
- (100): plane perpendicular to x-axis, intercepts at (a, ∞, ∞)
- (111): plane intercepts all three axes equally
- (110): plane parallel to z-axis

### Interplanar Spacing

**Derivation of Interplanar Spacing Formula:**

For a simple cubic crystal with lattice constant 'a' and plane (hkl):

The distance between successive (hkl) planes is:

**d_{hkl} = a/√(h² + k² + l²)**

**General case for non-cubic systems:**

For orthorhombic crystal:
**d_{hkl} = 1/√[(h/a)² + (k/b)² + (l/c)²]**

### Calculation: Interplanar Spacing between (111) and (121)

**Given:** Cubic system with lattice constant a (assume Si: a = 5.43 Å)

**For (111) planes:**
d₁₁₁ = a/√(1² + 1² + 1²) = a/√3

**For (121) planes:**
d₁₂₁ = a/√(1² + 2² + 1²) = a/√6

**Distance between (111) and (121) planes:**

The reciprocal lattice vectors are:
- G₁₁₁ = 2π(1, 1, 1)/a
- G₁₂₁ = 2π(1, 2, 1)/a

Angle between planes is given by:
cos θ = |G₁₁₁ · G₁₂₁|/(|G₁₁₁||G₁₂₁|)

G₁₁₁ · G₁₂₁ = 2π(1·1 + 1·2 + 1·1)/a² = 4(2π/a²)

|G₁₁₁| = 2π√3/a
|G₁₂₁| = 2π√6/a

cos θ = 4(2π/a²)/[(2π√3/a)(2π√6/a)] = 4/(√18) = 4/(3√2) = 2√2/3

**Direct spacing between the two families:**

Using Miller indices properties, the perpendicular distance between one (111) and one (121) plane in their respective families is:

**Δd = |d₁₁₁ - d₁₂₁| = a(1/√3 - 1/√6) = a(√2 - 1)/√6**

Numerically (for Si, a = 5.43 Å):
- d₁₁₁ = 5.43/√3 = 3.14 Å
- d₁₂₁ = 5.43/√6 = 2.22 Å
- Δd ≈ 0.92 Å

---

## Q3: Point Defects and Frenkel Defect Concentration

### Types of Point Defects

**1. Vacancy (Schottky defect):**
- One atom is missing from its lattice position
- Creates an empty site
- Neutral point defect in pure crystals

**2. Interstitial:**
- Extra atom sits between lattice positions
- Causes local distortion
- Creates strain in crystal

**3. Frenkel Defect:**
- A vacancy-interstitial pair
- One atom leaves its lattice site and occupies an interstitial position
- Net result: one empty lattice site + one atom in interstitial
- Common in ionic crystals (AgCl, AgBr)

**4. Substitutional Impurity:**
- Foreign atom replaces an atom at a lattice site
- If smaller: causes compression; if larger: causes tension
- Affects electrical and mechanical properties

### Frenkel Defect Concentration Derivation

**Assumptions:**
- Crystal has N lattice sites
- N_i interstitial positions
- Temperature T, Boltzmann constant k

**Energy Formation:**
- Energy to create one Frenkel pair = W_F (typically 0.5-1 eV)
- This energy creates one vacant lattice site and one interstitial atom

**Using Statistical Mechanics:**

At thermal equilibrium, number of Frenkel defects n_F is given by:

**n_F = √(N × N_i) × exp(-W_F/2kT)**

**Derivation:**

The probability of forming a Frenkel pair:
- Probability that a lattice site is vacant: p_v
- Probability that an interstitial site is occupied: p_i

For uncorrelated formation, but constrained by geometry:

The number of ways to form n_F defects from N sites and N_i interstitial positions:

Ω = (N!/(n_F!(N-n_F)!)) × (N_i!/(n_F!(N_i-n_F)!))

Helmholtz free energy: F = -kT ln Ω + n_F × W_F

Minimizing F with respect to n_F:

∂F/∂n_F = 0

-kT[ln(n_F/(N-n_F)) + ln(n_F/(N_i-n_F))] + W_F = 0

For n_F << N and n_F << N_i:

kT × 2ln(n_F) + kT[ln(N(N_i))] = W_F

ln(n_F²/(N·N_i)) = -W_F/kT

**n_F = √(N·N_i) × exp(-W_F/2kT)**

**Important Points:**
- Concentration ∝ √N (not linear like Schottky)
- Decreases exponentially with formation energy
- For Ge/Si at room temp: typically 10⁻¹² to 10⁻¹⁵ cm⁻³
- Increases strongly with temperature

**Numerical Example:**
For a crystal with N = 10²² lattice sites, N_i = N, W_F = 0.9 eV at T = 300K:

n_F = √(10²² × 10²²) × exp(-0.9/(2 × 8.617×10⁻⁵ × 300))
= 10²² × exp(-17.45)
≈ 10²² × 3.2 × 10⁻⁸
≈ 3.2 × 10¹⁴ cm⁻³

---

## Q4: Symmetry Operations and 5-fold Rotation

### What are Symmetry Operations?

Symmetry operations are transformations that leave a crystal structure unchanged. Main operations:

**1. Rotation (n-fold rotational axis):**
Rotation by angle 2π/n brings crystal back to equivalent position

**2. Reflection (mirror plane):**
Reflection across a plane maps structure to itself

**3. Inversion (inversion center):**
Point reflection through a center

**4. Rotoreflection (rotoreflection axis):**
Rotation followed by reflection

**5. Translation (lattice translation):**
Shift by lattice vector preserves structure

### Proof: 5-fold Rotation Symmetry is NOT Possible

**Mathematical Proof:**

Consider a crystal lattice with lattice constant 'a'. If a 5-fold rotational axis exists, then rotating the crystal by θ = 2π/5 = 72° should map it to an equivalent position.

Let's place two lattice points at:
- Point A at origin (0, 0)
- Point B at (a, 0)

**After rotation by 72° around some axis perpendicular to the plane:**
Point B moves to B' at position (a cos 72°, a sin 72°)

For the crystal to have 5-fold symmetry, B' must coincide with a lattice point. The lattice points near B' have coordinates (m·a, n·a) where m, n are integers.

This requires:
- a cos 72° = m·a  →  cos 72° = m (m integer)
- a sin 72° = n·a  →  sin 72° = n (n integer)

**Calculating cos 72°:**
cos 72° = cos(2π/5) = (√5 - 1)/4 ≈ 0.309

This is NOT an integer.

**Contradiction:** cos 72° cannot equal any integer. Therefore, 5-fold rotation is impossible in a crystal lattice.

### Allowed Rotational Symmetries

For periodic lattice, only n-fold rotations with n = 1, 2, 3, 4, 6 are possible.

**Proof for these values:**

For lattice periodicity, the rotation matrix R(2π/n) must map lattice vectors to lattice vectors.

The trace of rotation matrix: Tr(R) = 1 + 2cos(2π/n)

For lattice periodicity: Tr(R) must be integer (from crystallographic restriction).

1 + 2cos(2π/n) = integer

cos(2π/n) = 0, ±1/2, ±1

- cos(2π/n) = 1  →  n = 1 (identity)
- cos(2π/n) = 1/2  →  n = 6 (60°)
- cos(2π/n) = 0  →  n = 4 (90°)
- cos(2π/n) = -1/2  →  n = 3 (120°)
- cos(2π/n) = -1  →  n = 2 (180°)

**Conclusion:** 5-fold (and any n ≥ 7) rotational symmetry is forbidden in crystals due to lattice periodicity constraint.

---

## Q5: Bonding in Cu, NaCl, and Si

### Cu: Metallic Bonding

**Structure:** Face-centered cubic (FCC), coordination number = 12

**Nature of Bonding:**
- Valence electrons (3d¹⁰4s¹) are delocalized
- Form a "sea of electrons" that moves freely throughout the structure
- Positive Cu⁺ cores are held together by attractive force from electron sea
- Non-directional bonding

**Characteristics:**
- **High electrical conductivity:** Free electrons carry current
- **High thermal conductivity:** Electrons transport heat
- **Ductile and malleable:** Layers can slide without breaking bonds
- **High melting point:** Strong electron-ion interactions (≈1085°C)
- **Metallic luster:** Electrons absorb and re-emit light

**Bonding Energy:** Typically 3-4 eV per atom

### NaCl: Ionic Bonding

**Structure:** Rock salt (FCC), each Na⁺ surrounded by 6 Cl⁻ and vice versa

**Nature of Bonding:**
- Complete electron transfer from Na to Cl
- Na loses 3s¹ electron → Na⁺ (noble gas config)
- Cl gains electron → Cl⁻ (noble gas config)
- Electrostatic attraction between oppositely charged ions
- Non-directional, long-range Coulombic force

**Ionic Bond Energy:**

Born-Landé equation (simplified):
**U = -N_A·M·z⁺·z⁻·e²/(4πε₀r₀) × (1 - 1/n)**

Where:
- M = Madelung constant (≈1.748 for NaCl)
- z⁺, z⁻ = ionic charges
- r₀ = nearest neighbor distance
- n = Born exponent (≈9 for NaCl)

**Characteristics:**
- **Brittle:** Layers shifted bring like charges together → repulsion → fracture
- **High melting point:** Strong ionic attraction (≈801°C)
- **Insulators in solid state:** No free electrons
- **Conductors in molten state:** Ions mobile
- **Soluble in polar solvents:** Ions solvated by molecules
- **Directional bonding:** Between specific ions

**Typical Bonding Energy:** 7-8 eV per ion pair

### Si: Covalent Bonding

**Structure:** Diamond cubic, each Si bonded to 4 nearest neighbors tetrahedrally

**Nature of Bonding:**
- Valence electrons (3s²3p²) are shared between atoms
- Each Si atom shares its 4 valence electrons with 4 neighbors
- Forms strong directional bonds
- Electron pair shared between two atoms

**Covalent Bond Structure:**

Each Si-Si bond involves:
- Overlap of sp³ hybrid orbitals
- Electron density concentrated between nuclei
- Strong attractive force

**Characteristics:**
- **High melting point:** Very strong directional bonds (≈1414°C)
- **Hard and brittle:** Rigid structure; electrons localized
- **Poor electrical conductivity:** Electrons bound to atoms
- **Semiconducting:** Small energy gap (1.1 eV) allows thermal/optical excitation
- **Non-directional in appearance:** But bonds are directional
- **Insoluble in polar solvents:** No ions to solvate

**Bonding Energy:** Typically 4-5 eV per bond

### Comparative Bonding Properties Table

| Property | Cu (Metallic) | NaCl (Ionic) | Si (Covalent) |
|----------|---------------|-------------|---------------|
| Electrical Conductivity | Very High | Low (solid) | Moderate (semi) |
| Solubility | Insoluble | Soluble (polar) | Insoluble |
| Mechanical | Ductile | Brittle | Brittle |
| Melting Point | 1085°C | 801°C | 1414°C |
| Bonding Type | Delocalized | Electrostatic | Localized pair |
| Electron Mobility | High | Low | Medium |

---

## Q6: Unit Cell Structures of Li and Cu

### Lithium (Li)

**Crystal Structure:** Body-Centered Cubic (BCC)

**Unit Cell Parameters:**
- Lattice constant: a = 3.50 Å (at room temperature)
- Atoms per unit cell: 1 (center) + 8 × (1/8) (corners) = 2 atoms

**Atomic Positions:**
- Corner positions: (0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,1,1)
- Body center: (1/2, 1/2, 1/2)

**Coordination Number:** 8 (each atom touches 8 nearest neighbors)

**Nearest Neighbor Distance:**
- Along body diagonal: √3·a/2 = 1.732 × 3.50/2 = 3.03 Å

**Packing Fraction:**
- Volume of 2 atoms = 2 × (4/3)πr³
- From nearest neighbor: 2r = √3·a/2  →  r = √3·a/4
- Packing Fraction = (2 × 4π/3 × (√3·a/4)³)/a³ = (π√3)/8 ≈ 0.68 or 68%

**Density Calculation:**
- Atomic mass of Li = 6.94 g/mol
- Number of atoms per unit cell = 2
- Mass per unit cell = (2 × 6.94)/(6.022 × 10²³) g = 2.304 × 10⁻²³ g
- Volume per unit cell = (3.50 × 10⁻⁸)³ = 4.29 × 10⁻²³ cm³
- **Density ρ = 2.304 × 10⁻²³ / 4.29 × 10⁻²³ ≈ 0.537 g/cm³** (Experimental: 0.534 g/cm³)

### Copper (Cu)

**Crystal Structure:** Face-Centered Cubic (FCC)

**Unit Cell Parameters:**
- Lattice constant: a = 3.61 Å (at room temperature)
- Atoms per unit cell: 8 × (1/8) (corners) + 6 × (1/2) (faces) = 4 atoms

**Atomic Positions:**
- Corner positions: (0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1), (0,1,1), (1,1,1)
- Face centers: (1/2, 1/2, 0), (1/2, 0, 1/2), (0, 1/2, 1/2), and three more

**Coordination Number:** 12 (each atom touches 12 nearest neighbors)

**Nearest Neighbor Distance:**
- Along face diagonal: a/√2 = 3.61/√2 = 2.55 Å

**Packing Fraction:**
- Volume of 4 atoms = 4 × (4/3)πr³
- From nearest neighbor: 2r = a/√2  →  r = a/(2√2)
- Packing Fraction = (4 × 4π/3 × (a/(2√2))³)/a³ = (π√2)/6 ≈ 0.74 or 74%

**Density Calculation:**
- Atomic mass of Cu = 63.55 g/mol
- Number of atoms per unit cell = 4
- Mass per unit cell = (4 × 63.55)/(6.022 × 10²³) g = 4.22 × 10⁻²² g
- Volume per unit cell = (3.61 × 10⁻⁸)³ = 4.71 × 10⁻²³ cm³
- **Density ρ = 4.22 × 10⁻²² / 4.71 × 10⁻²³ ≈ 8.96 g/cm³** (Experimental: 8.96 g/cm³)

### Comparative Structure Properties

| Property | Li (BCC) | Cu (FCC) |
|----------|----------|----------|
| Crystal System | Cubic | Cubic |
| Lattice Constant (Å) | 3.50 | 3.61 |
| Atoms per Unit Cell | 2 | 4 |
| Coordination Number | 8 | 12 |
| Packing Fraction | 68% | 74% |
| Density (g/cm³) | 0.534 | 8.96 |
| Nearest Neighbor (Å) | 3.03 | 2.55 |

---

# UNIT 2: QUANTUM MECHANICS AND MATTER WAVES

## Q1: Wave-Particle Duality and Features of Matter Waves

### Historical Background and Concept

**De Broglie's Hypothesis (1924):**

Every material particle, like an electron or photon, has both wave and particle characteristics. This is the fundamental principle of wave-particle duality.

De Broglie proposed:
**λ = h/p = h/(mv)**

Where:
- λ = wavelength (de Broglie wavelength)
- h = Planck's constant = 6.626 × 10⁻³⁴ J·s
- p = momentum
- m = mass
- v = velocity

### Wave-Particle Nature Duality

**Particle Nature (Classical):**
- Has definite position and momentum
- Follows Newton's laws
- Has mass and kinetic energy
- Exhibits inertia
- Examples: electrons, protons, atoms

**Wave Nature (Quantum):**
- Has wavelength and frequency
- Exhibits diffraction and interference
- Can have phase and group velocity
- Waves overlap and superpose
- Examples: electromagnetic radiation, matter waves

**Both Coexist:** An electron behaves as particle in some experiments (collision, energy measurement) and as wave in others (diffraction, interference).

### Features of Matter Waves

**1. De Broglie Wavelength:**

λ = h/p = h/(mv)

**For electron with kinetic energy E:**
E = (1/2)mv² = p²/(2m)
p = √(2mE)

**λ = h/√(2mE)**

Numerically: λ (in Å) = 12.27/√E (E in eV)

Example: 100 eV electron has λ = 12.27/√100 = 1.23 Å

**2. Frequency Associated with Matter Wave:**

From wave equation: c = νλ (not applicable for matter)

Instead, use: E = hν

**ν = E/h**

For particle at rest, E = mc² (relativistic):
ν = mc²/h

For moving particle: E = KE + mc² = (1/2)mv² + mc²

**3. Phase Velocity and Group Velocity:**

Phase velocity: **v_p = ν·λ = (E/h) × (h/p) = E/p = E/(mv)**

For non-relativistic particle: E = (1/2)mv²
v_p = (1/2)mv²/(mv) = v/2

**v_p = v/2** (half the particle velocity)

Group velocity: **v_g = dω/dk**

For matter waves: ω = E/ℏ, k = p/ℏ = mv/ℏ

v_g = d(E/ℏ)/d(mv/ℏ) = (1/ℏ)·(dE/d(mv/ℏ)) = (m/ℏ)·(dE/dm)

For KE = (1/2)mv²: dE/dm·v = v

**v_g = v** (equals particle velocity)

**4. Wavelength in Different Media:**

When matter wave enters medium:
- Wavelength changes: λ' = λ/n (like light)
- Frequency remains constant: ν' = ν
- Velocity changes: v' = ν·λ' = v/n

**5. Diffraction of Matter Waves:**

When electrons pass through crystal with spacing d:
- Bragg's condition: 2d sinθ = nλ
- Produces diffraction pattern like X-rays
- Confirms wave nature experimentally

### Summary Table of Matter Wave Features

| Feature | Expression | Significance |
|---------|------------|--------------|
| Wavelength | λ = h/p | Decreases with momentum |
| Frequency | ν = E/h | Higher for more energetic particles |
| Phase Velocity | v_p = v/2 | Half of particle velocity |
| Group Velocity | v_g = v | Equals particle velocity |
| Wave number | k = 2π/λ = p/ℏ | Related to momentum |
| Angular frequency | ω = E/ℏ | Related to energy |

---

## Q2: Uncertainty Principle and Its 3 Applications

### Heisenberg's Uncertainty Principle Statement

**The uncertainty principle states:**

It is impossible to simultaneously determine both the position and momentum of a particle with arbitrary precision. The product of uncertainties satisfies:

**Δx · Δp ≥ h/(4π) = ℏ/2**

Where:
- Δx = uncertainty in position
- Δp = uncertainty in momentum
- h = Planck's constant
- ℏ = h/(2π) = 1.055 × 10⁻³⁴ J·s

**Alternative Form (Energy-Time):**
**Δt · ΔE ≥ ℏ/2**

### Physical Interpretation

1. **Fundamental Limitation:** Not due to experimental imperfection but inherent to nature
2. **Trade-off:** To precisely know position, momentum must be uncertain, and vice versa
3. **Quantum Nature:** Implies particles don't have well-defined trajectories
4. **Scale:** Negligible for macroscopic objects; significant for atomic/subatomic particles

### Application 1: Size of Hydrogen Atom

**Problem:** Find minimum size of electron orbit in hydrogen atom

**Solution:**

For electron confined to atom:
- Let minimum Bohr radius = r₀
- Position uncertainty: Δx ≈ r₀
- Momentum uncertainty: Δp ≈ ℏ/r₀ (from uncertainty principle)

For circular orbit:
- Kinetic energy: T = p²/(2m) ≥ (Δp)²/(2m) = ℏ²/(2mr₀²)

- Potential energy: V = -e²/(4πε₀r₀) = -ke²/r₀

- Total energy: E = ℏ²/(2mr₀²) - ke²/r₀

**Minimizing energy:**
dE/dr₀ = -ℏ²/(mr₀³) + ke²/r₀² = 0

r₀ = ℏ²/(mke²) = a₀

**This gives Bohr radius: a₀ ≈ 0.53 Å**

**Conclusion:** Uncertainty principle explains why electron doesn't collapse into nucleus—minimum orbital size is determined by Δx·Δp ≈ h.

### Application 2: Lifetime of Excited States

**Problem:** Excited atomic states have width/lifetime. What's the connection to uncertainty principle?

**Solution:**

From energy-time uncertainty: **ΔE · Δt ≥ ℏ/2**

When atom is in excited state for time τ (lifetime):
- Energy uncertainty: ΔE ≥ ℏ/(2τ)

For excited state decaying to ground state:
- Line width in spectrum: Γ = ΔE = ℏ/τ

**Natural Line Width:**
Γ = ℏ/τ

For state with lifetime τ = 10⁻⁸ s:
Γ = (1.055 × 10⁻³⁴)/(10⁻⁸) = 1.055 × 10⁻²⁶ J
= 6.6 × 10⁻⁸ eV

**Consequence:** Even perfect spectral lines have finite width due to finite lifetime of energy states.

### Application 3: Tunneling and Barrier Penetration

**Problem:** Classical particle cannot pass energy barrier E > kinetic energy. Quantum particle can tunnel through.

**Explanation:**

For particle with kinetic energy E approaching barrier of height V₀ > E:

Inside barrier:
- Position uncertainty: Δx ≈ penetration depth
- Energy shortage: ΔE = V₀ - E
- From uncertainty: Δt ~ ℏ/ΔE

From uncertainty principle:
Δx ≥ ℏ/(2Δp) = ℏ/(2√(2m·ΔE))

Δx ≥ ℏ/(2√(2m(V₀ - E)))

**Tunneling Probability:**

For particle tunneling through barrier of width a:

P ~ exp(-2κa)

Where: κ = √(2m(V₀ - E))/ℏ

**Physical Meaning:**
- Particle can temporarily violate energy conservation
- Energy-time uncertainty allows ΔE = V₀ - E for time Δt ~ ℏ/ΔE
- Particle penetrates distance Δx ~ ℏ/(2m(V₀ - E)) into barrier
- Can emerge on other side with certain probability
- Exponentially decreases with barrier height and width

**Practical Examples:**
- Alpha decay: α particles tunnel out of nucleus
- Tunnel diodes: electrons tunnel through p-n junction
- Scanning tunneling microscope (STM): electrons tunnel through vacuum gap

---

## Q3: Schrödinger Equations and Wave Function Significance

### Time-Dependent Schrödinger Equation

**Derivation:**

Start with plane wave representing particle with energy E and momentum p:

ψ(x,t) = A·exp[i(kx - ωt)]

Where: k = p/ℏ, ω = E/ℏ

Taking derivatives:
∂ψ/∂t = -iω·ψ = -iE/ℏ·ψ  →  **iℏ(∂ψ/∂t) = E·ψ** ... (1)

∂ψ/∂x = ik·ψ = ip/ℏ·ψ  →  -ℏ²(∂²ψ/∂x²) = p²·ψ** ... (2)

From kinetic energy: T = p²/(2m)
From total energy: E = T + V = p²/(2m) + V

Substituting in (1) and (2):

**iℏ(∂ψ/∂t) = [-ℏ²/(2m)·(∂²ψ/∂x²) + V(x,t)]·ψ**

**Three-dimensional form:**

**iℏ(∂ψ/∂t) = [-ℏ²/(2m)∇²ψ + V(r,t)]ψ**

Or:  **iℏ(∂ψ/∂t) = Ĥψ**

Where Ĥ is Hamiltonian operator.

### Time-Independent Schrödinger Equation

**Derivation (Separation of Variables):**

Assume solution of form: ψ(x,t) = φ(x)·g(t)

Substituting into time-dependent equation:

iℏ·φ(x)·(dg/dt) = [-ℏ²/(2m)·(d²φ/dx²) + V(x)]·g(t)·φ(x)

Dividing by φ·g:

iℏ·(1/g)·(dg/dt) = [-ℏ²/(2m)·(1/φ)·(d²φ/dx²) + V(x)]

Left side depends only on t; right side only on x.

For equation to hold, both must equal constant E:

**iℏ·(1/g)·(dg/dt) = E**  →  g(t) = e^(-iEt/ℏ)

**[-ℏ²/(2m)·(d²φ/dx²) + V(x)]·φ(x) = E·φ(x)**

**Time-Independent Schrödinger Equation (1-D):**

**-ℏ²/(2m)·(d²ψ/dx²) + V(x)·ψ = E·ψ**

**Three-dimensional form:**

**-ℏ²/(2m)∇²ψ + V(r)ψ = Eψ**

Or: **(Ĥ)ψ = Eψ** (Eigenvalue equation)

### Physical Significance of Wave Function ψ(x,t)

**1. Probability Interpretation (Born Interpretation):**

- ψ itself has no direct physical meaning
- |ψ(x,t)|² = ψ*·ψ = probability density
- |ψ(x,t)|²·dx = probability of finding particle between x and x+dx at time t

**Normalization Condition:**

∫_{-∞}^{∞} |ψ(x,t)|²·dx = 1

(Total probability = 1)

**2. Expectation Values:**

- **Position:** <x> = ∫ x·|ψ|²·dx (average position)
- **Momentum:** <p> = ∫ ψ*·(-iℏ·dψ/dx)·dx (average momentum)
- **Energy:** <E> = ∫ ψ*·Ĥ·ψ·dx (average energy)

**3. Properties:**

- Continuous at boundaries (for bound states)
- Single-valued (unlike classical fields)
- ψ→0 as x→±∞ (for bound states in finite potential)
- First derivative continuous (for finite potentials)

**4. Superposition:**

Any linear combination ψ = c₁ψ₁ + c₂ψ₂ is also valid solution.
|c₁|² and |c₂|² give probabilities of finding particle in state 1 or 2.

**5. Phase Information:**

- Phase of ψ is important for interference
- |ψ|² alone doesn't contain complete information
- Phase determines probability amplitudes

**6. Observable Quantities:**

Observable properties correspond to eigenvalues of operators:
- Position x̂: eigenvalue = x
- Momentum p̂ = -iℏ(d/dx): eigenvalue = p
- Energy Ĥ: eigenvalue = E

---

## Q4: Phase Velocity and Group Velocity (Detailed)

### Definitions

**Phase Velocity:**

The velocity at which a constant phase point moves:

**v_p = ω/k = λ·ν = ν·(h/p) = E/p**

Where:
- ω = angular frequency = 2πν
- k = wave number = 2π/λ

**Group Velocity:**

The velocity at which energy (and information) travels in a wave packet:

**v_g = dω/dk**

### For Matter Waves: Derivation

**Energy-momentum relation:**

From E = hν and p = h/λ:
- ω = E/ℏ
- k = p/ℏ

For non-relativistic particle: E = p²/(2m) + E₀

**ω = (1/ℏ)·[p²/(2m)] + constant**
= (1/ℏ)·(ℏ·k)²/(2m) + const
= (ℏ·k²)/(2m) + const

**Phase velocity:**

v_p = ω/k = (ℏ·k)/(2m) = (ℏ·k)/(2m) = p/(2m)

From p = mv:
**v_p = mv/(2m) = v/2**

(Half the particle velocity!)

**Group velocity:**

v_g = dω/dk = d/dk[(ℏ·k²)/(2m)]
= (2·ℏ·k)/(2m)
= (ℏ·k)/m
= p/m

**v_g = p/m = v**

(Equals particle velocity!)

### Verification of Properties

**Property 1: v_p + v_g ≠ constant**

For matter waves:
- v_p = v/2
- v_g = v
- v_p + v_g = 3v/2 (not constant)

**Property 2: Dispersion Relation**

For matter: ω = (ℏ/2m)·k²

This is **dispersive** (ω not linear in k)
→ v_p ≠ v_g (different velocities for different frequencies)

**Property 3: Information Transfer**

- Energy and information travel at v_g = v
- Particle moves at v
- Therefore information travels with particle ✓

**Property 4: Energy Transport**

Energy flux: I = (1/2)v_g·ω·A² (proportional to group velocity)

---

## Q5: Proof of v_g = v and v_p = v/2 for Non-relativistic Particles

### Part (i): Group Velocity = Particle Velocity

**Starting Point:**

De Broglie relations:
- λ = h/p
- ν = E/h

**Wave vector and frequency:**
- k = 2π/λ = 2π·p/h = p/ℏ
- ω = 2πν = 2π·E/h = E/ℏ

**For non-relativistic particle:**
E = (1/2)mv² = p²/(2m)

**Angular frequency:**
ω = E/ℏ = p²/(2mℏ) = (ℏ·k²)/(2m)

**Group velocity definition:**

v_g = dω/dk

Taking derivative with respect to k:

dω/dk = d/dk[(ℏ·k²)/(2m)]
= (ℏ·2k)/(2m)
= (ℏ·k)/m

Since p = ℏ·k:

**v_g = p/m = mv/m = v** ✓

**Physical Meaning:** The group velocity represents the speed of a wave packet (bundle of waves). For matter, this packet travels at the same speed as the particle itself, so energy and information move with the particle.

### Part (ii): Phase Velocity = v/2 for Non-relativistic Particles

**Definition of phase velocity:**

v_p = ω/k

From above:
- ω = p²/(2mℏ) = (ℏ·k²)/(2m)
- k = p/ℏ

**Substituting:**

v_p = ω/k = [(ℏ·k²)/(2m)] / k
= (ℏ·k)/(2m)
= p/(2m)

From p = mv:

**v_p = mv/(2m) = v/2** ✓

**Alternative derivation using phase velocity definition:**

v_p = λ·ν

From de Broglie:
- λ = h/p = h/(mv)
- ν = E/h = (1/2)mv²/h

**v_p = [h/(mv)] × [(1/2)mv²/h]**
= [(1/2)mv²]/(mv)
= v/2 ✓

**Important Relationship:**

v_g = 2·v_p

This means:
- Energy moves 2× faster than phase fronts
- Characteristic of dispersive media (matter waves are dispersive)

---

# UNIT 3: FREE ELECTRON THEORY AND BAND STRUCTURE

## Q1: Classical Free Electron Theory - Successes and Failures

### Basic Assumptions

The Drude-Lorentz model assumes:

1. Valence electrons move freely inside metal (free electron gas)
2. Electrons confined within metal by potential barrier
3. Elastic collisions every τ seconds (relaxation time)
4. Collision probability = 1/τ (independent of velocity)
5. Between collisions, electrons follow Newton's laws
6. Ions are immobile and don't interact with electrons

### Equation of Motion

For electron under external field **E**:

**m(dv/dt) = -e·E - (m·v)/τ**

First term: force due to electric field
Second term: frictional/collision term (damping)

**At steady state (dv/dt = 0):**

**v = (-e·τ/m)·E = -μ·E**

Where: **μ = eτ/m** (mobility)

**Current density:**

**J = n·e·v = n·e²·τ/m·E = σ·E**

Where: **σ = n·e²·τ/m** (conductivity)

### Successes of Classical Theory

**1. Ohm's Law:**
J = σ·E is satisfied with conductivity σ = n·e²τ/m

**2. Thermal Conductivity (Wiedemann-Franz Law):**
κ = (π²/3)·(k_B/e)²·σ·T

Relation between thermal and electrical conductivity matches experiments well

**3. Temperature Dependence of Resistance:**
ρ(T) = ρ₀[1 + α(T - T₀)]

Linear temperature dependence observed in metals

**4. Hall Effect:**
Correctly explains Hall coefficient R_H = -1/(n·e) for single carrier type

**5. Positive Charge Carriers:**
Correctly predicts Hall effect with positive carriers (later explained by hole concept in semiconductors)

### Failures of Classical Theory

**1. Specific Heat Capacity (Major Failure):**

**Classical prediction:**
Each free electron contributes (3/2)k_B to heat capacity (3 translational degrees of freedom)

For N electrons: C_V = (3/2)·N·k_B

**Experimental observation:**
Heat capacity of metals ∝ T (linear), much smaller than predicted

The electronic contribution is negligible at room temperature!

Contradiction: Where are the electrons predicted by classical theory? (≈10²² electrons/cm³)

**2. Magnetic Susceptibility:**

**Classical prediction:**
Diamagnetic susceptibility χ ~ -10⁻¹⁰ (temperature independent)

**Experiment:**
χ ~ -10⁻⁵ (much larger diamagnetism)

Also, paramagnetic behavior observed for some metals

**3. Paramagnetism:**

**Classical theory:**
Cannot explain paramagnetism of electrons (both diamagnetism and paramagnetism arise in quantum theory)

**4. Electrical Conductivity vs Temperature:**

**Classical prediction:**
σ ∝ 1/√T (since τ ∝ 1/√T from collision theory)

**Experiment:**
σ ∝ 1/T in many metals at high temperature

**5. Thermionic Emission:**

**Classical prediction:**
Electron-ion collisions should prevent emission; no electron can leave metal

**Experiment:**
Electrons easily escape hot metal surface (work function ~ 1-5 eV, but thermal energy k_BT ~ 0.05 eV at 500K)

Classical explanation: impossible!

**6. Photoelectric Effect:**

Maximum kinetic energy of ejected electrons: E_max = hν - W

W = work function (typically 2-5 eV)

Classical prediction: light intensity → electron energy → no threshold frequency (continuous absorption)

Experiment: Below threshold frequency ν₀ = W/h, no electrons regardless of intensity

**7. Ferromagnetism:**

Classical theory cannot explain alignment of atomic magnetic moments at high temperatures (should randomize)

### Failure Resolution

These failures require **quantum mechanics:**
- Pauli exclusion principle → only few electrons contribute to heat capacity
- Quantum statistics → different energy distribution
- Fermi-Dirac distribution → explains low heat capacity at room T
- Band structure → explains conductivity better than free electron gas

---

## Q2: Quantum Free Electron Theory - Successes and Failures

### Basic Assumptions

1. Electrons confined in 3D "box" potential (infinite square well)
2. Solutions from Schrödinger equation with V = 0 inside, V = ∞ outside
3. Applies Pauli exclusion principle (max 2 electrons per state)
4. Uses Fermi-Dirac distribution
5. Fills states from lowest energy up to Fermi energy E_F

### Energy Levels in 3D Box

**From time-independent Schrödinger equation:**

ψ(x,y,z) = A·sin(n_x·πx/L)·sin(n_y·πy/L)·sin(n_z·πz/L)

Where: n_x, n_y, n_z = 1, 2, 3, ... (quantum numbers)

**Energy eigenvalues:**

**E_{n_x,n_y,n_z} = (ℏ²π²/2mL²)·(n_x² + n_y² + n_z²)**

In k-space: **E = (ℏ²k²)/(2m)** where k = π(n²)/L

**Density of states:**

g(E) = (2V/π²)·(2m/ℏ²)^(3/2)·√E

### Fermi Energy at T = 0

At absolute zero, all states filled up to Fermi energy E_F:

**E_F = (ℏ²/2m)·(3π²n)^(2/3)**

Where n = N/V (electron density)

### Successes of Quantum Theory

**1. Specific Heat - Correctly Predicts Linear Dependence:**

At low T: C_V = (π²/2)·(N·k_B²·T)/E_F ∝ T

Explains why heat capacity is small despite 10²² electrons! (At T ~ 300K, only ~ k_BT electrons excited above ground state)

**2. Thermionic Emission:**

Some electrons near Fermi surface have energy to escape if thermal fluctuation provides Δ E ≈ W (work function)

Richardson equation: I ∝ T²·exp(-W/k_BT)

Matches experiments reasonably

**3. Magnetic Susceptibility (Pauli Paramagnetism):**

Predicts susceptibility χ_P ~ (N·μ_B²)/(E_F·V) (density of states at Fermi surface)

Explains both diamagnetic and paramagnetic contributions

**4. Conductivity Temperature Dependence:**

With collision time τ ∝ 1/T at high T:
σ ∝ 1/T (matches experiment better than classical)

**5. Photoelectric Effect - Can Be Explained:**

Only electrons within k_BT of Fermi surface can escape
Most electrons "stuck" below work function barrier
Threshold frequency arises naturally from work function concept

**6. Electrical Conductivity:**

σ = n·e²τ/m (same form but n_eff = density of states at E_F ≪ total electrons)

Explains why only fraction of electrons participate

### Failures and Limitations of Quantum Free Electron Theory

**1. Doesn't Account for Band Structure:**

Assumes all states available at any energy (continuous spectrum)

Reality: Energy bands separated by band gaps
- Causes distinction between conductors, semiconductors, insulators
- Single-band model fails for these materials

**2. Transition Metals - Cannot Explain:**

Fe, Co, Ni show ferromagnetism despite no obvious unpaired spin
Theory doesn't account for d-band effects and exchange interactions

**3. Semiconductors - Completely Fails:**

Theory predicts conductivity proportional to temperature (∝ T^(3/2))
But semiconductors show σ ∝ exp(E_g/2k_BT) - exponential dependence!

Gap E_g between bands not in free electron model

**4. Negative Charge Carriers Not Explained:**

Hall effect shows positive charge carriers in some metals (holes)
Free electron theory predicts only electrons

Actually requires band theory concept of electron-hole duality

**5. Mobility Cannot Be Properly Calculated:**

Theory uses collision time τ as parameter but doesn't derive it from first principles
Scattering mechanisms (phonon, impurity, etc.) need detailed band structure

**6. Superconductivity - No Explanation:**

Some materials show zero resistance below critical temperature
Free electron theory predicts ρ > 0 always

Requires electron-phonon coupling (not in theory)

**7. Quantum Oscillations (Landau Levels) - Missing:**

In strong magnetic field, free electron theory gives wrong oscillation periods
True periodic structure requires proper electron dynamics in band structure

**8. Effective Mass Concept Needs Band Structure:**

In metals, electrons near Fermi surface have different mass m* ≠ m
Arises from band curvature, absent in infinite box model

### Conclusion

Quantum free electron theory successfully explains:
- Low-temperature heat capacity
- Basic transport phenomena
- Fermi statistics applications

But fails for:
- Materials with band gaps (semiconductors, insulators)
- Ferromagnetism
- Many properties requiring band structure details

**Needed:** Full band theory (Bloch theory) - Topic of Q3

---

## Q3: Bloch Function and K·P Model

### Bloch Function - Definition and Derivation

**Bloch's Theorem:**

In a periodic potential V(r + R) = V(r) (R = lattice vector), the wave function solutions satisfy:

**ψ_{n,k}(r) = e^{i·k·r}·u_{n,k}(r)**

Where:
- n = band index
- k = wave vector in Brillouin zone
- u_{n,k}(r) = periodic part with same periodicity as lattice: u_{n,k}(r + R) = u_{n,k}(r)

**Proof (Outline):**

Time-independent Schrödinger with periodic potential:

Ĥψ = [-(ℏ²/2m)∇² + V(r)]ψ = Eψ

Key insight: Translation operator T̂_R: T̂_R f(r) = f(r + R)

Since V is periodic: [Ĥ, T̂_R] = 0 (Hamiltonian commutes with translation)

Therefore T̂_R and Ĥ have simultaneous eigenstates.

Eigenvalue of T̂_R: e^{ikR} (from periodicity)

This gives form: ψ(r + R) = e^{ikR}·ψ(r)

Let: ψ(r) = e^{i·k·r}·u(r)

Then: ψ(r + R) = e^{ik·(r+R)}·u(r + R) = e^{ikR}·e^{ikr}·u(r) ✓

And: u(r + R) = u(r) ✓

### Physical Interpretation of Bloch Function

- **e^{ik·r}**: Plane wave part (momentum-like character)
- **u_{n,k}(r)**: Periodic part (contains lattice periodicity)

Electron in periodic potential "carries" lattice periodicity while having momentum k

The periodic part u_{n,k} modulates the plane wave

### K·P Method (Perturbation Theory for Band Structure)

**Physical Idea:**

Near a special k-point (often k = 0 or Γ point), expand band structure as power series in k

**Formalism:**

**Perturbed Hamiltonian:**

Ĥ(k) = Ĥ(0) + (ℏ/m)k·p̂ + (ℏ²/2m)k²

Where:
- Ĥ(0) = -(ℏ²/2m)∇² + V(r) at k = 0
- p̂ = -iℏ∇ (momentum operator)
- k·p̂ term: linear in k (first order)
- k² term: quadratic in k (second order)

**For two-band model near band edge:**

The perturbation matrix elements are:

H = [E_c(k)     P·k        ]
    [P*·k    E_v(k)    ]

Where:
- E_c(k), E_v(k) = band energies at k
- P = (ℏ/m)<c|p|v> = momentum matrix element

**Band structure near k = 0:**

E_±(k) = (E_c + E_v)/2 ± √[((E_c - E_v)/2)² + |P·k|²]

For small k (linear approximation):

**E_c(k) ≈ E_c + (ℏ²k²)/(2m_c*)**
**E_v(k) ≈ E_v - (ℏ²k²)/(2m_v*)**

Where effective masses:

**m_e* ≈ m[1 + (2|P|²)/(m·E_g)]** (conduction band)
**m_h* ≈ m[1 - (2|P|²)/(m·E_g)]** (valence band)

**Physical Significance:**

1. Band structure near edge can be calculated from simpler parameters (P, E_g)
2. Effective masses differ from free electron mass
3. Band curvature determines transport properties
4. Works near special symmetry points (Γ, X, L, etc.)

**Validity:**

Good for: k in few % of Brillouin zone size
Requires: k·p matrix elements (calculated from first-principles or experiments)

---

## Q4: Temperature Effect on Fermi-Dirac Function and Average Energy

### Fermi-Dirac Distribution Function

The probability that an electron state with energy E is occupied at temperature T:

**f(E,T) = 1/[exp((E - E_F)/(k_BT)) + 1]**

Where:
- E_F = Fermi level (chemical potential)
- k_B = Boltzmann constant = 1.38 × 10⁻²³ J/K
- T = absolute temperature

**Properties at T = 0:**
- f(E) = 1 for E < E_F (all states filled)
- f(E) = 0 for E > E_F (all states empty)
- Sharp step at E = E_F

### Effect of Temperature on Fermi-Dirac Function

**At T = 0:**

f(E,0) = {1,  E < E_F
          {1/2, E = E_F
          {0,  E > E_F

**At T > 0:**

- f(E,T) remains 1 for E ≪ E_F
- f(E,T) approaches 1/2 at E = E_F
- f(E,T) becomes 0 for E ≫ E_F
- Smooth transition replaces sharp step
- Tail extends k_BT above and below E_F

**Key observations:**
1. Function remains symmetric around E_F: f(E_F + x) + f(E_F - x) = 1

2. At E = E_F ± k_BT:
   f(E_F ± k_BT) = 1/[exp(±1) + 1] ≈ 0.27 and 0.73

3. "Thermal broadening" = 2k_BT (width of transition region)

**T → 0:** Sharp Fermi surface; abrupt cutoff
**T → ∞:** f(E) → 1/2 everywhere (no preference)

### Derivation of Average Energy

**Average energy of electrons in metal:**

<E> = ∫₀^∞ E·g(E)·f(E,T)·dE / ∫₀^∞ g(E)·f(E,T)·dE

Where g(E) = density of states

**For free electron gas:**

g(E) = C√E, where C is constant

**Substituting:**

<E> = ∫₀^∞ E^(3/2)/(exp((E-E_F)/(k_BT))+1) dE / ∫₀^∞ E^(1/2)/(exp((E-E_F)/(k_BT))+1) dE

**At T = 0:**

Numerator: ∫₀^{E_F} E^(3/2) dE = (2/5)E_F^(5/2)
Denominator: ∫₀^{E_F} E^(1/2) dE = (2/3)E_F^(3/2)

**<E>|_{T=0} = (3/5)E_F**

This is the **average Fermi energy** - the mean energy of occupied states at T = 0

**At Low Temperatures (T ≪ T_F where T_F = E_F/k_B):**

**<E> ≈ (3/5)E_F + (π²/12)(k_BT)²/E_F + ...**

The temperature correction arises from thermal excitation of electrons near Fermi surface.

**Physical Interpretation:**

At T = 0: Average energy is 60% of E_F (not E_F itself!) because states are filled from bottom

As T increases: Electrons excite to higher states, <E> increases

The tail of f(E,T) above E_F allows electrons to occupy states with E > E_F

**Result for Total Energy:**

Total energy: U = N·<E> = (3/5)N·E_F + (π²/4)N·k_B²T²/E_F + ...

**Heat Capacity:**

C_V = (∂U/∂T)|_V = (π²/2)N·k_B·(k_BT)/E_F = γT

Where: γ = (π²/2)·(N·k_B)/E_F

**Conclusion:** Explains linear heat capacity C ∝ T observed in experiments!

---

## Q5: Total Wave Functions in Band = Number of Primitive Unit Cells

### Statement and Significance

**Theorem:**

The total number of possible wave functions (or allowed k-states) in an energy band equals the number of primitive unit cells in the crystal.

**Corollary:**

If N atoms in crystal (hence N primitive cells) → N allowed k-values in any band

Each k-state can hold 2 electrons (spin up/down) → maximum 2N electrons per band

### Proof

**Step 1: Bloch Wave Functions and Periodicity**

From Bloch theorem: ψ_{n,k}(r) = e^{i·k·r}·u_{n,k}(r)

With periodic boundary conditions (Born-von Karman):
ψ(r + L) = ψ(r) for linear size L

This gives: e^{ikL} = 1 → k = 2πm/L where m = 0, ±1, ±2, ...

**Step 2: Number of Independent k-states in First Brillouin Zone**

For crystal with N atoms in 1D, each separated by lattice constant 'a':
- Length: L = N·a
- Allowed k-values in first Brillouin zone: -π/a ≤ k ≤ π/a
- Spacing between k-values: Δk = 2π/L = 2π/(N·a)

**Number of k-values:**
N_k = (2π/a)/(2π/(N·a)) = (2π/a) × (N·a/2π) = N

**Step 3: Three-Dimensional Generalization**

For 3D crystal with volume V = N·V_cell (N primitive cells, each volume V_cell):

From periodic boundary conditions:
k_x = 2πn_x/L_x,   k_y = 2πn_y/L_y,   k_z = 2πn_z/L_z

First Brillouin zone volume: Ω = (2π)³/V_cell

Density of k-states: ρ_k = V/(2π)³ = (N·V_cell)/(2π)³

**Number of k-states in 1BZ:**
N_k = ρ_k × Ω = [V/(2π)³] × [(2π)³/V_cell]
    = V/V_cell
    = N

**Step 4: Quantum Mechanical Basis**

Each k-state corresponds to one eigenfunction (for given band n) of Schrödinger equation with periodic boundary conditions.

The constraint that wave function is continuous and single-valued determines allowed k-vectors.

In first Brillouin zone, exactly N independent states per band.

### Density of States Implication

**In one band:**
- N k-values possible
- Each k corresponds to one state (per spin: 2 states including spin)

If band contains N_cells primitive cells:
- **Total states in band = N_cells × (2 if including spin)**

**Density of states (for one band):**

g(E) = (number of states)/(bandwidth) ~ constant for nearly-free electron band near band edge

### Application: Maximum Electrons per Band

For N atoms in crystal:
- N primitive cells
- Each band: N allowed k-states
- Each k-state: 2 electrons (spin up/down)
- **Maximum electrons per band: 2N**

**Fully filled band:** Contains 2N electrons (completely filled)
**Half-filled band:** Contains N electrons
**Empty band:** 0 electrons

---

## Q6: Effective Electrons in Filled Band and Material Classification

### Effective Number of Free Electrons in Completely Filled Band

**Theorem:**

The effective number of free electrons in a completely filled band is **zero**.

**Meaning:** Even though filled band contains 2N electrons, it behaves electrically like it has no free carriers.

### Derivation

**Conductivity from k-space:**

Current = (charge) × (velocity)

For state |n,k⟩:

v_{n,k} = (1/ℏ)∇_k E_n(k)

Current from one electron: **J_{nk} = (-e)·v_{nk} = (-e/ℏ)·∇_k E_n(k)**

**Total current in band:**

J_total = Σ_k Σ_s (-e/ℏ)·∇_k E_n(k)·f(E_{nk})

Where f = 1 for filled band (all states occupied)

**For completely filled band:**

Each k-state has its "negative" k-state due to symmetry:
- State at k has energy E_n(k)
- State at -k has energy E_n(-k) = E_n(k) (band symmetry)

For all k in first Brillouin zone:
Σ_k ∇_k E_n(k) = 0

**Proof of summation:**

In first Brillouin zone:
∑_k [∇_k E_n(k)] = ∑_k d/dk[E_n(k)]

Due to periodic boundary conditions and symmetry, contributions from +k and -k cancel:

∫ ∇E(k) dk ~ ∫_{-k}^{+k} dE = E|_{-k}^{+k} = 0 (same energy at ±k)

**Therefore: J = 0** for filled band

### Physical Interpretation

**Why filled band doesn't conduct:**

1. For every electron moving in +x direction (velocity +v), there's another moving in -x direction (velocity -v)
2. Contributions cancel exactly
3. No net current despite huge number of electrons

**Analogy:** Like water in completely full glass - can't flow despite being full of molecules

### Material Classification Based on Band Filling

#### 1. **CONDUCTORS (METALS)**

- Topmost occupied band **partially filled**
- Free electrons exist that can move under electric field
- Examples: Cu, Al, Fe, Ag

**Characteristics:**
- Electrical conductivity σ >> 1 S/m (typically 10⁷ S/m)
- Finite density of states at Fermi surface
- Electrons readily excited by small applied field

#### 2. **INSULATORS**

- All bands are either **completely filled or completely empty**
- Energy gap E_g between filled valence band and empty conduction band
- **E_g large** (> 5 eV) - very difficult to excite electrons

**Examples:** Diamond (E_g = 5.5 eV), Quartz (SiO₂, E_g ≈ 9 eV), NaCl (E_g ≈ 9 eV)

**Characteristics:**
- Electrical conductivity σ << 1 S/m (typically < 10⁻¹¹ S/m)
- Almost no free carriers at room temperature
- Thermal excitation insufficient to create carriers

**At Fermi level:** E_F lies within band gap (no states available!)

#### 3. **SEMICONDUCTORS**

- Similar to insulators: filled valence band + empty conduction band
- **E_g small** (1-3 eV) - electrons can be excited thermally or optically
- At T = 0: insulator-like
- At T > 0: conductivity increases exponentially

**Examples:** Si (E_g = 1.1 eV), Ge (E_g = 0.67 eV), GaAs (E_g = 1.42 eV)

**Characteristics:**
- Conductivity strongly temperature dependent: σ ∝ exp(-E_g/2k_BT)
- Doping allows tuning of carrier concentration
- Intrinsic (undoped) conductivity low, but extrinsic (doped) can be high
- At room temperature: moderate conductivity (10⁻¹⁰ to 10⁻⁶ S/m)

### E_F Position and Density of States

**Conductors:**
- E_F within conduction band
- Non-zero density of states g(E_F)
- Electrons can easily accept energy and move

**Insulators:**
- E_F within band gap
- g(E_F) = 0 (no states at Fermi level)
- Energy gap prevents conduction

**Semiconductors:**
- E_F near middle of band gap (for intrinsic)
- g(E_F) = 0 (gap) but thermally available states above gap
- Exponential conductivity: σ ∝ exp(-E_g/2k_BT)

---

# UNIT 4: SEMICONDUCTOR DEVICES

## Q1: PN Junction Diode and Its Characteristics

### PN Junction Formation

A PN junction is formed by joining p-type and n-type semiconductors in intimate contact. [web:33][web:36]

**Initial Condition:**
- p-side: majority carriers = holes (h⁺), minority = electrons (e⁻)
- n-side: majority carriers = electrons (e⁻), minority = holes (h⁺)

**Diffusion Process (at t = 0 to t = equilibrium):**
1. Electrons diffuse from n→p (high concentration to low)
2. Holes diffuse from p→n (high concentration to low)
3. Recombination occurs in depletion region
4. Immobile positive ions (donors) left on n-side
5. Immobile negative ions (acceptors) left on p-side

**Formation of Space Charge Region:**
- Builds up electric field (pointing n→p)
- Creates potential barrier V_{bi} (built-in voltage)
- Stops further diffusion at equilibrium
- Width of depletion region: W ≈ √(2ε₀εᵣ(V_{bi} - V)/(e(N_A + N_D)/(N_A·N_D)))

### Equilibrium Condition

At equilibrium (no external voltage):
**Drift current = Diffusion current**

J_drift = J_diffusion

This gives built-in potential:

**V_{bi} = (k_BT/e)·ln(N_A·N_D/n_i²)**

Where n_i = intrinsic carrier concentration

For Si: V_{bi} ≈ 0.7-0.8 V (room temperature)
For Ge: V_{bi} ≈ 0.3-0.4 V

### Forward Bias Operation

**Applied voltage V_F > 0 (p-side positive):**

- Reduces effective barrier: V_{eff} = V_{bi} - V_F
- Depletion width reduces: W_F = √(2ε₀εᵣ(V_{bi} - V_F)/(e(N_A + N_D)/(N_A·N_D)))
- Majority carriers pushed across junction
- Forward current increases exponentially

**Ideal Diode Equation:**

**I = I_s·[exp(eV_F/(ηk_BT)) - 1]**

Where:
- I_s = saturation current (leakage current, very small)
- η = ideality factor (1-2, typically ≈ 1)
- V_F = forward voltage applied
- k_BT/e ≈ 26 mV at T = 300K

**Behavior regions:**
1. V_F < V_{knee}: I exponentially small (I << I_s)
2. V_F ≈ V_{knee} ≈ 0.6-0.7 V (Si): I rises sharply
3. V_F > V_{knee}: I increases with slope ~1/R_series (series resistance limits)

### Reverse Bias Operation

**Applied voltage V_R < 0 (p-side negative):**

- Increases effective barrier: V_{eff} = V_{bi} + |V_R|
- Depletion width increases: W_R increases
- Majority carriers pulled away from junction
- Only minority carriers cross junction (thermally generated)
- Reverse current ≈ I_s (nearly constant)

**Reverse Current:**

**I_R ≈ I_s = A·e·D·n_i²/(N_A·τ_A) + A·e·D·p_i²/(N_D·τ_D)**

Where:
- A = junction area
- D = diffusion coefficient
- n_i = intrinsic carrier concentration
- τ = minority carrier lifetime

**Breakdown (Reverse Saturation Increases):**

At high reverse bias, when |V_R| → V_{breakdown}:
- Electric field at junction ≈ E_c (critical field ~10⁵ V/cm)
- Two mechanisms:
  1. **Avalanche breakdown:** Minority carriers accelerated, create more carriers by impact ionization
  2. **Zener breakdown:** Electrons in valence band tunnel into conduction band across narrow band gap

Reverse current increases sharply (exponentially) near breakdown.

### VI Characteristic Curve Description

[For sketch: x-axis = Voltage V, y-axis = Current I]

**Forward Bias (right side, V > 0):**
- Negligible current until knee (≈0.6 V Si, ≈0.3 V Ge)
- Sharp exponential rise after knee
- Approximately linear at high forward voltage (series resistance dominates)
- Typical forward current: mA range

**Reverse Bias (left side, V < 0):**
- Nearly flat line at small negative saturation current I_s (μA to nA range)
- Extends over wide reverse voltage range
- Sharply rises (↓) at breakdown voltage V_{br} (typically 50-1000 V depending on doping)
- After breakdown: rapid current increase with slight voltage increase

**Key Parameters on Graph:**
1. Knee/threshold voltage ≈ 0.6 V (Si)
2. Reverse saturation current I_s ≈ μA
3. Breakdown voltage V_{br}
4. Dynamic resistance r_d = dV/dI (low forward, high reverse)

### Temperature Dependence

**Forward voltage decreases with T:**
- Knee voltage ≈ -2.2 mV/°C (typical)
- At higher T: I_s increases exponentially

**Reverse saturation increases with T:**
- I_s ∝ exp(-E_g/(k_BT))
- Roughly doubles every 5-10°C

---

## Q2: Common Base BJT - Working and Characteristics

### BJT Structure and Terminals

**n-p-n BJT (typical):**
- **Emitter (E):** Heavily doped n-region (n⁺)
- **Base (B):** Thin, lightly doped p-region
- **Collector (C):** Moderately doped n-region (larger area)

**Doping levels:** N_D(E) >> N_A(B) >> N_D(C)

### Common Base Configuration

**Connection:**
- Base terminal common to input and output circuits
- **Input:** Between Emitter-Base (EB)
- **Output:** Between Collector-Base (CB)

**Biasing for n-p-n:**
- EB junction: Forward biased (V_EB ≈ 0.7 V)
- CB junction: Reverse biased (V_CB > 0)

### Working Principle (Physical Operation)

**Step 1: Emitter Injection**

With EB forward biased:
- Electrons injected from emitter into base with high velocity
- Hole injection from base to emitter (small, n⁺ heavily doped)

**Step 2: Transport Through Base**

Thin base region (W_B ≈ 1 μm):
- Most injected electrons diffuse across without recombination
- Small fraction recombine with base holes (creates base current I_B)
- Electrons reach CB junction with E ≈ E_B

**Step 3: Collection at Collector**

At CB junction (reverse biased):
- Strong electric field sweeps electrons into collector
- Electrons accelerate, gaining additional kinetic energy
- Collected electrons become collector current I_C

**Amplification Mechanism:**

Small base current controls large collector current:
**I_C = α·I_E**

Where: α = I_C/I_E ≈ 0.95-0.99 (current gain, less than 1)

Also: I_E = I_B + I_C

Therefore: **I_C = (α/(1-α))·I_B = β·I_B**

Where: **β = α/(1-α)** = current amplification factor (50-300 typically)

### Input Characteristics (I_E vs V_EB at constant V_CB)

**Shape:**
- Forward diode-like behavior
- Negligible I_E until V_EB ≈ 0.6-0.7 V (cut-in)
- Rapid exponential rise for V_EB > 0.6 V
- Follows diode equation: I_E = I_{E0}·[exp(eV_EB/(k_BT)) - 1]

**Effect of V_CB:**

Higher V_CB (more reverse bias at collector):
- Depletion width increases at CB junction
- Base width W_B effectively decreases (Early effect)
- Slightly higher I_E for same V_EB (weak effect)
- Curves shift upward for increasing V_CB

**Input Resistance:**
r_{in} = ΔV_EB/ΔI_E ≈ (k_BT·V_T)/(I_E) ≈ 25-50 Ω (low input impedance)

Where V_T = k_BT/e ≈ 26 mV at 300K

### Output Characteristics (I_C vs V_CB at constant I_E)

**Three operating regions:**

**1. Cut-off Region (I_E ≈ 0):**
- Both EB and CB junctions reverse biased
- Only small leakage current I_CB0 flows (opposite direction)
- I_C ≈ I_CB0 ≈ few nA

**2. Active Region (I_E moderate, V_CB > 0.5 V):**
- EB forward biased, CB reverse biased (normal operation)
- I_C nearly independent of V_CB (approximately constant)
- Small slope = finite output resistance r_{out}
- For fixed I_E: I_C ≈ α·I_E (nearly constant)

**Output resistance:**
r_{out} = ΔV_CB/ΔI_C ≈ V_A/I_C (can be kΩ to MΩ)

Where V_A = Early voltage (≈50-200 V)

**3. Saturation Region (V_CB < 0.2 V):**
- CB junction becomes forward biased (lightly)
- Both junctions conduct in "wrong" direction
- Excess base current accumulates in base
- I_C no longer proportional to I_E
- V_CE(sat) ≈ 0.1-0.2 V (very small)
- Transistor acts like switch (ON state)

### Important CB Parameters

| Parameter | Value | Significance |
|-----------|-------|--------------|
| Current gain α | 0.95-0.99 | I_C/I_E ratio |
| Current gain β | 50-300 | I_C/I_B ratio |
| Input impedance r_{in} | 10-100 Ω | Low (good for current input) |
| Output impedance r_{out} | kΩ-MΩ | High (good for voltage output) |
| Voltage gain | (r_{out}/r_{in})·(ΔI_C/ΔI_E) ≈ high | Significant amplification |
| Frequency response | High f_{T} | Suitable for RF circuits |

### Advantages and Applications of CB Configuration

1. **High voltage gain:** Due to large (r_{out}/r_{in}) ratio
2. **Low input impedance:** Matches well with low-impedance sources
3. **High output impedance:** Matches well with high-impedance loads
4. **Good frequency response:** Less capacitive effects than CE
5. **Low noise:** Suitable for RF amplification
6. **High f_{max}:** Good for microwave frequencies

**Applications:**
- Low-noise RF preamplifiers
- High-frequency amplifiers
- Impedance matching circuits
- Cascode amplifiers

---

## Q3: Intrinsic Semiconductor - Derivation of Electron Concentration and Law of Mass Action

### Definition: Intrinsic Semiconductor

A pure semiconductor material with no impurities (extrinsic dopants).
- Electrons and holes generated only by thermal excitation across band gap
- n_e = n_h (charge neutrality in intrinsic material)
- Both electrons and holes serve as charge carriers

Examples: Pure Si, Ge at moderate temperatures (before heavy doping)

### Derivation of Electron Concentration in Conduction Band

**Step 1: Fermi-Dirac Distribution**

Probability that a state at energy E is occupied:

f(E) = 1/[exp((E - E_F)/(k_BT)) + 1]

For conduction band (above E_F by typically 0.3-0.5 eV at room T):
(E - E_F) >> k_BT

Therefore: exp((E - E_F)/(k_BT)) >> 1

**Boltzmann Approximation:**

f(E) ≈ exp(-(E - E_F)/(k_BT))

**Step 2: Density of States in Conduction Band**

For parabolic band near conduction band edge E_c:

E = E_c + (ℏ²k²)/(2m_e*)

g_c(E) = C_c·(E - E_c)^(1/2)

Where: **C_c = (2m_e*)^(3/2)/(π²ℏ³)** (constant for given band structure)

More explicitly: **g_c(E) = (2π)^(1/2)·(2m_e*)^(3/2)·(E - E_c)^(1/2)/h³**

**Step 3: Calculate Electron Concentration**

Number of electrons per unit volume at energy E to E+dE:

dn_e = g_c(E)·f(E)·dE

Total electron concentration:

**n_e = ∫_{E_c}^{∞} g_c(E)·f(E)·dE**

= ∫_{E_c}^{∞} C_c·(E - E_c)^(1/2)·exp(-(E - E_F)/(k_BT))·dE

**Step 4: Change of Variables**

Let: ε = E - E_c (energy above conduction band edge)

n_e = C_c·exp(-(E_c - E_F)/(k_BT))·∫_0^{∞} ε^(1/2)·exp(-ε/(k_BT))·dε

**Step 5: Evaluate Integral**

∫_0^{∞} ε^(1/2)·exp(-ε/(k_BT))·dε = (√π/2)·(k_BT)^(3/2)

Therefore:

n_e = C_c·(√π/2)·(k_BT)^(3/2)·exp(-(E_c - E_F)/(k_BT))

**Define effective density of states:**

**N_c = 2·(2πm_e*k_BT/h²)^(3/2)**

This is the "effective" conduction band edge density.

**Final Result:**

**n_e = N_c·exp(-(E_c - E_F)/(k_BT))**

Or equivalently:

**n_e = N_c·exp(-Δ_E/(k_BT))**

Where: Δ_E = E_c - E_F (energy of E_c below E_F)

**Numerically for Si (300K):**
- N_c ≈ 2.8 × 10¹⁹ cm⁻³

### Hole Concentration in Valence Band (Similar Derivation)

By parallel argument:

**p_h = N_v·exp(-(E_F - E_v)/(k_BT))**

Where:
- E_v = valence band edge
- N_v = 2·(2πm_h*k_BT/h²)^(3/2) (effective valence band density)
- m_h* = effective hole mass

**For Si (300K):**
- N_v ≈ 1.0 × 10¹⁹ cm⁻³

### Law of Mass Action (Product of Carrier Concentrations)

**Multiply electron and hole concentrations:**

n_e·p_h = N_c·N_v·exp(-(E_c - E_F)/(k_BT))·exp(-(E_F - E_v)/(k_BT))

= N_c·N_v·exp(-(E_c - E_v)/(k_BT))

= N_c·N_v·exp(-E_g/(k_BT))

**Where: E_g = E_c - E_v (band gap energy)**

**Final Law of Mass Action:**

**n_i² = N_c·N_v·exp(-E_g/(k_BT))**

Or: **n_i = √(N_c·N_v)·exp(-E_g/(2k_BT))**

Where: **n_i = intrinsic carrier concentration** (electrons and holes are equal in intrinsic material)

### Physical Significance of Mass Action Law

**1. Independence from E_F:**

The product n_i² depends ONLY on:
- Band gap E_g (material property)
- Temperature T
- Effective masses (material property)

Does NOT depend on Fermi level position! (Unlike individual n or p)

**2. Temperature Dependence:**

n_i² ∝ exp(-E_g/(k_BT))

- Lower E_g → higher n_i (Ge > Si > GaAs)
- Higher T → higher n_i (exponentially)
- Roughly doubles every 5-10°C for Si

**3. At Equilibrium in Intrinsic Material:**

n_e = n_h = n_i

And Fermi level position:

E_F = (E_c + E_v)/2 + (k_BT/2)·ln(N_v/N_c)

For Si: m_e* ≈ m_h*, so E_F ≈ (E_c + E_v)/2 (mid-gap)

**4. Numerical Example for Si at 300K:**

- E_g = 1.12 eV
- N_c ≈ 2.8 × 10¹⁹ cm⁻³
- N_v ≈ 1.0 × 10¹⁹ cm⁻³
- k_BT ≈ 0.026 eV

n_i² = (2.8 × 10¹⁹)(1.0 × 10¹⁹)·exp(-1.12/0.026)
    = 2.8 × 10³⁸·exp(-43.08)
    = 2.8 × 10³⁸ × 1.3 × 10⁻¹⁹
    ≈ 3.6 × 10¹⁹

**n_i ≈ 1.5 × 10¹⁰ cm⁻³**

(Experimental value: 1.5 × 10¹⁰ cm⁻³ ✓)

---

## Q4: Extrinsic Semiconductors - Electron Concentration and Donor Doping

### Definition: Extrinsic Semiconductor

Semiconductor doped with impurities (dopants):
- **Donor dopants:** Group V elements in Si/Ge lattice (P, As, Sb)
  - 5 valence electrons → 4 bonded in lattice, 1 extra
  - Easily donates electron to conduction band
  - Creates n-type semiconductor

- **Acceptor dopants:** Group III elements in Si/Ge lattice (B, Al, Ga)
  - 3 valence electrons → 4 needed for bonding
  - Accepts electron from valence band
  - Creates p-type semiconductor

Focus here: **n-type with donor doping**

### Energy Levels: Donor States

When donor atom (e.g., P in Si) substitutes for host atom:
- 4 electrons participate in Si-P bonds
- 5th electron weakly bound by Coulomb attraction

**Binding energy of donor electron:**

E_d = 13.6 eV × (m_e*/m_e) × (1/εᵣ²)

For Si:
- m_e*/m_e ≈ 0.26 (effective mass ratio)
- εᵣ ≈ 11.7 (dielectric constant)

E_d ≈ 13.6 × 0.26 / (11.7)² ≈ 0.026 eV ≈ 2.6 meV

**Donor energy level:**
E_D ≈ E_c - E_d

(About 26 meV below conduction band edge for Si)

**Ionization:** At room temperature (k_BT ≈ 26 meV), donor atoms easily ionized!

### Derivation of Electron Concentration in n-Type

**Step 1: Charge Neutrality Condition**

Total positive charge = Total negative charge

**N_D^+ + p = n + N_A^-**

Where:
- N_D^+ = ionized donor concentration
- p = hole concentration
- n = electron concentration
- N_A^- = ionized acceptor concentration (for compensated doping, usually ≈ 0)

**For uncompensated n-type (N_A = 0):**

**n + N_D^+ = p + 0**

**Step 2: Ionization Equilibrium**

Donor ionization:
D ⇌ D^+ + e^-

With ionization constant:

K_D = [D^+][e^-]/[D] ∝ exp(-E_d/(k_BT))

Number of ionized donors vs. total donors:

**N_D^+ = N_D/(1 + 2·exp(-(E_c - E_D)/(k_BT)))**

For T > 100K (most conditions at room T and above):
Nearly all donors ionized: **N_D^+ ≈ N_D**

**Step 3: Temperature Regimes**

**(A) High temperature (T >> T_D = E_d/k_B ≈ 300K for Si):**

All donors ionized: N_D^+ = N_D

Intrinsic region: n ≈ n_i >> N_D (dominated by thermal generation)

From law of mass action: n·p = n_i²

With n > p: n·p ≈ n_i² → **n ≈ n_i**

**(B) Intermediate temperature (T ≈ 100-300K):**

Extrinsic region: n ≈ N_D (donor-dominated conduction)

Most donors ionized, but not all intrinsic carriers

From charge neutrality: n ≈ N_D^+ ≈ N_D

From mass action: p = n_i²/n ≈ n_i²/N_D (holes from intrinsic generation)

**Total: n ≈ N_D + n_i (if N_D > n_i)**

**(C) Low temperature (T << T_D):**

Few donors ionized: n < N_D

**Freeze-out region:** n decreases with decreasing T

From partial ionization:
n² ≈ N_c·N_D·exp(-E_d/(k_BT))

**Step 4: Final Expression in Extrinsic Region

At room temperature (T ≈ 300K), if donors mostly ionized:

From charge neutrality: **n = N_D + p**

From mass action: **n·p = n_i²**

Combining (for N_D >> n_i):

p = n_i²/n

n - n_i²/n = N_D

n² - N_D·n - n_i² = 0

**Solving quadratic:**

n = [N_D + √(N_D² + 4n_i²)]/2

For N_D >> n_i:

**n ≈ N_D**

(Electron concentration approximately equals donor concentration)

### Proof: Electron Concentration ∝ √N_D (at specific conditions)

**This proof applies in a specific low-doping/low-temperature regime:**

**Assumptions:**
1. Low doping: N_D small
2. Partial ionization: Not all donors ionized
3. Intrinsic carriers negligible: n >> n_i initially

**From partial ionization (freeze-out region):**

Ionized donor density: n = N_D^+

From ionization equilibrium and Fermi statistics:

**n·(N_D - n) ∝ N_c·exp(-E_d/(k_BT))**

Rearranging:

n² - N_D·n + N_c·N_D·exp(-E_d/(k_BT)) = 0

**For low doping where E_F is far below E_c:**

Using Boltzmann approximation more carefully:

**n = √[N_D·N_c·2·exp(-E_d/(k_BT))]**

This shows: **n ∝ √N_D·√N_c·√[exp(-E_d/(k_BT))]**

For fixed T and material: **n ∝ √N_D** ✓

### Numerical Example

**Si doped with N_D = 10¹⁶ cm⁻³ at 300K:**

- Intrinsic n_i ≈ 1.5 × 10¹⁰ cm⁻³
- N_D >> n_i, so extrinsic region

**Electron concentration:**
n ≈ N_D = 10¹⁶ cm⁻³ (room temperature)

**Hole concentration:**
p = n_i²/n = (1.5 × 10¹⁰)²/(10¹⁶)
  = 2.25 × 10²⁰/(10¹⁶)
  = 2.25 × 10⁴ cm⁻³

**Ratio:**
n/p = 10¹⁶/(2.25 × 10⁴) ≈ 4.4 × 10¹¹

(Electrons dominate; negligible holes!)

---

## Q5: Half Wave and Full Wave Rectifiers

### Half Wave Rectifier (HWR)

#### Circuit Description

**Components:**
- AC source: V_s(t) = V_m·sin(ωt)
- Transformer (optional): steps voltage up/down
- Diode D: conducts when anode +, blocks when −
- Load resistor: R_L
- No filter capacitor (for analysis)

**Circuit topology:**
```
AC source → Diode D → Load R_L
            (Cathode)
```

#### Working Principle

**Positive half-cycle (V_s > 0):**
- Anode is positive, cathode is at lower potential
- Diode is forward biased (ON)
- Current flows through R_L: i = V_s/R_L
- Voltage across R_L: V_o = V_s

**Negative half-cycle (V_s < 0):**
- Anode is negative, cathode is positive
- Diode is reverse biased (OFF)
- No current flows: i = 0
- Voltage across R_L: V_o = 0 (or small reverse leakage)

**Result:** Output is only positive half-waves (0 to V_m)

#### Output Waveform

- Positive half-sine waves: V_m·sin(ωt) for 0 ≤ t ≤ π/ω
- Zero for π/ω ≤ t ≤ 2π/ω
- Repeats at 2ω frequency
- Average (DC) value: **V_avg = V_m/π ≈ 0.318·V_m**

#### Mathematical Analysis

**Average output voltage (DC):**

V_avg = (1/T)∫_0^T V_o dt

where T = 2π/ω (period)

V_avg = (1/(2π))∫_0^π V_m·sin(ωt)·d(ωt)

= (V_m/(2π))·[-cos(ωt)]_0^π

= (V_m/(2π))·[-cos(π) + cos(0)]

= (V_m/(2π))·[1 + 1]

**V_avg = V_m/π**

**RMS voltage:**

V_rms² = (1/T)∫_0^T V_o² dt = (1/(2π))∫_0^π V_m²·sin²(ωt)·d(ωt)

Using sin²θ = (1 - cos(2θ))/2:

V_rms² = (V_m²/(2π))·(π/2) = V_m²/4

**V_rms = V_m/2**

**Form Factor:**

F = V_rms/V_avg = (V_m/2)/(V_m/π) = π/2 ≈ **1.57**

(Higher FF means more ripple)

**Ripple Factor:**

**r = √(V_rms² - V_avg²)/V_avg = √[(V_m/2)² - (V_m/π)²]/(V_m/π)**

= √[π²/4 - 1]·π/π² ≈ **1.21** (or 121%)

**Peak Inverse Voltage (PIV):**

Maximum reverse voltage across diode (when blocking):
**PIV = V_m**

(Important for diode selection)

#### Efficiency

**Rectification efficiency:**

η = (P_DC)/(P_AC) = (V_avg²/R_L)/(V_rms²/R_L) = (V_avg/V_rms)²

= (V_m/π)² / (V_m/2)² = (2/π)² ≈ **0.406 or 40.6%**

(Less than 50% of input power reaches load as DC)

### Full Wave Rectifier (FWR)

#### Center-Tap Full Wave Rectifier (CT-FWR)

**Circuit Components:**
- AC source with center tap (2 secondary coils)
- Two diodes D1, D2
- Load resistor R_L
- Center tap grounded (reference)

**Circuit Topology:**
```
Upper coil → D1 ↘
              ├ R_L (to ground)
              ↙ D2
Lower coil → /
```

#### Working Principle

**Positive half-cycle (upper terminal positive):**
- D1 forward biased (ON)
- D2 reverse biased (OFF)
- Current through D1 → R_L → D2
- Output voltage: V_o = V_s (forward drop neglected)

**Negative half-cycle (lower terminal positive):**
- D1 reverse biased (OFF)
- D2 forward biased (ON)
- Current through D1 ← R_L ← D2 (opposite direction)
- Output voltage: V_o = -V_s (but appears positive across R_L)

**Result:** Both half-cycles produce positive output

#### Output Waveform (CT-FWR)

- Full-wave rectified sine: |V_m·sin(ωt)| for entire period
- Positive voltage at frequency 2f (twice AC frequency)
- Average DC: **V_avg = 2V_m/π ≈ 0.637·V_m**

#### Mathematical Analysis (CT-FWR)

**Average voltage:**

V_avg = (1/T)∫_0^T |V_o| dt = (1/T)∫_0^{T/2} V_m·sin(ωt)·d(ωt) + (1/T)∫_{T/2}^T V_m·sin(ωt)·d(ωt)

= 2 × (V_m/(2π))∫_0^π V_m·sin(θ)·dθ

= (2V_m/π)·[-cos(θ)]_0^π

**V_avg = 2V_m/π**

**RMS voltage:**

V_rms = √[(1/T)∫_0^T V_m²·sin²(ωt)·dt]

For full wave: Integration over full period with full-wave rectified sine

V_rms² = V_m²·(1/T)∫_0^{T/2} sin²(ωt)·dt

= V_m²·(1/T)·(T/4)

**V_rms = V_m/√2**

**Ripple Factor:**

r = √[V_rms² - V_avg²]/V_avg

= √[(V_m/√2)² - (2V_m/π)²]/(2V_m/π)

= √[π²/2 - 4]·π/(2π²)

= √[4.93 - 4]/2 ≈ **0.48** (or 48%)

(Much lower than HWR!)

**Rectification Efficiency:**

η = (V_avg/V_rms)² = (2V_m/π)²/(V_m/√2)²

= (4/π²)/(1/2) = 8/π² ≈ **0.812 or 81.2%**

(Much higher than HWR's 40.6%!)

**Peak Inverse Voltage (PIV):**

When one diode blocks (other conducts), the full input voltage appears across the blocking diode:
**PIV = 2V_m**

(Twice that of HWR - requires higher-rated diodes)

#### Bridge Full Wave Rectifier (Bridge FWR)

**Alternative circuit:** Uses 4 diodes in bridge, single-phase input (no center tap needed)

**Working:**
- 2 diodes conduct during +half cycle
- Other 2 conduct during −half cycle
- Output always positive

**Advantages:**
- No center-tap transformer needed
- Same performance as CT-FWR
- More efficient transformer use
- PIV = V_m (lower than CT)

**Output same as CT:** V_avg = 2V_m/π, r ≈ 0.48

### Comparative Table: HWR vs FWR

| Parameter | Half Wave | Full Wave (CT) | Full Wave (Bridge) |
|-----------|-----------|----------------|-------------------|
| V_avg | V_m/π | 2V_m/π | 2V_m/π |
| V_rms | V_m/2 | V_m/√2 | V_m/√2 |
| Ripple Factor r | 1.21 (121%) | 0.48 (48%) | 0.48 (48%) |
| Efficiency η | 40.6% | 81.2% | 81.2% |
| PIV | V_m | 2V_m | V_m |
| Output Freq | f | 2f | 2f |
| Diodes | 1 | 2 | 4 |
| Transformer | Any | Center-tap | Standard |
| DC Output | Unidirectional | Unidirectional | Unidirectional |

### Practical Considerations

**With Filter Capacitor:**

In real circuits, smoothing capacitor C is connected across load:

**Filter action:**
- Capacitor charges to V_m during conducting interval
- Holds charge and discharges through R_L when diode blocks
- Ripple voltage (peak-to-peak): V_r ≈ V_m/(f·RC) for HWR
- For FWR: V_r ≈ V_m/(2f·RC) (less ripple due to 2f ripple frequency)

**Peak voltage stress:**
- Diode must block at least PIV
- Forward current surge when capacitor charges
- Series resistance limits this surge

**Transformer Utilization Factor (TUF):**

Only portion of transformer secondary used efficiently

For CT-FWR: TUF ≈ 2/π ≈ 0.637

For Bridge FWR: TUF ≈ 4/π ≈ 1.27 (better utilization)

---

## Summary

This completes all solutions for Unit 4:
1. **PN Junction Diode** - Formation, biasing modes, VI characteristics, temperature dependence
2. **Common Base BJT** - Structure, biasing, working principle, input/output characteristics, parameters
3. **Intrinsic Semiconductors** - Electron concentration derivation, mass action law with proof
4. **Extrinsic Semiconductors** - Donor doping, electron concentration, √N_D relationship
5. **Rectifiers** - HWR and FWR with complete derivations of efficiency, ripple, PIV

All derivations include step-by-step mathematics, physical interpretation, and numerical examples suitable for B.Tech 1st year Kurukshetra University examination.

---

# END OF DOCUMENT

**Note:** For Units 1, 2, and 3 with similar detailed coverage, please request separately. Each unit contains comprehensive theory, derivations, diagrams descriptions, and exam-oriented solutions.
