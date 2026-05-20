# PPV Attenuation Analysis

Empirical blast vibration analysis using the **Holmberg-Persson** attenuation model, fitted to 200 field monitoring points.

## Results

| Parameter | Value |
|-----------|-------|
| K | 2950.93 |
| α (alpha) | 1.369 |
| R² | 0.75 |

Formula: **PPV = K × SD^(−α)** where SD = D / √W

---

## Project structure

```
ppv-attenuation/
├── blast_monitoring_data.csv   ← raw field data (200 points)
├── ppv_analysis.py             ← Day 1 + Day 2 analysis script
├── requirements.txt
├── outputs/
│   ├── master_ppv_clean.csv    ← cleaned dataset with log columns
│   ├── hp_params.json          ← K, alpha, R² for UDEC input
│   ├── ppv_scatter.png         ← log-log attenuation scatter
│   ├── hp_fit.png              ← H-P model fit vs measured
│   └── ppv_components.png      ← V / L / T component breakdown
└── README.md
```

---

## Setup

### 1. Clone the repo
```bash
git clone https://github.com/YOUR_USERNAME/ppv-attenuation.git
cd ppv-attenuation
```

### 2. Create a virtual environment (recommended)
```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

### 4. Run the analysis
```bash
python ppv_analysis.py
```

Outputs are saved in the current directory. Move them to `outputs/` if you want to keep the root clean.

---

## Column reference

| Column | Description |
|--------|-------------|
| `Point` | Monitoring point ID |
| `D_m` | Distance from blast (m) |
| `W_kg` | Charge per delay (kg) |
| `SD` | Scaled distance = D / √W (m/kg⁰·⁵) |
| `PPV_mms` | Resultant PPV (mm/s) |
| `PPV_V` | Vertical PPV component |
| `PPV_L` | Longitudinal PPV component |
| `PPV_T` | Transverse PPV component |

---

## Using K and alpha in UDEC

Load `hp_params.json` in your UDEC pre-processing script:

```python
import json
with open("hp_params.json") as f:
    params = json.load(f)
K, alpha = params["K"], params["alpha"]
# Generate PPV input values: PPV = K * SD**(-alpha)
```
