In our study of synaptic plasticity, we use the following notation:

- **+1** → LTP (Long-Term Potentiation)  
- **-1** → LTD (Long-Term Depression)  
- **0** → Basal (No Change)

## STDP Pattern Notation

Following the convention in **Graupner & Brunel (2012)** (*Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location*. PNAS, **109**(10), 3993. *Figure 2.C*), we define synaptic plasticity patterns as:

| **Pattern** | **Sequence** | **Description** |
|------------|-------------|----------------|
| **DP**   | `0, -1, 1, 0` | desired result |
| **DPD**  | `0, -1, 1, -1, 0` | Basal-Depression-Potentiation-Depression-Basal |
| **DPD'** | `-1, 1, -1` | Depression-Potentiation-Depression |
| **D**    | `0, -1, 0` | Basal-Depression-Basal |
| **D'**   | `-1` | All LTD |
| **P**    | `0, 1, 0` | Basal-Potentiation-Basal |
| **P'**   | `1` | All LTP |
| **All Basal** | `0, 0, 0` | All Basal |
| **PD**   | `0, 1, -1, 0` | Basal-Potentiation-Depression-Basal |


## notation on parameter tuning
If you want to implement the code and do some parameter analysis, the first thing you need to import **state_analysis.py**. Then you need import data from the txt file generated from the code. 

