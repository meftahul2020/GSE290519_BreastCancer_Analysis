# GSE290519_BreastCancer_Analysis
## **RNA-Seq Analysis of Breast Cancer (GSE290519) using DESeq2 🧬**  

### **Overview**  
This project analyzes RNA-Seq data from the **GSE290519** dataset to identify differentially expressed genes in breast cancer cell lines under **Dexamethasone (DEX) treatment**. The analysis was conducted using **R** and the **DESeq2** package for differential gene expression analysis.  

---

### **Dataset Information 📊**  
- **Dataset:** [GSE290519](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290519) (Breast Cancer Cell Line RNA-Seq)  
- **Treatment Groups:** Vehicle vs. DEX-treated  
- **Objective:** Identify genes significantly affected by **Dexamethasone (DEX)** treatment  

---

### **Methods & Tools 🛠**  
- **Programming Language:** R  
- **Packages Used:**  
  - `GEOquery` (Retrieve GEO dataset)  
  - `DESeq2` (Differential gene expression analysis)  
  - `ggplot2` (Data visualization)  
  - `pheatmap` (Heatmap visualization)  
  - `tidyverse` (Data manipulation)  

---

### **Project Workflow 🚀**  
1. **Data Retrieval**  
   - Downloaded RNA-Seq count data from **GEO (GSE290519)** using `GEOquery`  
   - Extracted and processed metadata for sample grouping  

2. **Data Preprocessing**  
   - Transformed count data into a DESeq2-compatible format  
   - Ensured metadata consistency with expression data  

3. **Differential Expression Analysis (DESeq2)**  
   - Performed normalization and variance stabilization  
   - Identified differentially expressed genes (DEX vs. Vehicle)  
   - Filtered significant genes (**p-adjusted < 0.05 & |log2FC| > 1**)  

4. **Visualization & Insights**  
   - **MA Plot**: Overall gene expression changes  
   - **Volcano Plot**: Highlighting significantly different genes  
   - **Heatmap**: Top 50 differentially expressed genes  
   - **PCA Plot**: Sample clustering based on expression profiles  

---

### **Results & Key Findings 📌**  
- Several genes showed significant expression changes after **DEX treatment**  
- **Top 50 genes** were identified and visualized in a heatmap  
- **Principal Component Analysis (PCA)** confirmed clear clustering between treatment groups  

---

### **How to Run the Analysis 🔧**  
#### **1️⃣ Clone this repository**  
```bash
git clone https://github.com/your-username/GSE290519_BreastCancer_Analysis.git
cd GSE290519_BreastCancer_Analysis
```
#### **2️⃣ Install Required Packages in R**  
```r
install.packages(c("dplyr", "readr", "tidyverse", "ggplot2", "pheatmap"))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2"))
```
#### **3️⃣ Run the R Script**  
```r
source("GSE290519_Analysis.R")
```
#### **4️⃣ View Results**  
- Check `DESeq2_results.csv` for full analysis  
- Open `plots/` folder for visualization outputs  

---

### **Future Work & Improvements 🔬**  
✅ **Expand dataset scope** to include more treatment conditions  
✅ **Functional enrichment analysis** (GO, KEGG pathways) for biological interpretation  
✅ **Integration with machine learning** for biomarker discovery  

---

### **Contributions & Contact 📩**  
If you're interested in collaborating or have suggestions, feel free to open an issue or connect with me on [LinkedIn](https://www.linkedin.com/in/meftahul-islam/)!  

🔗 **GitHub Repository:** https://github.com/meftahul2020/GSE290519_BreastCancer_Analysis

📧 **Email:** mefta.kuet.bme.18@gmail.com  
